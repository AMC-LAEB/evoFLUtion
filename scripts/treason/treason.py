#!/usr/bin/env python3 

import argparse, os, sys, snakemake
import pandas as pd
import numpy as np
from Bio import SeqIO
from time import gmtime, strftime
from utils import check_max_ambig, check_sequence_length
from labels import m2n, metcols

cwd = os.getcwd()
filedir = os.path.abspath(os.path.dirname(__file__))
refdir = os.path.join(filedir,"..", "..","data","reference")
dropdir = os.path.join(filedir, "..", "..","data","to_drop") #gisaid ID of isolates that should not be included

def ArgumentParser():
    """Argument parser"""

    parser = argparse.ArgumentParser(prog = "treason.py", #changing this later
        formatter_class = argparse.RawTextHelpFormatter,
        description = "Filter clinical A/H3N2 samples from GISAID and construct phylogenetic tree per season (and previous sseason)" ) 
    
    parser.add_argument('-d','--data-folder', required=True, action="store", type=str, help="data folder where raw GISAID sequence and metadata are stored")
    parser.add_argument('-o','--output', required=False, action="store", type=str, help="output directory to store all generated output file (default: ./)")
    parser.add_argument('-s','--segment', required=False, action="store", type=str, default="HA", choices=["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"], help="segment abbreviations that will be analyzed")
    parser.add_argument('-st','--subtype', required=False, action="store", type=str, default="H3N2", choices=["H3N2", "H1N1pdm"], help="subtype to be analyzed (default: A/H3N2)")

    parser.add_argument('-p','--protein', required=False, action="store_true", 
                        help="if sequences need be translated into protein sequences (coding region only) and if mutations in final LBI need to be translated")
    parser.add_argument('-on','--only-nonsyn', required=False, action="store_true", help="if '-p' flag is specified, only report non-synonymous mutations in tree files")
    
    parser.add_argument('-mxa','--max-ambig', required=False, action="store", type=float, default=0.01, help="maximum percentage of ambiguous nucleotides allowed (default: 0.01)")
    parser.add_argument('-ml','--min-length', required=False, action="store", type=float, default=0.95, help="miminum percentage of length w.r.t. the reference segment")
   
    #parser.add_argument() #seperate hemispheres
    parser.add_argument('-sh','--separate-hemispheres', required=False, action="store_true", help="perform analysis for northern and southern hemisphere separately")
    parser.add_argument('-i','--interval', required=False, action="store", type=str, default='2015-2019', 
                        help="interval period in years form which individual season analyses need to be made (default: 2015-2019)")
    parser.add_argument('-is','--interval-start', required=False, action="store", type=str, default="Southern",
                        choices=["Northern", "Southern", "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"],
                        help="specify the desired start of the time frame (Northern hemisphere flu season starts in October and southern Hemisphere seasons starts in May) (ignored when '-sh') (default:Southern)")
    parser.add_argument('-ip','--interval-period', required=False, action="store", type=int, default=12, help="number of months the period takes. max=12 (default=12)")
    #parser.add_argument('-sm','--start-month', required=False, action="store", type=str, default="may", help="name of the month from which the season should start (default: may)")

    parser.add_argument('-ss','--sub-sample',required=False, nargs="?", const=10, action="store",type=int, 
                        help="If season needs to be down sample for a max number of sequences/country/month (if specified default:10)")
    
    parser.add_argument('-co', '--clinical-only', required=False, action="store_true", help="whether only (direct) clinical sample should be included in the analysis")

    parser.add_argument('-c','--seed', required=False, action="store", type=int, default=29, help="seed to use for random sampling (default:29)")
    parser.add_argument('-nc','--no-seed', required=False, action="store_true", help="If flag is specified no seed will be used")

    parser.add_argument('-r','--redo', required=False, action="store_true", help="Redo all steps except the initial merging step if data already exists")
    parser.add_argument('-ra','--redo-all', required=False, action="store_true", help="Redo ALL steps if data already exists")
    parser.add_argument('-t','--threads',required=False,action='store',type=int,help="Number of threads (default 1)")
    parser.add_argument('-f','--force-rule', required=False,action='store',choices=["FilterClinical","GetSeason","MSA","PhyloTree","TreeTime","FilterMolecularClockOutliers",
                                                                                    "RedoMSA","RedoPhyloTree","RedoTreeTime", "LBI", "TranslateMutations", "Translate"],
                        help="Force execution of specific snakemake rule")

    parser.add_argument('-v', '--verbose', required=False,action="store_false",help="print a bunch of stuff")

    #check if arguments were entered 
    if len(sys.argv[2:])<1: 
        parser.print_help()
        sys.exit(-1)
    else:
        return parser.parse_args()

def main():
    args = ArgumentParser()

    #get threads
    threads = args.threads if args.threads else 1
    print("Number of threads is: ", threads) 

    #get first year in time frame and last year in time frame 
    if "-" in args.interval and len(args.interval)==9:
        tfs, tfe = int(args.interval.split("-")[0]), int(args.interval.split("-")[-1])
    else:
        sys.stderr.write(f"Error: interval flag is not properly specified. Specify as: yyyy-yyyy with years being at least one year apart.")
        sys.exit(-1)

    #get interval start 
    ins = args.interval_start
    ins = 10 if ins == "Northern" else 5 if ins == "Southern" else m2n[args.interval_start]#interval stars
    
    #get interval period
    ip = args.interval_period
    if ip >12 or ip <=0:
        sys.stderr.write(f"Error: interval period is not properly specified. min=1 month, max=12 months")
        sys.exit(-1)

    if args.separate_hemispheres:
        periods = []
        for h, m in {"sh":5, "nh":10}.items():
            pys, pms = tfs, m
            pye, pme = pys if +ip <= 12 else pys+1, pms+ip if (pms+ip) <= 12 else (pms+ip)-12
            while pye <= tfe:
                periods.append(f"{h}-{'0' +str(pms) if pms<10 else str(pms)}{str(pys)[-2:]}-{'0' +str(pme-1) if pme-1<10 else str(pme-1)}{str(pye)[-2:]}")
                pys, pms = pye, pme
                pye, pme = pys if ins+ip <= 12 else pys+1, pms+ip if (pms+ip)<= 12 else (pms+ip)-12
    else:
        #if not per hemisphere
        pys, pms = tfs, ins #period year start, period month start
        pye, pme = pys if ins+ip <= 12 else pys+1, pms+ip if (pms+ip) <= 12 else (pms+ip)-12
        periods = []
        while pye <= tfe:
            periods.append(f"{'0' +str(pms) if pms<10 else str(pms)}{str(pys)[-2:]}-{'0' +str(pme-1) if pme-1<10 else str(pme-1)}{str(pye)[-2:]}")
            pys, pms = pye, pme
            pye, pme = pys if ins+ip <= 12 else pys+1, pms+ip if (pms+ip)<= 12 else (pms+ip)-12

    #get segments of interest
    segment = args.segment
    #get subtype
    subtype = args.subtype

    #get data folder
    datafolder = os.path.join(cwd, args.data_folder)
    if not os.path.isdir (datafolder):
        sys.stderr.write(f"Error: cannot find{datafolder}. Check if the directory is properly specified.")
        sys.exit(-1)

    #get to drop files
    for f in os.listdir(dropdir):
        if f.split("_")[1] == segment and f.split("_")[0] == subtype:
            to_drop_file = os.path.join(dropdir,f)
    
    #get the output directory
    if args.output:
        output = os.path.join(cwd, args.output)
        os.makedirs(output, exist_ok=True)
    else:
        output = cwd

    #for metadata only getting the columns of interest
    
    if not os.path.isfile(f"{output}/sequences/{subtype}_metadata_gisaid.csv") or args.redo_all:
        print (f"Merging Raw GISAID sequences into a single file per segment: {strftime('%Y-%m-%d %H:%M:%S', gmtime())}")
        segment_records= {}
        for f in os.listdir(datafolder):
            if os.path.isdir(os.path.join(datafolder,f)) and f == segment:
                segment_records = [record for f2 in os.listdir(os.path.join(datafolder, f)) if f2.endswith(".fasta") \
                                  for record in SeqIO.parse(os.path.join(os.path.join(datafolder, f), f2), "fasta")]
            
            #load sequences
            elif os.path.isfile(os.path.join(datafolder,f)):
                if f.endswith(".fasta") and segment in f:
                    segment_records = list(SeqIO.parse(os.path.join(datafolder,f),"fasta"))
                
                #read metadata
                elif f.endswith(".xls") or f.endswith(".xlsx") or f.endswith(".csv"):
                    try:
                        if f.endswith(".xls") or f.endswith(".xlsx"):
                            metadf = pd.concat([metadf, pd.read_excel(os.path.join(datafolder,f),usecols=metcols)])
                        else: 
                            metadf = pd.concat([metadf, pd.read_csv(os.path.join(datafolder, f))], usecols=metcols)
                    except:
                        if f.endswith(".xls") or f.endswith(".xlsx"):
                            metadf = pd.read_excel(os.path.join(datafolder,f),usecols=metcols) 
                        else: 
                            metadf =  pd.read_csv(os.path.join(datafolder,f),usecols=metcols) 
    

        #filter out all sequences for which there is no complete collection date
        metadf["Collection_Date"] = pd.to_datetime(metadf["Collection_Date"], errors='coerce')
        metadf = metadf.dropna(subset=["Collection_Date"]).reset_index(drop=True)

        #also filter out sequences that have no country specification 
        metadf = metadf.drop(metadf[metadf["Location"].str.count("/") < 1].index).reset_index(drop=True)
        #creating sequence merge files of all segments in the output folder
        seqdir = os.path.join(output, "sequences")
        if not os.path.isdir(seqdir):
            os.mkdir(seqdir)

        #write merge fasta files 
        itr = set(pd.read_csv(to_drop_file)["Isolate_Id"]) #ids to remove
        print (f"creating merge file for sequences: {strftime('%Y-%m-%d %H:%M:%S', gmtime())}")
            
        #remove duplicate records and low quality records
        seen = set()
        unique = []

        for record in segment_records:
            record_id = record.id.split("|")[0]
            strain_name = record.id.split("|")[2]
                
            #check for duplicates
            if record_id not in itr and strain_name not in seen and record_id in metadf["Isolate_Id"].tolist() and record_id not in seen:
                    
                #check for low quality 
                if check_sequence_length(segment, str(record.seq), args.min_length) and check_max_ambig(str(record.seq), args.max_ambig): 
                    seen.update([record_id, strain_name])
                    unique.append(record)

        #write output 
        merge_file = os.path.join(seqdir, f"{subtype}_{segment}_gisaid.fasta")
        with open(merge_file, 'w') as fw:
            SeqIO.write(unique, fw, "fasta")
        print (f"finished: {strftime('%Y-%m-%d %H:%M:%S', gmtime())}")

        #write merge metadata file to rawdata folder
        print (f"creating merge file for metadata: {strftime('%Y-%m-%d %H:%M:%S', gmtime())}")
        try:
            metadf.to_csv(os.path.join(seqdir, f"{subtype}_metadata_gisaid.csv"), index=False)
        except:
            sys.stderr.write(f"Could not find metadata file(s) in {datafolder}. Check if metadata is in correct directory and in correct format")
            sys.exit(-1)
        print (f"finished: {strftime('%Y-%m-%d %H:%M:%S', gmtime())}")

    #get snakefile 
    snakefile = os.path.join(filedir, "treason.smk")
    if not os.path.isfile(snakefile):
        sys.stderr.write(f"Error: cannot find {snakefile}.")#\nCheck if installed properly")
        sys.exit(-1)

    subsample = args.sub_sample if args.sub_sample is not None else False

    #setup snakemake config
    config = {
        "segment":segment,
        "subtype":subtype,
        "output":output,
        "refdir":refdir,
        "translate":args.protein,
        "o_nonsyn":args.only_nonsyn,
        "periods": periods,
        "subsample":subsample,
        "clinical_only": args.clinical_only,
        "seed": args.seed,
        "no_seed":args.no_seed,
    }

    quite_mode = args.verbose
    
    redo = True if args.redo or args.redo_all else False
    
    force_rule= args.force_rule 

    #start snakemake
    if force_rule is not None and len(force_rule) > 0:
        snakemake.snakemake(snakefile, printshellcmds=True, forceall=redo, forcerun=force_rule, config=config, cores=threads, lock=False, latency_wait=15,
                            quiet=quite_mode)
    else:
        snakemake.snakemake(snakefile, printshellcmds=True, forceall=redo, config=config, cores=threads, lock=False, latency_wait=15,
                            quiet=quite_mode)


if __name__ == "__main__":
    main()

