#!/usr/bin/env python3 

import argparse, os, sys, snakemake
import pandas as pd
import numpy as np
from Bio import SeqIO
from time import gmtime, strftime
from utils import get_time_interval

cwd = os.getcwd()
filedir = os.path.abspath(os.path.dirname(__file__))
refdir = os.path.join(filedir,"..","data","reference")
dropdir = os.path.join(filedir, "..","data","to_drop") #gisaid ID of isolates that should not be included

def ArgumentParser():
    """Argument parser"""

    parser = argparse.ArgumentParser(prog = "treason.py", #changing this later
        formatter_class = argparse.RawTextHelpFormatter,
        description = "Filter clinical A/H3N2 samples from GISAID and construct phylogenetic tree per season (and previous sseason)" ) 
    
    parser.add_argument('-d','--data-folder', required=True, action="store", type=str, help="data folder where raw GISAID sequence and metadata are stored")
    parser.add_argument('-o','--output', required=False, action="store", type=str, help="output directory to store all generated output file (default: ./)")
    parser.add_argument('-s','--segments', required=False, nargs='+', action="store", type=str, help="segment abbreviations that will be analyzed")
    parser.add_argument('-st','--subtype', required=False, action="store", type=str, default="H3N2", choices=["H3N2", "H1N1pdm"], help="subtype to be analyzed (default: A/H3N2)")

    parser.add_argument('-p','--protein', required=False, action="store_true", help="if sequences need be translated into protein sequences (coding region only) and if mutations in final LBI need to be translated")
    parser.add_argument('-on','--only-nonsyn', required=False, action="store_true", help="if '-p' flag is specified, only report non-synonymous mutations in tree files")
    
    parser.add_argument('-mxa','--max-ambig', required=False, action="store", type=float, default=0.01, help="maximum percentage of ambiguous nucleotides allowed (default: 0.01)")
    parser.add_argument('-ml', '--min-length', required=False, action="store", type=float, default=0.95, help="miminum percentage of length w.r.t. the reference segment")
   
    parser.add_argument('-i','--interval', required=False, action="store", type=str, default='2015-2019', 
                        help="interval period in years form which individual season analyses need to be made (default: 2015-2019)")
    parser.add_argument('-sm','--start-month', required=False, action="store", type=str, default="may", help="name of the month from which the season should start (default: may)")
    parser.add_argument('-ss','--sub-sample',required=False, action="store_true", help="If season needs to be down sample for a max number of sequences/country/month (specify '-pcm')")
    parser.add_argument('-pcm', '--per-country-month', required=False, action="store", type=int, default=10, help="maximum number of sequence per country per month (default: 10)")
    parser.add_argument('-ips','--include-previous-season', required=False, action="store_true", help="if the sequence of the previous season should be included in the next season")
    
    parser.add_argument('-icb','--include-cell-based', required=False, action="store_true", help="include cell based isolates there are less than 500 (clinical) sequences per season")
    parser.add_argument('-ieb','--include-egg-based', required=False, action="store_true", 
                        help="Include cell based isolates there are less than 500 (clinical) sequences per season (if icb is specified cell-based isolates will be prioritized)")
    parser.add_argument('-cs', '--complete-subset', required=False, action="store_true", 
                        help="included sequences in the order of clinical, cell-based, egg-based, remaining if the number sequences in the season subset are less than 500")

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
    if args.threads:
        threads = args.threads
    else:
        threads = 1
    print("Number of threads is: ", threads) 

    #get interval of interest
    if "-" in args.interval and len(args.interval)==9:
        timeframe = args.interval
        full_intervals = get_time_interval(timeframe, args.start_month)
        time_intervals = []
        for i in full_intervals:
            time_intervals.append(f"{''.join(i[0].split('-')[::-1])}_{''.join(i[-1].split('-')[::-1])}")
    else:
        sys.stderr.write(f"Error: interval flag is not properly specified. Specify as: yyyy-yyyy with years being at least one year apart.")#\nCheck if installed properly")
        sys.exit(-1)
    
    #get season from time_interval
    seasons = []
    for i in time_intervals:
        season = ""
        for y in i.split("_"):
            season += "".join(list(y)[-2:])
        seasons.append(season)

    #if previous season need to be included add icp for every season but the first
    if args.include_previous_season:
        seasons.sort()
        full_seasons = []
        for i, s in enumerate(seasons):
            full_seasons.append(s)
            if i != 0:
                seasons[i-1]
                full_seasons.append("".join(list(seasons[i-1])[:2]) + s)
    else: #full_season equals seasons > this is to avoid snakemake issues
        full_seasons = seasons.copy()

    #get segments of interest
    all_segments  = ['PB2','PB1','PA','HA','NP','NA','M','NS']
    if args.segments:
        segments = [seg for seg in args.segments if seg in all_segments] 
    else:
        segments = all_segments
    if len(segments) == 0:
        sys.stderr.write(f"Error: specified segments are not recognized. Choose from {', '.join(all_segments)}")#\nCheck if installed properly")
        sys.exit(-1)

    #get subtype
    subtype = args.subtype

    #get data folder
    datafolder = os.path.join(cwd, args.data_folder)
    if not os.path.isdir (datafolder):
        sys.stderr.write(f"Error: cannot find{datafolder}. Check if the directory is properly specified.")
        sys.exit(-1)

    #get to drop files
    to_drop_files = []
    for f in os.listdir(dropdir):
        if f.split("_")[0] in segments:
            to_drop_files.append(os.path.join(dropdir,f))

    #get the output directory
    if args.output:
        output = os.path.join(cwd,args.output)
        if not os.path.isdir(output):
            try:
                os.mkdir(output)
            except:
                try:
                    os.mkdir("/".join(output.split("/")[:-1]))
                    os.mkdir(output)
                except:
                    sys.stderr.write(f"Error: cannot create {output}. Try creating the output directory manually and run again")
                    sys.exit(-1)
    else:
        output = cwd

    #find metadata file(s) > combined for all files -- and see if this needs to be merged
    #and find fasta file(s) per segments -- also see if these need to be merged

    #for metadata only getting the columns of interest
    metcols = ["Isolate_Id", "PB2 Segment_Id", "PB1 Segment_Id", "PA Segment_Id", "HA Segment_Id", "NP Segment_Id",
               "NA Segment_Id", "MP Segment_Id", "NS Segment_Id", "Isolate_Name", "Passage_History", "Location", "Collection_Date"]

    if not os.path.isfile(f"{output}/raw/{subtype}_metadata_gisaid_raw.csv") or args.redo_all:
        print ("Merging Raw GISAID sequences into a single file per segment")
        print (strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        segment_records= {}
        for f in os.listdir(datafolder):
            if os.path.isdir(os.path.join(datafolder,f)):
                if f in segments:
                    segment_records[f] = []
                    for f2 in os.listdir(os.path.join(datafolder,f)):
                        if f2.endswith(".fasta"):
                            segment_records[f].extend(list(SeqIO.parse(os.path.join(datafolder,f,f2),"fasta")))
            elif os.path.isfile(os.path.join(datafolder,f)) and f.endswith(".fasta"):
                segid = list(set(f.split("_")).intersection(segments))[0]
                if len(segid) > 0:
                    if segid not in segment_records.keys():
                        segment_records[segid] = []
                    segment_records[segid].extend(list(SeqIO.parse(os.path.join(datafolder,f),"fasta")))
            elif os.path.isfile(os.path.join(datafolder,f)) and (f.endswith(".xls") or f.endswith(".xlsx")):
                try:
                    metadf = pd.concat([metadf, pd.read_excel(os.path.join(datafolder,f),usecols=metcols)])
                except:
                    metadf = pd.read_excel(os.path.join(datafolder,f),usecols=metcols)  
    
        #filter out segments for which no data is present > otherwise errors will occurs
        segments = list(set(segments).intersection(segment_records.keys()))  

        #creating sequence merge files of all segments in the output folder
        rawdir = os.path.join(output, "raw")
        if not os.path.isdir(rawdir):
            os.mkdir(rawdir)

        #write merge fasta files 
        for segment, records in segment_records.items():
            print (f"creating merge file for {segment} sequences")
            print (strftime("%Y-%m-%d %H:%M:%S", gmtime()))
            
            #remove duplicate records
            seen = set()
            unique = []
            for record in records:
                if record.id not in seen:
                    seen.add(record.id)
                    unique.append(record)

            #write output 
            merge_file = os.path.join(rawdir, f"{subtype}_{segment}_gisaid_raw.fasta")
            with open(merge_file, 'w') as fw:
                SeqIO.write(unique, fw, "fasta")
            print (f"finished: {strftime('%Y-%m-%d %H:%M:%S', gmtime())}")

        #write merge metadata file to rawdata folder
        print ("creating merge file for metadata")
        print (strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        try:
            metadf.to_csv(os.path.join(rawdir, f"{subtype}_metadata_gisaid_raw.csv"), index=False)
        except:
            sys.stderr.write(f"Could not find metadata file(s) in {datafolder}. Check if metadata is in correct directory and in correct format")
            sys.exit(-1)
        print (f"finished: {strftime('%Y-%m-%d %H:%M:%S', gmtime())}")

    #get snakefile 
    snakefile = os.path.join(filedir, "treason.smk")
    if not os.path.isfile(snakefile):
        sys.stderr.write(f"Error: cannot find {snakefile}.")#\nCheck if installed properly")
        sys.exit(-1)

    
    #setup snakemake config
    config = {
        "segments":segments,
        "subtype":subtype,
        "output":output,
        "refdir":refdir,
        "translate":args.protein,
        "o_nonsyn":args.only_nonsyn,
        "mxa": args.max_ambig,
        "mnlp":args.min_length,
        "timeframe": timeframe,
        "to_drop_files":to_drop_files,
        "time_intervals":time_intervals,
        "seasons":seasons,
        "full_seasons":full_seasons,
        "start_month": args.start_month,
        "subsample":args.sub_sample,
        "pcm": args.per_country_month,
        "icp":args.include_previous_season,
        'icb':args.include_cell_based,
        "ieb":args.include_egg_based,
        "complete_subset":args.complete_subset,
        "seed": args.seed,
        "no_seed":args.no_seed,
    }

    quite_mode = args.verbose
    
    redo = True if args.redo or args.redo_all else False

    treachery_rules = [ "RedoMSA","RedoPhyloTree","RedoTreeTime"]
    if args.force_rule is not None:
        if type(args.force_rule) == str:
            args.force_rule = [args.force_rule]
        force_rule_now = [rule for rule in args.force_rule if rule not in treachery_rules]
        force_rule_later = [rule for rule in args.force_rule if rule in treachery_rules]
    else:
        force_rule_now = []
        force_rule_later = []
    config["force_rule_later"] = force_rule_later


    #start snakemake
    if len(force_rule_now) > 0:
        snakemake.snakemake(snakefile, printshellcmds=True, forceall=redo, forcerun=force_rule_now, config=config, cores=threads, lock=False, latency_wait=15,
                            quiet=quite_mode)
    else:
        snakemake.snakemake(snakefile, printshellcmds=True, forceall=redo, config=config, cores=threads, lock=False, latency_wait=15,
                            quiet=quite_mode)


if __name__ == "__main__":
    main()

