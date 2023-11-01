#!/usr/bin/env python3 

import argparse, os, sys, snakemake
import pandas as pd
import numpy as np
from Bio import SeqIO
from time import gmtime, strftime
from datetime import datetime
from utils import get_time_interval,get_dates,redo_fasta_header
import random
import statistics as st

cwd = os.getcwd()
filedir = os.path.abspath(os.path.dirname(__file__))
refdir = os.path.join(filedir,"..","data","reference")
dropdir = os.path.join(filedir, "..","data","to_drop") #gisaid ID of isolates that should not be included

#get season definitions > doing this manually
#script works with seasons from 2000-2001 till 2018-2019
seasons = {"0001":[datetime(2000,5,1), datetime(2001,4,30)], "0102":[datetime(2001,5,1), datetime(2002,4,30)], "0203":[datetime(2002,5,1), datetime(2003,4,30)],
           "0304":[datetime(2003,5,1), datetime(2004,4,30)], "0405":[datetime(2004,5,1), datetime(2005,4,30)], "0506":[datetime(2005,5,1), datetime(2006,4,30)],
           "0607":[datetime(2006,5,1), datetime(2007,4,30)], "0708":[datetime(2007,5,1), datetime(2008,4,30)], "0809":[datetime(2008,5,1), datetime(2009,4,30)],
           "0910":[datetime(2009,5,1), datetime(2010,4,30)], "1011":[datetime(2010,5,1), datetime(2011,4,30)], "1112":[datetime(2011,5,1), datetime(2012,4,30)],
           "1213":[datetime(2012,5,1), datetime(2013,4,30)], "1314":[datetime(2013,5,1), datetime(2014,4,30)], "1415":[datetime(2014,5,1), datetime(2015,4,30)], 
           "1516":[datetime(2015,5,1), datetime(2016,4,30)], "1617":[datetime(2016,5,1), datetime(2017,4,30)], "1718":[datetime(2017,5,1), datetime(2018,4,30)],
           "1819":[datetime(2018,5,1), datetime(2019,4,30)]}

def determine_season(d):
    """
    determine the flu seasons based on the date
    """
    ss = list(seasons.keys())
    season_starts = [v[0] for v  in seasons.values()]
    season_ends = [v[-1] for v in seasons.values()]
 
    if d > season_ends[-1] or d < season_starts[0]:
        return None

    for i, se in enumerate(season_ends):
        if d < se:
            return ss[i]

def ArgumentParser():
    """Argument parser"""

    parser = argparse.ArgumentParser(prog = "mol_clock_est.py", #changing this later
        formatter_class = argparse.RawTextHelpFormatter,
        description = "establish molecular clock rate for influenza segment of interest with certain season bound" ) 
    
    parser.add_argument('-d','--data-folder', required=True, action="store", type=str, help="data folder where raw GISAID sequence and metadata are stored")
    parser.add_argument('-o','--output', required=False, action="store", type=str, help="output directory to store all generated output file (default: ./)")
    parser.add_argument('-s','--segment', required=False, action="store", type=str, help="segment abbreviations that will be analyzed")
    parser.add_argument('-st','--subtype', required=False, action="store", type=str, default="H3N2", choices=["H3N2", "H1N1pdm"], help="subtype to be analyzed (default: A/H3N2)")

    parser.add_argument('-sr', '--season-range', required=False, action="store",type=str, default="2015-2019", help="season range of interest from which trees can be included in random subset (default:2015-2019)")
    parser.add_argument('-nr', '--replicate-number', required=False, action="store", type=int, default=1, help="number of replicates (random sample) for which molecular clock rates needs to be established (default:1)")
    parser.add_argument('-ss', '--subset-size', required=False, action="store", type=int, default=1000, help="number of sequences in every replicate (default: 1000)" )
    
    parser.add_argument('-r', '--redo', required=False, action="store_true", help="specify if analysis needs to be redone")
    parser.add_argument('-c','--seed', required=False, action="store", type=int, help="seed to use for random sampling (default:29)")
    parser.add_argument('-t','--threads',required=False,action='store',type=int,help="Number of threads (default 1)")

    parser.add_argument('-v', '--verbose', required=False,action="store_true",help="print a bunch of stuff")


    #check if arguments were entered 
    if len(sys.argv[2:])<1: 
        parser.print_help()
        sys.exit(-1)
    else:
        return parser.parse_args()

def main():
    args = ArgumentParser()

    segment = args.segment
    subtype = args.subtype
    n_subsets = args.replicate_number
    subset_size = args.subset_size

    #get threads
    if args.threads:
        threads = args.threads
    else:
        threads = 1
    print("Number of threads is: ", threads) 


    #gather fasta and metadata file(s) from data directory
    datafolder = os.path.join(cwd, args.data_folder)
    if not os.path.isdir (datafolder):
        sys.stderr.write(f"Error: cannot find{datafolder}. Check if the directory is properly specified.")
        sys.exit(-1)

    fasta_files = []
    metadata_files = []
    for f in os.listdir(datafolder):
        if f.endswith(".fasta"):
            fasta_files.append(os.path.join(datafolder,f))
        elif f.endswith(".csv") or f.endswith(".xlsx") or f.endswith(".xls"):
            metadata_files.append(os.path.join(datafolder,f))

    #exit if no fasta or metadata files are found
    if len(fasta_files) == 0:
        sys.stderr.write(f"Error: cannot find fasta files in {datafolder}. Check if the directory is properly specified and if sequence file(s) is in the correct format.")
        sys.exit(-1)
    if len(metadata_files) == 0:
        sys.stderr.write(f"Error: cannot find metadata files in {datafolder}. Check if the directory is properly specified and if the metadata file(s) is in the correct format.")
        sys.exit(-1)

    #get output dir
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

    #make sequence dir > both replicate files and their corresponding dates files
    replicate_dir = os.path.join(output, "replicates")
    if not os.path.isdir(replicate_dir):
        os.mkdir(replicate_dir)

    #read fasta file(s) and remove outliers 
    to_drop = pd.read_csv(os.path.join(dropdir, f"{subtype}_{segment}_gitr.csv"))["Isolate_Id"].to_list()
    records = {}
    for f in fasta_files:
        for r in SeqIO.parse(f, "fasta"):
            isolate = r.id.split("|")[0]
            if isolate not in to_drop:
                r.id = redo_fasta_header(str(r.id), segment)
                r.id = r.id.split("|")[0] + "|EPI" + "|".join(r.id.split("|")[1:])
                r.name = r.id
                r.description = r.id
                records[isolate] = r

    #get seasons of interest and create equal sampling per season
    try:
        start_season, end_season = args.season_range.split("-")
    except:
        sys.stderr.write(f"Error: season range is not properly specified. see help manual and try again.")
        sys.exit(-1)
    #instead of years get actual season specification
    if start_season.startswith("20"):
        try:
            start_season = f"{start_season.lstrip('20')}{str(int(start_season.lstrip('20'))+1)}"
        except:
            start_season = f"0{list(start_season)[-1]}0{int(list(start_season)[-1])+1}"
    if end_season.startswith("20"):
        end_season = f"{str(int(end_season.lstrip('20'))-1)}{end_season.lstrip('20')}"

    try:
        season_range = [s for s in seasons.keys() if list(seasons.keys()).index(s)>=list(seasons.keys()).index(start_season) and list(seasons.keys()).index(s)<=list(seasons.keys()).index(end_season)]
    except:
        sys.stderr.write(f"Error: season range is not properly specified. see help manual and try again.")
        sys.exit(-1)
    season_size = int(subset_size/len(season_range))
    subset_ids = [i for i in range(1,n_subsets+1)]

    #if redo is specified 
    if args.redo or len(os.listdir(replicate_dir)) < (n_subsets*2):
        #read metadata file(s) and remove outliers
        for f in metadata_files:
            df = pd.read_csv(f) if f.endswith(".csv") else pd.read_excel(f)
            try:
                metadata = pd.concat([metadata, df])
            except:
                metadata = df

        metadata = metadata[~metadata["Isolate_Id"].isin(to_drop)]
        metadata["Collection_Date"] = pd.to_datetime(metadata["Collection_Date"],format="%Y-%m-%d" )

        season_seqs = {s:[] for s in seasons.keys()}
        #determine season of each isolate
        for isolate in records.keys():
            d = str(np.datetime_as_string(metadata[metadata["Isolate_Id"]==isolate]["Collection_Date"].values[0]))
            d = datetime(int(d.split("-")[0]), int(d.split("-")[1]), int(d.split("-")[2].split("T")[0]))
            season = determine_season(d)
            if season is not None:
                season_seqs[season].append(isolate)

        for i in subset_ids:
            #generate random sample
            random_sample = []
            for season in season_range:
                if len(season_seqs[season])> season_size:
                    random_sample.extend([random.sample(season_seqs[season],season_size)])
                else:
                    random_sample.extend([season_seqs[season]])

            random_sample = [s for sl in random_sample for s in sl]
            selected_records = []
            for s in random_sample:
                selected_records.append(records[s])

            #write output file
            with open(os.path.join(replicate_dir, f"{subtype}_{segment}_mce_seqs_{i}.fasta"), "w") as fw:
                SeqIO.write(selected_records, fw, "fasta")
                
            #create date file 
            for f in metadata_files:
                df = pd.read_csv(f) if f.endswith(".csv") else pd.read_excel(f)
                try:
                    md = pd.concat([md, df])
                except:
                    md = df
            md = md[~md["Isolate_Id"].isin(to_drop)]
            sub_df = md[md["Isolate_Id"].isin(random_sample)]
            dates_file = os.path.join(replicate_dir, f"{subtype}_{segment}_mce_dates_{i}.csv")
            get_dates(segment,sub_df,dates_file)
    
    #get snakefile 
    snakefile = os.path.join(filedir, "mol_clock_est.smk")
    if not os.path.isfile(snakefile):
        sys.stderr.write(f"Error: cannot find {snakefile}.")#\nCheck if installed properly")
        sys.exit(-1)
    
    #setup snakemake config
    config = {
        "segment":segment,
        "subtype":subtype,
        "subset_ids":subset_ids,
        "output":output,
        "refdir":refdir,
        "seed": args.seed if type(args.seed)==int else False,
    }

    #quite_mode = args.verbose
    redo = True if args.redo else False

    #start snakemake
    snakemake.snakemake(snakefile, printshellcmds=True, forceall=redo, config=config, cores=threads, lock=False, )#quiet=quite_mode)

    #print results
    rates = []
    for f in os.listdir(os.path.join(output, "treetime")):
        if f.endswith("_treetime.txt"):
            with open(os.path.join(output, "treetime", f),"r") as fr:
                for l in fr:
                    if "--rate" in l:
                        rate = float(l.split("\t")[-1].strip("\n"))
                        rates.append(rate)
    if len(rates) > 1:
        print (f"estimated molecular clock rate is: {st.mean(rates)} for {n_subsets}")

if __name__ == "__main__":
    main()

