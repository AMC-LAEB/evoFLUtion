import os, random, sys, copy, dendropy,re
from Bio import SeqIO, Seq, SeqRecord
import multiprocessing as mp 
import pandas as pd
import numpy as np
from time import gmtime, strftime
from datetime import datetime
from dateutil.relativedelta import relativedelta
from utils import *
from labels import *
import dendropy

############### VARS ####################
#molecular clock rates used for treetime determined for sequences from 2000 to 2019
clock_rates = {"H3N2":{'PB2':0.00227, 'PB1':0.00202, 'PA':0.00192, 'HA':0.00415, 'NP':0.00177, 'NA':0.00266, 'M':0.00190, 'NS':0.00167},
                "H1N1pdm":{'PB2':0.00277, 'PB1':.00230, 'PA':0.00279, 'HA':0.00357, 'NP':0.00221, 'NA':0.00326, }}

############### MAIN RULE ####################
if config["translate"]:
    rule all:
        input:
            prot = expand(f"{config['output']}/protein/{config['subtype']}_{config['segment']}_{{period}}_proteins.fasta",  period=config['periods']),
            final = expand(f"{config['output']}/protein/{config['subtype']}_{config['segment']}_{{period}}.nexus", period=config['periods']),
else:
    rule all:
        input:
            LBI = expand(f"{config['output']}/lbi/{config['subtype']}_{config['segment']}_{{period}}_LBI.nexus",  period=config['periods']),

############### RULES ####################
###### Nucleotides ######
rule FilterClinical:
    #also editing the header of the clinical sequence here for se
    input:
        fastas = f"{config['output']}/sequences/{config['subtype']}_{config['segment']}_gisaid.fasta", 
        metadata = f"{config['output']}/sequences/{config['subtype']}_metadata_gisaid.csv",
    params:
        subtype = config["subtype"],  
        segment = config['segment'],
        output_dir =config['output'],
        co = config["clinical_only"],

        #specifying these files as parameters as work around > these are output files if clinical only is not specified
        cell_based_sequences = f"{config['output']}/sequences/{config['subtype']}_{config['segment']}_gisaid_cell_based.fasta",
        cell_based = f"{config['output']}/sequences/{config['subtype']}_metadata_gisaid_cell_based.csv",
        egg_based_sequences = f"{config['output']}/sequences/{config['subtype']}_{config['segment']}_gisaid_egg_based.fasta",
        egg_based = f"{config['output']}/sequences/{config['subtype']}_metadata_gisaid_egg_based.csv",
        remaining_sequences = f"{config['output']}/sequences/{config['subtype']}_{config['segment']}_gisaid_remaining.fasta",
        remaining = f"{config['output']}/sequences/{config['subtype']}_metadata_gisaid_remaining.csv",
        
    output:
        sequences = f"{config['output']}/sequences/{config['subtype']}_{config['segment']}_gisaid_clinical.fasta",
        clinical = f"{config['output']}/sequences/{config['subtype']}_metadata_gisaid_clinical.csv",
    threads: 1
    message: "removing non-clinical isolates and previously reported outliers"
    priority: 1
    run:
        #read metadata and remove potential duplicates
        metadata = pd.read_csv(str(input.metadata)).drop_duplicates()

        #filter clinical sequences
        clinical = metadata[metadata["Passage_History"].isin(cpl)]
        
       
        print (f"Filtering clinical sequences: {strftime('%Y-%m-%d %H:%M:%S', gmtime())}")
        records = [record for record in SeqIO.parse(f"{params.output_dir}/sequences/{params.subtype}_{params.segment}_gisaid.fasta", "fasta") if record.id.split("|")[0] in set(clinical["Isolate_Id"])]
        if len(records) > 0:
            with open(f"{params.output_dir}/sequences/{params.subtype}_{params.segment}_gisaid_clinical.fasta", "w") as fw:
                SeqIO.write(records, fw, "fasta")
        print (f"Finished: {strftime('%Y-%m-%d %H:%M:%S', gmtime())}")
            
        #write output clinical metadata file
        clinical.to_csv(str(output.clinical), index=False)

        if not params.co: #if we're not filtering on clinical data
            #filter cell based
            cell_based = metadata[metadata["Passage_History"].isin(cbpl)].reset_index(drop=True) 
            print (f"Filtering cell-based sequences: {strftime('%Y-%m-%d %H:%M:%S', gmtime())}")
            records = [record for record in SeqIO.parse(f"{params.output_dir}/sequences/{params.subtype}_{params.segment}_gisaid.fasta", "fasta") if record.id.split("|")[0] in set(cell_based["Isolate_Id"])]
            if len(records) > 0:
                with open(f"{params.output_dir}/sequences/{params.subtype}_{params.segment}_gisaid_cell_based.fasta", "w") as fw:
                    SeqIO.write(records, fw, "fasta")
            print (f"Finished: {strftime('%Y-%m-%d %H:%M:%S', gmtime())}")
                
            cell_based.to_csv(params.cell_based, index=False)

            #filter egg based 
            egg_based = metadata[metadata["Passage_History"].isin(epl)].reset_index(drop=True) 
            print (f"Filtering egg-based sequences: {strftime('%Y-%m-%d %H:%M:%S', gmtime())}")
            records = [record for record in SeqIO.parse(f"{params.output_dir}/sequences/{params.subtype}_{params.segment}_gisaid.fasta", "fasta") if record.id.split("|")[0] in set(egg_based["Isolate_Id"])]
            if len(records) > 0:
                with open(f"{params.output_dir}/sequences/{params.subtype}_{params.segment}_gisaid_egg_based.fasta", "w") as fw:
                    SeqIO.write(records, fw, "fasta")
            print (f"Finished: {strftime('%Y-%m-%d %H:%M:%S', gmtime())}")
                
            egg_based.to_csv(params.egg_based, index=False)
            
            #filter remaining 
            remaining = metadata[~metadata["Passage_History"].isin(cpl+cbpl+epl)].reset_index(drop=True)
            print (f"Filtering remaining sequences: {strftime('%Y-%m-%d %H:%M:%S', gmtime())}")
            records = [record for record in SeqIO.parse(f"{params.output_dir}/sequences/{params.subtype}_{params.segment}_gisaid.fasta", "fasta") if record.id.split("|")[0] in set(remaining["Isolate_Id"])]
            if len(records) > 0:
                with open(f"{params.output_dir}/sequences/{params.subtype}_{params.segment}_gisaid_egg_based.fasta", "w") as fw:
                    SeqIO.write(records, fw, "fasta")
            print (f"Finished: {strftime('%Y-%m-%d %H:%M:%S', gmtime())}")
            remaining.to_csv(params.remaining, index=False)

rule Getperiods:
    input:
        #sequences = expand(f"{config['output']}/sequences/{config['subtype']}_{config['segment']}_gisaid_clinical.fasta", ),
        sequences = f"{config['output']}/sequences/{config['subtype']}_{config['segment']}_gisaid.fasta",
        clinical = f"{config['output']}/sequences/{config['subtype']}_metadata_gisaid_clinical.csv",
        metadata = f"{config['output']}/sequences/{config['subtype']}_metadata_gisaid.csv",

    params:
        output_dir = f"{config['output']}/sequences",
        seed = config['seed'], #in case seed is used
        no_seed = config['no_seed'], #whether or not seed needs to be used
        subsample = config['subsample'],
        periods = config["periods"],
        subtype = config["subtype"],
        segment = config["segment"],
        co = config["clinical_only"],

        #specifying these files as parameters as work around > these are output files if clinical only is not specified
        #cell_based_sequences = expand(f"{config['output']}/sequences/{config['subtype']}_{config['segment']}_gisaid_cell_based.fasta", ),
        cell_based = f"{config['output']}/sequences/{config['subtype']}_metadata_gisaid_cell_based.csv",
        #egg_based_sequences = expand(f"{config['output']}/sequences/{config['subtype']}_{config['segment']}_gisaid_egg_based.fasta", ),
        egg_based = f"{config['output']}/sequences/{config['subtype']}_metadata_gisaid_egg_based.csv",
        #remaining_sequences = expand(f"{config['output']}/sequences/{config['subtype']}_{config['segment']}_gisaid_remaining.fasta", ),
        remaining = f"{config['output']}/sequences/{config['subtype']}_metadata_gisaid_remaining.csv",
    output:
        sequences = expand(f"{config['output']}/sequences/{config['subtype']}_{config['segment']}_{{period}}.fasta", period=config['periods']),
        metadata = expand(f"{config['output']}/sequences/{config['subtype']}_metadata_{{period}}.csv", period=config['periods'])
    run:
        #set seed if required
        if not params.no_seed:
            random.seed(params.seed)

        #read metadata files
        metadata = pd.read_csv(input.metadata)
        clinical = pd.read_csv(input.clinical, parse_dates=["Collection_Date"])
        clinical["Country"] = pd.Series([i.split(" / ")[1] for i in clinical["Location"]]).replace(cs)
        dfs = [clinical]
        if not params.co:
            cell_based = pd.read_csv(params.cell_based, parse_dates=["Collection_Date"])
            cell_based["Country"] = pd.Series([i.split(" / ")[1] for i in cell_based["Location"]]).replace(cs)
            dfs.append(cell_based)
            egg_based = pd.read_csv(params.egg_based, parse_dates=["Collection_Date"])
            egg_based["Country"] = pd.Series([i.split(" / ")[1] for i in egg_based["Location"]]).replace(cs)
            dfs.append(egg_based)
            remaining = pd.read_csv(params.remaining, parse_dates=["Collection_Date"])
            remaining["Country"] = pd.Series([i.split(" / ")[1] for i in remaining["Location"]]).replace(cs)
            dfs.append(remaining)

        #get sequences as dict
        seqs = {r.id.split("|")[0]:r for r in SeqIO.parse(input.sequences, "fasta")}
        
        #for each period select sequences
        for period in params.periods:
            i2s = [] #ids to select

            #get dates from period ASSUMING 
            start_date = period.split("-")[-2].lstrip("0")
            start_date = datetime(int(f'20{start_date[-2:]}'), int(start_date[:-2]),1)
            end_date = period.split("-")[-1].lstrip("0")
            end_date = datetime(int(f'20{end_date[-2:]}'), int(end_date[:-2]), ldm[int(end_date[:-2])])

            spc = {} #sequences per country
            for df in dfs:
                
                if len(period.split("-")) > 2:#if per hemisphere
                    if period.split("-") == "nh": #northern hemisphere
                        subdf = df[df["Country"].isin(nhc)]
                    else: #southern hemisphere
                        subdf = df[df["Country"].isin(shc)]
                else:
                    subdf = df
                
                #filter on period data
                subdf = subdf[(subdf["Collection_Date"]>=start_date)&(subdf["Collection_Date"]<=end_date)]
    
                # if subsampling is required
                if params.subsample != False:
                    for country in subdf["Country"].unique():
                        if country not in spc.keys():
                            spc[country] = {}

                        m = start_date
                        while m < end_date:
                            subsubdf = subdf[(subdf["Collection_Date"]>=m)&(subdf["Collection_Date"]<=m+relativedelta(months=1))]
                            if m not in spc[country].keys():
                                spc[country][m] = []
    
                            if len (spc[country][m]) < params.subsample:
                                cids =subsubdf[subsubdf["Country"]==country]["Isolate_Id"].to_list()
                                ts = params.subsample - len(spc[country][m])#target size 
                                if ts > 0:
                                    try:
                                        spc[country][m].extend(random.sample(cids, ts))
                                    except:
                                        if len(cids) > 0:
                                            spc[country][m].extend(cids)

                            m += relativedelta(months=1)
                else: #not subsampling so select all ids
                    i2s.extend(subdf["Isolate_Id"].to_list())

            #add subsampled ids to i2s
            if len(spc) > 0:
                #i2s.extend([i for l in spc.values() for m,i in l])
                i2s.extend([i for d in spc.values() for l in d.values() for i in l])

            #write fasta 
            with open(f'{params.output_dir}/{params.subtype}_{params.segment}_{period}.fasta', 'w') as fw:
                SeqIO.write([r for rid, r in seqs.items() if rid in i2s], fw, "fasta")

            #write metadata
            md = metadata[metadata["Isolate_Id"].isin(i2s)]
            md.to_csv(f"{params.output_dir}/{params.subtype}_metadata_{period}.csv", index=False)

rule MSA:
    input:
        sequences =  f"{config['output']}/sequences/{config['subtype']}_{config['segment']}_{{period}}.fasta"
    params:
        refdir = config["refdir"],
        segment = config['segment'],
        subtype = config['subtype']
    output:
        msa = f"{config['output']}/alignment/{config['subtype']}_{config['segment']}_{{period}}_MSA.fasta"
    threads:  workflow.cores if workflow.cores < 3 else 3
    message: "Performing multiple sequence alignment for {wildcards.period} sequences"
    run: 
        #get reference
        for f in os.listdir(params.refdir):
            if params.segment in f and params.subtype in f:
                reference = os.path.join(params.refdir, f)

        #run mafft
        cmd = ['mafft', '--auto', '--thread', str(threads), '--keeplength', '--addfragments', input.sequences, reference,'>', output.msa]
        shell(" ".join(cmd))

        #remove ref > want ref for initial alignment but not relevant else
        records = list(SeqIO.parse(output.msa, "fasta"))
        with open(output.msa,"w")as fw:
            SeqIO.write(records[1:],fw,"fasta")

rule PhyloTree:
    input:
        msa = rules.MSA.output.msa
    params:
        model = "GTR",
        seed = config['seed'],
        no_seed = config['no_seed'],
        treefile = temp(f"{config['output']}/alignment/{config['subtype']}_{config['segment']}_{{period}}_MSA.fasta.treefile"),
    output:
        tree = f"{config['output']}/tree/{config['subtype']}_{config['segment']}_{{period}}.treefile",
    message: "Constructing phylogenetic maximum likelihood tree for {wildcards.period}"
    threads: workflow.cores if workflow.cores < 3 else 3
    run:
        #run IQTREE 
        if not params.no_seed:
            iqtree_command = ['iqtree', '-s', input.msa, '-B', '1000', '-alrt', '1000', '-m', params.model, 
                              '-nt', str(threads), '--seed', str(params.seed), '-redo']
        else:
            iqtree_command = ['iqtree', '-s', input.msa, '-B', '1000', '-alrt', '1000', '-m', params.model,
                              '-nt', str(threads), '-redo']
        
        shell(" ".join(iqtree_command))

        #move tree file
        shell (f'mv {params.treefile} {output.tree}')
        
        #remove all IQtree files but the actual treefile
        files_ex = [".bionj",".ckp.gz",".iqtree",".mldist",".model.gz",".contree",".log",".splits.nex",".uniqueseq.phy",".model", ".parstree"]
        for ex in files_ex:
            if os.path.isfile(input.msa + ex):
                os.remove(input.msa+ ex)

rule TreeTime:
    input:
        metadata = expand(f"{config['output']}/sequences/{config['subtype']}_metadata_{{period}}.csv",period=config['periods']),
        msa = rules.MSA.output.msa, 
        tree = rules.PhyloTree.output.tree,
    params: 
        molclock_rates = clock_rates,
        segment = config['segment'],
        subtype = config['subtype'],
        treetime_folder = temp(f"{config['output']}/treetime/{config['subtype']}_{config['segment']}_{{period}}_tt"),
        timetree = temp(f"{config['output']}/treetime/{config['subtype']}_{config['segment']}_{{period}}_tt/timetree.nexus"),
        divergence_tree = temp(f"{config['output']}/treetime/{config['subtype']}_{config['segment']}_{{period}}_tt/divergence_tree.nexus"),
        mol_clock = temp(f"{config['output']}/treetime/{config['subtype']}_{config['segment']}_{{period}}_tt/molecular_clock.txt"),
    output:
        temp_msa = temp(f"{config['output']}/treetime/{config['subtype']}_{config['segment']}_{{period}}.fasta"),
        dates = f"{config['output']}/treetime/{config['subtype']}_{config['segment']}_{{period}}_inputdates.csv",
        timetree = f"{config['output']}/treetime/{config['subtype']}_{config['segment']}_{{period}}_timetree.nexus",
        divergence_tree = f"{config['output']}/treetime/{config['subtype']}_{config['segment']}_{{period}}_divergence.nexus",
        treelog = f"{config['output']}/treetime/{config['subtype']}_{config['segment']}_{{period}}_treetime.txt",
        mol_clock = temp(f"{config['output']}/treetime/{config['subtype']}_{config['segment']}_{{period}}_molecular_clock.txt"),
        #dates_estimate = f"{config['output']}/treetime/{config['segment']}/dates.tsv",
    message: "Constructing time tree for {wildcards.period}"
    threads: 1 #no multi threading possible for treetime
    run: 
        #get the correct metadata file(s) and load metadata as pandas df
        for f in input.metadata:
            if wildcards.period in f:
                metadata = pd.read_csv(f, parse_dates=["Collection_Date"], index_col="Isolate_Id")
        
        #iqtree fucks up the ids in the msa so fix that
        ids = [i.taxon.label.replace(" ", "_") for i in dendropy.Tree.get(path=input.tree, schema="newick").leaf_node_iter()]
        ids = {i.split("|")[0]:i for i in ids}

        records = [SeqRecord.SeqRecord(r.seq, id=ids[r.id.split("|")[0]]) for r in SeqIO.parse(input.msa, "fasta")] 
        with open(output.temp_msa, "w") as fw:
            SeqIO.write(records, fw, "fasta")
        #print (ids["EPI_ISL_363969"] =="EPI_ISL_363969|1487224|A/Cote_D_Ivoire/14/2019|A_/_H3N2|HA")
        
        #generate dates file
        header_date = []
        for rid, header in ids.items():
            d = metadata.loc[rid, "Collection_Date"].date()
            fd = d.year + (d.timetuple().tm_yday/365)
            header_date.append([header, fd])
        
        dates = pd.DataFrame.from_records(header_date, columns=["accession", "date"])
        dates.to_csv(output.dates, index=False)   
        #get clock rates 
        clock_rate = params.molclock_rates[params.subtype][params.segment]
        clock_stdev = clock_rate*0.2

        #run treetime 
        cmd = ['treetime', '--tree', input.tree, '--aln', output.temp_msa, '--dates', output.dates, '--outdir', params.treetime_folder,
               '--clock-rate', str(clock_rate), '--clock-std-dev', str(clock_stdev)]
        #we need command line output > os tee as I also want to see what's happening
        cmd.extend(['|', 'tee', output.treelog])
        shell (" ".join(cmd))

        #mv treetime files of interest and then remove treetime folder > don't want to consume memory
        shell(f"mv {params.timetree} {output.timetree}")
        shell(f"mv {params.divergence_tree} {output.divergence_tree}")
        shell(f"mv {params.mol_clock} {output.mol_clock}")
        shell(f"rm -r -d {params.treetime_folder}")       
  
rule LBI:
    input:
        timetree = f"{config['output']}/treetime/{config['subtype']}_{config['segment']}_{{period}}_timetree.nexus"
    params:
        treeformat = "nexus",
        outputformat = "nexus",
        tau = 0.3, #might make these command line options
        normalize = True, #might make these command line options
    output:
        lbi_tree = f"{config['output']}/lbi/{config['subtype']}_{config['segment']}_{{period}}_LBI.nexus"
    message: "constructing LBI tree for {wildcards.period}"
    run:
        #read tree
        tree = dendropy.Tree.get(path=input.timetree, schema=params.treeformat)

        #prep tree
        tree = prep_tree_for_lbi(tree)

        #calculate LBI
        tree = calculate_LBI(tree, params.tau, normalize=params.normalize)

        #write output tree
        tree.write(path=output.lbi_tree, schema=params.outputformat)

###### Proteins ######
rule Translate:
    input:
        sequences = f"{config['output']}/sequences/{config['subtype']}_{config['segment']}_{{period}}.fasta",
        temp_msa = temp(f"{config['output']}/treetime/{config['subtype']}_{config['segment']}_{{period}}.fasta"),
    params:
        refdir = config["refdir"],
        segment = config['segment'],
    output:
        proteins = f"{config['output']}/protein/{config['subtype']}_{config['segment']}_{{period}}_proteins.fasta",
    message: "translating mcc sequences for {wildcards.period} into proteins"
    run:
        #get length of reference sequence
        for f in os.listdir(params.refdir):
            if params.segment in f:
                reference = os.path.join(params.refdir, f)
        ref_length = len(list(SeqIO.parse(reference,"fasta"))[0].seq)

        #get msa records 
        msa_recs = {}
        for record in SeqIO.parse(input.temp_msa, "fasta"):
            msa_recs[record.id] = record.seq
        
        proteins = []
        for record in SeqIO.parse(input.sequences, "fasta"):
            msaseq = msa_recs[record.id] 
            #get coding region from MSA sequence
            seq = Seq.Seq(str(msaseq).replace("-","n")).translate()
            #make protein record and add to  protein dict        
            protrec = SeqRecord.SeqRecord(seq, id=record.id, name=record.name, description=record.description)
            proteins.append(protrec)
        
        #write output file
        with open(output.proteins,"w") as fw:
            SeqIO.write(proteins, fw, "fasta")


rule TranslateMutations:
    input:
        #sequences = f"{config['output']}/clinical/{config['segment']}_{{period}}_mcc.fasta",
        temp_msa = temp(f"{config['output']}/treetime/{config['subtype']}_{config['segment']}_{{period}}.fasta"),
        tree = rules.LBI.output.lbi_tree,
    params:
        seed = config['seed'],
        no_seed = config['no_seed'],
        o_nonsyn = config['o_nonsyn'], 
        treeformat = "nexus", #lbi output format should not change
    output:
        tree = f"{config['output']}/protein/{config['subtype']}_{config['segment']}_{{period}}.nexus",
    message: "translating the mutations within the phylogenetic tree for {wildcards.period}"
    run:
        #making a dict of the isolates with the sequences > using MSA 
        isolates = {}
        for record in SeqIO.parse(input.temp_msa,"fasta"):
            if record.id not in isolates.keys():
                #only appending from start codon onwards as tree was build this way > important for numbering 
                isolates[record.id] = str(record.seq)
        
        #load the tree
        tree = dendropy.Tree.get(path=input.tree, schema=params.treeformat)

        #translate mutations
        if params.no_seed:
            tree = translate_tree_mutations(isolates,tree,params.o_nonsyn)
        else:
            tree = translate_tree_mutations(isolates,tree,params.o_nonsyn,seed=params.seed)
    
        #write tree to output file
        tree.write(path=output.tree, schema=params.treeformat)




