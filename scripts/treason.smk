import os, random, sys, copy, dendropy,re
from Bio import SeqIO, Seq, SeqRecord
import pandas as pd
import numpy as np
from time import gmtime, strftime
from utils import *

############### VARS ####################
#molecular clock rates used for treetime determined for sequences from 2015/01/01 to 2022/11/28
clock_rates = {'PB2':0.00219, 'PB1':0.00202, 'PA':0.00115, 'HA':0.00449, 'NP':0.00177, 'NA':0.00244, 'M':0.00190, 'NS':0.00167}

############### WILDCARDS ####################
wildcard_constraints:
    segment = "[A-Za-z0-9]{1,3}"


############### MAIN RULE ####################
if config["translate"]:
    rule all:
        input:
            prot = expand(f"{config['output']}/protein/{{segment}}_{{period}}_proteins.fasta", segment=config['segments'], period=config['full_seasons']),
            final = expand(f"{config['output']}/protein/{{segment}}_{{period}}.nexus", segment=config['segments'], period=config['full_seasons'])
else:
    rule all:
        input:
            LBI = expand(f"{config['output']}/lbi/{{segment}}/{{segment}}_{{period}}_LBI.nexus", segment=config['segments'], period=config['full_seasons'])

############### RULES ####################
###### Nucleotides ######
rule FilterClinical:
    #also editing the header of the clinical sequence here for se
    input:
        fastas = expand(f"{config['output']}/raw/{{segment}}_gisaid_raw.fasta", segment=config['segments']),
        metadata = f"{config['output']}/raw/metadata_gisaid_raw.csv"
    params:
        icb = config["icb"],
        cell_based_sequences = expand(f"{config['output']}/sequences/{{segment}}_gisaid_cell_based.fasta", segment=config['segments']),
        cell_based = f"{config['output']}/sequences/metadata_gisaid_cell_based.csv",
        ieb = config["ieb"],
        complete_subset = config["complete_subset"],
        egg_based_sequences = expand(f"{config['output']}/sequences/{{segment}}_gisaid_egg_based.fasta", segment=config['segments']),
        egg_based = f"{config['output']}/sequences/metadata_gisaid_egg_based.csv",
        remaining_sequences = expand(f"{config['output']}/sequences/{{segment}}_gisaid_remaining.fasta", segment=config['segments']),
        remaining = f"{config['output']}/sequences/metadata_gisaid_remaining.csv",
        segments = config['segments'],
        mnlp = config['mnlp'], #min length of the reference
        mxa = config['mxa'], #max % of ambiguous nucleotides allowed
        output_dir = f"{config['output']}/sequences",
        to_drop_files = config['to_drop_files'], #file for each segment with GISAID IDs that need to be removed 
    output:
        sequences = expand(f"{config['output']}/sequences/{{segment}}_gisaid_clinical.fasta", segment=config['segments']),
        clinical = f"{config['output']}/sequences/metadata_gisaid_clinical.csv"
    threads: workflow.cores if workflow.cores < len(config['segments']) else len(config['segments'])
    message: "removing non-clinical isolates and previously reported outliers"
    priority: 1
    run:
        #read metadata and remove potential duplicates
        metadata = pd.read_csv(input.metadata).drop_duplicates()

        #filter clinical sequences
        clinical, clinical_sequences = filter_clinical(metadata, input.fastas, ncpu=threads, rh=True, adth=False, mxa=params.mxa, mnlp=params.mnlp)

        #read to drop files and remove isolates that need to be dropped
        to_drop = []
        for f in params.to_drop_files:
            dropf = pd.read_csv(f)
            to_drop.extend(list(dropf["Isolate_Id"]))

        for segment, records in clinical_sequences.items():
            temp_records = []
            for record in records:
                if record.id.split("|")[0] not in to_drop:
                    temp_records.append(record)
            clinical_sequences[segment] = temp_records
        
        clinical = clinical[~clinical["Isolate_Id"].isin(to_drop)]

        #write output fastas for all clinical sequences
        for segment, records in clinical_sequences.items():
            with open(f"{params.output_dir}/{segment}_gisaid_clinical.fasta", "w") as fw:
                SeqIO.write(records, fw, "fasta")
    
        #write output clinical metadata file
        clinical.to_csv(output.clinical, index=False)

        #if cell based is requested repeat the process for cell based isolates
        if params.icb or params.complete_subset:
            cell_based, cell_based_sequences = filter_cell_based(metadata, input.fastas, ncpu=threads, rh=True, adth=False, mxa=params.mxa, mnlp=params.mnlp) 
        
            for segment, records in cell_based_sequences.items():
                temp_records = []
                for record in records:
                    if record.id.split("|")[0] not in to_drop:
                        temp_records.append(record)
                cell_based_sequences[segment] = temp_records
        
            cell_based = cell_based[~cell_based["Isolate_Id"].isin(to_drop)]

            #write output fastas for all cell based sequences
            for segment, records in cell_based_sequences.items():
                fname = f"{params.output_dir}/{segment}_gisaid_cell_based.fasta"
                with open(fname, "w") as fw:
                    SeqIO.write(records, fw, "fasta")
    
            #write output clinical metadata file
            cell_based.to_csv(params.cell_based, index=False)

        #if egg based is requested repeat the process for egg based isolates
        if params.ieb or params.complete_subset:
            egg_based, egg_based_sequences = filter_egg_based(metadata, input.fastas, ncpu=threads, rh=True, adth=False, mxa=params.mxa, mnlp=params.mnlp)

            for segment, records in egg_based_sequences.items():
                temp_records = []
                for record in records:
                    if record.id.split("|")[0] not in to_drop:
                        temp_records.append(record)
                egg_based_sequences[segment] = temp_records

            egg_based = egg_based[~egg_based["Isolate_Id"].isin(to_drop)]
            
            #write output fastsas for all egg based sequences
            for segment, records in egg_based_sequences.items():
                fname = f"{params.output_dir}/{segment}_gisaid_egg_based.fasta"
                with open(fname, "w") as fw:
                    SeqIO.write(records, fw, "fasta")
    
            #write output clinical metadata file
            egg_based.to_csv(params.egg_based, index=False)
        
        if params.complete_subset:
            remaining, remaining_sequences = filter_remaining(metadata, input.fastas, ncpu=threads, rh=True, adth=False, mxa=params.mxa, mnlp=params.mnlp)

            for segment, records in remaining_sequences.items():
                temp_records = []
                for record in records:
                    if record.id.split("|")[0] not in to_drop:
                        temp_records.append(record)
                remaining_sequences[segment] = temp_records

            remaining = remaining[~remaining["Isolate_Id"].isin(to_drop)]
            
            #write output fastsas for all egg based sequences
            for segment, records in remaining_sequences.items():
                fname = f"{params.output_dir}/{segment}_gisaid_remaining.fasta"
                with open(fname, "w") as fw:
                    SeqIO.write(records, fw, "fasta")
    
            #write output clinical metadata file
            remaining.to_csv(params.remaining, index=False)

rule GetSeason:
    input:
        sequences = expand(f"{config['output']}/sequences/{{segment}}_gisaid_clinical.fasta", segment=config['segments']),
        metadata = f"{config['output']}/sequences/metadata_gisaid_clinical.csv"
    params:
        output_dir = f"{config['output']}/sequences",
        seed = config['seed'], #in case seed is used
        no_seed = config['no_seed'], #whether or not seed needs to be used
        subsample = config['subsample'],
        timeframe = config['timeframe'],
        icp = config['icp'],
        sm = config['start_month'],
        spcpm = config['pcm'],
        icb = config['icb'],
        ieb = config['ieb'],
        complete_subset = config['complete_subset'],
        icb_cutoff = 500, 
        cell_based_sequences = expand(f"{config['output']}/sequences/{{segment}}_gisaid_cell_based.fasta", segment=config['segments']),
        cell_based = f"{config['output']}/sequences/metadata_gisaid_cell_based.csv",
        egg_based_sequences = expand(f"{config['output']}/sequences/{{segment}}_gisaid_egg_based.fasta", segment=config['segments']),
        egg_based = f"{config['output']}/sequences/metadata_gisaid_egg_based.csv",
        remaining_sequences = expand(f"{config['output']}/sequences/{{segment}}_gisaid_remaining.fasta", segment=config['segments']),
        remaining = f"{config['output']}/sequences/metadata_gisaid_remaining.csv",
    output:
        sequences = temp(expand(f"{config['output']}/sequences/{{segment}}_{{period}}.fasta", segment=config['segments'], period=config['full_seasons'])),
        metadata = temp(expand(f"{config['output']}/sequences/metadata_{{season}}.csv", season=config['seasons']))
    run:
        #read clinical data
        clinical = pd.read_csv(input.metadata)

        #get the full time intervals
        full_time_intervals = get_time_interval(params.timeframe, params.sm)
        
        #check if subsample is specified
        if params.subsample:
            if params.no_seed: 
                to_select_per_interval = select_per_year_per_month(clinical, full_time_intervals, params.spcpm, seed=None)
            else:
                to_select_per_interval = select_per_year_per_month(clinical, full_time_intervals, params.spcpm, seed=params.seed)
        else:
            to_select_per_interval = select_per_year(clinical, full_time_intervals)

        #check if include cell based is specified > if so check if select per interval is smaller than the cufoff
        if params.icb or params.complete_subset:
            #read cell based metadata
            cell_based = pd.read_csv(params.cell_based)
            for i, to_select in to_select_per_interval.items():
                idx = list(to_select_per_interval.keys()).index(i)
                #if len(to_select) < params.icb_cutoff:
                if params.no_seed: 
                    to_select.extend(extend_to_select_cell(clinical, cell_based, full_time_intervals[idx], to_select, params.spcpm, seed=None))
                else:
                    to_select.extend(extend_to_select_cell(clinical, cell_based, full_time_intervals[idx], to_select, params.spcpm, seed=params.seed))
                to_select_per_interval[i] = to_select
        
        #check if include egg based is specified > if so check if select per interal is smaller than the cutoff 
        if params.ieb or params.complete_subset:
            egg_based = pd.read_csv(params.egg_based)
            for i, to_select in to_select_per_interval.items():
                idx = list(to_select_per_interval.keys()).index(i)
                #check wheter or not include cell based is specified
                
                #if len(to_select) < params.icb_cutoff:
                if params.icb or params.complete_subset:  
                    if params.no_seed:
                        to_select.extend(extend_to_select_egg(clinical, egg_based, full_time_intervals[idx], cell_based, to_select, params.spcpm, seed=None))
                    else:
                        to_select.extend(extend_to_select_egg(clinical, egg_based, full_time_intervals[idx], cell_based, to_select, params.spcpm, seed=params.seed))
                else:
                    if params.no_seed:
                        to_select.extend(extend_to_select_egg(clinical, egg_based, full_time_intervals[idx], None, to_select, params.spcpm, seed=None))
                    else:
                        to_select.extend(extend_to_select_egg(clinical, egg_based, full_time_intervals[idx], None, to_select, params.spcpm, seed=params.seed))  
                to_select_per_interval[i] = to_select

        #check if subset needs to be completed 
        if params.complete_subset:
            remaining = pd.read_csv(params.remaining)
            for i, to_select in to_select_per_interval.items():
                idx = list(to_select_per_interval.keys()).index(i)

                #if len(to_select) < params.icb_cutoff:
                if params.no_seed:
                    to_select.extend(extend_to_select_remaining(clinical, cell_based, egg_based, remaining, full_time_intervals[idx],to_select, params.spcpm, seed=None))
                else:
                    to_select.extend(extend_to_select_remaining(clinical, cell_based, egg_based, remaining, full_time_intervals[idx],to_select, params.spcpm, seed=params.seed))
                
                to_select_per_interval[i] = to_select
        
        for i, to_select in to_select_per_interval.items():
            fe = "_".join(["".join(i.split("_")[0].split("-")[::-1]), "".join(i.split("_")[-1].split("-")[::-1])])
            season = ""
            for y in fe.split("_"):
                season += "".join(list(y)[-2:])

            #include previous season if requested 
            if params.icp:
                indx = list(to_select_per_interval.keys()).index(i)
                if indx > 0:
                    icp_to_select = to_select.copy()
                    previ = indx-1
                    icp_to_select.extend(list(to_select_per_interval.values())[previ])
                else:
                    #store first season as no icp tree will be made for this season
                    first_season = season
                    icp_to_select = []
            else: #to avoid problems
                icp_to_select = []

            for seqfile in input.sequences:
                segment = seqfile.split("/")[-1].split("_")[0]
                selected_records = []
                icp_records = []
                for record in SeqIO.parse(seqfile, "fasta"):
                    if record.id.split("|")[0] in to_select:
                        selected_records.append(record)
                    if record.id.split("|")[0] in icp_to_select:
                        icp_records.append(record)

                if params.icb or params.complete_subset:
                    cell_based_f = seqfile.replace("clinical", "cell_based")
                    for record in SeqIO.parse(cell_based_f, "fasta"):
                        if record.id.split("|")[0] in to_select:
                            selected_records.append(record)
                        if record.id.split("|")[0] in icp_to_select:
                            icp_records.append(record)

                if params.ieb or params.complete_subset:
                    egg_based_f = seqfile.replace("clinical", "egg_based")
                    for record in SeqIO.parse(egg_based_f, "fasta"):
                        if record.id.split("|")[0] in to_select:
                            selected_records.append(record)
                        if record.id.split("|")[0] in icp_to_select:
                            icp_records.append(record)
                
                if params.complete_subset:
                    remaining_f = seqfile.replace("clinical", "remaining")
                    for record in SeqIO.parse(remaining_f, "fasta"):
                        if record.id.split("|")[0] in to_select:
                            selected_records.append(record)
                        if record.id.split("|")[0] in icp_to_select:
                            icp_records.append(record)
                    
                if params.icp:
                    if season != first_season:
                        prev = "".join(list(list(to_select_per_interval.keys())[indx-1].split("-")[0])[2:])
                        with open(f"{params.output_dir}/{segment}_{prev+season}.fasta", "w")as fw:
                            SeqIO.write(icp_records, fw, "fasta")
                    with open(f"{params.output_dir}/{segment}_{season}.fasta", "w")as fw:
                        SeqIO.write(selected_records, fw, "fasta")
                else:
                    with open(f"{params.output_dir}/{segment}_{season}.fasta", "w")as fw:
                        SeqIO.write(selected_records, fw, "fasta")
                    
            subset = clinical[clinical["Isolate_Id"].isin(to_select)].reset_index(drop=True)
            if params.icb or params.complete_subset:
                cell_subset = cell_based[cell_based["Isolate_Id"].isin(to_select)].reset_index(drop=True)
                subset = pd.concat([subset,cell_subset])
            if params.ieb or params.complete_subset:
                egg_subset = egg_based[egg_based["Isolate_Id"].isin(to_select)].reset_index(drop=True)
                subset = pd.concat([subset, egg_subset])
            if params.complete_subset:
                remaining_subset = remaining[remaining["Isolate_Id"].isin(to_select)].reset_index(drop=True)
            subset.to_csv(f"{params.output_dir}/metadata_{season}.csv", index=False)

rule MSA:
    input:
        sequences =  f"{config['output']}/sequences/{{segment}}_{{period}}.fasta"
    params:
        refdir = config["refdir"]
    output:
        msa = temp(f"{config['output']}/alignment/{{segment}}/{{segment}}_{{period}}_MSA.fasta")
    threads:  workflow.cores if workflow.cores < 3 else 3
    message: "Performing multiple sequence alignment for {wildcards.segment} {wildcards.period} sequences"
    run: 
        #get reference
        for f in os.listdir(params.refdir):
            if wildcards.segment in f:
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
        treefile = temp(f"{config['output']}/alignment/{{segment}}/{{segment}}_{{period}}_MSA.fasta.treefile"),
    output:
        tree = temp(f"{config['output']}/tree/{{segment}}/{{segment}}_{{period}}.treefile"),
    message: "Constructing phylogenetic maximum likelihood tree for {wildcards.segment} {wildcards.period}"
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
        metadata = expand(f"{config['output']}/sequences/metadata_{{season}}.csv", season=config['seasons']),
        msa = rules.MSA.output.msa, 
        tree = rules.PhyloTree.output.tree,
    params: 
        molclock_rates = clock_rates,
        treetime_folder = temp(f"{config['output']}/treetime/{{segment}}_{{period}}_tt"),
        timetree = temp(f"{config['output']}/treetime/{{segment}}_{{period}}_tt/timetree.nexus"),
        divergence_tree = temp(f"{config['output']}/treetime/{{segment}}_{{period}}_tt/divergence_tree.nexus"),
        mol_clock = temp(f"{config['output']}/treetime/{{segment}}_{{period}}_tt/molecular_clock.txt"),
    output:
        dates = temp(f"{config['output']}/treetime/{{segment}}/{{segment}}_{{period}}_inputdates.csv"),
        timetree = temp(f"{config['output']}/treetime/{{segment}}/{{segment}}_{{period}}_timetree.nexus"),
        divergence_tree = temp(f"{config['output']}/treetime/{{segment}}/{{segment}}_{{period}}_divergence.nexus"),
        treelog = temp(f"{config['output']}/treetime/{{segment}}/{{segment}}_{{period}}_treetime.txt"),
        mol_clock = temp(f"{config['output']}/treetime/{{segment}}/{{segment}}_{{period}}_molecular_clock.txt"),
        #dates_estimate = f"{config['output']}/treetime/{{segment}}/dates.tsv",
    message: "Constructing time tree for {wildcards.segment} {wildcards.period}"
    threads: workflow.cores if workflow.cores < 3 else 3
    run: 
        #get the correct metadata file(s) and load metadata as pandas df
        #print (wildcards.period)
        if len(wildcards.period) > 4:
            season = "".join(list(wildcards.period)[2:])
            prevs = "".join(list(wildcards.period)[:4])
            for f in input.metadata:
                if season in f or prevs in f:
                    try: 
                        metadata = pd.concat([metadata, pd.read_csv(f)])
                    except:
                        metadata = pd.read_csv(f)
        else:
            for f in input.metadata:
                if wildcards.period in f:
                    metadata = pd.read_csv(f)
        
        #generate dates file
        get_dates(wildcards.segment, metadata, output.dates) 

        #get clock rates 
        clock_rate = params.molclock_rates[wildcards.segment]
        clock_stdev = clock_rate*0.2

        #run treetime 
        cmd = ['treetime', '--tree', input.tree, '--aln', input.msa, '--dates', output.dates, '--outdir', params.treetime_folder,
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
        timetree = f"{config['output']}/treetime/{{segment}}/{{segment}}_{{period}}_timetree.nexus"
    params:
        treeformat = "nexus",
        outputformat = "nexus",
        tau = 0.3, #might make these command line options
        normalize = True, #might make these command line options
    output:
        lbi_tree = f"{config['output']}/lbi/{{segment}}/{{segment}}_{{period}}_LBI.nexus"
    message: "constructing LBI tree for {wildcards.segment} {wildcards.period}"
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
        sequences = f"{config['output']}/sequences/{{segment}}_{{period}}.fasta",
        msa = f"{config['output']}/alignment/{{segment}}/{{segment}}_{{period}}_MSA.fasta",
    params:
        refdir = config["refdir"],
    output:
        proteins = f"{config['output']}/protein/{{segment}}_{{period}}_proteins.fasta",
    message: "translating mcc sequences for {wildcards.segment} {wildcards.period} into proteins"
    run:
        #get length of reference sequence
        for f in os.listdir(params.refdir):
            if wildcards.segment in f:
                reference = os.path.join(params.refdir, f)
        ref_length = len(list(SeqIO.parse(reference,"fasta"))[0].seq)

        #get msa records 
        msa_recs = {}
        for record in SeqIO.parse(input.msa, "fasta"):
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
        #sequences = f"{config['output']}/clinical/{{segment}}_{{period}}_mcc.fasta",
        msa = f"{config['output']}/alignment/{{segment}}/{{segment}}_{{period}}_MSA.fasta",
        tree = rules.LBI.output.lbi_tree,
    params:
        seed = config['seed'],
        no_seed = config['no_seed'],
        o_nonsyn = config['o_nonsyn'], 
        treeformat = "nexus", #lbi output format should not change
    output:
        tree = f"{config['output']}/protein/{{segment}}_{{period}}.nexus",
    message: "translating the mutations within the phylogenetic tree for {wildcards.segment} {wildcards.period}"
    run:
        #making a dict of the isolates with the sequences > using MSA 
        isolates = {}
        for record in SeqIO.parse(input.msa,"fasta"):
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




