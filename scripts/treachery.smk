import os, random, sys, copy, dendropy,re
from Bio import SeqIO, Seq, SeqRecord
import pandas as pd
import numpy as np
from time import gmtime, strftime
from utils import *

############## VARS ####################
#molecular clock rates used for treetime determined for sequences from 2015/01/01 to 2022/11/28
clock_rates = {'PB2':0.00219, 'PB1':0.00202, 'PA':0.00115, 'HA':0.00449, 'NP':0.00177, 'NA':0.00244, 'M':0.00190, 'NS':0.00167}

############### WILDCARDS ####################
wildcard_constraints:
    segment = "[A-Za-z0-9]{1,3}"

#get seasons and periods correctly
if type(config['periods']) == tuple:
    periods = list(config['periods'])
elif type(config['periods']) == str:
    periods = str(config['periods']).split(",") if "," in str(config['periods']) else "0"+str(config["periods"]) if len(str(config["periods"])) < 4 else str(config['periods'])
else: #type is list
    periods = config['periods']

if type(config['seasons']) == tuple:
    seasons = list(config['seasons'])
elif type(config['seasons']) == str:
    seasons = str(config['seasons']).split(",") if "," in str(config['seasons']) else "0"+str(config["seasons"]) if len(str(config["seasons"])) < 4 else str(config['seasons'])
else: #type is list
    seasons = config['seasons']


############### MAIN RULE ####################
rule all:
    input:
        #msa = expand(f"{config['output']}/alignment/{{segment}}/{{segment}}_{{period}}_MSA.fasta", segment=config['segments'], period=config['periods']),
        #tree = expand(f"{config['output']}/tree/{{segment}}/{{segment}}_{{period}}.treefile", segment=config['segments'], period=config['periods']),
        timetree = expand(f"{config['output']}/treetime/{{segment}}/{{segment}}_{{period}}_timetree.nexus", segment=config['segments'].split(","), period=periods),


rule RedoMSA:
    input:
        sequences = f"{config['output']}/sequences/{{segment}}_{{period}}.fasta"
        #seqeunces = rules.Sample.output.sequences
    params:
        refdir = config["refdir"]
    output:
        msa = f"{config['output']}/alignment/{{segment}}/{{segment}}_{{period}}_MSA.fasta"
    threads:  workflow.cores if workflow.cores < 3 else 3
    message: "Redoing multiple sequence alignment for {wildcards.segment} {wildcards.period} sequences"
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

rule RedoPhyloTree:
    input:
        msa = rules.RedoMSA.output.msa
    params:
        model = "GTR",
        seed = config['seed'],
        no_seed = config['no_seed'],
        treefile = temp(f"{config['output']}/alignment/{{segment}}/{{segment}}_{{period}}_MSA.fasta.treefile"),
    output:
        tree = f"{config['output']}/tree/{{segment}}/{{segment}}_{{period}}.treefile",
    message: "Reconstructing phylogenetic maximum likelihood tree for {wildcards.segment} {wildcards.period}"
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

rule RedoTreeTime:
    input:
        metadata = expand(f"{config['output']}/sequences/metadata_{{season}}.csv",season=seasons),
        msa = rules.RedoMSA.output.msa, 
        tree = rules.RedoPhyloTree.output.tree,
    params: 
        molclock_rates = clock_rates,
        treetime_folder = temp(f"{config['output']}/treetime/{{segment}}_{{period}}"),
        timetree = temp(f"{config['output']}/treetime/{{segment}}_{{period}}/timetree.nexus"),
        divergence_tree = temp(f"{config['output']}/treetime/{{segment}}_{{period}}/divergence_tree.nexus"),
    output:
        dates = f"{config['output']}/treetime/{{segment}}/{{segment}}_{{period}}_inputdates.csv",
        timetree = f"{config['output']}/treetime/{{segment}}/{{segment}}_{{period}}_timetree.nexus",
        divergence_tree = f"{config['output']}/treetime/{{segment}}/{{segment}}_{{period}}_divergence.nexus",
    message: "Reconstructing time tree for {wildcards.segment} {wildcards.period}"
    run: 
        #get the correct metadata file(s) and load metadata as pandas df
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
        shell (" ".join(cmd))

        shell(f"mv {params.timetree} {output.timetree}")
        shell(f"mv {params.divergence_tree} {output.divergence_tree}")
        shell(f"rm -r -d {params.treetime_folder}")