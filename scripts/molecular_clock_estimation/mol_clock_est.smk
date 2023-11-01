import os, random, sys, copy, dendropy,re
from Bio import SeqIO, Seq, SeqRecord
import pandas as pd
import numpy as np
from time import gmtime, strftime
from utils import *

rule all:
    input:
        treelog = expand(f"{config['output']}/treetime/{{subtype}}_{{segment}}_mce_{{i}}_treetime.txt",  subtype=config['subtype'], segment=config['segment'], i=config['subset_ids']),

rule MSA:
    input:
        sequences =  f"{config['output']}/replicates/{{subtype}}_{{segment}}_mce_seqs_{{i}}.fasta"
    params:
        refdir = config["refdir"],
        segment = config["segment"],
        subtype = config["subtype"]
    output:
        msa = f"{config['output']}/alignment/{{subtype}}_{{segment}}_mce_{{i}}_MSA.fasta"
    threads:  workflow.cores if workflow.cores < 3 else 3
    #message: f"Performing multiple sequence alignment for {{subtype}} {{segment}} replicate {{wildcards.i}}"
    run: 
        #get reference
        for f in os.listdir(params.refdir):
            if wildcards.segment in f and wildcards.subtype in f:
                reference = os.path.join(params.refdir, f)
                print ("testie bestie")

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
        treefile = temp(f"{config['output']}/alignment/{{subtype}}_{{segment}}_mce_{{i}}_MSA.fasta.treefile"),
    output:
        tree = f"{config['output']}/tree/{{subtype}}_{{segment}}_mce_{{i}}.treefile",
    #message: "Constructing phylogenetic maximum likelihood tree for  {{subtype}} {{segment}} replicate {wildcards.i}"
    threads: workflow.cores if workflow.cores < 3 else 3
    run:
        #run IQTREE 
        if params.seed is not False:
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
        dates =  f"{config['output']}/replicates/{{subtype}}_{{segment}}_mce_dates_{{i}}.csv",
        msa = rules.MSA.output.msa, 
        tree = rules.PhyloTree.output.tree,
    params: 
        treetime_folder = temp(f"{config['output']}/treetime/{{subtype}}_{{segment}}_mce_{{i}}_tt"),
        regression_pdf = temp(f"{config['output']}/treetime/{{subtype}}_{{segment}}_mce_{{i}}_tt/root_to_tip_regression.pdf"),
        #timetree = temp(f"{config['output']}/treetime/{{subtype}}_{{segment}}_mce_{{i}}_tt/timetree.nexus"),
        #divergence_tree = temp(f"{config['output']}/treetime/{{subtype}}_{{segment}}_mce_{{i}}_tt/divergence_tree.nexus"),
        #mol_clock = temp(f"{config['output']}/treetime/{{subtype}}_{{segment}}_mce_{{i}}_tt/molecular_clock.txt"),
    output:
        regression_pdf = f"{config['output']}/treetime/{{subtype}}_{{segment}}_mce_{{i}}_tt.pdf",
        #timetree = f"{config['output']}/treetime/{{subtype}}_{{segment}}_mce_{{i}}_timetree.nexus",
        #divergence_tree = f"{config['output']}/treetime/{{subtype}}_{{segment}}_mce_{{i}}_divergence.nexus",
        treelog = f"{config['output']}/treetime/{{subtype}}_{{segment}}_mce_{{i}}_treetime.txt",
        #mol_clock = f"{config['output']}/treetime/{{subtype}}_{{segment}}_mce_{{i}}_molecular_clock.txt",
        #dates_estimate = f"{config['output']}/treetime/{{segment}}/dates.tsv",
    #message: "estimating molecular clock rates from {{subtype}} {{segment}} replicate {wildcards.i}"
    threads: workflow.cores if workflow.cores < 3 else 3
    run: 

        #run treetime 
        cmd = ['treetime', 'clock','--tree', input.tree, '--aln', input.msa, '--dates', input.dates, '--outdir', params.treetime_folder]
        #we need command line output > os tee as I also want to see what's happening
        cmd.extend(['|', 'tee', output.treelog])
        shell (" ".join(cmd))

        #mv treetime files of interest and then remove treetime folder > don't want to consume memory
        shell(f"mv {params.regression_pdf} {output.regression_pdf}")
        #shell(f"mv {params.divergence_tree} {output.divergence_tree}")
        #shell(f"mv {params.mol_clock} {output.mol_clock}")
        shell(f"rm -r -d {params.treetime_folder}")  