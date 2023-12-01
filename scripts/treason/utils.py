#!/usr/bin/env python3

import os, sys, random, copy, dendropy
import pandas as pd
import numpy as np
from Bio import SeqIO, Seq, SeqRecord
from datetime import datetime, date
import multiprocessing as mp 
from labels import *

#months to numbers 
m2n = {"january":1, "jan":1, "february":2, "feb":2, "march":3, "mar":3, "april":4, "apr":4, "may":5, "june":6, "jun":6,"july":7, "jul":7,
       "august":8, "aug":8, "september":9, "sep":9, "october":10, "oct":10, "november":11, "nov":11, "december":12, "dec":12}

def check_sequence_length(segment, sequence, mnlp=0.95):
    """
    check if sequence is smaller than {mnlp}% of the reference sequence for the segment
    input:
        - segment: segment of interest (str)
        - sequence: DNA sequence of strain that has to be checked (str)
        - mnlp: minimum length percentage (float)
    output:
        bool: True if sequence is equal of longer than min length, else False
    """
    #for reference_1968 
    reference_lengths = {'PB2':2280,'PB1':2274,'PA':2151,'HA':1701,'NP':1497,'NA':1410,'M':890,'NS':800}
    if len(sequence) < np.round(reference_lengths[segment]*mnlp):
        return False 
    else:
        return True

def check_max_ambig(sequence, mxa=0.01):
    """
    check if sequence has ambiguous nucleotide % greater than the mxa
    input:
        - sequence: DNA sequence of strain that has to be checked (str)
        - mxa: maximum ambiguous percentage (float)
    output:
        bool: True if sequence has ambiguous content lower than mxa, else False
    """
    #also checking if there aren't any illegal characters present
    valid_nucs = ['A','T','C','G','R','Y','B','D','K','M','H','V','S','W','N']
    if not all(i in valid_nucs for i in sequence.upper()):
        return False

    N_ambig = sequence.upper().count("N")
    
    if round(N_ambig/len(sequence),2) > mxa:
        return False
    else:
        return True

def prep_tree_for_lbi(tree):
    """
    time tree annotations don't work completely fine dendropy so fixing that here
    input:
        - tree: dendropy Tree object
    return:
        - tree: same dendropy Tree object only with correct annotations
    """
    for node in tree.preorder_node_iter():
        mutations = " " #not every node was assigned a mutations
        for k,v in node.annotations.values_as_dict().items(): #mutations and date keys are expected
            if k == "date" or k.endswith("date"):
                date = float(v)
            elif k == "mutations":
                if "date" in node.annotations.values_as_dict().keys():
                    mutations = [v] if len(v) > 0 else " "
                else:
                    v = v.lstrip('"') if v.startswith('"') else v
                    mutations = [v]
                    for i in list(node.annotations.values_as_dict().keys())[list(node.annotations.values_as_dict().keys()).index(k)+1].split('"')[0].split(","):
                        mutations.append(i)
        
            #dropping the old value
            node.annotations.drop(name=k)  
        
        #change mutations from list to string for convenience
        if mutations is not None and mutations!=" ":
            mutations = ";".join(mutations)

        #assigned the new attributes 
        setattr(node, "date", date)
        node.annotations.add_bound_attribute("date")    
        setattr(node, "mutations", mutations)
        node.annotations.add_bound_attribute("mutations") 
    return tree

def calculate_LBI(tree, tau=0.3, normalize=False):
    """
    calculate the local branching index for the isolates in the given time tree
    input:
        - tree: dendropy tree object of treetime tree for which lbi's need to be calculated
        - tau: tau value for LBI calculation
    output: 
        - tree: dendropy tree object with annotated LBI in the nodes
    """
    #calculate molecular clock length for each node 
    tree.seed_node.clock_length = 0.0
    for node in tree.postorder_node_iter():
        for child in node.child_nodes():
            child.clock_length = child.annotations.get_value("date") - node.annotations.get_value("date")

    #calculate the 'up' message from node i to its parent
    #up-message = tau(1-exp(-(b/tau))) + exp(-(b/tau)) < sum runs over the childern ij of node i
    for node in tree.postorder_node_iter():
        node.down_polarizer = 0
        node.up_polarizer = 0
        for child in node.child_nodes():
            node.up_polarizer += child.up_polarizer
        bl = node.clock_length / tau
        node.up_polarizer *= np.exp(-bl)
        node.up_polarizer += tau*(1-np.exp(-bl))

    #calculate the 'down' message from the parent (internal node) to its children
    #down-message of child ij = tau(1-exp(-(b/tau))) + 
    # exp(-(b/tau))(down-message of parent i + sum of the up messages of the other child nodes of parent i )
    for node in tree.preorder_internal_node_iter():
        for child1 in node.child_nodes():
            child1.down_polarizer = node.down_polarizer
            for child2 in node.child_nodes():
                if child1 != child2:
                    child1.down_polarizer += child.up_polarizer
            
            bl = child1.clock_length /tau
            child1.down_polarizer *= np.exp(-bl)
            child1.down_polarizer += tau*(1-np.exp(-bl))

    #calculate the LBI 
    #LBI = 
    max_LBI = 0.0
    for node in tree.postorder_node_iter():
        lbi = node.down_polarizer
        for child in node.child_nodes():
            lbi += child.up_polarizer

        node.lbi = lbi
        if node.lbi > max_LBI:
            max_LBI = node.lbi

    #normalize LBI if requested
    for node in tree.preorder_node_iter():
        if normalize:
            node.lbi /= max_LBI
        setattr(node, "lbi", node.lbi)
        node.annotations.add_bound_attribute("lbi") 
    
    return tree

def get_codon(sequence, pos):
    """
    get codon, codon number and position of nucleotide in codon based on sequence and nucleotide position
    '!' assuming that sequence starts with a start codon
    input:
        - sequence: sequence of interest from which codon needs to be extracted (str)
        - pos: position of nucleotide in sequence for which codon needs to be determined (int)
    output:
        - codon: codon (str)
        - nucleotide position in the codon
    """
    for i in range(0, len(sequence),3):
        if pos in range(i,i+3):
            for k, l in enumerate(range(i,i+3)):
                if l==pos:
                    return sequence[i:i+3], k
       
def get_leaf_nodes(node):
    """
    get all leaf nodes for a given node
    input:
        - node: node for which leaf nodes need to be generated (dendropy node object)
    output:
        - leaf: list of all edge node for the input node (list)
    """
    leafs = []
    for child in node.postorder_iter():
        if child.is_leaf():
            leafs.append(child)
    return leafs

def translate_tree_mutations(sequence_dict,tree,o_nonsyn, seed=None):
    """
    translate mutations in giving input tree
    input:
        - sequence_dict: dict with {id:sequence} for all isolates present in input tree
        - tree: dendropy tree object of tree for which the mutations need to be translated
        - o_nonsyn: if only non synonymous mutations or all mutations need to be translated 
        - seed: seed for random selecting of leaf nodes
    output:
        - tree: the same dendropy tree object as the input tree only with mutations annotated
    """
    if seed is not None:
        random.seed(seed)

    for node in tree.preorder_node_iter():
        aa_mutations = []
        #if there are no mutations drop this mutations attribute
        if len(node.annotations.get_value("mutations")) == 0:
            node.annotations.drop(name="mutations")
    
        else: #mutations to translate
            if ";" in node.annotations.get_value("mutations"):
                mutations = node.annotations.get_value("mutations").split(";")
            else:
                mutations = [node.annotations.get_value("mutations")]
            
            for mut in mutations:
                nuc_pos = int(mut[1:len(mut)-1])-1 #minus 1 to account for indexing
                #codon_num is index > starts at 0
                codon_num = 0 if nuc_pos < 0 else int(nuc_pos/3) if nuc_pos %3==0 else int((nuc_pos-1)/3) if (nuc_pos-1)%3==0 else int((nuc_pos-2)/3) #if (nuc_pos-2)%3==0 
                codon_num +=1 #correction 
                if node.is_leaf():
                    #check if mutations is with coding region > not relevant else
                    if nuc_pos > len(sequence_dict[node.taxon.label.replace(" ", "_")]):
                        continue
                    #sequence contains the mutations so codon obtain from sequence will be the mutated codon
                    new_codon, pos = get_codon(sequence_dict[node.taxon.label.replace(" ", "_")], nuc_pos)
                    if "-" in new_codon: #only complete codons
                        continue
                    codon = new_codon[:pos] + mut[0].lower() + new_codon[pos +1:]
                    if "-" in codon:
                        continue

                    #check if only non synonymous mutations need to be translated
                    if o_nonsyn:
                        if Seq.Seq(codon).translate() !=  Seq.Seq(new_codon).translate():
                            aa_mut = str(Seq.Seq(codon).translate()) + str(codon_num) + str(Seq.Seq(new_codon).translate())
                            aa_mutations.append(aa_mut)
                    else:
                        aa_mut = str(Seq.Seq(codon).translate()) + str(codon_num) + str(Seq.Seq(new_codon).translate())
                        aa_mutations.append(aa_mut)
                else:
                    #in theory all child nodes of the the internal carry the mutation of the node > should there be able to select a
                    #random child to get the AA mutation
                    leafs = [leaf.taxon.label.replace(" ","_") for leaf in get_leaf_nodes(node)]
                    rleaf = random.choice(leafs)
                    
                    #check if mutations is with coding region > not relevant else
                    if nuc_pos > len(sequence_dict[rleaf]):
                        continue
                    
                    new_codon, pos = get_codon(sequence_dict[rleaf], nuc_pos)
                    if "-" in new_codon: #only complete codons
                        continue
                    codon = new_codon[:pos] + mut[0].lower() + new_codon[pos +1:]
                    if "-" in codon:
                        continue
                
                    #check if only non synonymous mutations need to be translated
                    if o_nonsyn:
                        if Seq.Seq(codon).translate() !=  Seq.Seq(new_codon).translate():
                            aa_mut = str(Seq.Seq(codon).translate()) + str(codon_num) + str(Seq.Seq(new_codon).translate())
                            aa_mutations.append(aa_mut)
                    else:
                        aa_mut = str(Seq.Seq(codon).translate()) + str(codon_num) + str(Seq.Seq(new_codon).translate())
                        aa_mutations.append(aa_mut)
                        
        #annotate the amino acid mutations at the node
        if len(aa_mutations) >0:
            aa_mutations =  ";".join(aa_mutations)
            setattr(node, "aa_mutations",aa_mutations)
            node.annotations.add_bound_attribute("aa_mutations")

    return tree