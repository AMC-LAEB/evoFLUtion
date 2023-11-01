#!usr/bin/env python3
import os, dendropy
import pandas as pd, numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from collections import OrderedDict
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime, date

def get_sequence_files(file_dir):
    """
    get protein sequence (treason output) files from input directory and return dictionary with opened sequence files 
    sequence file will be returned as list with SeqIO records per season
    input:
        - file_dir: path to treason protein output director (str)
    output:
        - sequences: dictionary with sequences per season as list with SeqIO records
    """
    sequences = {}
    for f in os.listdir(file_dir):
        if f.endswith(".fasta"):
            #assuming season in on index 1
            season = f.split(".")[0].split("_")[2]
            sequences[season] = list(SeqIO.parse(os.path.join(file_dir,f),"fasta"))
    return sequences

def get_tree_files_per_season(file_dir):
    """
    get treefiles (treason output) from input directory and return dictionary with opened tree file as dendropy object
    input:
        - file_dir: path to treason protein output director (str)
    output:
        - trees: dictionary with dendropy trees object per season
    """
    trees = {}
    for f in os.listdir(file_dir):
        #assuming season in on index 1
        if f.endswith(".nexus"):
            season = f.split(".")[0].split("_")[2]
            trees[season] = dendropy.Tree.get(path=os.path.join(file_dir,f), schema="nexus")
    return trees

def get_mature_sequences(sequences,start_codon=17):
    """
    get the mature sequence for each sequence in the input sequences 
    assuming that seqeunces are complete AA CDS 
    input:
        -sequences: list SeqIO sequences records 
        -start_codon: start_codon of the mature sequences
    output:
        -mature_sequences: list SeqIO sequences records with mature proteins only
    """
    mature_sequences = sequences.copy()
    for rec in mature_sequences:
        rec.seq = Seq("".join(rec.seq[start_codon-1:])) #-1 to account for indexing
    return mature_sequences

def get_time_frames(seasons, time_frame_size=5, step=1):
    """
    get time_frames of size n with step m given list of input seasons
    input:
        - seasons: list of seasons to generate time frames for
        - time_frame_size: size of time frame (default 5)
        - step: step of season to move for next time frame (default 1)
    output:
        - timeframe: list  of time frames (time frame is a list of seasons)
    """
    seasons = sorted(seasons) #just in case

    if len(seasons) <= time_frame_size:
        return [seasons]
    
    timeframes = []
    for i in range(0,len(seasons)-1,step):
        tf = [seasons[i]]
        if i+time_frame_size < len(seasons): #if full time can be made
            for j in range(1,time_frame_size):
                tf.append(seasons[i+j])
            timeframes.append(tf)
    return timeframes

def determine_season(dates, start_month="05"):
    """
    determine the season of a sample given the input date
    ASSUMING SEASON STARTS IN MAY
    incomplete date > only year will be label to the start year of the season (2015 > season: 1516)
    input:
        - dates: list of date string > yyyy-mm-d
    output:

        - seasons: pd dataframe columns with seasons > same length as input column
    """

    if type(dates) == np.ndarray:
        dates = [d for dl in dates.tolist() for d in dl]
        array = True

        
    seasons = []
    for date in dates:

        year, month = date.split("-")[0:2] if "-" in date else [date, None]
        s1 =  "".join(list(year[2:]))
        if month is None:
            s2 = "0"+str(int(s1)+1) if len(str(int(s1)+1))==1 else str(int(s1)+1)
        else:
            if int(month) < int(start_month):
                s2 = s1
                s1 = "0"+str(int(s1)-1) if len(str(int(s1)-1))==1 else str(int(s1)-1)
                if s1 == "-1":
                    s1 = str(100+int(s1))
            else:
                s2 = "0"+str(int(s1)+1) if len(str(int(s1)+1))==1 else str(int(s1)+1)
        seasons.append(s1+s2)
    
    seasons = pd.DataFrame.from_dict({"seasons":seasons})
    return seasons

def seqs_per_season_from_metadata(metadata_file, seasons, seq_type=""):
    """
    count the number of sequence per season from GISAID metadata file
    input:
        - metadata_file: file path to metadata file of interest
        - seasons: list of season of interest
        - seq_type: type of sequences from column naming, e.g. 'clinical' (default: '')
    output:
        - seq_count: pandas dataframe with the number of sequences per season
    """
    #load metadata
    metadata = pd.read_csv(metadata_file)
    if "index" in metadata.columns:
        metadata = metadata.drop(["index"], axis=1)

    #infer season
    metadata[["season"]] = determine_season(metadata[["Collection_Date"]].values)

    seq_count = pd.DataFrame.from_dict({season:len(metadata[metadata["season"]==season]) for season in seasons}, orient="index")
    seq_count = seq_count.reset_index(level=0).rename(columns={"index":"season", 0:f"number of {seq_type} sequences".replace("  ", " ")})

    return seq_count

def annotate_aa_position(sequence, numbering_dict=None):
    """
    annotate each amino acid in the sequence with its position and return list of mutation/positions
    input:
        - seqeunce: amino acid sequence to be annotated
        - numbering_dict: dictionary with sequence index: mature numbering res number
    output:
        - pos_aa: list with position and aa combos 
    """
    pos_aa = []
    for i, aa in enumerate(str(sequence)):
        if numbering_dict is None:
            pos_aa.append(f"{str(i+1)}{aa}")
        else:
            pos_aa.append(f"{numbering_dict[i]}{aa}")
    return pos_aa

def get_variants(sequences):
    """
    get all the variants for each position given the input sequences > two or more variants per positions 
    input:
        - sequences: dict with list of sequences per season {season:[sequence list]}
    output:
        - variants: list of all variants per position, sorted by ascending position
    """
    #get observed amino acids 
    observed_aa = set()
    for season, seqs in sequences.items():
        for seq in seqs:
            variants = annotate_aa_position(seq.seq)
            observed_aa.update(variants)

    #sort the variants per position
    #ignore the variants that contain the amino acid "X" as this an ambiguous amino acid
    sort_variants = {}
    for var in observed_aa:
        if "X" not in var:
            pos = int("".join(var[:-1]))
            try:
                sort_variants[pos].append(var)
            except:
                sort_variants[pos] = [var]
    sort_variants = dict(OrderedDict(sorted(sort_variants.items())))

    variants = []
    for varlist in sort_variants.values():
        #if len(varlist) > 1:
        variants.extend(varlist)
    return variants

def determine_var_frequency(var, sequences,numbering_dict=None):
    """
    count the occurence of the input mutation in a certain sequence set
    input:
        - mutation: mutations of which the occurence needs to be counted (str)
        - sequences: sequence set for which mutation occurence needs to be counted (list with SeqIO sequence records)
        - numbering_dict: dictionary with sequence index: mature numbering res number
    output:
        - count: occurence of the mutations in the given sequence set
    """
    total = len(sequences)
    count = 0 
    if numbering_dict is None:
        pos = int(var[:-1])  #get the position of the mutation; 
    else:
        pos = list(numbering_dict.keys())[list(numbering_dict.values()).index(int(var[:-1]))]


    for record in sequences:
        if len(record.seq) >= pos-1:
            aa = str(record.seq)[pos-1] #account for indexing
            if aa == var[-1]:
                count +=1 
            elif aa == "X":
                total -= 1
        else:
            total -= 1
    
    frequency = (count/total) 
    #if original is observed at least once return frequency else return 0
    return frequency#return 0.0

def determine_frequency_pattern(freqs, minf=0.05, maxf=0.95):
    """
    determine frequency pattern for every mutation
    input:
        -freqs: dictionary containing frequency per mutation per season
    output:
        -freqs: same dictionary with two added 'columns': position, and status (frequency pattern)
    """
    #minor variant is frequencies remain lower than 5%
    #variant is dominant if frequency in last season > 95%
    #variant appears if frequency >5% in later season 
    #variant disappears if frequency in later seasons drops to < 5%
    #variant reappears if frequency is initially > 5% drops < 5% and than > 5%
    ###variant declines if frequency drops by at least 20% between start and end of period of interest
    ###variant increases if frequency increases by at least 20%  
    ####competition: variant is fixed at some point and drops in frequency resulting in a frequency below 95% 
    #if variant appears and disappears within timeframe, variant demises
    for mut, sf in freqs.items():
        slist = list(sf.keys())
        f = [float(fi) for fi in list (sf.values())]
        dominant,fixed, minor, reappears, disappears = False, False, False, False, False
        appears, decreases, increases, competition = False, False, False, False
        demises= False
        
        #if start is below minf 
        if f[0] < minf: #appearing, demising, (and later appearing again) or minor variant
            if f[-1] > maxf: #crazy scenario but it could be the case
                fixed = True
            else:
                for freq in f[1:]:
                    if freq > minf:
                        if not appears:
                            appears = True
                        elif appears and demises:
                            reappears = True
                        if freq > maxf:
                            competition = True
                    else:
                        if appears:
                            demises = True
                        
                if not appears and not demises:
                    minor = True
        elif f[0] >= maxf: #fixed and remains fixed or competition
            if f[-1] >= maxf:
                dominant = True
            else:
                competition = True
        else:
            if f[-1] >= maxf:
                fixed = True
            elif f[-1] > minf:
                for freq in f[1:-1]:
                    if freq < minf:
                        disappears = True
                    if freq >= maxf:
                        competition = True
                if disappears:
                    reappears = True
            else:
                disappears = True
            if f[0] < f[-1] and abs(f[0]-f[-1]) > 0.2:
                increases = True
            elif f[0] > f[-1] and abs(f[0]-f[-1]) > 0.2:
                decreases = True


        #set position 
        freqs[mut]["position"] = int(mut[:-1])
        #set "status" 
        if dominant:
            freqs[mut]["status"] = "dominant"
        elif fixed:
            freqs[mut]["status"] = "fixed"
        elif competition and f[-1] < minf:
            freqs[mut]["status"] = "disappearing"
        elif competition and f[-1] < minf:
            freqs[mut]["status"] = "demising"
        elif competition:
            freqs[mut]["status"] = "competing"
        elif minor:
            freqs[mut]["status"] = "minor"
        elif reappears and f[-1] > minf:
            freqs[mut]["status"] = "reappearing"
        elif reappears: #last frequency below cutoff so demising
            freqs[mut]["status"] = "demising"
        elif demises:
            freqs[mut]["status"] = "demising"
        elif disappears:
            freqs[mut]["status"] = "disappearing"
        elif appears:
            freqs[mut]["status"] = "appearing"
        elif increases:
            freqs[mut]["status"] = "increasing"
        elif decreases:
            freqs[mut]["status"] = "decreasing"
        else:
            freqs[mut]["status"] = "normal"             
    return freqs

def determine_frequency_pattern_general(var_freqs, fixed_freq=0.8, transient_freq=0.4, minor_freq=0.05):
    """
    determine frequency pattern of mutations but a more general labeling
    frequency reaches above the fixed_freq cutoff --> fixed
    frequency remaining between fixed_freq and the transient_freq cutoff --> competing
    frequency reaches below transient_freq cutoff --> transient
    frequency below minor cutoff over all time frames --> minor 
    input:
        - var_freqs: dict in dict that looks like {variant:{season:frequency}}
    """
    freq_pattern = {}
    for var, sf in var_freqs.items():
        seasons = list(sf.keys())
        freqs = [float(fi) for fi in list (sf.values())]
        
        #start by checking if the one of the frequencies reaches above the minor freq cutoff
        minor = True if not any(f > minor_freq for f in freqs) else False
        
        #labeling is kept simple as in basically based on the last observed frequency 
        fixed = True if freqs[-1] >= fixed_freq else False
        transient = True if freqs[-1] < transient_freq else False
        competing = True if freqs[-1] < fixed_freq and freqs[-1] >= transient_freq else False

        if minor:
            freq_pattern[var] = "minor"
        elif transient:
            freq_pattern[var] = "transient"
        elif competing:
            freq_pattern[var] = "competing"
        elif fixed:
            freq_pattern[var] = "fixed"
    
    return freq_pattern

def get_sequence_ids_per_variant(variant,sequences):
    """
    get the sequence IDs of the sequences that carry the AAP (certain amino acid at certain position)
    input:
        - aap: amino acid and its position (example: 311H)
        - seqeunces: list of SeqIO records
    output:
        - ids: list of the sequence IDs of the sequence that carry the AAP
    """
    pos = int("".join(variant[:-1]))
    aa= variant[-1]

    ids = []
    for rec in sequences:
        if rec.seq[pos-1] == aa:#-1 to account for indexing
            ids.append(rec.id)
    return ids

def get_background(sequence, positions):
    """
    get the mutation background of input sequences for the positions of interest
    input:
        - sequence: amino acid sequences (str/Seq)
        - positions: list of positions of interest (list)
    output:
        - background: list of AAP specific for the sequence
    """
    background = []
    for pos in positions:
        aa = sequence[pos-1] #-1 to correct for indexing
        background.append(str(pos)+aa)
    return background

def sort_mutations_by_position(mutations):
    """
    given an input list of mutations sort the mutation by position
    input:
        - mutations: list of mutations to be sorted
    output:
        - mutations: same list of mutation as input only sorted on position
    """
    pos_mut = {int("".join(mut[1:-1])):mut for mut in mutations}
    pos_mut_sorted = dict(OrderedDict(sorted(pos_mut.items())))
    return [mut for mut in pos_mut_sorted.values()]

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

def get_lbi_from_sequence_node(tree, id_list):
    """
    get the LBI from node in the input tree if there is a node that only carries leafs with ids in the id_list
    input:
        - tree: dendropy tree object of tree that needs to be checked 
        - id_list: list of ids of interest
    output:
        - lbi: lbi of the node if node with only the ids as leafs is found else None
    """
    if len(id_list) == 1:
        for node in tree.leaf_node_iter():
            if node.taxon.label.replace(" ", "_") == id_list[0]:
                return round(float(node.annotations.get_value("lbi")),3)
        return None #id not present in tree > no lbi
    else:
        for node in tree.postorder_internal_node_iter():
            leafs = get_leaf_nodes(node)
            labels = [leaf.taxon.label.replace(" ","_") for leaf in leafs]
            if sorted(id_list) == sorted(labels):
                return round(float(node.annotations.get_value("lbi")),3)

        return None #no node that matches

def find_cluster_nodes(tree,seqs,cluster_definitions):
    """
    find the internal nodes that make up clusters as defines in the cluster definitions in the tree
    that are ancestral to the sequences defined in seqs
    for all sequences for which no cluster node can be found, the leaf node itself will be returned
    input:
        - tree: dendropy tree object that needs to be searched
        - seqs: list of taxon/sequence IDs of interest
        - cluster_definitions: dictionary with max cluster sizes and there min percentage sequences
                               in seqs that are leafs nodes in this cluster
    output:
        - cluster_nodes: dict with cluster nodes and the sequences IDs of seqs that those nodes are 
                         parent of 
        - singleton: list of nodes to considered a singleton
    """
    #get the max cluster size based on the number of sequences
    max_clus_size = [cs for cs in sorted(list(cluster_definitions.keys())) if cs>len(seqs)][0]
    min_clus_size = min(cluster_definitions.keys())
    clus_per = cluster_definitions[max_clus_size] #get percentage for cluster

    singletons = []
    cluster_nodes = {}

    #check if min cluster size can be reached with the number sequences
    if len(seqs) < min_clus_size:
        for leaf_node in tree.leaf_node_iter():
            if leaf_node.taxon.label.replace(" ", "_") in seqs:
                singletons.append(leaf_node)
        return cluster_nodes, singletons
    
    seen_seqs = set()#keep track of the sequences that were assigned a cluster

    #get all nodes that meet the cluster requirement
    ancestral_nodes = {}
    for node in tree.preorder_internal_node_iter():
        all_leafs = [l.taxon.label.replace(" ","_") for l in get_leaf_nodes(node)]
        if len(all_leafs) <= max_clus_size:
            mut_leafs = [l for l in all_leafs if l in seqs]
            if len(mut_leafs)/len(all_leafs) >= clus_per:
                ancestral_nodes[node] = mut_leafs

    
    #get the biggest nodes > clusters carrying the most leafs
    seen_index = []
    for n in sorted(ancestral_nodes, key=lambda n: len(ancestral_nodes[n]), reverse=True):
        mut_leafs = ancestral_nodes[n]
        if len(set(mut_leafs).intersection(seen_seqs)) == 0:
            cluster_nodes[n] = mut_leafs
            seen_seqs.update(mut_leafs)

    #find the sequences that aren't part of a cluster
    cluster_seqs  = [v for l in ancestral_nodes.values() for v in l]
    single_leafs = [seq for seq in set(seqs).difference(seen_seqs) if seq not in cluster_seqs]
    #add the single nodes to the cluster nodes (we are using this LBI) (name of the dict is irrelevant)
    for node in tree.leaf_node_iter():
        if node.taxon.label.replace(" ","_") in single_leafs:
            #cluster_nodes[node] = [node.taxon.label.replace(" ","_")]
            singletons.append(node)
            seen_seqs.add(node.taxon.label.replace(" ","_"))

    return cluster_nodes, singletons

def calculate_weighted_lbi(cluster_nodes, singletons):
    """
    calculate the weighted average LBI based on the cluster nodes and the number of leaves for each cluster node
    input:
        - cluster_nodes: dict with cluster node and a list of leaves
        - singletons: list of singleton nodes
    output:
        - weighted_LBI: the weighted average LBI
    """
    lbi_weights = {}
    #add lbi for cluster nodes to lbi weights
    for node, children in cluster_nodes.items():
        lbi = round(float(node.annotations.get_value("lbi")),3)
        try:
            lbi_weights[lbi] += len(children)
        except:
            lbi_weights[lbi] = len(children)

    #add lbi for singletons to lbi weights
    for node in singletons:
        lbi = round(float(node.annotations.get_value("lbi")),3)
        try:
            lbi_weights[lbi] += 1
        except:
            lbi_weights[lbi] =  1
        
    return round(sum([k*v for k,v in lbi_weights.items()])/sum(list(lbi_weights.values())),3)

def get_branch_length(parent, child, tree):
    """
    get the branch length between the parent node and the child node > child must be a direct descendant
    input:
        - parent: ancestral node of branch 
        - child: 'leaf' node of branch 
        - tree: dendropy tree object
    output:
        - branch length between the parent and child node
        """
    all_edges =tree.edges()
    for edge in all_edges:
        if edge._get_tail_node() == parent and edge._get_head_node() == child:
            return edge.length

def get_total_branch_length(mrca, leaf, tree):
    """
    get the total branch length from MRCA to child nodes 
    input:
        - MRCA: ancestral node from which the total branch length needs to be determined 
        - leaf: tip/internal node to which the total branch length needs to be determined
    output:
        - total_branch_length: the total branch length between MRCA and leaf
    """
    total_branch_length = 0.00
    parent = leaf.parent_node
    while leaf !=mrca: #need to get the total length from MRCA to tip
        total_branch_length += get_branch_length(parent, leaf, tree)
        leaf = parent
        parent= leaf.parent_node
    return total_branch_length

def get_nodes_from_labels(isolates, tree):
    """
    get the node that belongs to isolate in the input isolates from the tree
    input:
        - isolates: list of isolates for which the tree nodes need to be found
        - tree: dendropy object of the phylogenetic tree in which the isolates are present
    output:
        - nodes: list of the nodes from the isolates
    """
    nodes = [""]*len(isolates)
    for node in tree.preorder_node_iter():
        if node.taxon is not None:
            if node.taxon.label.replace(" ", "_") in isolates:
                nodes[isolates.index(node.taxon.label.replace(" ", "_"))] = node
    return nodes

def get_mrca(nodes, tree):
    """
    get the Most Recent Common Ancestor of the two input nodes given the tree
    input:
        - nodes: list of dendropy node objects for which 
        - tree: dendropy tree object in which the nodes are present
    output:
        - mrca: mrca dendropy node object of mrca of node1 and node2.
    """
    for internal_node in tree.postorder_node_iter():
        children = get_leaf_nodes(internal_node)
        if all([node in children for node in nodes]):
            return internal_node

def get_depths(nodes, tree):
    """
    get the depths from the nodes in the input list to their mrca
    input:
        - nodes: list of dendropy node objects
        - tree: dendropy tree object of tree in which the nodes are located
    output:
        - depths: list of depths for each node in nodes to the mrca of all nodes 
    """
    # get mrca 
    mrca = get_mrca(nodes, tree)
    
    depths = []
    for node in nodes:
        depths.append(get_total_branch_length(mrca,node,tree))
    return depths

def get_node_depth(node):
    """
    determine depth of node in tree > number of splits
    input:
        - node: node for which the depth need to be determine
    output:
        - depths: max number of splits from leaf node to input node
    """
    if node.is_leaf():
        return 0
    leaf_nodes = get_leaf_nodes(node)

    max_depth = 0
    for leaf in leaf_nodes:
        depth = 1
        parent = leaf.parent_node
        while parent != node:
            depth +=1
            parent = parent.parent_node
        if depth > max_depth:
            max_depth = depth
    return max_depth

def count_branches(parent, child):
    """
    count the number of branch between parent node and child node
    input:
        - parent: parent node 
        - child: child node
    output:
        - n_branches: number of branches
    """ 
    n_branches = 1
    p = child.parent_node
    while p != parent:
        n_branches +=1 
        p = p.parent_node
    return n_branches
 
def calculate_wiener_index(distance_matrix):
    """
    Calculate the wiener index from the input matrix
    Wiener index is the total sum of paths between all isolates in the input matrix
    (note:)
    input:
        - distance_matrix: sub matrix of distance matrix containing only the isolates for which
            the wiener index needs to be determined
    output:
        - wiener: calculated wiener index (int)
    """
    #check matrix is numpy array or matrix and not pd dataframe
    if type(distance_matrix) == pd.DataFrame:
        distance_matrix = distance_matrix.to_numpy()
    
    #every path is counted twice > A -- B and B -- A, so divide by two
    wiener = distance_matrix.sum()
    if wiener > 0:
        wiener = wiener/2

    return wiener

def get_pairwise_distance_matrix(isolates, tree):
    """
    get the pair wise distance between all isolates in the input list from the phylogenetic tree 
    and return the pair wise distance matrix
    input:
        - isolates: list of isolates for which the pairwise distance matrix needs to be constructed
        - tree: dendropy object of the input phylogenetic tree
    """
    #create empty matrix
    pwd_matrix= np.zeros((len(isolates), len(isolates)))
    targets = get_nodes_from_labels(isolates,tree) #get nodes instead of the labels
    
    for i in range(0,len(targets)):
        t1 = targets[i] #target1 
        for j in range(i, len(targets)):
            t2 = targets[j] #targets

            #check if it's not the same node > distance will be 0 in that case
            if i==j:
                pwd_matrix[i,j] = 0
            else:
                #get mrca between the two nodes
                mrca = get_mrca(t1, t2, tree)

                #get distance from each node to the mrca 
                mrca_t1 = get_total_branch_length(mrca, t1, tree)
                mrca_t2 = get_total_branch_length(mrca, t2, tree)
                #get distance from t1 to t2
                t1_t2 = mrca_t1 + mrca_t2

                #assign distance in matrix
                pwd_matrix[i,j] = t1_t2
                pwd_matrix[j,i] = t1_t2
    return pwd_matrix

def determine_n_cherries(nodes): 
    """
    determine the number of cherry nodes 
    input:
        - nodes: list of dendropy node objects for which the number cherries need to be checked
    output:
        - n_cherries: number of cherry nodes found
    """
    n_cherries = 0
    if len(nodes) >= 2:
        for i in range(0, len(nodes)):
            for j in range(i+1, len(nodes)):
                if nodes[i].parent_node == nodes[j].parent_node:
                    #check if parent node only consists of two leaves > cherry shape
                    if len(nodes[i].parent_node.child_nodes()) <=2:
                        n_cherries +=1   
    return n_cherries

def determine_n_pitchforks(nodes, pitchfork_size=3):
    """
    determine the number of pitchforks found the phylogenetic tree
    input:
        - nodes: list of dendropy node objects for which the number cherries need to be checked
        - pitchfork_size: maximum number of nodes to label as pitchfork (default=3 > min cluster size)
    output:
        - n_pitchforks: number of cherry nodes found
    """
    n_pitchforks = 0
    seen_nodes = [] #node can only be counted once to pitchfork
    for i in range(0, len(nodes)):
        if nodes[i] not in seen_nodes:
            
            pitchfork_nodes = [nodes[i]]
            parent = nodes[i].parent_node #get parent node 

            while (len(get_leaf_nodes(parent))) <= pitchfork_size:
                pitchfork_nodes.extend(list(set(get_leaf_nodes(parent)).difference(pitchfork_nodes)))
                parent = parent.parent_node
  
            if len(set(nodes).difference(pitchfork_nodes)) == 0:
                seen_nodes.extend(pitchfork_nodes)
                n_pitchforks += 1

    return  n_pitchforks

def add_vline(ax, xpos, ypos, lent=1):
    line = plt.Line2D([ypos,ypos], [xpos+0.02, (xpos+lent)-0.02], color='#808080', transform=ax.transAxes)
    line.set_clip_on(False)
    ax.add_line(line)

def add_hline(ax,xpos, ypos,lent=1):
    line = plt.Line2D([xpos+0.01, (xpos+lent)-0.01],[ypos,ypos], color='#808080', transform=ax.transAxes)
    line.set_clip_on(False)
    ax.add_line(line)

def edit_timelabels_for_plot(dataframe):
    if "time point" in dataframe.columns:
        for i, row in dataframe.iterrows():
            if row["time point"].startswith("t"):
                tp = row['time point'].replace("t","$t_{") + "}$"
                if "f" in tp:
                    tp = tp.replace("f","f{").replace("}$", "}}$")
                    if "01" in tp:
                        tp=tp.replace("01", "0-1")
                dataframe.iloc[i, dataframe.columns.get_loc("time point")] = tp
    else:
        print ("no 'time point' column found in dataframe > dataframe stays unchanged")
    return dataframe

def generate_p_value_annotation(dataframe):
    """change pvalue to string for input dataframe and return a copy"""
    df_annot = dataframe.copy()
    if "p-value" in df_annot.columns:
        for i, row in df_annot.iterrows():
            npv = "≤0.05" if row["p-value"] <= 0.05 else str(round(row["p-value"] ,2))
            df_annot.iloc[i, df_annot.columns.get_loc("p-value")] = npv
    else:
        print ("no 'p-value' column found dataframe remains unchanged")
    return df_annot

def change_width(ax, new_value):
    for patch in ax.patches:
        current_width = patch.get_width()
        diff = current_width - new_value
        patch.set_width(new_value)
        patch.set_x(patch.get_x() + diff * .5)

def date_to_float_date(d):
    """convert date string (yyyy-mm-dd) to float date"""
    d = "-".join(d.split("-")[::-1])

    year = int(d.split("-")[-1])
    yd = datetime.strptime(d,'%d-%m-%Y').timetuple().tm_yday
    
    if ((year % 400 == 0) and (year % 100 == 0)) or ((year % 4 ==0) and (year % 100 != 0)): #leap year
        dfloat = yd/366
    else:
        dfloat = yd/365
    return year + dfloat