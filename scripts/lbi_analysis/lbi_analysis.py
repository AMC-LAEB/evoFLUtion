#!usr/bin/env python3
import os, sys, itertools,argparse, random
import pandas as pd, numpy as np, statistics as st
import scipy.stats as stats
from scipy.optimize import curve_fit
from collections import OrderedDict
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
import matplotlib.ticker as mtick
import matplotlib.colorbar as mcp
from utils import *
import multiprocessing as mp 
from tqdm import tqdm

sns.set_theme(style="white", font="Arial")
sns.set_context({'font.size': 10.0, 'axes.labelsize': 'medium', 'axes.titlesize': 'medium',
                 'xtick.labelsize': 'small', 'ytick.labelsize': 'small', 'legend.fontsize': 'small',
                 'legend.title_fontsize': None, 'axes.linewidth': 0.8, 'grid.linewidth': 0.8,
                 'lines.linewidth': 1.5, 'lines.markersize': 6.0, 'patch.linewidth': 1.0,
                 'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.minor.width': 0.3,
                 'ytick.minor.width': 0.3, 'xtick.major.size': 3.5, 'ytick.major.size': 3.5,
                 'xtick.minor.size': 2.0, 'ytick.minor.size': 2.0,})

#color palettes
base_color = "#B6B8B9"
cpal = [base_color,sns.color_palette("PiYG")[0],sns.color_palette("PiYG")[4]]
cpal2 = ["black"]#,"#8b8c89"]
group_order = ["A", "B", "C"]

#generate matplotlib color palette for lbi_heatmap 
colors = []
for i, c in enumerate(sns.color_palette("PiYG", n_colors=34)[4:9] + list(sns.color_palette("PiYG", n_colors=34)[24:29])):
    colors.extend([c,c])
pvalue_colmap = LinearSegmentedColormap.from_list(name="lbi heatmap",
        colors=list(zip([0.0,0.01,0.01,0.02,0.02,0.03,0.03,0.04,0.04,0.05,0.05,0.20,0.20,0.40,0.40,0.6,0.6,0.8,0.8,1.0],
                    colors)))

#solid vars
cwd = os.getcwd()
filedir = os.path.abspath(os.path.dirname(__file__))

all_epitope_sites = {"H3N2":{"HA":{"A":[122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146], 
                                   "B":[155,156,157,158,159,160,186,187,188,189,190,191,192,193,194,195,196,197,198], 
                                   "C":[44,45,46,47,48,49,50,51,52,53,54,274,275,276,277,278,279,280,281], 
                                   "D":[201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220], 
                                   "E":[78,79,80,81,82,83,84,85,86,87,89,90,91,92,93,94,261,262,263,264,265,266]},
                             "NA":{"":[118,119,151,152,165,178,179,198,222,224,227,274,276,277,292,294,371,406,425]}},
                     "H1N1pdm":{"HA":{"Sa":[124,125,153,154,155,156,157,159,160,161,162,163,164],
                                      "Sb":[184,185,186,187,188,189,190,191,192,193,194,195], 
                                      "Ca":[137,138,139,140,141,142,166,167,168,169,170,203,204,205,221,222,235,236,237],
                                      "Cb":[70,71,72,73,74,75]},
                                "NA":{"":[118,119,151,152,165,178,179,198,222,224,227,274,276,277,292,294,371,406,425]}},
                    }

all_egg_adaptive_mutations = {"H3N2":{"HA":{"major": ["H156Q", "H156R", "H183L", "G186V", "L194P", "S219Y","S219F", "N246H", "N246S"],
                                            "minor":["T160K", "H183L", "A196T", "H183L", "A196T"]}},
                              "H1N1pdm":{"HA":{"major":[], "minor":[]}}}

#starting positions of MATURE in the ORF > for HA skipping the signal protein
starting_positions = {"H3N2":{"PB2":1, "PB1":1, "PA":1, "HA":17, "NA":1},
                      "H1N1pdm":{"PB2":1, "PB1":1, "PA":1, "HA":18, "NA":1}}


def ArgumentParser():
    """Argument parser"""

    parser = argparse.ArgumentParser(prog = "LBI_analysis.py", #changing this later
        formatter_class = argparse.RawTextHelpFormatter,
        description = "Analyses and comparison of mutation specific LBI distribution from treason (protein) output" ) 
    
    parser.add_argument('-d','--data-folder', required=True, action="store", type=str, help="data folder where treason protein outputs are stored")
    parser.add_argument('-o','--output', required=False, action="store", type=str, help="output directory to store all generated output file (default: ./)")
    parser.add_argument('-s','--segment', required=False, action="store", type=str, default="HA", choices=['PB2','PB1','PA','HA','NA'], help="segment abbreviation that will be analyzed (default: HA)")
    parser.add_argument('-st','--subtype', required=False, action="store", type=str, default="H3N2", choices=["H3N2", "H1N1pdm"], help="subtype to be analyzed (default: A/H3N2)")
    
    parser.add_argument('-mc','--minor-cutoff', required=False, action="store", type=float, default=0.05, help="minimum frequency cut-off for macroscopic detection of mutations (default:0.05)")
    parser.add_argument('-fc','--fixed-cutoff', required=False, action="store", type=float, default=0.95, help="frequency cut-off for fixation of mutations")

    parser.add_argument('-p', '--plot', required=False, action="store_true", help="make plots")

    parser.add_argument('-e', '--epitope', required=False, action="store_true", help="Analyse epitope sites (HA and NA only)")

    parser.add_argument('-r','--redo', required=False, action="store_true", help="Redo all steps except the initial merging step if data already exists")
    parser.add_argument('-t','--threads',required=False,action='store',type=int,help="Number of threads (default 1)")
    #parser.add_argument('-v', '--verbose', required=False,action="store_false",help="print a bunch of stuff")

    #check if arguments were entered 
    if len(sys.argv[2:])<1: 
        parser.print_help()
        sys.exit(-1)
    else:
        return parser.parse_args()
    
def get_frequency(fdir, analysis, minor_cutoff, start_pos):
    seasons, mature_seqs, trees = get_sequences_trees_seasons(fdir, start_pos)
    mut_freqs = get_mutation_frequency_trajectories(mature_seqs, minor_cutoff, analysis)
    return mut_freqs
    
#analysis function for pooling
def run_lbi_analysis(analysis, fdir, mutation_freqs, start_pos, minor_cutoff=0.05, fixed_cutoff=0.95, epitope_positions=None, timepoints=["t_0", "t_1"]):
    seasons, mature_seqs, trees = get_sequences_trees_seasons(fdir, start_pos)
    mutation_classification = classify_trajectories(mutation_freqs, seasons, minor_cutoff, fixed_cutoff)
    mutation_season_seqs = get_mutation_seqs_season(mutation_freqs, mature_seqs)
    mutation_nodes = get_mutation_nodes(mutation_classification, mutation_season_seqs, trees, timepoints)

    if epitope_positions is not None:
        mutation_info = get_additional_info(mutation_freqs, epitope_positions)
        mutation_lbis = collect_mutation_lbis(mutation_nodes, mutation_classification ,mutation_info)
    else:
        mutation_lbis = collect_mutation_lbis(mutation_nodes, mutation_classification)

    #add columns with analysis number
    mutation_lbis["analysis"] = analysis
    return mutation_lbis
  
def statistical_lbi_analysis(a, lbi_frame, epitope=False, timepoints=["t_0", "t_onset"]):
    """
    statistical mann-whitney U analysis of the different mutation groups at t_0 and t_onset
    """
    if epitope:
        mutation_types = ["epitope", "non-epitope", "all"]
    else:
        mutation_types = ["all"]
    stat_results = []
    #for a in lbi_frame["analysis"].unique():
    for t in timepoints:
        tp = "$" + t.split("_")[0] +"_{"+ t.split("_")[-1]+ "}$"
        df = lbi_frame.query(f"season=={t}")
        for mt in mutation_types:
            subdf = df.copy() if mt=="all" else df[df["epitope site"]==True] if mt=="epitope" else df[df["epitope site"]==False]
            
            if t =="t_0": #a vs b& c and c vs a&b only at t0
                try:
                    uvalue, pvalue = stats.mannwhitneyu(subdf[subdf["group"]=="A"]["LBI"], subdf[subdf["group"]!="A"]["LBI"])
                    stat_results.append([a, mt,tp, "A vs. \nB & C", pvalue])
                except:
                    stat_results.append([a, mt,tp, "A vs. \nB & C", np.NaN])
                try:
                    uvalue, pvalue = stats.mannwhitneyu(subdf[subdf["group"]=="C"]["LBI"], subdf[subdf["group"]!="C"]["LBI"])
                    stat_results.append([a, mt,tp, "C vs. \nA & B", pvalue])
                except:
                    stat_results.append([a, mt,tp, "C vs. \nA & B", np.NaN])
            
            try:
                uvalue, pvalue = stats.mannwhitneyu(subdf[subdf["group"]=="A"]["LBI"], subdf[subdf["group"]=="B"]["LBI"])
                stat_results.append([a, mt,tp, "A vs. B", pvalue])
            except:
                stat_results.append([a, mt,tp, "A vs. B", np.NaN])
                    
    stat_results = pd.DataFrame.from_records(stat_results, columns=["analysis", "mutation type", "time point", "groups compared", "p-value"])
    return stat_results

def frequency_trajectory_plot(mutation_freqs, mutation_classification, seasons, minor_cutoff=0.05):

    #plot mutation trajectory > some mutations have two trajectories and we'd want to color those appropriately 
    mfs = mutation_classification.copy()
    mfs[seasons] = np.NaN
    for i, row in mfs.iterrows():
        if seasons.index(row["t_0"]) != 0:
            subseasons = seasons[seasons.index(row["t_0"])-1:seasons.index(row["t_d"])+1]
        else:
            subseasons = seasons[seasons.index(row["t_0"]):seasons.index(row["t_d"])+1]
        if seasons.index(row["t_d"]) != len(seasons)-1:
            subseasons.append(seasons[seasons.index(row["t_d"])+1])
        for s in subseasons:
            mfs.loc[i,s] = mutation_freqs[mutation_freqs["mutation"]==row["mutation"]][s].values[0]

    #initialize figure
            
    fig, ax = plt.subplots(1,1, figsize=(8,3))
    sns.lineplot(data=pd.melt(mfs, id_vars=["analysis","mutation","group", "t_0", "t_1", "t_d"], var_name="season", value_name="frequency"), ax=ax, x="season",
                y="frequency", units="mutation", estimator=None, hue="group", hue_order=["A", "B", "C"], palette=cpal, lw=1)
    #cosmetics
    ax.spines['top'].set_visible(False)
    ax.spines["right"].set_visible(False) 
    ax.axhline(y=minor_cutoff,c="black",linestyle=":", linewidth=1)
    ax.text(x="1920", y=minor_cutoff, s="macroscopic\ndetection", fontsize=8, va="center")
    #ax.axhline(y=fixed_cutoff,c="black",linestyle=":", linewidth=1)

    #redo tick labels
    new_labels = [f"20{''.join(season[:2])}-\n20{''.join(season[2:])}" for season in seasons]
    ax.set_xticklabels(new_labels, size=8)
    ax.legend (loc="center left", bbox_to_anchor=(1.0,0.5), title="Group")


    return fig

def calculate_support_values(stat_results):
    #calculate support for each time point, group, mutation type, subgroup combo
    support_values = []
    for tp in stat_results["time point"].unique():
        df = stat_results[stat_results["time point"]==tp]
        for g in df["groups compared"].unique():
            subdf = df[df["groups compared"]==g]
            for mt in subdf["mutation type"].unique():
                subsubdf = subdf[subdf["mutation type"]==mt]
                #for sg in subsubdf["subgroup"].unique():
                #    subsubsubdf = subsubdf[subsubdf["subgroup"]==sg]
                support = sum(p < 0.05 for p in subsubdf["p-value"].values)
                support_values.append([tp, g, mt, support])
    support_values = pd.DataFrame.from_records(support_values, columns=["time point", "groups compared", "mutation type", "support"])  
    return support_values

def epitope_LBI_heatmap_plot(stat_plot, support_values):

    def label_columns(ax, df):
        xpos = 1.005
        scale = 1./df.columns.size
        for level in reversed(range(df.columns.nlevels)):
            pos = df.columns.size
            labels = [(k, sum(1 for i in g)) for k,g in itertools.groupby(df.columns.get_level_values(level))]

            for label, rpos in reversed(labels):
                pos -= rpos
                lypos = (pos + rpos/2)*scale 
                
                if level == 0: 
                    if label.startswith("t"):
                        label = label.replace("t","$t_{") + "}$"
                        
                    ax.text(lypos,xpos, label, ha='center', transform=ax.transAxes, size=10)

                elif level==1:
                    ax.text(lypos,xpos, label, ha='center', transform=ax.transAxes, size=10)
                elif level ==2:
                    if rpos > 1:
                        for i in range(1, rpos+1):
                            lypos = (pos + i)*scale -0.02
                            ax.text(lypos,xpos, label, ha='right', transform=ax.transAxes, size=10)
                    else:
                        ax.text(lypos,xpos, label, ha='center', transform=ax.transAxes, size=10)

                if level == 0:
                    add_hline(ax, pos*scale, xpos*scale+1.06,rpos*scale)
                elif level == 1:
                    add_hline(ax, pos*scale, xpos*scale+0.99,rpos*scale)

            xpos += .09

    def add_support_values(ax, df):
        ypos = -0.05
        xpos = 0.045
        for i, row in df.iterrows():
            support = f'{row["support"]}%'
            ax.text(xpos,ypos, support, ha='center', transform=ax.transAxes, size=10)
            xpos += .083

    #intiliazile figure 
    fig = sns.clustermap(stat_plot, row_cluster=False,col_cluster=False,vmin=0, vmax=1.0,cmap=pvalue_colmap,
                     linewidths=0.00000001,figsize=(14,7))
    
    #cosmetics 
    #restyle row labeling 
    fig.ax_heatmap.yaxis.set_ticks_position("left")
    plt.setp(fig.ax_heatmap.xaxis.get_majorticklabels(), fontsize=10)
    fig.ax_heatmap.set_yticks([])

    #restyle column labeling
    fig.ax_heatmap.xaxis.set_ticks_position("top")
    plt.setp(fig.ax_heatmap.xaxis.get_majorticklabels(), fontsize=10)
    fig.ax_heatmap.set_xticks([])
    label_columns(fig.ax_heatmap, stat_plot)

    add_support_values(fig.ax_heatmap, support_values)

    fig.ax_row_dendrogram.set_visible(False)
    fig.ax_col_dendrogram.set_visible(False)
    fig.ax_heatmap.set(xlabel=None, ylabel=None)
    plt.subplots_adjust(left=0.2)

    space_before = [('$t_{0}$', 'epitope', 'A vs. \nB & C'), ('$t_{0}$', 'non-epitope', 'A vs. \nB & C'),('$t_{onset}$', 'all', 'A vs. B'),
                    ('$t_{1}$', 'all', 'A vs. B')]
    for i,col in enumerate(stat_plot.columns):

        if col in space_before:
            fig.ax_heatmap.axvline(i, color='white', lw=2)
            if '$t_{f_{0}}$' in col:
                fig.ax_heatmap.axvline(i, color='white', lw=3)

    bounds = [0.0,0.01,0.02,0.03,0.04,0.05,0.2,0.4,0.6,0.8,1.0]#[0.0,0.05,1.0]
    norm = BoundaryNorm(bounds, pvalue_colmap.N)
    cb = mcp.ColorbarBase(fig.cax, cmap=pvalue_colmap, boundaries=bounds, ticks=bounds)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=9)

    fig.ax_cbar.set_position((0.975, .43, .01, .25))
    fig.ax_cbar.set_title("P value", size=10, loc="left", weight="bold")

    fig.ax_heatmap.text(12.2,-18, "Time point", ha='left', size=10, weight="bold")
    fig.ax_heatmap.text(12.2,-11, "Substitutions", ha='left', size=10, weight="bold")
    fig.ax_heatmap.text(12.2,-3, "Comparison", ha='left', size=10, weight="bold")

    fig.ax_heatmap.text(12.2,105, "Support*", ha='left', size=10, weight="bold")
    fig.ax_heatmap.text(12.2,110, "*(P value<0.05)", ha='left', size=10)

    return fig 

def simple_LBI_heatmap_plot(stat_plot, support_values):

    def label_columns(ax, df):
        xpos = 1.005
        scale = 1./df.columns.size
        for level in reversed(range(df.columns.nlevels)):
            pos = df.columns.size
            labels = [(k, sum(1 for i in g)) for k,g in itertools.groupby(df.columns.get_level_values(level))]

            for label, rpos in reversed(labels):
                pos -= rpos
                lypos = (pos + rpos/2)*scale 
                
                if level == 0:
                    if label.startswith("t"):
                        label = label.replace("t","$t_{") + "}$"
                        
                        if "f" in label:
                            label = label.replace("}$", "}}$").replace("f", "f_{")
                    ax.text(lypos,xpos, label, ha='center', transform=ax.transAxes, size=10)

                elif level ==2:
                    if rpos > 1:
                        for i in range(1, rpos+1):
                            lypos = (pos + i)*scale -0.02
                            ax.text(lypos,xpos, label, ha='right', transform=ax.transAxes, size=8)
                    else:
                        ax.text(lypos,xpos, label, ha='center', transform=ax.transAxes, size=8)

                if level == 0:
                    add_hline(ax, pos*scale, xpos*scale+.8,rpos*scale)
               
            if level != 1:
                xpos += .09

    def add_support_values(ax, df):
        ypos = -0.05
        xpos = 0.15
        for i, row in df.iterrows():
            support = f'{row["support"]}%'
            ax.text(xpos,ypos, support, ha='center', transform=ax.transAxes, size=8)
            xpos += .24

    #intiliazile figure 
    fig = sns.clustermap(stat_plot, row_cluster=False,col_cluster=False,vmin=0, vmax=1.0,cmap=pvalue_colmap,
                     linewidths=0.00000000000001,figsize=(4.5,6.5))

    #cosmetics
    #restyle row labeling 
    fig.ax_heatmap.yaxis.set_ticks_position("left")
    plt.setp(fig.ax_heatmap.xaxis.get_majorticklabels(), fontsize=10)
    fig.ax_heatmap.set_yticks([])

    #restyle column labeling
    fig.ax_heatmap.xaxis.set_ticks_position("top")
    plt.setp(fig.ax_heatmap.xaxis.get_majorticklabels(), fontsize=10)
    fig.ax_heatmap.set_xticks([])
    label_columns(fig.ax_heatmap, stat_plot)

    add_support_values(fig.ax_heatmap, support_values)

    fig.ax_row_dendrogram.set_visible(False)
    fig.ax_col_dendrogram.set_visible(False)
    fig.ax_heatmap.set(xlabel=None, ylabel=None)
    plt.subplots_adjust(left=0.2)

    space_before = [('$t_{0}$', 'epitope', 'A vs. \nB & C'), ('$t_{0}$', 'non-epitope', 'A vs. \nB & C'),('$t_{onset}$', 'all', 'A vs. B'),
                    ('$t_{1}$', 'all', 'A vs. B')]
    for i,col in enumerate(stat_plot.columns):

        if col in space_before:
            fig.ax_heatmap.axvline(i, color='white', lw=2)
            if '$t_{f_{0}}$' in col:
                fig.ax_heatmap.axvline(i, color='white', lw=3)

    bounds = [0.0,0.01,0.02,0.03,0.04,0.05,0.2,0.4,0.6,0.8,1.0]#[0.0,0.05,1.0]
    norm = BoundaryNorm(bounds, pvalue_colmap.N)
    cb = mcp.ColorbarBase(fig.cax, cmap=pvalue_colmap, boundaries=bounds, ticks=bounds)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=8)

    fig.ax_cbar.set_position((0.9, .43, .03, .25))
    fig.ax_cbar.set_title("P value", size=9, loc="left", weight="bold")

    fig.ax_heatmap.text(4.2,-11, "Time point", ha='left', size=9, weight="bold")
    fig.ax_heatmap.text(4.2,-3, "Comparison", ha='left', size=9, weight="bold")

    fig.ax_heatmap.text(4.2,105, "Support*", ha='left', size=9, weight="bold")
    fig.ax_heatmap.text(4.2,108, "*(P value<0.05)", ha='left', size=6)

    return fig

def main():
    args = ArgumentParser()

    #get threads:
    threads = args.threads if args.threads else 1
    print("Number of threads is: ", threads) 

    #get data folder
    datafolder = os.path.join(cwd, args.data_folder)
    if not os.path.isdir (datafolder):
        sys.stderr.write(f"Error: cannot find{datafolder}. Check if the directory is properly specified.")
        sys.exit(-1)
    
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

    #get segment and subtype
    segment = args.segment
    subtype = args.subtype

    #get cutoff
    minor_cutoff = args.minor_cutoff
    fixed_cutoff = args.fixed_cutoff

    #get epitope site if requested
    epitope_analysis = args.epitope
    if epitope_analysis and segment in ["HA", "NA"]:
        epitope_sites = all_epitope_sites[subtype][segment]
        epitope_positions =[site for site in epitope_sites.values()]
        epitope_positions = sorted([pos for pos_list in epitope_positions for pos in pos_list])

    #get analysis
    analyses_dirs = {}
    for f in os.listdir(datafolder):
        if os.path.isdir(os.path.join(datafolder,f)):
            try:
                analysis = int(f.split(".")[0].split("_")[-1])
            except:
                sys.stderr.write(f"Error: analyses must be numbered at the end: e.g. H3N2_HA_1, protein_29")
                sys.exit(-1)
            if f.startswith("protein"):
                analyses_dirs[analysis] = os.path.join(datafolder, f)
            else: #expecting direct treason output
                for f2 in os.listdir(os.path.join(datafolder, f)):
                    if f2.startswith("protein"):
                        analyses_dirs[analysis] = os.path.join(datafolder, f,f2)
                        break
            if analysis not in analyses_dirs.keys():
                sys.stderr.write(f"Error: could not find protein directory for {f}. check and try again")
                sys.exit(-1)
    analyses = dict(OrderedDict(sorted(analyses_dirs.items())))
                
    #to do: sequence origin analysis here later on
    freq_file = os.path.join(output, "mutation_frequencies.csv")
    lbi_file = os.path.join(output, "mutation_lbi.csv")
    stat_file = os.path.join(output, "lbi_stats.csv")
    
    #get mature_sequences and trees per analysis
    start_pos = starting_positions[subtype][segment]
    
    #get mutation frequencies
    mutation_freqs = []
    if not os.path.isfile(freq_file) or args.redo:
        print (f"determining frequency trajectories for {len(analyses)} replicates")
        with mp.Pool(processes=threads) as p:
            for result in tqdm([p.apply_async(get_frequency, args=(fdir, analysis, minor_cutoff, start_pos)) for analysis, fdir in analyses.items()]):
                try:
                    mutation_freqs = pd.concat([mutation_freqs, result.get()])
                except:
                    mutation_freqs = result.get()
        mutation_freqs.to_csv(freq_file)
        print ("finished determining frequency trajectories")
    else: 
        mutation_freqs = pd.read_csv(freq_file).drop(["Unnamed: 0"], axis=1)
        #quick check if load file matches input file
        if len([a for a in mutation_freqs["analysis"].unique() if a not in analyses.keys()]) > 0 or \
            len([a for a in analyses.keys() if a not in mutation_freqs["analysis"].unique()]) >0 :
            
                sys.stderr.write(f"Error: mutation_frequencies.csv found in output directory does not match the analysis found in the input directory. Please specify '--redo")
                sys.exit(-1)

    #remove egg adaptive mutations > HA
    if segment == "HA" and subtype=="H3N2":
        egg_adaptive_mutations = [i for t,l in all_egg_adaptive_mutations[subtype][segment].items() for i in l]
        mutation_freqs = mutation_freqs[~mutation_freqs["mutation"].isin(egg_adaptive_mutations)]

    #get LBI distributions per mutation per analysis
    if not os.path.isfile(lbi_file) or args.redo:
        print (f"running LBI analysis for {len(analyses)} replicates")
        with mp.Pool(processes=threads) as p:
            if epitope_analysis:
                for result in tqdm([p.apply_async(run_lbi_analysis, args=(a,fdir, mutation_freqs[mutation_freqs["analysis"]==a].reset_index(), start_pos, minor_cutoff, fixed_cutoff, epitope_positions)) for a,fdir in analyses.items()]):
                    try:
                        mutation_lbis = pd.concat([mutation_lbis, result.get()])
                    except:
                        mutation_lbis = result.get()
            else:
                for result in tqdm([p.apply_async(run_lbi_analysis, args=(a,fdir, mutation_freqs[mutation_freqs["analysis"]==a].reset_index(), start_pos, minor_cutoff, fixed_cutoff )) for a,fdir in analyses.items()]):
                    try:
                        mutation_lbis = pd.concat([mutation_lbis, result.get()])
                    except:
                        mutation_lbis = result.get()
        mutation_lbis.to_csv(lbi_file)
        print ("finished LBI analysis")
    else:
        mutation_lbis = pd.read_csv(lbi_file,dtype={'season': str, 't_0':str, 't_1':str, 't_d':str}).drop(["Unnamed: 0"], axis=1)
        #quick check if load file matches input file
        if len([a for a in mutation_lbis["analysis"].unique() if a not in analyses.keys()]) > 0 or \
            len([a for a in analyses.keys() if a not in mutation_lbis["analysis"].unique()]) >0 :
            
                sys.stderr.write(f"Error: mutation_lbis.csv found in output directory does not match the analysis found in the input directory. Please specify '--redo")
                sys.exit(-1)

    
    #get LBI stats
    if not os.path.isfile(stat_file) or args.redo:
        print  (f"running statistical analysis for {len(analyses)} replicates")
        with mp.Pool(processes=threads) as p:
            for result in tqdm([p.apply_async(statistical_lbi_analysis, args=(a, mutation_lbis[mutation_lbis["analysis"]==a], epitope_analysis, ["t_0", "t_1"])) for a in analyses.keys()]):
                try:
                    stat_results = pd.concat([stat_results, result.get()])
                except:
                    stat_results = result.get()
        stat_results.to_csv(stat_file, index=False)
    else:
        stat_results = pd.read_csv(stat_file)
        if epitope_analysis is True and "epitope" not in list(stat_results["mutation type"].unique()):
            print ("redoing statistical analysis requested output does not match found results")
            with mp.Pool(processes=threads) as p:
                for result in tqdm([p.apply_async(statistical_lbi_analysis, args=(a, mutation_lbis[mutation_lbis["analysis"]==a], epitope_analysis, ["t_0", "t_1"])) for a in analyses.keys()]):
                    try:
                        stat_results = pd.concat([stat_results, results.get()])
                    except:
                        stat_results = result.get()
            stat_results.to_csv(stat_file, index=False)
        elif epitope_analysis is False and "epitope" in list(stat_results["mutation type"].unique()):
            stat_results = stat_results[stat_results["mutation type"]=="all"]

    #check if plots need to be made:
    if args.plot:

        ####### random frequency trajectory plot (1/100)
        a =random.choice(list(analyses.keys()))
        mut_freqs = mutation_freqs[mutation_freqs["analysis"]==a]
        seasons = [c for c in mut_freqs.columns if c not in ["mutation", "analysis"]]
        mut_classification = classify_trajectories(mut_freqs, seasons, minor_cutoff, fixed_cutoff)
        plot_file = os.path.join(output,f"frequency_trajectory_{a}.png")
        fig = frequency_trajectory_plot(mut_freqs, mut_classification, seasons, minor_cutoff)
        fig.savefig(plot_file,  dpi=600)


        ###### LBI heat map
        #get support values
        support_values = calculate_support_values(stat_results)

        #get support values ready for plotting
        support_values["time point"] = pd.Categorical(support_values["time point"], categories=["$t_{0}$", "$t_{1}$"])
        support_values["mutation type"] = pd.Categorical(support_values["mutation type"], categories=["all", "epitope", "non-epitope"])
        support_values["groups compared"] = pd.Categorical(support_values["groups compared"], categories=["A vs. \nB & C", "B vs. \nA & C", "C vs. \nA & B", "A vs. B"])
        support_values = support_values.sort_values(["time point", "mutation type", "groups compared"]).set_index(["time point", "mutation type", "groups compared"])

        #get stat_results ready for plotting 
        stat_plot = stat_results.copy()
        stat_plot["time point"] = pd.Categorical(stat_plot["time point"], categories=["$t_{0}$", "$t_{1}$"])
        stat_plot["groups compared"] = pd.Categorical(stat_plot["groups compared"], categories=["A vs. \nB & C", "B vs. \nA & C",  "C vs. \nA & B", "A vs. B"])
        stat_plot = stat_plot.pivot(index=["analysis"], columns=["time point", "mutation type", "groups compared"], values="p-value").sort_values(by=["time point", "mutation type", "groups compared" ], axis=1)

        ## make plot
        plot_file  = os.path.join(output,"LBI_analysis.png")
        if epitope_analysis:
            fig = epitope_LBI_heatmap_plot(stat_plot, support_values)
            fig.savefig(plot_file,  dpi=600)
        else:
            fig = simple_LBI_heatmap_plot(stat_plot, support_values)
            fig.savefig(plot_file,  dpi=600)

if __name__ == "__main__":
    main()



