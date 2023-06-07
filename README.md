<!-- omit in toc -->
# evoFLUtion - Predicting the Evolutionary Fate of Hemagglutinin Amino Acid Substitutions in Seasonal Influenza A/H3N2 Viruses

Scripts and data for the analysis and influenza A/H3N2 (HA) amino acid substitutions.  
Analysis is comprised out of three parts:

1. Sequence download from GISAID (which needs to be performed manually)
2. Automated molecular clock phylogeny construction and LBI calculation per season via Treason
3. Results analysis and visualization via jupyter Notebook

- [Requirements](#requirements)
- [Sequence download](#sequence-download)
- [Treason - Molecular clock phylogenetic trees per season](#treason---molecular-clock-phylogenetic-trees-per-season)
- [Substitution frequency trajectories and LBI distributions](#substitution-frequency-trajectories-and-lbi-distributions)
- [Isolates used in the original analysis](#isolates-used-in-the-original-analysis)
- [References](#references)

## Requirements

All tools and package needed for analysis are included in [conda enviroment setup file](environment.yml), so some version of conda is required.  
The conda environment can be setup via ``conda env create -f environment.yml``.  
The conda environment can be activated via ``conda activate flu``. 
The conda environment can be updated via ``conda env update -f environment.yml``.

## Sequence download

For the analysis performed in the manuscript we downloaded all sequences and metadata that resulted from the following search in the GISAID EpiFlu Database: Type=A, H=3, N=2, Host=Human, Location=-all-, Collection date From=2000-01-01, Collection date To=2019-12-31, segment=HA. All other settings were left at default.  
Since GISAID only allows a maximum download of 20000 sequences, we downloaded the sequences in batches.  

The strains/isolates/sequences used for this study are listed under [isolate_analysis.csv](data/isolates_analysis.csv)

## Treason - Molecular clock phylogenetic trees per season

To construct molecular clock phylogenetic trees per season Treason can be run via `python scripts/treason.py`.  
We ran treason with the following options `python scripts/treason.py -d {gisaid data dir} -o {output dir} -s HA -p -on -ss -cs -nc`

```.
usage: treason.py [-h] -d DATA_FOLDER [-o OUTPUT] [-s SEGMENTS [SEGMENTS ...]] [-p] [-on] [-mxa MAX_AMBIG] [-ml MIN_LENGTH] [-i INTERVAL] [-sm START_MONTH] [-ss] [-pcm PER_COUNTRY_MONTH] [-ips] [-icb]
                  [-ieb] [-cs] [-c SEED] [-nc] [-r] [-ra] [-t THREADS]
                  [-f {FilterClinical,GetSeason,MSA,PhyloTree,TreeTime,FilterMolecularClockOutliers,RedoMSA,RedoPhyloTree,RedoTreeTime,LBI,TranslateMutations,Translate}] [-v]

Filter clinical A/H3N2 samples from GISAID and construct phylogenetic tree per season (and previous season)

options:
  -h, --help            show this help message and exit
  -d DATA_FOLDER, --data-folder DATA_FOLDER
                        data folder where raw GISAID sequence and metadata are stored
  -o OUTPUT, --output OUTPUT
                        output directory to store all generated output file (default: ./)
  -s SEGMENTS [SEGMENTS ...], --segments SEGMENTS [SEGMENTS ...]
                        segment abbreviations that will be analyzed
  -p, --protein         if sequences need be translated into protein sequences (coding region only) and if mutations in final LBI need to be translated
  -on, --only-nonsyn    if '-p' flag is specified, only report non-synonymous mutations in tree files
  -mxa MAX_AMBIG, --max-ambig MAX_AMBIG
                        maximum percentage of ambiguous nucleotides allowed (default: 0.01)
  -ml MIN_LENGTH, --min-length MIN_LENGTH
                        miminum percentage of length w.r.t. the reference segment
  -i INTERVAL, --interval INTERVAL
                        interval period in years form which individual season analyses need to be made (default: 2015-2019)
  -sm START_MONTH, --start-month START_MONTH
                        name of the month from which the season should start (default: may)
  -ss, --sub-sample     If season needs to be down sample for a max number of sequences/country/month (specify '-pcm')
  -pcm PER_COUNTRY_MONTH, --per-country-month PER_COUNTRY_MONTH
                        maximum number of sequence per country per month (default: 10)
  -ips, --include-previous-season
                        if the sequence of the previous season should be included in the next season
  -icb, --include-cell-based
                        include cell based isolates there are less than 500 (clinical) sequences per season
  -ieb, --include-egg-based
                        Include cell based isolates there are less than 500 (clinical) sequences per season (if icb is specified cell-based isolates will be prioritized)
  -cs, --complete-subset
                        included sequences in the order of clinical, cell-based, egg-based, remaining if the number sequences in the season subset are less than 500
  -c SEED, --seed SEED  seed to use for random sampling (default:29)
  -nc, --no-seed        If flag is specified no seed will be used
  -r, --redo            Redo all steps except the initial merging step if data already exists
  -ra, --redo-all       Redo ALL steps if data already exists
  -t THREADS, --threads THREADS
                        Number of threads (default 1)
  -f {FilterClinical,GetSeason,MSA,PhyloTree,TreeTime,FilterMolecularClockOutliers,RedoMSA,RedoPhyloTree,RedoTreeTime,LBI,TranslateMutations,Translate}, --force-rule {FilterClinical,GetSeason,MSA,PhyloTree,TreeTime,FilterMolecularClockOutliers,RedoMSA,RedoPhyloTree,RedoTreeTime,LBI,TranslateMutations,Translate}
                        Force execution of specific snakemake rule
  -v, --verbose         print a bunch of stuff
```

## Substitution frequency trajectories and LBI distributions

To computed substitutions as well as their trajectories and LBI distributions from the Treason output, the [LBI_analysis jupyter notebook](notebooks/LBI_analysis.ipynb) can be used. All plots from the study can be produced using this notebook.

## Isolates used in the original analysis

All GISAID ids per bootstrap replicate per season for the entire analysis are listed in [isolates_analysis.csv](data/isolates_analysis.csv).

## References

|Name |Publication|Website/link|
|:---|:---|:---|
|Biopython|Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., ... & De Hoon, M. J. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422-1423.|<https://biopython.org/>|
|DendroPY|Sukumaran, J. and Mark T. Holder. 2010. DendroPy: A Python library for phylogenetic computing. Bioinformatics 26: 1569-1571.|<https://dendropy.org/>|
|Conda|NA|<https://conda.io/>|
|Git|NA|<https://git-scm.com/>|
|IPython|Perez, F., Granger, B. (2007). IPython: A System for Interactive Scientific Computing|<https://ipython.org/>|
|IQ-TREE|Nguyen, L.-T., Schmidt, H. A., von Haeseler, A., Minh, B. Q. (2015). IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies. Molecular Biology and Evolution, 32, 268-274. |<http://www.iqtree.org/#download>|
|Jupyter Notebooks|Kluyver, T., Ragan-Kelley, B., Granger, F., Bussonnier, M., Frederic, J., Kelley, K. ... & Willing, C. (2016). Jupyter Notebooks - a publishing format for reproducible computational workflows|<https://jupyter.org/>|
|MAFFT| Katoh, K., Toh, H. (2008). Recent developments in the MAFFT multiple sequence alignment program. Briefings in bioinformatics 2008, 9(4), 286-298.|<https://mafft.cbrc.jp/alignment/software/>|
|Matplotlib|Hunter, J. D. (2007). Matplotib: A 2D Graphics Environment|<https://matplotlib.org/>|
|Python|G. van Rossum, Python tutorial, Technical Report CS-R9526, Centrum voor Wiskunde en Informatica (CWI), Amsterdam, May 1995.|<https://www.python.org/>|
|Seaborn|Waskom, M. L. (2021). seaborn: statistical data visualization |<https://seaborn.pydata.org/>|
|SH-aLRT|Guindon, Stéphane, et al. New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. Systematic biology 59.3 (2010): 307-321.|<https://doi.org/10.1093/sysbio/syq010>|
|Snakemake|Köster, J. and S.J.B. Rahmann, Snakemake—a scalable bioinformatics workflow engine. 2012. 28(19): p. 2520-2522.|<https://snakemake.readthedocs.io/en/stable/>|
|UFBoot|Hoang, D.T., Chernomor, O., von Haeseler, A., Minh, B. Q., Vinh, L. S. (2018). UFBoot2: Improving the ultrafast bootstrap approximation. Molecular Biology and Evolution, 35, 518–522.|<https://doi.org/10.1093/molbev/msx281>|