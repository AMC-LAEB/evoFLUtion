#!/usr/bin/env python3

import os, sys, random, copy, dendropy
import pandas as pd
import numpy as np
from Bio import SeqIO, Seq, SeqRecord
from datetime import datetime, date
import multiprocessing as mp 
from passage_labels import *

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

def redo_fasta_header(header,segment):
    """
    change fasta header from >Isolate_ID|Segment_ID|Strain_name|type|segment to 
    >isolate_ID|Segment_ID|strain_name|segment (no _ in the other individual fields)
    date will be as float and strain name won't contain any special character
    input:
        - header: sequence header that needs to be changed (str)
    """
    if len(header.split("|")) == 5:
        isolate,segid,strain,t,segment = header.split("|") #get the individual items 
    elif len(header.split("|")) ==3:
        isolate,segid,strain =  header.split("|")
    else: #to avoid problems
        return header

    
    #check if there are special characters in the strain 
    if "(" in strain: #something in parenthesis strain name 
        strain = strain.split("(")[0] #only taking the first part of the strain name
    if any(not c.isalnum() for c in strain):
        strain = ''.join(c for c in strain if c.isalnum() or c=="/")
    
    return "|".join([isolate,segid,strain,segment])

def add_date_to_header(header, d, sep="|"):
    """
    input:
        - d: string of data in format yyyy-mm-dd (str)
    """
    #get year
    if "-" in d:
        year = int(d.split("-")[0])
    else:
        year = int(d) #assuming that there is only a year present in this case
        return f"{header}|{str(year)}"
    
    #convert date to float
    if len(d) > 7: #date
        ym = datetime.strptime(d,'%Y-%m-%d').timetuple().tm_yday
    else:
        ym = datetime.strptime(d,'%Y-%m').timetuple().tm_yday
    if ((year % 400 == 0) and (year % 100 == 0)) or ((year % 4 ==0) and (year % 100 != 0)): #leap year
        dfloat = ym/366
    else:
        dfloat = ym/365
    
    floatd = year + dfloat
    
    return f"{header}|{str(floatd)}"

def filter_fasta(clinical,fasta,rh=True, adth=True, mxa=0.01, mnlp=0.95, message="clinical"):
    """
    filter fasta files per segment 
    input:
        - clinical: dataframe containing metadata for all clinical gisaid sequences (pandas dataframe)
        - fasta: fasta file containing sequences that need to be filtered
        for other parameters see filter_clinical
    output:
        - to_drop: isolate_Id of isolates that need to be dropped from further analyses (set)
        for other output see filter_clinical
    """
    seen = set() #records that have been observed to avoid duplicates
    to_drop = set() #isolates that did not pass QC
        
    #get segment 
    segment = fasta.split("/")[-1].split("_")[0]
    print (f"filtering {message} sequences for {segment}")
        
    clinical_records = []
        
    #read fasta
    for record in SeqIO.parse(fasta,"fasta"):
        isolate = record.id.split("|")[0]
        segid = record.id.split("|")[1]
            

        if isolate in list(clinical["Isolate_Id"]) and (isolate not in to_drop) and (segid not in seen):
            rdf = clinical[clinical["Isolate_Id"]==isolate]
            #correct version on the segment (there could be multiple uploads)
            if rdf[f"{segment} Segment_Id"].iloc[0].split("|")[0].lstrip("EPI")==segid: 

                #check if sequence fullfills the minimum/maximum requirements
                if check_sequence_length(segment, str(record.seq), mnlp) and check_max_ambig(str(record.seq), mxa): 
                        
                    #renaming segment id to match with metadata
                    record.id = record.id.replace(segid,f"EPI{segid}")
                    record.name = record.name.replace(segid,f"EPI{segid}")
                    record.description = record.description.replace(segid,f"EPI{segid}")
                    
                    #edit record if needed 
                    if rh:
                        record.id = redo_fasta_header(str(record.id), segment)
                        record.name = record.id
                        record.description = record.id
                    if adth:
                        for i, v in enumerate(rdf["Collection_Date"]):
                            if i == 0:
                                date = v
                        #record.id = add_date_to_header(str(record.id), date)
                        record.id = f"{str(record.id)}|{date}"

                    clinical_records.append(record)
                    
                else: #dropping if requirements aren't met
                    to_drop.add(isolate)
        seen.add(segid)
    print (f"finished {message} sequences for {segment}")
    return segment, clinical_records, to_drop

def remove_failed(segment, clinical, records):
    """
    remove failed sequences from > isolate will have failed for another segment and has not been removed yet
    input:
        - segment: segment for which failed sequences need to be removed (str)
        - clinical: list of isolate_Ids of sequence that passed the quality control (list)
        - records: records that need to be filtered (list)
    output:
        - segment: segment in question (for pooling)
        - filtered: filtered records (list)
    """
    filtered = []
    for record in records:
        if record.id.split("|")[0] in clinical:
            filtered.append(record)
    return segment, filtered

def filter_clinical(metadata,sequence_files, ncpu=1, rh=True, adth=True, mxa=0.01, mnlp=0.95):
    """
    Filter clinical records from GISAID metadata file per Isolate, if one of the segments doesn't meet 
    the requirements entire isolate will be removed

    '!' assuming to FASTA header looks like >Isolate_ID|Segment_ID|Strain_name|type|segment
    '!' assuming that sequence file name start with segment indication 
    input:
        - metadata: GISAID metadata file containing passage history for all samples (pandas dataframe)
        - sequence_files: list with GISAID raw FASTA files from which the clinical sequences need to be filtered (list)
        - rh: redo header, if change header must changed with redo_header_func (Bool)
        - adth: if sequencing date should be add to the header (date will be a float)
        - mnlp: minimum length percentage for the sequence(float)
        - mxa: maximum ambiguous percentage for the sequence(float)
    output:
        - clinical: dataframe containing GISAID metadata of all clinical isolate that fullfilled the criteria for all segmenets
        - clinical_records: dictionary with the clinical sequence record per segments
    """
    #get clincical samples from the metadata 
    clinical = metadata[metadata["Passage_History"].isin(cpl)].reset_index()

    #also checking if collection date consist of at least year and month
    to_drop = set()
    for i, v in enumerate(clinical["Collection_Date"]):
        try:
            datetime.strptime(v, '%Y-%m-%d')
        except:
            to_drop.add(i)
    clinical = clinical.drop(to_drop).reset_index()

    #and location at least consists of country and continent
    to_drop = set()
    for i, v in enumerate(clinical["Location"]):
        if len(v.split(" / ")) < 2:
            to_drop.add(i)
    clinical = clinical.drop(to_drop)#.reset_index()
    

    #rename MP segment to M, if this hasn't been done
    if "MP Segment_Id" in list(clinical.columns):
        clinical = clinical.rename(columns={"MP Segment_Id": "M Segment_Id"})

    #there are no clinical records found return none 
    if len(clinical) == 0: 
       return None 
    
    clinical_records = {}
    to_drop = set() #isolates that do not match the requirements 

    #filter clinical sequences 
    pool = mp.Pool(processes=ncpu) 
    results = [pool.apply_async(filter_fasta, args=(clinical, fasta, rh, adth, mxa, mnlp)) for fasta in sequence_files]
    output = [p.get() for p in results]
    pool.close()
    pool.join()

    for output_res in output:
        segment, clin, drop  = output_res
        clinical_records[segment] = clin
        to_drop.update(drop)

    #filter clinical with records that need to be removed
    clinical = clinical[~clinical["Isolate_Id"].isin(to_drop)]
    
    #filter clinical record from to drop
    pool = mp.Pool(processes=ncpu)
    results = [pool.apply_async(remove_failed, args=(segment, list(clinical["Isolate_Id"]), records)) for segment, records in clinical_records.items()]
    output = [p.get() for p in results]
    pool.close()
    pool.join()

    for output_res in output:
        segment, filtered = output_res
        clinical_records[segment] = filtered
        
    return clinical, clinical_records

def filter_cell_based(metadata,sequence_files, ncpu=1, rh=True, adth=True, mxa=0.01, mnlp=0.95):
    """
    Filter cell_based records from GISAID metadata file per Isolate, if one of the segments doesn't meet 
    the requirements entire isolate will be removed

    '!' assuming to FASTA header looks like >Isolate_ID|Segment_ID|Strain_name|type|segment
    '!' assuming that sequence file name start with segment indication 
    input:
        - metadata: GISAID metadata file containing passage history for all samples (pandas dataframe)
        - sequence_files: list with GISAID raw FASTA files from which the cell based sequences need to be filtered (list)
        - rh: redo header, if change header must changed with redo_header_func (Bool)
        - adth: if sequencing date should be add to the header (date will be a float)
        - mnlp: minimum length percentage for the sequence(float)
        - mxa: maximum ambiguous percentage for the sequence(float)
    output:
        -cell_based: dataframe containing GISAID metadata of all cell based isolates that fullfilled the criteria for all segmenets
        -cell_based_records: dictionary with the cellbased sequence record per segments
    """
    #get cell based samples from the metadata 
    cell_based = metadata[metadata["Passage_History"].isin(cbpl)].reset_index(drop=True) 

    #also checking if collection date consist of at least year and month
    to_drop = set()
    for i, v in enumerate(cell_based["Collection_Date"]):
        try:
            datetime.strptime(v, '%Y-%m-%d')
        except:
            to_drop.add(i)
    cell_based = cell_based.drop(to_drop).reset_index()

    #and location at least consists of country and continent
    to_drop = set()
    for i, v in enumerate(cell_based["Location"]):
        if len(v.split(" / ")) < 2:
            to_drop.add(i)
    cell_based = cell_based.drop(to_drop)#.reset_index()

    #rename MP segment to M, if this hasn't been done
    if "MP Segment_Id" in list(cell_based.columns):
        cell_based = cell_based.rename(columns={"MP Segment_Id": "M Segment_Id"})

    #there are no cell_based records found return none 
    if len(cell_based) == 0: 
       return None 
    
    cell_based_records = {}
    to_drop = set() #isolates that do not match the requirements 

    #filter cell_based sequences 
    pool = mp.Pool(processes=ncpu) 
    results = [pool.apply_async(filter_fasta, args=(cell_based, fasta, rh, adth, mxa, mnlp, "cell based")) for fasta in sequence_files]
    output = [p.get() for p in results]
    pool.close()
    pool.join()

    for output_res in output:
        segment, cbased, drop  = output_res
        cell_based_records[segment] = cbased
        to_drop.update(drop)

    #filter cell_based with records that need to be removed
    cell_based = cell_based[~cell_based["Isolate_Id"].isin(to_drop)]
    
    #filter cell_based record from to drop
    pool = mp.Pool(processes=ncpu)
    results = [pool.apply_async(remove_failed, args=(segment, list(cell_based["Isolate_Id"]), records)) for segment, records in cell_based_records.items()]
    output = [p.get() for p in results]
    pool.close()
    pool.join()

    for output_res in output:
        segment, filtered = output_res
        cell_based_records[segment] = filtered
        
    return cell_based, cell_based_records

def filter_egg_based(metadata,sequence_files, ncpu=1, rh=True, adth=True, mxa=0.01, mnlp=0.95):
    """
    Filter egg_based records from GISAID metadata file per Isolate, if one of the segments doesn't meet 
    the requirements entire isolate will be removed

    '!' assuming to FASTA header looks like >Isolate_ID|Segment_ID|Strain_name|type|segment
    '!' assuming that sequence file name start with segment indication 
    input:
        - metadata: GISAID metadata file containing passage history for all samples (pandas dataframe)
        - sequence_files: list with GISAID raw FASTA files from which the cell based sequences need to be filtered (list)
        - rh: redo header, if change header must changed with redo_header_func (Bool)
        - adth: if sequencing date should be add to the header (date will be a float)
        - mnlp: minimum length percentage for the sequence(float)
        - mxa: maximum ambiguous percentage for the sequence(float)
    output:
        -egg_based: dataframe containing GISAID metadata of all egg-based isolates that fullfilled the criteria for all segmenets
        -egg_based_records: dictionary with the egg-based sequence record per segments
    """
    #get egg based samples from the metadata 
    egg_based = metadata[metadata["Passage_History"].isin(epl)].reset_index(drop=True) 

    #also checking if collection date consist of at least year and month
    to_drop = set()
    for i, v in enumerate(egg_based["Collection_Date"]):
        try:
            datetime.strptime(v, '%Y-%m-%d')
        except:
            to_drop.add(i)
    egg_based = egg_based.drop(to_drop).reset_index()

    #and location at least consists of country and continent
    to_drop = set()
    for i, v in enumerate(egg_based["Location"]):
        if len(v.split(" / ")) < 2:
            to_drop.add(i)
    egg_based = egg_based.drop(to_drop)#.reset_index()

    #rename MP segment to M, if this hasn't been done
    if "MP Segment_Id" in list(egg_based.columns):
        egg_based = egg_based.rename(columns={"MP Segment_Id": "M Segment_Id"})

    #there are no cell_based records found return none 
    if len(egg_based) == 0: 
       return None 
    
    egg_based_records = {}
    to_drop = set() #isolates that do not match the requirements 

    #filter cell_based sequences 
    pool = mp.Pool(processes=ncpu) 
    results = [pool.apply_async(filter_fasta, args=(egg_based, fasta, rh, adth, mxa, mnlp, "egg based")) for fasta in sequence_files]
    output = [p.get() for p in results]
    pool.close()
    pool.join()

    for output_res in output:
        segment, ebased, drop  = output_res
        egg_based_records[segment] = ebased
        to_drop.update(drop)

    #filter cell_based with records that need to be removed
    egg_based = egg_based[~egg_based["Isolate_Id"].isin(to_drop)]
    
    #filter cell_based record from to drop
    pool = mp.Pool(processes=ncpu)
    results = [pool.apply_async(remove_failed, args=(segment, list(egg_based["Isolate_Id"]), records)) for segment, records in egg_based_records.items()]
    output = [p.get() for p in results]
    pool.close()
    pool.join()

    for output_res in output:
        segment, filtered = output_res
        egg_based_records[segment] = filtered
        
    return egg_based, egg_based_records

def filter_remaining(metadata,sequence_files, ncpu=1, rh=True, adth=True, mxa=0.01, mnlp=0.95):
    """
    Filter non-clinical, non-cell-based, and non-egg-based records from GISAID metadata file per Isolate, if one of the segments doesn't meet 
    the requirements entire isolate will be removed

    '!' assuming to FASTA header looks like >Isolate_ID|Segment_ID|Strain_name|type|segment
    '!' assuming that sequence file name start with segment indication 
    input:
        - metadata: GISAID metadata file containing passage history for all samples (pandas dataframe)
        - sequence_files: list with GISAID raw FASTA files from which the cell based sequences need to be filtered (list)
        - rh: redo header, if change header must changed with redo_header_func (Bool)
        - adth: if sequencing date should be add to the header (date will be a float)
        - mnlp: minimum length percentage for the sequence(float)
        - mxa: maximum ambiguous percentage for the sequence(float)
    output:
        - remaining: dataframe containing GISAID metadata of all non-clinical, non-cell-based, and non-egg-based isolates that fullfilled the criteria for all segmenets
        - remaining_records: dictionary with thenon-clinical, non-cell-based, and non-egg-based sequence records per segments
    """
    #get non-clinical, non-cell-based, and non-egg-based samples from the metadata 
    remaining = metadata[~metadata["Passage_History"].isin(cpl)].reset_index(drop=True) #pandas doesn't seem to do and statement? odd???
    remaining = remaining[~remaining["Passage_History"].isin(cbpl)].reset_index(drop=True) 
    remaining = remaining[~remaining["Passage_History"].isin(epl)].reset_index(drop=True) 

    #also checking if collection date consist of at least year and month
    to_drop = set()
    for i, v in enumerate(remaining["Collection_Date"]):
        try:
            datetime.strptime(v, '%Y-%m-%d')
        except:
            to_drop.add(i)
    remaining = remaining.drop(to_drop).reset_index()

    #and location at least consists of country and continent
    to_drop = set()
    for i, v in enumerate(remaining["Location"]):
        if len(v.split(" / ")) < 2:
            to_drop.add(i)
    remaining = remaining.drop(to_drop)#.reset_index()

    #rename MP segment to M, if this hasn't been done
    if "MP Segment_Id" in list(remaining.columns):
        remaining = remaining.rename(columns={"MP Segment_Id": "M Segment_Id"})

    #there are no cell_based records found return none 
    if len(remaining) == 0: 
       return None 
    
    remaining_records = {}
    to_drop = set() #isolates that do not match the requirements 

    #filter cell_based sequences 
    pool = mp.Pool(processes=ncpu) 
    results = [pool.apply_async(filter_fasta, args=(remaining, fasta, rh, adth, mxa, mnlp, "remaining")) for fasta in sequence_files]
    output = [p.get() for p in results]
    pool.close()
    pool.join()

    for output_res in output:
        segment, rem, drop  = output_res
        remaining_records[segment] = rem
        to_drop.update(drop)

    #filter cell_based with records that need to be removed
    remaining = remaining[~remaining["Isolate_Id"].isin(to_drop)]
    
    #filter cell_based record from to drop
    pool = mp.Pool(processes=ncpu)
    results = [pool.apply_async(remove_failed, args=(segment, list(remaining["Isolate_Id"]), records)) for segment, records in remaining_records.items()]
    output = [p.get() for p in results]
    pool.close()
    pool.join()

    for output_res in output:
        segment, filtered = output_res
        remaining_records[segment] = filtered
        
    return remaining, remaining_records

def get_time_interval(timeframe, sm="may"):
    """
    get the full time interval per year per month
    input:
        - timeframe: timeframe to select from (str) should look like ("2015-2019")
        - sm: starting month from which year bin should start with (str) (default may > flu season)
    output:
        - full_time_intervals: list in list of time interval per moth per year from starting month
    """
    #get time frame and intervals   
    starty, endy  = timeframe.split("-")
        
    full_time_intervals = [] #list in list
    nexti = int(starty)
    while nexti < int(endy):
        full_int = []
        month = m2n[sm.lower()]
        while month <= 12: #start the interval 
            if len(str(month)) ==1:
                full_int.append(f"{str(nexti)}-0{str(month)}")
            else:
                full_int.append(f"{str(nexti)}-{str(month)}")
            month +=1 
        month = 1 #reset month for year
        nexti += 1 #and go to the next year
        while month < m2n[sm]: #complete the interval
            if len(str(month)) ==1:
                full_int.append(f"{str(nexti)}-0{str(month)}")
            else:
                full_int.append(f"{str(nexti)}-{str(month)}")
            month +=1

        full_time_intervals.append(full_int)
    return full_time_intervals

def select_per_year_per_month(metadata, time_intervals, spmpc=10, seed=None):
    """
    select GISAID sequences per year for given time frame start for the specified month (sm)
    and select at max spypc (sequence per year per country) 
    input: 
        - metadata: pandas data frame of the metadata file of the GISAID sequences to select from
        - time_intervals: list of time interval per year per month of interest (list) (gen with get_time_interval)
        - spypc: max sequences per country per year to select to avoid bias (int)
        - seed: seed for random sampling 
    """
    #set seed 
    if seed is not None:
        random.seed(seed)

    #subset metadata > only interest in Isolate_Id, Collection_Date, and Location
    subset = metadata[["Isolate_Id", "Collection_Date", "Location"]]
    to_drop = [] #used to drop 'incomplete' isolates, co is specified

    #split date into year month and day
    subset[["year", "month", "day"]] = subset["Collection_Date"].str.split("-", expand=True)
    subset["year-month"] = subset["year"] + "-" + subset["month"]
   
    #split Location in country
    subset["Country"] = "" #adding country column for later purposes
    for i, l in enumerate(subset["Location"]):
        if " / " not in l: #only region no country so can't use those
            to_drop.append(subset.iloc[i,subset.columns.get_loc("Isolate_Id")])
        elif " / " in l:
            #subset.iloc[i,subset.columns.get_loc("Location")] = l.split(" / ")[0]
            subset.iloc[i,subset.columns.get_loc("Country")] = l.split(" / ")[1]
            
    subset = subset[~subset["Isolate_Id"].isin(to_drop)]

    to_select_py = {}
    for interval in time_intervals:
        group = f"{interval[0]}_{interval[-1]}"
        to_select_py[group] = []

        for i in interval:
            intdata = subset[subset["year-month"]==i]

            for country in intdata["Country"].unique():
                isolates = list(intdata[intdata["Country"]==country]["Isolate_Id"])
                if len(isolates) > spmpc:
                    
                    to_select_py[group].extend(list(random.sample(isolates,spmpc)))
                else:
                    to_select_py[group].extend(isolates)

    return to_select_py

def get_year_month_location_count(selected_isolates, year_month_location_tracker={}):
    """
    get the year month location combination count for the isolatd in the inputted selected_isolates df
    input:
        - selected_isolates: pandas dataframe with isolates, collection_dates and location
        - year_month_location_tracker: year month location count dict to which observations needs to be added (empty dict by default)
    output:
        - yml:  year month location count dict update with the selected isolates
    """
    if len(selected_isolates) > 0:

        selected_isolates[["year", "month", "day"]] = selected_isolates["Collection_Date"].str.split("-", expand=True)
        selected_isolates["year-month"] = selected_isolates["year"] + "-" + selected_isolates["month"]

        for i, row in selected_isolates.iterrows():
            location = row["Location"].split(" / ")[1]
            yml = row["year-month"] +"-" + location
            try:
                year_month_location_tracker[yml] +=1
            except:
                year_month_location_tracker[yml] = 1
    return year_month_location_tracker

def extend_to_select_cell(clinical, cell_based, time_interval, selected=[], spmpc=10, seed=None):
    """
    select GISAID sequences per year for given time frame start for the specified month (sm)
    and select at max spypc (sequence per year per country) 
    input: 
        - clinical: pandas data frame of the clinical metadata file of the GISAID sequences to select from
        - cell_based: pandas data frame of the cell based metadata file of the GISAID sequences to select from
        - time_interval: list of time interval per month for the season of interest (list)
        - selected: list of the already selected sequences per season
        - spypc: max sequences per country per year to select to avoid bias (int)
        - seed: seed for random sampling 
    """
    #set seed 
    if seed is not None:
        random.seed(seed)

    #subset metadata > only interest in Isolate_Id, Collection_Date, and Location
    clinical_subset = clinical[["Isolate_Id", "Collection_Date", "Location"]]
    #get the selected samples from clincical
    selected_isolates = clinical[clinical["Isolate_Id"].isin(selected)]

    year_month_location_tracker = get_year_month_location_count(selected_isolates)

    cell_subset = cell_based[["Isolate_Id", "Collection_Date", "Location"]]

    #split date into year month and day
    cell_subset[["year", "month", "day"]] = cell_subset["Collection_Date"].str.split("-", expand=True)
    cell_subset["year-month"] = cell_subset["year"] + "-" + cell_subset["month"]
   
    to_drop = [] #used to drop 'incomplete' isolates, co is specified

    #split Location in country
    cell_subset["Country"] = "" #adding country column for later purposes
    for i, l in enumerate(cell_subset["Location"]):
        if " / " not in l: #only region no country so can't use those
            to_drop.append(cell_subset.iloc[i,cell_subset.columns.get_loc("Isolate_Id")])
        elif " / " in l:
            #subset.iloc[i,subset.columns.get_loc("Location")] = l.split(" / ")[0]
            cell_subset.iloc[i,cell_subset.columns.get_loc("Country")] = l.split(" / ")[1]
            
    cell_subset = cell_subset[~cell_subset["Isolate_Id"].isin(to_drop)]

    to_select = []
    for i in time_interval:
        intdata = cell_subset[cell_subset["year-month"]==i]

        for country in intdata["Country"].unique():
            try:
                already_selected = year_month_location_tracker[i +"-"+country]
            except:
                already_selected = 0

            isolates = list(intdata[intdata["Country"]==country]["Isolate_Id"])
            if already_selected < spmpc:
                remaining = spmpc - already_selected
                if len(isolates) > (remaining):
                    
                    to_select.extend(list(random.sample(isolates,remaining)))
                else:
                    to_select.extend(isolates)

    return to_select

def extend_to_select_egg(clinical, egg_based, time_interval, cell_based=None, selected=[], spmpc=10, seed=None):
    """
    select GISAID sequences per year for given time frame start for the specified month (sm)
    and select at max spypc (sequence per year per country) 
    input: 
        - clinical: pandas data frame of the clinical metadata file of the GISAID sequences to select from
        - egg_based: pandas data frame of the egg based metadata file of the GISAID sequences to select from
        - time_interval: list of time interval per month for the season of interest (list)
        - cell_based: pandas data frame of the cell based metadata file of the GISAID sequences to select from if None, not used
        - selected: list of the already selected sequences per season
        - spypc: max sequences per country per year to select to avoid bias (int)
        - seed: seed for random sampling 
    """
    #set seed 
    if seed is not None:
        random.seed(seed)

    #subset metadata > only interest in Isolate_Id, Collection_Date, and Location
    clinical_subset = clinical[["Isolate_Id", "Collection_Date", "Location"]]
    #get the selected samples from clincical
    selected_isolates = clinical[clinical["Isolate_Id"].isin(selected)]

    year_month_location_tracker = get_year_month_location_count(selected_isolates)

    #repeat this process with the cell based sequences if cell based was specified
    if cell_based is not None:
        cell_subset = cell_based[["Isolate_Id", "Collection_Date", "Location"]]
        selected_isolates = cell_subset[cell_subset["Isolate_Id"].isin(selected)]

        year_month_location_tracker = get_year_month_location_count(selected_isolates, year_month_location_tracker)

    egg_subset = egg_based[["Isolate_Id", "Collection_Date", "Location"]]

    #split date into year month and day
    egg_subset[["year", "month", "day"]] = egg_subset["Collection_Date"].str.split("-", expand=True)
    egg_subset["year-month"] = egg_subset["year"] + "-" + egg_subset["month"]
   
    to_drop = [] #used to drop 'incomplete' isolates, co is specified

    #split Location in country
    egg_subset["Country"] = "" #adding country column for later purposes
    for i, l in enumerate(egg_subset["Location"]):
        if " / " not in l: #only region no country so can't use those
            to_drop.append(egg_subset.iloc[i,egg_subset.columns.get_loc("Isolate_Id")])
        elif " / " in l:
            #subset.iloc[i,subset.columns.get_loc("Location")] = l.split(" / ")[0]
            egg_subset.iloc[i,egg_subset.columns.get_loc("Country")] = l.split(" / ")[1]
            
    egg_subset = egg_subset[~egg_subset["Isolate_Id"].isin(to_drop)]

    to_select = []
    for i in time_interval:
        intdata = egg_subset[egg_subset["year-month"]==i]

        for country in intdata["Country"].unique():
            try:
                already_selected = year_month_location_tracker[i +"-"+country]
            except:
                already_selected = 0

            isolates = list(intdata[intdata["Country"]==country]["Isolate_Id"])
            if already_selected < spmpc:
                remaining = spmpc - already_selected
                if len(isolates) > (remaining):
                    
                    to_select.extend(list(random.sample(isolates,remaining)))
                else:
                    to_select.extend(isolates)

    return to_select

def extend_to_select_remaining(clinical, cell_based, egg_based, remaining, time_interval, selected=[], spmpc=10, seed=None):
    """
    select GISAID sequences per year for given time frame start for the specified month (sm)
    and select at max spypc (sequence per year per country) 
    input: 
        - clinical: pandas data frame of the clinical metadata file of the GISAID sequences to select from
        - egg_based: pandas data frame of the egg based metadata file of the GISAID sequences to select from
        - time_interval: list of time interval per month for the season of interest (list)
        - cell_based: pandas data frame of the cell based metadata file of the GISAID sequences to select from if None, not used
        - selected: list of the already selected sequences per season
        - spypc: max sequences per country per year to select to avoid bias (int)
        - seed: seed for random sampling 
    """
    #set seed 
    if seed is not None:
        random.seed(seed)

    #subset metadata > only interest in Isolate_Id, Collection_Date, and Location
    clinical_subset = clinical[["Isolate_Id", "Collection_Date", "Location"]]
    #get the selected samples from clincical
    selected_isolates = clinical_subset[clinical_subset["Isolate_Id"].isin(selected)]

    year_month_location_tracker = get_year_month_location_count(selected_isolates)

    #repeat this process with the cell based sequences if cell based was specified
    cell_subset = cell_based[["Isolate_Id", "Collection_Date", "Location"]]
    selected_isolates = cell_subset[cell_subset["Isolate_Id"].isin(selected)]
    year_month_location_tracker = get_year_month_location_count(selected_isolates, year_month_location_tracker)
    #and for egg-based sequences
    egg_subset = egg_based[["Isolate_Id", "Collection_Date", "Location"]]
    selected_isolates = egg_subset[egg_subset["Isolate_Id"].isin(selected)]
    year_month_location_tracker = get_year_month_location_count(selected_isolates, year_month_location_tracker)

    #split date into year month and day
    remaining_subset = remaining[["Isolate_Id", "Collection_Date", "Location"]]
    remaining_subset[["year", "month", "day"]] = remaining_subset["Collection_Date"].str.split("-", expand=True)
    remaining_subset["year-month"] = remaining_subset["year"] + "-" + remaining_subset["month"]
   
    to_drop = [] #used to drop 'incomplete' isolates, co is specified

    #split Location in country
    remaining_subset["Country"] = "" #adding country column for later purposes
    for i, l in enumerate(remaining_subset["Location"]):
        if " / " not in l: #only region no country so can't use those
            to_drop.append(remaining_subset.iloc[i,remaining_subset.columns.get_loc("Isolate_Id")])
        elif " / " in l:
            #subset.iloc[i,subset.columns.get_loc("Location")] = l.split(" / ")[0]
            remaining_subset.iloc[i,remaining_subset.columns.get_loc("Country")] = l.split(" / ")[1]
            
    remaining_subset = remaining_subset[~remaining_subset["Isolate_Id"].isin(to_drop)]

    to_select = []
    for i in time_interval:
        intdata = remaining_subset[remaining_subset["year-month"]==i]

        for country in intdata["Country"].unique():
            try:
                already_selected = year_month_location_tracker[i +"-"+country]
            except:
                already_selected = 0

            isolates = list(intdata[intdata["Country"]==country]["Isolate_Id"])
            if already_selected < spmpc:
                rem = spmpc - already_selected
                if len(isolates) > (rem):
                    
                    to_select.extend(list(random.sample(isolates,rem)))
                else:
                    to_select.extend(isolates)

    return to_select

def select_per_year(metadata, time_intervals, ):
    """
    select GISAID sequences per year for given time frame start for the specified month (sm)
    and select at max spypc (sequence per year per country) 
    input: 
        - metadata: pandas data frame of the metadata file of the GISAID sequences to select from
        - time_intervals: list of time interval per year per month of interest (list) (gen with get_time_interval)
    """
    #this function can be a lot shorter but i wrote select per year per month first and changing two lines was faster
    #than writing an entire new function

    #subset metadata > only interest in Isolate_Id, Collection_Date, and Location
    subset = metadata[["Isolate_Id", "Collection_Date", "Location"]]
    to_drop = [] #used to drop 'incomplete' isolates, co is specified

    #split date into year month and day
    subset[["year", "month", "day"]] = subset["Collection_Date"].str.split("-", expand=True)
    subset["year-month"] = subset["year"] + "-" + subset["month"]
   
    #split Location in country
    subset["Country"] = "" #adding country column for later purposes
    for i, l in enumerate(subset["Location"]):
        if " / " not in l: #only region no country so can't use those
            to_drop.append(subset.iloc[i,subset.columns.get_loc("Isolate_Id")])
        elif " / " in l:
            #subset.iloc[i,subset.columns.get_loc("Location")] = l.split(" / ")[0]
            subset.iloc[i,subset.columns.get_loc("Country")] = l.split(" / ")[1]
            
    subset = subset[~subset["Isolate_Id"].isin(to_drop)]

    to_select_py = {}
    for interval in time_intervals:
        group = f"{interval[0]}_{interval[-1]}"
        to_select_py[group] = []

        for i in interval:
            intdata = subset[subset["year-month"]==i]

            for country in intdata["Country"].unique():
                isolates = list(intdata[intdata["Country"]==country]["Isolate_Id"])
                to_select_py[group].extend(isolates)

    return to_select_py

def get_dates(segment, metadata, outfile):
    """
    generates dates.tsv file for each segment contain the sequence header and date as a float
    input:
        - segment: segment of interest (str)
        - metadata: metadata of pandas dataframe from which information will be extracted (pd dataframe)
        - outfile: file to which dates needs to be written to (str) (.csv format)
    '!' assuming sequence header looks like Isolate_ID|Segment_ID|Isolate_Name|segment  
    """
    #add header to metadata
    for i, row in metadata.iterrows():
        header = "|".join([row["Isolate_Id"],row[f"{segment} Segment_Id"].split("|")[0], row["Isolate_Name"],segment]).replace(" ","").replace("-","")
        header = f"{header}P" if segment=="M" else header #GISAID refers to M as MP so MP is still in sequence HEADER
        metadata.at[i,'header'] = header

    #get Header and Collection_Date columns from metadata
    segrec = metadata[["header", "Collection_Date"]]

    #translate to collection date to float
    for i,v in enumerate(segrec["Collection_Date"]):
        v = "-".join(v.split("-")[::-1])

        year = int(v.split("-")[-1])
        yd = datetime.strptime(v,'%d-%m-%Y').timetuple().tm_yday
    
        if ((year % 400 == 0) and (year % 100 == 0)) or ((year % 4 ==0) and (year % 100 != 0)): #leap year
            dfloat = yd/366
        else:
            dfloat = yd/365
        segrec.iloc[i, segrec.columns.get_loc("Collection_Date")] = year + dfloat

    #rename columns 
    segrec = segrec.rename(columns={"header": "accession", "Collection_Date":"date"})

    #write output file
    segrec.to_csv(outfile,index=False)

def get_treetime_outliers(files):
    """
    get isolates that were reported as outliers from treetime screen output
    input:
        - files: list of treetime screen output files for segments to be analysed (list)
    output:
        - bad_isolates: set of the bad_isolate from all files (set)
    """
    bad_isolates = {}
    for f in files:
        segment = f.split("/")[-1].split("_")[0]
        bad_isolates[segment] = set()
        with open(f,"r") as fr:
            lines = fr.readlines()
        for l in lines:
            l = l.rstrip("\n")
            if len(l) >0 and "\t" in l:
                if l.split("\t")[1].startswith("EPI_ISL"):
                    isolate = l.split("\t")[1].split("|")[0] #get isolate ID
                    bad_isolates[segment].add(isolate)
    return bad_isolates 

def filter_mc_fasta(fasta,bad_isolates, outfile):
    """
    remove bad_isolate from FASTA file and write a new output file
    input:
        - fasta: fasta file that needs to be filtered (str)
        - bad_isolates: set/list of isolates that need to be removed (set)
        - outfile: name of output fasta file (str)
    """
    #segment 
    segment = fasta.split("/")[-1].split("_")[0]
    print(f"removing failed molecular clock isolates for {segment}")
    
    records = [] #filter out bad_isolates
    for record in SeqIO.parse(fasta, "fasta"):
        if record.id.split("|")[0] not in bad_isolates:
            records.append(record)

    #write outfile
    with open(outfile,"w") as fw:
        SeqIO.write(records,fw,"fasta")

def filter_molecular_clock_sequences(sequence_files, bad_isolates, outdir, ncpu=1):
    """
    remove all isolates that failed molecular clock estimation and write new files
    input:
        - sequence_files: sequence files from which bad isolates need to be removed (list)
        - bad_isolates: set/list of isolates that need to be removed (set)
        - output: directory to which new files need to be written (str)
        - ncpu: number of cpu to be used for pooling (int)
    """
    sequence_outfiles = {}
    for fl in sequence_files:
        outfile = os.path.join(outdir, f"{'_'.join(fl.split('/')[-1].split('.')[0].split('_')[:-1])}.fasta")
        sequence_outfiles[fl] = outfile
    
    pool = mp.Pool(processes=ncpu)
    results = [pool.apply_async(filter_mc_fasta, args=(fasta, bad_isolates, sequence_outfiles[fasta])) for fasta in sequence_files]
    #output = [p.get() for p in results]
    pool.close()
    pool.join()

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