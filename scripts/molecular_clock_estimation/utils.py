#!/usr/bin/env python3
import os, sys, random, copy, dendropy
import pandas as pd
import numpy as np
from Bio import SeqIO, Seq, SeqRecord
from datetime import datetime, date
import multiprocessing as mp 

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