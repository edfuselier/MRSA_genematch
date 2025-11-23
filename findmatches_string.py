#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 13:48:59 2024

@author: efuselier

TThe script begins by parsing the STRING database files
(`93061.protein.*`) to extract locus tags, protein sequences, gene
names, and product information for all genes annotated in
Staphylococcus aureus subsp. aureus NCTC 8325. These data are saved
into two output files, `stringsequences.xlsx` and `stringinfo.xlsx`.

It then processes GenBank files downloaded from the NIH database
(e.g., `ATCC_43300.gb`, `USA100_N315.gb`, `USA300_FPR3757.gb`,
`VRSA1.gb`), extracting sequences, gene names, and other metadata.
For each gene, the script searches for the closest match within the
STRING database. NIH-provided gene names are used by default;
however, if a STRING entry achieves a Levenshtein similarity ratio
greater than 0.8, the corresponding STRING identifier is used as the
gene name. Locus tags, STRING matches, gene names from each source,
and additional annotations are compiled into a separate Excel file
for each strain.

Genes that do not have associated sequences in the NIH GenBank files
are tracked in separate Excel files, one per strain, with filenames
ending in `_skipped.xlsx`.                                                                            
"""

import genedata_processes as gdp
import pandas as pd
from multiprocessing import Pool
import time
from functools import partial

# extract sequences, gene names, and descriptions from string data
# creates "stringsequences.xlsx" and "stringinfo.xlsx"
gdp.stringdata2excel()

genbankfiles = ["ATCC_43300.gb","USA100_N315.gb","USA300_FPR3757.gb","VRSA1.gb"]

for gb_file in genbankfiles:

    gdp.genbank2excel(gb_file)
    #stringdata2excel(): 
    
    df = pd.read_excel(gb_file.removesuffix('.gb')+'.xlsx')
    ids = df.values[:,0];
    names = df.values[:,1];
    seq = df.values[:,2];
    
    stringseq = gdp.getstringseq()
    
    n = len(seq)
    stringmatch = ['']*n
    conf = [-1]*n
    
    t1 = time.time()
    with Pool() as executor:
    #    stringmatch, conf = executor.map(partial(gdp.findmatch,proteinsequences=stringseq),seq[0:5])
        X = executor.map(partial(gdp.findmatch,proteinsequences=stringseq),seq)
    
    stringmatch = [x[0] for x in X]
    conf = [x[1] for x in X]
    
    t2 = time.time()
    print('total runtime to compute matches:',t2-t1)
    
    seq2stringid, stringid2stringname = gdp.makestringdictionaries()
    
    stringids = ['']*len(seq)
    stringnames = ['']*len(seq)
    for i in range(n):
    #    print(txt.format(i,n)) 
        stringids[i] = seq2stringid[stringmatch[i]]
        stringnames[i] = stringid2stringname[stringids[i]]   
    #    stringids[i],conf[i] = getstringid(seq[i])
    #    stringnames[i] = stringid2name[stringids[i]]
    
    gb_file = gb_file.removesuffix('.gb')
    
    chosen_names = names.copy()
    chosen_names_nounknowns = names.copy()
    name_choice = ['NIH']*len(chosen_names)
    
    for i in range(len(chosen_names)):
        x = str(names[i])
        if conf[i]>.8:
            chosen_names[i]=stringnames[i]
            name_choice[i]='string'        
        
        if x=='nan':  
            if conf[i]>.8:
                chosen_names[i]=stringnames[i]
                name_choice[i]='string'
            else:
                chosen_names[i]='unknown'
                
        if chosen_names[i]=='unknown':
            chosen_names_nounknowns[i] = ids[i]
        else:
            chosen_names_nounknowns[i] = chosen_names[i]        
                
    # insert the skipped ones
    
    strain2strdf = pd.DataFrame(list(zip(ids,stringids,chosen_names,chosen_names_nounknowns,names,stringnames,conf,name_choice)),columns=['gene_id','string_id','chosen_name','chosen_name2','dbname','stringdbname','confidence','name_origin'])
    
    df = pd.read_excel(gb_file.removesuffix('.gb')+'_skipped.xlsx')
    skipped_ids = list(df.loc[:,'locus_tag'])
    skipped_names = list(df.loc[:,'name'])
    for i in range(len(skipped_names)):
        x = str(skipped_names[i])
        if x=='nan':  
            skipped_names[i]='unknown'
    
    df_new = pd.DataFrame({'gene_id':skipped_ids,'chosen_name':skipped_names,'name_origin':['NIH']*len(skipped_names),'confidence':['No sequence in NIH genbank']*len(skipped_names),'dbname':skipped_names})
    strain2strdf = pd.concat([strain2strdf, df_new], ignore_index=True)
    
    strain2strdf.to_excel(gb_file+'_stringmatches.xlsx',index=False)

