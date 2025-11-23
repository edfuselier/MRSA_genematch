#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 08:34:46 2024

@author: efuselier

This files contains helper functions used in findmatches_string.py and
findmatches_all.py

"""
import pandas as pd
from Bio import SeqIO
import numpy as np
import Levenshtein as lv

def matchratio(str1,str2):
    return lv.ratio(str1,str2)

def findmatch(str1,proteinsequences):
    mrat = np.array([matchratio(str1,i) for i in proteinsequences])
    return proteinsequences[np.argmax(mrat)], max(mrat)

def genbank2start(genbankfile):
    for gb_record in SeqIO.parse(open(genbankfile,"r"), "genbank") :
        print('genbank parsing successful')
    
    features = gb_record.features
    
    starts = []
    geneids = []
    for feature in features:
        if feature.type == 'CDS':
            geneids.append(feature.qualifiers['locus_tag'][0])
            starts.append(str(feature.location.start))
        
    return dict(zip(starts,geneids))    

def genbank2excel(genbankfile):
    
    '''    for gb_record in SeqIO.parse(open(genbankfile,"r"), "genbank") :
        print('genbank parsing successful')

    # list of features of various types
    # features[0].type is "source" - info about the organism
    # features[n>0].type is either "gene" (basic name, aliases or gene)
    # or "CDS" (more info about gene, including its 'sequence')
    features = gb_record.features
    '''
        
    gbdict = SeqIO.to_dict(SeqIO.parse(open(genbankfile,"r"),"genbank"))
    
    features = []
    for x in gbdict.keys():
        features = features + gbdict[x].features
            
    
    # get indices of all features of type 'CDS'
    # NOTE: The sequences (i.e. "translation") for some genes may be missing.
    # If so, these are tracked in the "_skipped.xlsx" spreadsheet.
    y = []
    yskipped = []
    for feature in features:
        if feature.type == 'CDS':
            # create empty default values
            geneval = ['']
            locus_tag = ['']
            transval = ['']
            product = ['']
            if 'gene' in feature.qualifiers:
                geneval = feature.qualifiers['gene']
            if 'translation' in feature.qualifiers:
                transval = feature.qualifiers['translation']
            if 'locus_tag' in feature.qualifiers:
                locus_tag = feature.qualifiers['locus_tag']
            if 'product' in feature.qualifiers:
                product = feature.qualifiers['product']    
                
    #        x = {'gene':geneval,'locus_tag':tag,'translation':transval}
            if len(transval[0])>0:
                y.append([locus_tag[0],geneval[0],transval[0],product[0]])
            else:
                yskipped.append([locus_tag[0],geneval[0],transval[0],product[0]])
                
    df = pd.DataFrame(y,columns=['locus_tag','name','sequence','product'])
    filename = genbankfile.removesuffix('.gb')
    df.to_excel(filename+'.xlsx',index=False)
    
    df = pd.DataFrame(yskipped,columns=['locus_tag','name','sequence','product'])
    filename = genbankfile.removesuffix('.gb')
    df.to_excel(filename+'_skipped.xlsx',index=False)
       
def stringdata2excel():   
    stringseq_file = "93061.protein.sequences.v12.0.fa"
    stringinfo_file = "93061.protein.info.v12.0.txt"

    seq2stringid = {}
    y = []
    stringproteinsequences = []
    for seq_record in SeqIO.parse(open(stringseq_file,"r"), "fasta") :
        record_id = seq_record.id
        sequence = seq_record.seq
        y.append([record_id,sequence])
        stringproteinsequences.append(sequence)
        seq2stringid[sequence] = record_id
        
    df = pd.DataFrame(y,columns=['string_id','sequence'])
    df.to_excel('stringsequences.xlsx',index=False)

    # column 0 is string protein id 
    # (NOTE string id seems same as species # + genbank locus_tag)
    # column 1 is "preferred name"
    dfinfo = pd.read_table(stringinfo_file)
    dfinfo.to_excel('stringinfo.xlsx',index=False)


def makestringdictionaries():
    stringseq_file = "stringsequences.xlsx"
    stringinfo_file = "stringinfo.xlsx"
    df = pd.read_excel(stringinfo_file)
    #dfinfo.to_excel('stringinfo.xlsx',index=False)
    stringids = df.values[:,0]
    name = df.values[:,1]
    stringid2stringname = dict(zip(stringids,name))
    df = pd.read_excel(stringseq_file)
    stringids = df.values[:,0]
    seq = df.values[:,1]
    seq2stringid = dict(zip(seq,stringids))
    return seq2stringid, stringid2stringname

def getstringseq():
    dfseq = pd.read_excel('stringsequences.xlsx')
    return dfseq.values[:,1]
#    dfinfo = pd.read_excel('stringinfo.xlsx')
    
    
    
    
    
    
    

    
    
    
    
    