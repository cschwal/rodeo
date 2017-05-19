#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 17 19:48:13 2017

@author: bryce
"""

from Bio import Entrez

class Orf(object):
    def __init(self):
        self.species = ""
        self.genus = ""
        self.sequence = ""
        self.start = -1
        self.length = -1
        self.end = -1
        self.accession = ""

def get_tag_contents(line, tag):
    """
    Get contents contained inside a tag. Only works when both theopen and 
    closing tag are on the same line
    """
    if "</" not in line:
        print("Closing tag not on same line for line:")
        print(line)
    line = line.strip().split("</" + tag + ">")
    line = line[0].split("<" + tag + ">")[1]
    return line


#TODO is this all neighboring genes or all genes or what... not sure yet
#different how Parth uses efetch, he uses elink -db nuccore to feed to 
#efetch, whereas I don't use elink at all
#TODO find out why Parth uses the 'nuccore' db for elink
#TODO handle idtype args other than "acc"
#TODO handle batch sizes > 1
#TODO find out what other things we need from the efetch results

def get_efetch_results(query, batchSize=1, retmax=10**9, db="protein", idtype="acc"):
    """
    Get efetch results from query accession sequence
    Currently only handles accession sequence queries
    """
    record = Entrez.read(Entrez.esearch("protein",term=query, idtype=idtype))
    total_count = record["Count"]
    if int(total_count) < 1:
        print("ERROR: Esearch returns no results for query " + query)
        return
    handle = Entrez.esearch( db=db,term=query,retmax=retmax )#TODO not sure about the db arg...
    giList = Entrez.read(handle)['IdList']
    
    #post NCBI query
    search_handle     = Entrez.epost(db=db, id=",".join(giList))
    search_results    = Entrez.read(search_handle)
    webenv,query_key  = search_results["WebEnv"], search_results["QueryKey"] 
    
    #fetch all results in batch of batchSize entries at once
    for start in range( 0,len(giList),batchSize ):
        handle = Entrez.efetch(db=db, rettype="native", retmode="xml", 
                               retstart=start, retmax=batchSize, webenv=webenv, 
                               query_key=query_key)
        all_orfs = []
        line_no = 0
        #I notice that sometimes there are multiple close tags
        #So i want to make sure we keep track of open seq entry and close seq-entry tags
        in_entry = False 
        
        for line in handle.readlines():
            line_no += 1
            if "<Seq-entry>" in line:
                current_orf = Orf()
                in_entry = True
            
            if "taxname" in line:
                tag_contents = get_tag_contents(line, "Org-ref_taxname")
                tag_contents = tag_contents.split(' ')
                current_orf.genus = tag_contents[0]
                current_orf.species = tag_contents[1]
                
            if "<Textseq-id_accession>" in line:
                current_orf.accession = get_tag_contents(line=line, 
                                                         tag="Textseq-id_accession")
                
            if "<Textseq-id_version>" in line:
                version = get_tag_contents(line=line, tag="Textseq-id_version")
                if version != '':
                    current_orf.accession = current_orf.accession + "." + version
           
            if "IUPACaa" in line:
                pp_seq = get_tag_contents(line=line, tag="IUPACaa")
                current_orf.sequence = pp_seq
                
                
            if "</Seq-entry>" in line and in_entry:
                #TODO this is a sloppy way of checking whether or not this seq-entry
                #was actually an ORF rather than a label for a BGC or something else
                in_entry = False
                if hasattr(current_orf, "sequence"): 
                    all_orfs.append(current_orf)
                    
    return all_orfs
            
    
