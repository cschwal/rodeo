#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 17 19:48:13 2017

@author: bryce
"""

from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

class Intergenic_seq(object):
    def __init__(self, start_idx, end_idx, nn_seq):
        self.start_idx = start_idx
        self.end_idx = self.end_idx
        
    
    def translate(self):
        pass
    
class Orf(object):
    def __init__(self):
        self.species = ""
        self.genus = ""
        self.sequence = ""
        self.start = -1
        self.length = -1
        self.end = -1
        self.accession = ""
        self.direction = "None"
        self.id_gi = -1
        
class Efetch_result(object):
    def __init__(self, query, query_type="acc"):
        self.id = query
        self.idtype = query_type
        self.orfs = []
        self.intergenic_seqs = []
        self.cluster_accession = ""
        self.cluster_sequence = ""
        self.cluster_length = ""
        self.genus = ""
        self.species = ""
        
    
    def trim_to_n_neighbors(self, N):
        #make sure theyre sorted..
        self.orfs.sort(key=lambda orf: orf.start)
        for i in range(len(self.orfs)):
            orf = self.orfs[i]
            if orf.accession == self.id or orf.id_gi == self.id:
                #we found the query node, now get rid of all other orfs outside the 
                #neighborhood. 
                print(i)
                if i == 0 or i == len(self.orfs):
                    self.orfs = self.orfs[max(0,i-N-1):min(i+N+1,len(self.orfs))]
                else:
                    self.orfs = self.orfs[max(0,i-N):min(i+N+1,len(self.orfs))]
                break
        return
                
    def trim_to_n_nucleotides(self, N):
        #make sure theyre sorted..
        self.orfs.sort(key=lambda orf: orf.start)
        start_index = -1 #unset till we find it
        query_pos_left = 0
        query_pos_right = 0
        end_index = -1 #unset till we find it
        for i in range(len(self.orfs)):
            orf = self.orfs[i]
            if orf.accession == self.id or orf.id_gi == self.id:
                query_pos_left = min(orf.start,orf.end)
                query_pos_right = max(orf.start,orf.end)
                
        for i in range(len(self.orfs)):
            orf = self.orfs[i]
            if orf.start > query_pos_left - N and start_index == -1:
                start_index = i
            elif orf.start > query_pos_right + N:
                end_index = i
                break
            
        self.orfs = self.orfs[start_index:end_index]
        return
    
    def get_intergenic_seqs(self):
        for orf in self.orfs:
            seq_start = orf.start
            seq_end = orf.end
            nn_subseq = self.cluster_sequence[seq_start:seq_end+1]
            self.intergenic_seqs.append(Intergenic_seq(nn_subseq))
    
    def print_info(self):
        print("="*50)
        counter = 0
        for orf in self.orfs:
            print(counter)
            counter+=1
            if orf.accession == self.id or orf.id_gi == self.id:
                print("QUERY Accession:  " + orf.accession)
            else:
                print("Accession:  " + orf.accession)
            print("Coords:  " + str(orf.start) + " to " +  str(orf.end))
        print("="*50)

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

def get_efetch_results(query, batchSize=1, retmax=10**9, idtype="acc"):
    """
    Get efetch results from query accession sequence
    Currently only handles accession sequence queries
    """
    record = Entrez.read(Entrez.esearch("protein",term=query))
    total_count = record["Count"]
    if int(total_count) < 1:
        print("ERROR: Esearch returns no results for query " + query)
        return
    IdList = record["IdList"]#TODO may be multiple, need an outer for loop
    link_records = Entrez.read(Entrez.elink(dbfrom="protein",db="nuccore",id=IdList))
    nuccore_ids=[]
    for record in link_records:
        nuccore_ids.append(record['LinkSetDb'][0]['Link'][0]['Id'])  #TODO Not sure when [0] should or shouldn't be needed...
    
    #post NCBI query
    search_handle     = Entrez.epost(db="nuccore", id=",".join(nuccore_ids))
    search_results    = Entrez.read(search_handle)
    webenv,query_key  = search_results["WebEnv"], search_results["QueryKey"] 
    
    #fetch all results in batch of batchSize entries at once
    for start in range( 0,len(nuccore_ids),batchSize ):
        handle = Entrez.efetch(db="nuccore", dbfrom="protein", rettype="native",
                               retmode="xml", retstart=start, retmax=batchSize, 
                               webenv=webenv, query_key=query_key)
#        temp_file = open(query + "_nuccore_fetch.txt", "w+")
#        temp_file.write(handle.read())
        record = Entrez.read(handle)
        #Sometimes its Seq-entry_set, and others its Seq-entry_seq
        #Seq-entry_set is for regular acc, and Seq-entry_seq seems to be for 
        #WP_* accession sequences (non-unique)?
        try:
            record['Bioseq-set_seq-set'][0]['Seq-entry_set']
            key_type = "set"
        except Exception as e:
            print(e)
            try:
                record['Bioseq-set_seq-set'][0]['Seq-entry_seq']
                key_type = "seq"
            except Exception as e:
                print(e)
                print("ERROR:\t Not sure how to handle this accession...")
        
        if key_type == "set":
            #Get information from record.
            return_query = Efetch_result(query)
            #Backslashes look cleaner than one long string
#            for orf in record['Bioseq-set_seq-set']\
#                               [0]['Seq-entry_set']\
#                               ['Bioseq-set']['Bioseq-set_seq-set']:
            for i in range(len(record['Bioseq-set_seq-set']\
                               [0]['Seq-entry_set']\
                               ['Bioseq-set']['Bioseq-set_seq-set'])):
                print(i)
                #if not an aa sequence, then this is the orf entry describing
                #the entire cluster. so, don't create orf, just get info for query
                
                if 'aa' not in record['Bioseq-set_seq-set'][0]['Seq-entry_set']\
                                    ['Bioseq-set']['Bioseq-set_seq-set'][i]\
                                    ['Seq-entry_seq']['Bioseq']['Bioseq_inst']\
                                    ['Seq-inst']['Seq-inst_mol'].attributes['value']:
                                        
                   #Get cluster accession just in case i need it later                     
                    try:
                        cluster_accession = record['Bioseq-set_seq-set']\
                                                   [0]['Seq-entry_set']\
                                                   ['Bioseq-set']['Bioseq-set_seq-set']\
                                                   [i]['Seq-entry_seq']['Bioseq']\
                                                   ['Bioseq_id'][0]['Seq-id_ddbj']\
                                                   ['Textseq-id']['Textseq-id_accession']
                        return_query.cluster_accession = cluster_accession
                        version = record['Bioseq-set_seq-set']\
                                                   [0]['Seq-entry_set']\
                                                   ['Bioseq-set']['Bioseq-set_seq-set']\
                                                   [i]['Seq-entry_seq']['Bioseq']\
                                                   ['Bioseq_id'][0]['Seq-id_ddbj']\
                                                   ['Textseq-id']['Textseq-id_version']
                        return_query.cluster_accession += "."  + str(version)
                    
                    except Exception as e:
                        print(e)
                    
                    #Get DNA sequence
                    try:
                        sequence  = record['Bioseq-set_seq-set'][0]\
                                        ['Seq-entry_set']['Bioseq-set']\
                                        ['Bioseq-set_seq-set'][i]\
                                        ['Seq-entry_seq']['Bioseq']\
                                        ['Bioseq_inst']['Seq-inst']\
                                        ['Seq-inst_seq-data']['Seq-data']\
                                        ['Seq-data_iupacna']['IUPACna']
                        return_query.cluster_sequence = sequence
                        return_query.cluster_length = len(sequence)
                    
                    except Exception as e:
                        print(e)
                    
                    #Get genus/species
                    try:
                        tax_info = record['Bioseq-set_seq-set'][0]['Seq-entry_set']\
                                        ['Bioseq-set']['Bioseq-set_descr']\
                                        ['Seq-descr'][0]['Seqdesc_source']\
                                        ['BioSource']['BioSource_org']['Org-ref']\
                                        ['Org-ref_orgname']['OrgName']['OrgName_name']\
                                        ['OrgName_name_binomial']['BinomialOrgName']
                        return_query.genus = tax_info['BinomialOrgName_genus']
                        return_query.species = tax_info['BinomialOrgName_species']
                    
                    except Exception as e:
                        print(e)
                        
                                        
                #Otherwise, it is an ORF,                   
                else:  
                    current_orf = Orf()
                    orf_info = record['Bioseq-set_seq-set']\
                               [0]['Seq-entry_set']\
                               ['Bioseq-set']['Bioseq-set_seq-set'][i]
                    #Get id_gi
                    try:
                        id_gi = orf_info['Seq-entry_seq']['Bioseq']\
                                        ['Bioseq_id'][1]['Seq-id_gi']
                        current_orf.id_gi = id_gi
                    except Exception as e:
                        print("ERROR:\t Exception when getting id_gi")
                        print(e)
                      
                    #Get accession_ids
                    #TODO find where the db is stored and just use that
                    for db in ["ddbj", "embl", "ncbi", "gb"]:
                        try:
                            success = False
                            accession_id = orf_info['Seq-entry_seq']['Bioseq']['Bioseq_id']\
                                                [0]['Seq-id_' + db]['Textseq-id']\
                                                ['Textseq-id_accession']
                            current_orf.accession = accession_id
                            
                            version = orf_info['Seq-entry_seq']['Bioseq']['Bioseq_id']\
                                            [0]['Seq-id_' + db]['Textseq-id']\
                                            ['Textseq-id_version']
                            current_orf.accession += "." + str(version)
                            success = True
                        except Exception as e:
                            print("ERROR:\t Exception when getting accession with db = " + db)
                            print(e)
                        if success:
                            break
                        
                    #Get sequence
                    try:
                        sequence = orf_info['Seq-entry_seq']['Bioseq']\
                        ['Bioseq_inst']['Seq-inst']['Seq-inst_seq-data']['Seq-data']\
                        ['Seq-data_iupacaa']['IUPACaa']
                        current_orf.sequence = sequence
                    except Exception as e:
                        print("ERROR:\t Exception when getting sequence")
                        print(e)
                    
                    #Get sequence coordinates
                    try:
                        for potential_orf in record['Bioseq-set_seq-set'][0]\
                                                ['Seq-entry_set']['Bioseq-set']\
                                                ['Bioseq-set_annot'][0]\
                                                ['Seq-annot_data']\
                                                ['Seq-annot_data_ftable']:
                            gid = potential_orf['Seq-feat_product']['Seq-loc']\
                                                ['Seq-loc_whole']['Seq-id']\
                                                ['Seq-id_gi']
                            if gid == current_orf.id_gi:
                                seq_location_info = potential_orf['Seq-feat_location']\
                                                                    ['Seq-loc']\
                                                                    ['Seq-loc_int']\
                                                                    ['Seq-interval']
                                                        
                                seq_end = seq_location_info['Seq-interval_to']
                                seq_start = seq_location_info['Seq-interval_from']
                                current_orf.start = int(seq_start)
                                current_orf.end = int(seq_end) 
                                direction = "fwd"
                                if seq_start > seq_end:
                                    direction = "rev"
                                current_orf.direction = direction
                    except Exception as e:
                        print("ERROR:\t Exception when getting coordinates")
                        print(e)
                        
                    return_query.orfs.append(current_orf)
                        
                
    return return_query


    

            
    
