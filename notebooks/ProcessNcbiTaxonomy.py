#!/usr/bin/env python
# coding: utf-8
__all__ = ["ProcessNcbiTaxonomy"]

import extract_data as ex
import os
import pandas as pd
import tarfile


class ProcessNcbiTaxonomy(ex.ExtractData):
    """
    Extract node and edge data from NCBI Taxonomy files.
    """
    def __init__(self, filename, fileformat, compression):
        self.filename = filename
        self.fileformat = fileformat
        self.compression = compression
        self.node_name = "Organism"
        self.edge = "Organism-IS_A-Organims"
        
        # metadata: property, type, description, example
        node_properties = [
            ["id", "string", "NCBI taxonomy id", "taxonomy:9606"],
            ["name", "string", "organism name", "Homo sapiens"],
            ["synonyms", "string", "alternative organism names", ""],
            ["division", "string", "division in taxonomy tree", ""],
            ["rank", "string", "rank in taxonomy tree", ""],
        ]
        edge_properties = [
            ["from", "string", "NCBI taxonomy id", "taxonomy: 9605"],
            ["to", "string", "NCBI taxonomy id", "taxonomy: 9606"],
            
        self.node_metadata = {self.node_name, node_properties}
        self.edge_metadata = {self.edge_name, edge_properties}
            
    def get_node_metadata(self):
        return self.node_metadata
            
    def get_edge_metadata(self):
        return self.edge_metadata

    def extract_data(self):
        with tarfile.open(self.filename) as tf:
            names_obj = tf.extractfile('names.dmp')
            nodes_obj = tf.extractfile('nodes.dmp')
            divisions_obj = tf.extractfile('division.dmp')
            
            # method to strip tabs present in the dmp files
            trim = lambda x: x.strip()
  
            # import NCBI Taxonomy Names
            name_columns = ['id', 'name', 'nameCategory']
            names = pd.read_csv(names_obj, sep='|', header=None,
                                usecols=[0,1,3], names=name_columns,
                                converters={0: trim, 1: trim, 3: trim})
        
            name_categories = ['scientific name', 'synonym', 'blast name', 'genbank common name', 'equivalent name']
            names = names[names['nameCategory'].isin(name_categories)]

            # extract scientific names and add as a new column to names dataframe
            sci_name = names.query("nameCategory == 'scientific name'").copy()
            sci_name = sci_name.rename(columns={'name': 'scientificName'})
            sci_name = sci_name[['id', 'scientificName']]

            names = names.merge(sci_name, on='id', how='left')    
            names = names.groupby(['id', 'scientificName'])['name'].apply(lambda x: '|'.join(x)).reset_index(name='synonyms')
            names.rename(columns={'scientificName': 'name'}, inplace=True)
     
            # import NCBI Taxonomy Nodes
            node_columns = ['id', 'parentId', 'rank', 'divisionId']
            nodes = pd.read_csv(nodes_obj, sep='|', header=None, 
                                usecols=[0,1,2,4], names=node_columns,
                                converters={0: trim, 1: trim, 2: trim, 4: trim})

            # import NCBI Taxonomy Divisions
            division_columns = ['divisionId', 'division']
            divisions = pd.read_csv(divisions_obj, sep='|', header=None,
                                    usecols=[0,2], names=division_columns, 
                                    converters={0: trim, 2: trim})
        
            # add division and name data to the nodes
            nodes = nodes.merge(divisions, on='divisionId')
            nodes = nodes.merge(names, on='id')
        
            # add prefix to identifiers (see identfiers.org)
            # Wimalaratne, S., Juty, N., Kunze, J. et al. Uniform resolution of compact identifiers for biomedical data.
            # Sci Data 5, 180029 (2018). https://doi.org/10.1038/sdata.2018.29
            # https://n2t.net/e/cdl_ebi_prefixes.yaml
            nodes['id'] = 'taxonomy:' + nodes['id']
            nodes['parentId'] = 'taxonomy:' + nodes['parentId']

            # save nodes and edges
            self.save_nodes(nodes[['id','name', 'synonyms', 'division', 'rank']], self.node_name)
  
            edges = nodes[['id', 'parentId']].copy()
            #edges.drop_duplicates(inplace=True)
            edges.rename(columns={'id': 'from', 'parentId': 'to'}, inplace=True)
            edges = edges[(edges['from'] != 'taxonomy:1')]
            self.save_edges(edges[['from', 'to']], self.edge_name)