#!/usr/bin/env python
# coding: utf-8
__all__ = ["ExtractUniProtData"]

import extract_data as ex
import os
import gzip
import time
import dask.dataframe as dd
import pandas as pd
import re
from lxml import etree
import tarfile


class ExtractUniProt(ex.ExtractData):
    """
    Extract node and relationship data from UniProt xml files.
    """

    def __init__(self, filename, fileformat, compression):
        self.filename = filename
        self.fileformat = fileformat
        self.compression = compression
        self.namespace = ".//{http://uniprot.org/uniprot}"
        self.df_protein = None
        self.df_host = None
        self.df_cleavage = None

    def extract_data(self):
        print(self.filename)
        print(self.namespace)
        print(self.get_postfix())

        with gzip.open(self.filename, "rb") as file:
            self.__parse_uniprot(file)
            
        self.save_nodes(self.df_protein, "Protein")

        if self.df_host.shape[0] > 0:
            self.save_relationships(self.df_host, "Organism-HOSTS-Organism")

        if self.df_cleavage.shape[0] > 0:
            self.save_relationships(self.df_cleavage, "Protein-CLEAVED_TO-Protein")

    def __parse_uniprot(self, file):
        proteins = []
        hosts = []
        chains = []

        for event, entry in etree.iterparse(
            file, tag="{http://uniprot.org/uniprot}entry", huge_tree=True
        ):
            accessions = self.__get_accessions(entry)
            entry_name = self.__get_entry_name(entry)
            names = self.__get_protein_names(entry)
            sequence = self.__get_sequence(entry)
            gene_id = self.__get_gene_id(entry)
            gene_name = self.__get_gene_names(entry)
            taxid = self.__get_taxid(entry)

            # collect protein data
            proteins.append(
                {
                "id": accessions[0],
                "name": names[0],
                "accessions": accessions,
                "entryName": entry_name,
                "synonymes": names,
                "geneName": gene_name,
                "gene_id": gene_id,
                "taxid": taxid,
                "sequence": sequence,
                "fullLength": True
                }
            )

            # collect host data
            hosts += self.get_host_taxids(entry, taxid)

            # collect chain (cleaved peptide) data
            chains += self.get_chains(entry, accessions[0], sequence)

            entry.clear()

        self.df_protein = pd.DataFrame(proteins + chains)
        self.df_host = pd.DataFrame(hosts)
        df_chains = pd.DataFrame(chains)
        self.df_cleavage = df_chains[['id', 'parentId', 'begin', 'end']].copy()
        self.df_cleavage.rename(columns={'id': 'from', 'parentId': 'to'})

    def __get_accessions(self, entry):
        return [
            "uniprot:" + accession.text
            for accession in entry.findall(self.namespace + "accession")
        ]

    def __get_entry_name(self, entry):
        return entry.find(self.namespace + "name").text

    def __get_sequence(self, entry):
        # TODO matches <sequence> in isoform records, too. Why???. Use existance of checksum attribute to find the correct sequence.?
        sequences = [
            seq.text
            for seq in entry.findall(self.namespace + "sequence")
            if seq.get("checksum")
        ]
        return sequences[0]

    def __get_protein_names(self, entry):
        return [
            name.text
            for name in entry.find(self.namespace + "protein").findall(
                self.namespace + "fullName"
            )
        ] + [
            name.text
            for name in entry.find(self.namespace + "protein").findall(
                self.namespace + "shortName"
            )
        ]

    def __get_gene_names(self, entry):
        try:
            return [
                name.text
                for name in entry.find(self.namespace + "gene").findall(
                    self.namespace + "name"
                )
            ]
        except Exception:
            return []

    def __get_gene_id(self, entry):
        # <dbReference type="GeneID" id="4156347"/>
        try:
            return [
                "ncbigene:" + dbref.get("id")
                for dbref in entry.findall(self.namespace + "dbReference")
                if dbref.get("type") == "GeneID"
            ][0]
        except Exception:
            return ""

    def __get_taxid(self, entry):
        #  <organism> <dbReference type="NCBI Taxonomy" id="9606"/> </organism>
        taxids = [
            "taxonomy:" + dbref.get("id")
            for dbref in entry.find(self.namespace + "organism").findall(
                self.namespace + "dbReference"
            )
            if dbref.get("type") == "NCBI Taxonomy"
        ]
        return taxids[0]

    def get_host_taxids(self, entry, taxid):
        try:
            rows = []
            host_taxids = [
                "taxononmy:" + dbref.get("id")
                for dbref in entry.find(self.namespace + "organismHost").findall(
                    self.namespace + "dbReference"
                )
                if dbref.get("type") == "NCBI Taxonomy"
            ]
            rows = [
                {"from": taxid, "to": host_taxid} for host_taxid in host_taxids
            ]
            return rows
        except Exception:
            return []

    def get_chains(self, entry, accession, sequence):
        # <feature type="chain" id="PRO_0000377950" description="Uncharacterized protein 037L">
        #    <location>
        #       <begin position="1"/>
        #       <end position="253"/>
        # </location>

        rows = []
        for feature in entry.findall(self.namespace + "feature"):
            if feature is not None and feature.get("type") == "chain":
                chain_id = 'uniprot.chain:' + feature.get("id")
                name = feature.get("description")
                location = feature.find(self.namespace + "location")
                
                begin = location.find(self.namespace + "begin")
                pos1 = ""
                if begin is not None:
                    pos1 = begin.get("position")
                    
                end = location.find(self.namespace + "end")
                pos2 = ""
                if end is not None:
                    pos2 = end.get("position")
                
                # cleave sequence if both cleavage sites are specified
                if pos1 is not None and pos1.isdigit() and pos2 is not None and pos2.isdigit():
                    seq = sequence[int(pos1):int(pos2)]
                else:
                    seq = ''
                    
                rows.append(
                    {
                        "id": chain_id,
                        "name": name,
                        "sequence": seq,
                        "parentId": accession,
                        "begin": pos1,
                        "end": pos2,
                    }
                )

        return rows

    
class ExtractNcbiTaxonomyPd(ex.ExtractData):
    """
    Extract node and relationship data from NCBI Taxonomy files.
    """
    def __init__(self, filename, fileformat, compression):
        self.filename = filename
        self.fileformat = fileformat
        self.compression = compression

    def extract_data(self):
        tar = tarfile.open(self.filename)
        path = self.get_path()
        tar.extract('names.dmp',path)
        tar.extract('nodes.dmp', path)
        tar.extract('division.dmp', path)
        tar.close()

        #### Import NCBI Taxonomy Names  
        
        # method to strip tabs present in the dmp files
        trim = lambda x: x.strip()
  
        # 145 sec.
        name_columns = ['id', 'name', 'nameCategory']

        names = pd.read_csv(os.path.join(path, 'names.dmp'), sep='|', header=None,
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
        nodes = pd.read_csv(os.path.join(path, 'nodes.dmp'), sep='|', header=None, 
                            usecols=[0,1,2,4], names=node_columns,
                            converters={0: trim, 1: trim, 2: trim, 4: trim})

        division_columns = ['divisionId', 'division']
        divisions = pd.read_csv(os.path.join(path,'division.dmp'), sep='|', header=None,
                                usecols=[0,2], names=division_columns, 
                                converters={0: trim, 2: trim})
        
        # add division and name data to the nodes
        nodes = nodes.merge(divisions, on='divisionId')
        nodes = nodes.merge(names, on='id')
        
        # add prefix to identifiers (see identfiers.org)
        nodes['id'] = 'taxonomy:' + nodes['id']
        nodes['parentId'] = 'taxonomy:' + nodes['parentId']

        # save nodes and relationships
        self.save_nodes(nodes[['id','name', 'synonyms', 'division', 'rank']], "Organism")

        relationships = nodes[['id', 'parentId']].copy()
        #relationships.drop_duplicates(inplace=True)
        relationships.rename(columns={'id': 'from', 'parentId': 'to'}, inplace=True)
        relationships = relationships[(relationships['from'] != 'taxonomy:1')]
        self.save_relationships(relationships[['from', 'to']], 'Organism-IS_A-Organism')

        
class ExtractNcbiTaxonomyDd(ex.ExtractData):
    """
    Extract node and relationship data from NCBI Taxonomy files.
    """
    def __init__(self, filename, fileformat, compression):
        self.filename = filename
        self.fileformat = fileformat
        self.compression = compression

    def extract_data(self):
        tar = tarfile.open(self.filename)
        path = self.get_path()
        tar.extractall(path)
        tar.close()

        #### Import NCBI Taxonomy Names  
        columns = ['identifier', 'name', 'nameCategory']
        name_categories = ['scientific name', 'synonym', 'blast name', 'genbank common name', 'equivalent name']
        
        trim = lambda x: x.strip()
  
        names = dd.read_csv(os.path.join(path, 'names.dmp'), sep='|', header=None,
                            usecols=[0,1,3], names=columns,
                            converters={0: trim, 1: trim, 3: trim})
        names = names[names['nameCategory'].isin(name_categories)]
        # print(names.describe())
        # print(names.head(25))
        
        #names['identifier'] = names['identifier'].astype(str)

        sci_name = names.query("nameCategory == 'scientific name'").copy()
        sci_name = sci_name.rename(columns={'name': 'scientificName'})
        sci_name = sci_name[['identifier', 'scientificName']]

        names = names.merge(sci_name, on='identifier', how='left')
        print("With sci names")
        print(names.describe())
        print(names.head(25))
  
        names = names.rename(columns={'name': 'synonyms'})
        names = names.groupby(['identifier', 'scientificName'])['synonyms'].apply(lambda x: '|'.join(x), meta=('synonyms', 'object')).reset_index()
        
        print("---- Names")
        print(names.describe())
        print(names.head(25))
     
        #### Import NCBI Taxonomy Nodes
        node_columns = ['identifier', 'parentId', 'rank', 'divisionId']
        nodes = dd.read_csv(os.path.join(path, 'nodes.dmp'), sep='|', header=None, 
                            usecols=[0,1,2,4], names=node_columns,
                            converters={0: trim, 1: trim, 2: trim, 4: trim})
 
        #nodes['divisionId'] = nodes['divisionId'].astype(str)
        print("------ Nodes")
        print(nodes.info())
        print(nodes.head())
    
#         division_columns = ['divisionId', 'division']
#         # use pandas here since this is a very small dataset
#         divisions = pd.read_csv(os.path.join(path,'division.dmp'), sep='|', header=None,
#                                 usecols=[0,2], names=division_columns, 
#                                 converters={0: trim, 2: trim})
#         #divisions['divisionId'] = divisions['divisionId'].astype(str)
        
#         print("------- Divisions")
#         print(divisions.describe())
#         print(divisions.head())
        
              
#         nodes = nodes.merge(divisions, on='divisionId')
#         print("Nodes-Division:", nodes.shape[0])
#         print(nodes.describe())
#         print(nodes.head(10))
        
        nodes1 = dd.merge(nodes, names, on='identifier')
        
        #nodes['identifier'] = 'taxonomy:' + nodes['identfier']
        #nodes['parentId'] = 'taxonomy:' + nodes['parentId']

        
        print("Nodes-Division-Names:", nodes1.shape[0])
        print(nodes1.describe())
        print(nodes1.head(25))
        
        return

        self.save_nodes(nodes1[['identifier','name', 'synonyms', 'division', 'rank']], "Organism")

        relationships = nodes1[['parentId', 'id']].copy()
        relationships = relationships.rename(columns={'id': 'from', 'parentId': 'to'})
        relationships.drop_duplicates(inplace=True)
        relationships = relationships[(relationships['from'] != 'taxonomy:1')]

        self.save_relationships(relationships, 'Organism-IS_A-Organism')


class ExtractNcbiGenesPd(ex.ExtractData):
    """
    Extract node and relationship data from NCBI Gene files.
    """
    def __init__(self, filename, fileformat, compression):
        self.filename = filename
        self.fileformat = fileformat
        self.compression = compression

    def extract_data(self):
        print(self.filename)
        print(self.get_postfix())
        
        column_names = ['#tax_id', 'GeneID', 'Symbol', 'Synonyms', 'description', 'type_of_gene', 'Symbol_from_nomenclature_authority']   

        genes = pd.read_csv(self.filename, usecols=column_names, dtype=str, sep='\t')
        genes.rename(columns={'#tax_id': 'taxonomyId', 'GeneID': 'id', 'Symbol': 'name', 'Synonyms': 'synonyms', 'type_of_gene': 'geneType', 'Symbol_from_nomenclature_authority': 'officialSymbol'}, inplace=True)

        genes.replace('-', '', inplace=True)
        genes['id'] = 'ncbigene:' + genes['id']
        genes['taxonomyId'] = 'taxonomy:' + genes['taxonomyId']

        # genes[['id', 'name', 'synonyms', 'description', 'officialSymbol', 'geneType', 'taxonomyId']].to_parquet(os.path.join(self.get_path(), 'Gene' + self.get_postfix() + '.parquet'), compression='brotli', index=False)

        self.save_nodes(genes, 'Gene')
   
         #genes.drop(columns=['name', 'synonyms', 'description', 'officialSymbol', 'geneType'], inplace=True)

        relationships = genes[['taxonomyId', 'id']].copy()
        relationships.drop_duplicates(inplace=True)
        relationships.rename(columns={'taxonomyId': 'from', 'id': 'to'}, inplace=True)

        self.save_relationships(relationships, 'Organism-ENCODES-Gene')
        

class ExtractNcbiGenes(ex.ExtractData):
    """
    Extract node and relationship data from NCBI Gene files.
    """
    def __init__(self, filename, fileformat, compression):
        self.filename = filename
        self.fileformat = fileformat
        self.compression = compression

    def extract_data(self):
        print(self.filename)
        print(self.get_postfix())
        

        filename, file_extension = os.path.splitext(self.filename)

        import shutil
        with gzip.open(self.filename, 'rb') as f_in:
            with open(filename, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        # import gzip
        # with gzip.open(self.filename, 'rb') as f:
        #     with open(ufilename, 'wb') as f_out:
        #         f_out.write(f.read())

        column_names = ['#tax_id', 'GeneID', 'Symbol', 'Synonyms', 'description', 'type_of_gene', 'Symbol_from_nomenclature_authority']   

        genes = dd.read_csv(filename, usecols=column_names, dtype=str, sep='\t', blocksize='0.1 GB')
        print('Gene partitions:', genes.npartitions)
        
        genes = genes.rename(columns={'#tax_id': 'taxonomyId', 'GeneID': 'id', 'Symbol': 'name', 'Synonyms': 'synonyms', 'type_of_gene': 'geneType', 'Symbol_from_nomenclature_authority': 'officialSymbol'})
        genes = genes.replace('-', '')
        genes['id'] = 'ncbigene:' + genes['id']
        genes['taxonomyId'] = 'taxonomy:' + genes['taxonomyId']

        #genes.drop(columns=['name', 'synonyms', 'description', 'officialSymbol', 'geneType'], inplace=True)
        #genes = genes.drop(columns=['name', 'synonyms', 'description', 'officialSymbol', 'geneType'])
        self.save_nodes(genes, 'Gene')

        #relationships = genes[['taxonomyId', 'id']].copy()
        relationships = genes[['taxonomyId', 'id']]
        relationships = relationships.rename(columns={'taxonomyId': 'from', 'id': 'to'})

        self.save_relationships(relationships, 'Organism-ENCODES-Gene')