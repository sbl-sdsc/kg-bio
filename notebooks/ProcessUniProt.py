#!/usr/bin/env python
# coding: utf-8
__all__ = ["ProcessUniProtXL"]

import extract_data as ex
import os
import gzip
import pandas as pd
import dask.dataframe as dd
import dask.bag as db
from lxml import etree
import csv


class ProcessUniProtXL(ex.ExtractData):
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
        self.node_name = "Gene"
        self.edge = "Organism-ENCODES-Gene"
        
        # metadata: property, type, description, example
        node_properties = [
            ["id", "string", "protein identifier", "ncbigene:59272"],
            ["name", "string", "protein name", "ACE2"],
            ["synonyms", "string", "alternative protein names", "ACEH"],
            ["ids", "string", "alternative identifiers", "angiotensin converting enzyme 2"],
            ["geneName", "string", "gene name", ""],
            ["geneId", "string", "gene identifier", "ncbigene:59272"],
            ["taxonomyId", "string", "NCBI taxonomy id", "taxonomy:9606"],
            ["sequence", "string", "protein sequence", "taxonomy:9606"],
            ["fullLength", "boolean", "True if full length protein", "True"],
        ]
        edge_properties = [
            ["from", "string", "NCBI taxonomy id", "taxonomy:9606"],
            ["to", "string", "gene identifier", "ncbigene:59272"],
            
        self.node_metadata = {self.node_name, node_properties}
        self.edge_metadata = {self.edge_name, edge_properties}

    def extract_data(self):
        with gzip.open(self.filename, "rb") as file:
            self.__parse_uniprot(file)

    def __parse_uniprot(self, file):
        protein_columns = [
            "id",
            "name",
            "synonymes",
            "ids",
            "geneName",
            "geneId",
            "taxid",
            "sequence",
            "fullLength",
        ]
        host_columns = ["from", "to"]
        chain_columns = ["id", "name", "sequence", "parentId", "begin", "end"]

        proteins = []
        hosts = []
        chains = []

        protein_count = 0
        host_count = 0
        chain_count = 0

        path = self.get_path()
        protein_filename = os.path.join(path, "Protein.parquet")
        os.makedirs(protein_filename, exist_ok=True)
        host_filename = os.path.join(path, "Organisms-HOSTS-Organism.parquet")
        os.makedirs(host_filename, exist_ok=True)
        chain_filename = os.path.join(path, "Protein-CLEAVED_TO-Protein.parquet")
        os.makedirs(chain_filename, exist_ok=True)

        for event, entry in etree.iterparse(
            file, tag="{http://uniprot.org/uniprot}entry"
        ):

            accessions = self.__get_accessions(entry)
            entry_name = self.__get_entry_name(entry)
            names = self.__get_protein_names(entry)
            taxid = self.__get_taxid(entry)
            sequence = self.__get_sequence(entry)
            
            synonymes = [entry_name]
            if len(names > 1)
                synonymes = synonymes + names[1:]
            synonymes = '|'.join(synonymes)
            
            proteins.append(
                [
                    accessions[0],
                    names[0],
                    synonymes,
                    accessions,
                    self.__get_gene_names(entry),
                    self.__get_gene_id(entry),
                    taxid,
                    sequence,
                    True,
                ]
            )

            # save data in batches
            if len(proteins) == 100000:
                self.save_df(protein_filename, proteins, protein_columns, protein_count)
                proteins = []
                protein_count += 1

            hosts.extend(self.get_host_taxids(entry, taxid))

            chains.extend(self.get_chains(entry, accessions[0], sequence))
            if len(chains) > 100000:
                self.save_df(chain_filename, chains, chain_columns, chain_count)
                chains = []
                chain_count += 1

            entry.clear()

        # save any remaining records
        if len(proteins) > 0:
            self.save_df(protein_filename, proteins, protein_columns, protein_count)

        if len(hosts) > 0:
            self.save_df(host_filename, hosts, host_columns, host_count, True)

        if len(chains) > 0:
            self.save_df(chain_filename, chains, chain_columns, chain_count)

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
            rows = [[taxid, host_taxid] for host_taxid in host_taxids]
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
                chain_id = "uniprot.chain:" + feature.get("id")
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
                if (
                    pos1 is not None
                    and pos1.isdigit()
                    and pos2 is not None
                    and pos2.isdigit()
                ):
                    seq = sequence[int(pos1) : int(pos2)]
                else:
                    seq = ""

                rows.append([chain_id, name, seq, accession, pos1, pos2])

        return rows

    def save_df(self, filename, records, column_names, count, drop_duplicates=False):
        df = pd.DataFrame(records, columns=column_names)
        if drop_duplicates:
            df.drop_duplicates(inplace=True)
        df.to_parquet(os.path.join(filename, f"part.{count}.parquet"))
