#!/usr/bin/env python
# coding: utf-8
__all__ = ["ProcessNcbiGenes"]

import extract_data as ex
import os
import gzip
import time
import dask.dataframe as dd
import pandas as pd
import re
import shutil


class ProcessNcbiGenes(ex.ExtractData):
    """
    Extract node and edge data from NCBI Gene files.
    """

    def __init__(self, filename, fileformat, compression):
        self.filename = filename
        self.fileformat = fileformat
        self.compression = compression
        self.node_name = "Gene"
        self.edge = "Organism-ENCODES-Gene"
        
        # metadata: property, type, description, example
        node_properties = [
            ["id", "string", "gene identifier", "ncbigene:59272"],
            ["name", "string", "gene symbol", "ACE2"],
            ["synonyms", "string", "alternative gene symbols", "ACEH"],
            ["description", "string", "gene name", "angiotensin converting enzyme 2"],
            ["geneType", "string", "gene type", "protein-coding"],
            ["taxonomyId", "string", "NCBI taxonomy id", "taxonomy:9606"],
        ]
        edge_properties = [
            ["from", "string", "NCBI taxonomy id", "taxonomy:9606"],
            ["to", "string", "gene identifier", "ncbigene:59272"],
        ]
            
        self.node_metadata = {self.node_name, node_properties}
        self.edge_metadata = {self.edge_name, edge_properties}

    def extract_data(self):
        filename, file_extension = os.path.splitext(self.filename)

        with gzip.open(self.filename, "rb") as f_in:
            with open(filename, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        column_names = ["GeneID", "Symbol", "Synonyms", "description", "type_of_gene", "#tax_id"]
        new_names = [x[0] for x in self.node_properties]

        genes = dd.read_csv(filename, usecols=column_names, dtype=str, sep="\t", blocksize="0.1 GB")
        genes = genes.rename(columns=dict(zip(column_names, new_names)))

        genes = genes.replace("-", "")
        genes["id"] = "ncbigene:" + genes["id"]
        genes["taxonomyId"] = "taxonomy:" + genes["taxonomyId"]

        self.save_nodes(genes, self.node_name)

        edges = genes[["taxonomyId", "id"]]
        edges = edges.rename(columns={"taxonomyId": "from", "id": "to"})

        self.save_edges(edges, self.edge)
