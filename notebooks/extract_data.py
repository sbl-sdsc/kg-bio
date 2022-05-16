#!/usr/bin/env python
# coding: utf-8
__all__ = ["ExtractData"]

import os
import gzip
import time
import pandas as pd
import re
from lxml import etree


class ExtractData:
    """
    Generic Object to extract node and edge data.

    It is the parent class of all more specific data extraction classes.
    """

    def __init__(self, filename, fileformat='parquet', compression='brotli'):
        self.filename = filename
        self.fileformat = fileformat
        self.compression = compression
        self.node_metadata = {"node_name": ["property"]}
        self.edge_metadata = {"edge_name": ["property"]}

    def get_node_metadata(self):
        return self.node_metadata
            
    def get_edge_metadata(self):
        return self.edge_metadata
    
    def get_postfix(self):
        """
        Return postfix for output file names
        """
        prefix = re.search("cache/(.+?)/", self.filename).group(1)
        datetime = re.search(f"{prefix}/(.+?)/", self.filename).group(1)
        return f"_{prefix}_{datetime}"

    def get_path(self):
        return self.filename.rsplit("/", maxsplit=1)[0]
    
    def save_nodes(self, dataframe, node_name):
        self.__save_dataframe(dataframe, 'nodes', node_name)
        
    def save_edges(self, dataframe, edge_name):
        self.__save_dataframe(dataframe, 'edges', edge_name)

    def __save_dataframe(self, dataframe, subdirname, name):
        path = os.path.join(self.get_path(), subdirname)
        os.makedirs(os.path.join(path), exist_ok=True)
        
        if self.fileformat == 'parquet':
            print(self.compression)
            dataframe.to_parquet(
            #os.path.join(path, name + self.get_postfix() + '.parquet'), compression=self.compression, index=False)
            os.path.join(path, name + self.get_postfix() + '.parquet'), compression=self.compression)
        elif self.fileformat == 'csv':
            dataframe.to_csv(
            os.path.join(path, name + self.get_postfix() + '.csv'), compression=self.compression, index=False)

        # self.node_metadata = {self.node_name, gene_metadata}
        # self.edge_metadata = {self.edge_name, organism_encodes_gene_metadata}
            
