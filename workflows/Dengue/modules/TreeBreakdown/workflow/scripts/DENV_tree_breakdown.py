#!/usr/bin/env python
# coding: utf-8

# Load in modules
import sys
import pandas as pd
import math
import datetime as dt
import dendropy
from tqdm import tqdm
import json

# Process date formats
def get_date_str(dec_date):
    date = dt.datetime(int(dec_date), 1, 1) + dt.timedelta(days=(dec_date % 1) * 365)
    return dt.datetime.strftime(date, '%Y-%m-%d')

if __name__ == "__main__":
    # Command line arguments for input and output files
    input_tree_file = sys.argv[1]
    output_csv_file = sys.argv[2]

    # Read in tree file
    tree = dendropy.Tree.get(path=input_tree_file, schema='nexus')

    # Name internal nodes
    count = 0
    for node in tree.preorder_node_iter():
        count += 1
        node.label = str(count)  # Ensure the label is a string
    
    # Store information from tree
    edge_info = []
    for edge in tqdm(tree.postorder_edge_iter()):
        if edge.tail_node:
            head_node = edge.tail_node  # tail_node in edge indicates older node
            head_dec_date = float(head_node.annotations['num_date'].value)
            head_date = get_date_str(head_dec_date)
            
            tail_node = edge.head_node  # head_node in edge indicates younger node
            tail_dec_date = float(tail_node.annotations['num_date'].value)
            tail_date = get_date_str(tail_dec_date)
                            
            head_country = head_node.annotations['Country'].value
            tail_country = tail_node.annotations['Country'].value
                    
            edge_info.append(
                {
                    'head_node': head_node.label,
                    'tail_node': tail_node.label,
                    'length': edge.length,
                    'head_country': head_country,
                    'tail_country': tail_country,
                    'head_date': head_date,
                    'head_dec_date': head_dec_date,
                    'tail_date': tail_date,
                    'tail_dec_date': tail_dec_date
                }
            )

    tree_edge_df = pd.DataFrame(edge_info)

    # Writing the DataFrame to a TSV file correctly indented
    tree_edge_df.to_csv(output_csv_file, sep='\t', index=False)
