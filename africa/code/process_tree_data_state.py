import argparse
import pandas as pd
import math
import datetime as dt
import dendropy
import tqdm
import json

# Constants
DATE_FORMAT = '%Y-%m-%d'

def get_date_str(dec_date):
    """
    Converts a decimal date to a string in YYYY-MM-DD format.
    """
    date = dt.datetime(int(dec_date), 1, 1) + dt.timedelta(days=(dec_date % 1) * 365)
    return dt.datetime.strftime(date, DATE_FORMAT)

# Setting up the argument parser
parser = argparse.ArgumentParser(description='Annotate and count the number of state changes from the results of TreeTime migration')
parser.add_argument('-t', '--tree', metavar='INPUT FILE', action='store', type=str, required=True, help='Unprocessed output file from TreeTime migration (with country as an attribute)')
parser.add_argument('-ot', '--outfile_tsv', metavar='OUTPUT_TSV', required=True, help='Output TSV file to be written to')
parser.add_argument('-on', '--outfile_newick', metavar='OUTPUT_NEWICK', required=True, help='Output Newick tree file to be written to')

args = parser.parse_args()

# Processing the tree
tree = dendropy.Tree.get(path=args.tree, schema='newick')

## Name internal nodes
count = 0
for node in tree.preorder_node_iter():
    count += 1
    node.label = str(count)  # Ensure labels are strings

## Store information from tree
edge_info = []
for edge in tqdm.tqdm(tree.postorder_edge_iter(), desc="Processing edges"):
    if edge.tail_node:
        head_node = edge.tail_node  # Tail node in edge indicates older node
        head_dec_date = float(head_node.annotations['num_date'].value)
        head_date = get_date_str(head_dec_date)
        
        tail_node = edge.head_node  # Head node in edge indicates younger node
        tail_dec_date = float(tail_node.annotations['num_date'].value)
        tail_date = get_date_str(tail_dec_date)
                        
        head_country = head_node.annotations.get_value('Country', default=None)
        tail_country = tail_node.annotations.get_value('Country', default=None)
        
        # Initialize location fields
        head_location = head_country
        tail_location = tail_country
        
        # Check if the country is Brazil and retrieve state if available
        if head_country == 'Brazil':
            head_state = head_node.annotations.get_value('State', default='Unknown')
            head_location = f"Brazil - {head_state}"
        
        if tail_country == 'Brazil':
            tail_state = tail_node.annotations.get_value('State', default='Unknown')
            tail_location = f"Brazil - {tail_state}"
        
        # Append edge info including head_location and tail_location
        edge_info.append(
            {
                'head_node': head_node.label,
                'tail_node': tail_node.label,
                'length': edge.length,
                'head_location': head_location,
                'tail_location': tail_location,
                'head_date': head_date,
                'head_dec_date': head_dec_date,
                'tail_date': tail_date,
                'tail_dec_date': tail_dec_date
            }
        )

# Writing results to TSV file
tree_edge_df = pd.DataFrame(edge_info)

# Optionally, you can reorder columns or select specific columns
# For example:
# tree_edge_df = tree_edge_df[['head_node', 'tail_node', 'length', 'head_location', 'tail_location', 'head_date', 'tail_date']]

tree_edge_df.to_csv(args.outfile_tsv, sep='\t', index=False)

# Writing the tree to a Newick file
tree.write(path=args.outfile_newick, schema='newick')
