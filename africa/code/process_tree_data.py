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
parser = argparse.ArgumentParser(description='Annotate and count the number of state changes from the results of TreeTime mugration')
parser.add_argument('-i', '--tree', metavar='INPUT FILE', action='store', type=str, required=True, help='Unprocessed output file from TreeTime mugration (with country as an attribute)')
parser.add_argument('-ot', '--outfile_tsv', metavar='OUTPUT_TSV', required=True, help='Output TSV file to be written to')
parser.add_argument('-on', '--outfile_newick', metavar='OUTPUT_NEWICK', required=True, help='Output Newick tree file to be written to')

args = parser.parse_args()

# Processing the tree

tree = dendropy.Tree.get(path=args.tree, schema='newick')

## name internal nodes
count = 0
for node in tree.preorder_node_iter():
    count += 1
    node.label = count
    
## store information from tree
edge_info = []
for edge in tqdm.tqdm(tree.postorder_edge_iter()):
    if edge.tail_node:
        head_node = edge.tail_node ## tail_node in edge indicates older node
        head_dec_date = float(head_node.annotations['num_date'].value)
        head_date = get_date_str(head_dec_date)
        
        tail_node = edge.head_node ## head_node in edge indicates younger node
        tail_dec_date = float(tail_node.annotations['num_date'].value)
        tail_date = get_date_str(tail_dec_date)
                        
        head_country = head_node.annotations['country'].value
        tail_country = tail_node.annotations['country'].value
                
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

# Writing results to TSV file
tree_edge_df = pd.DataFrame(edge_info)
tree_edge_df.to_csv(args.outfile_tsv, sep='\t', index=False)

# Writing the tree to a Newick file
tree.write(path=args.outfile_newick, schema='newick')
