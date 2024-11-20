#!/usr/bin/env python

import argparse
from Bio import Phylo
from treetime import TreeAnc

def main():
    parser = argparse.ArgumentParser(
        description="Run TreeAnc on a refined phylogenetic tree.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--input-tree', type=str, required=True, help="Input refined tree (Newick format)")
    parser.add_argument('--alignment', type=str, required=True, help="Sequence alignment file (FASTA format)")
    parser.add_argument('--output-tree', type=str, required=True, help="Output tree with TreeAnc applied (Newick format)")
    parser.add_argument('--ignore-missing-alns', action='store_true', help="Ignore missing alignments")
    
    args = parser.parse_args()
    
    # Read the input tree
    tree = Phylo.read(args.input_tree, 'newick')
    
    # Initialize TreeAnc with the specified parameters
    tt = TreeAnc(tree=tree, aln=args.alignment, ignore_missing_alns=args.ignore_missing_alns)
    
    # Perform optimization or any other TreeAnc operations as needed
    tt.optimize_tree(prune_short=True)
    
    # Save the processed tree to the output file
    Phylo.write(tt.tree, args.output_tree, 'newick')
    
    print(f"TreeAnc processing complete. Output saved to {args.output_tree}")

if __name__ == "__main__":
    main()
