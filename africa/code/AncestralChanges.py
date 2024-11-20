#!/usr/bin/env python
"""
Created on Tue 9 April 2019
Authors: David Rasmussen
Edited: Eduan Wilkinson
Refactored to handle dynamic traits and dates: Rhys Inward
"""

import baltic as bt
import pandas as pd
from optparse import OptionParser

def main():
    usage = "Usage: %prog [options]"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-i", "--input", dest="tree_file", default="annotated_tree.nexus",
                      help="Input tree file in Newick or Nexus format. [default: %default]")
    parser.add_option("-m", "--metadata", dest="metadata_file", default=None,
                      help="Metadata CSV file with a 'date' column to determine absolute time. [default: %default]")
    parser.add_option("-o", "--output", dest="out_file", default="annotated_tree_events.csv",
                      help="Output CSV file for annotated tree events. [default: %default]")
    parser.add_option("-t", "--trait", dest="trait", default="country",
                      help="The trait (e.g., 'country' or 'State') to analyze for changes. [default: %default]")
    
    (options, args) = parser.parse_args()
    
    # Load the phylogenetic tree from the specified file
    myTree = bt.loadNewick(options.tree_file, absoluteTime=False)

    # Determine the absolute time for the last sampled tip:
    # If a metadata file is provided, use the maximum date from it
    if options.metadata_file:
        metadata_df = pd.read_csv(options.metadata_file)
        if 'date' in metadata_df.columns:
            # Assume the date column contains numerical or parseable values
            absolute_time = metadata_df['date'].max()
        else:
            raise ValueError("The metadata file does not have a 'date' column.")
    else:
        # If no metadata is provided, raise an error as we cannot determine absolute time
        raise ValueError("A metadata file must be provided to determine the absolute time.")
    
    myTree.setAbsoluteTime(absolute_time)  # Set this to the time of the last sampled tip

    myTree.traverse_tree()  # Required to set heights
    myTree.treeStats()  # Report stats about the tree

    changes = 0
    times = []
    origins = []
    destinations = []

    # Iterate over a flat list of branches in the tree
    for k in myTree.Objects:
        # Retrieve the trait from the current node
        if options.trait in k.traits:
            node_trait = k.traits[options.trait]
        elif options.trait.lower() in k.traits:
            # Handle the case where the trait name might be lower-cased in the data
            node_trait = k.traits[options.trait.lower()]
        else:
            node_trait = 'UNKNOWN'
            k.traits[options.trait] = node_trait

        # Retrieve the trait from the parent node
        if k.parent and k.parent.traits:
            if options.trait in k.parent.traits:
                parent_trait = k.parent.traits[options.trait]
            elif options.trait.lower() in k.parent.traits:
                parent_trait = k.parent.traits[options.trait.lower()]
            else:
                parent_trait = 'UNKNOWN'
        else:
            parent_trait = 'UNKNOWN'
        
        # Check if there's a change in the trait from the parent node to the current node
        if node_trait != parent_trait:
            changes += 1
            times.append(k.absoluteTime)
            origins.append(parent_trait)
            destinations.append(node_trait)
            
    print("Total number of state changes: " + str(changes))

    # Create a DataFrame with the recorded event data
    df = pd.DataFrame({
        'EventTime': times,
        'Origin': origins,
        'Destination': destinations
    })

    # Save the DataFrame to the specified output CSV file
    df.to_csv(options.out_file, index=False)
    print(f"Events written to {options.out_file}")

if __name__ == "__main__":
    main()
