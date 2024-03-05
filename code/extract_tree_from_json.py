import sys
import json
from ete3 import Tree

def main(json_file_path, output_file_path):
    # read in json file
    with open(json_file_path, 'r') as infile:
        data = json.load(infile)

    # extract tree
    data_tree = data['tree']

    # construct tree
    tree = Tree(name=data_tree['name'])
    tree.dist = 0
    tree.add_features(Country=data_tree['node_attrs']['Country']['value'],
                      num_date=data_tree['node_attrs']['num_date']['value'])

    # recursive function to add children
    def add_children(node_obj, node):
        for child_obj in node_obj['children']:
            child_node = node.add_child(name=child_obj['name'])
            child_node.dist = child_obj['node_attrs']['div']
            child_node.add_features(Country=child_obj['node_attrs']['Country']['value'],
                                    num_date=child_obj['node_attrs']['num_date']['value'])
            if 'children' in child_obj and len(child_obj['children']) > 0:
                add_children(child_obj, child_node)

    # build
    add_children(data_tree, tree)

    # export tree as annotated nexus (with annotations)
    pre_output_filename = output_file_path
    tree.write(outfile=pre_output_filename, format=1, features=['Country', 'num_date'], format_root_node=True)

    # fix formatting issues in exported tree
    fix_formatting(pre_output_filename)

def fix_formatting(pre_output_filename):
    # read in exported tree
    with open(pre_output_filename, 'r') as infile:
        tree_str = infile.read().strip()

    # replace &&NHX with &
    tree_str = tree_str.replace('&&NHX:', '&')
    # replace :num_date with ,num_date
    tree_str = tree_str.replace(':num_date', ',num_date')
    # add #NEXUS\nbegin trees;\ntree one = at the beginning of the file
    tree_str = '#NEXUS\nbegin trees;\ntree one = ' + tree_str
    # add end; at the end of the file
    tree_str += ' end;'

    # export fixed tree
    with open(pre_output_filename, 'w') as outfile:
        outfile.write(tree_str)

if __name__ == "__main__":
    json_file_path = sys.argv[1]
    output_file_path = sys.argv[2]
    main(json_file_path, output_file_path)
