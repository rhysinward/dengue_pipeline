import json
from augur.utils import json_to_tree

# Open and load the JSON file
json_file_path = "~/rhys/Brazil_phylo/auspice/dengue_Dengue_3.json"  # Replace with your actual file path
with open(json_file_path, "r") as json_fh:
    json_dict = json.load(json_fh)

# Convert JSON to tree
tree = json_to_tree(json_dict, root=True)
Phylo.write(tree, "output_tree.nexus", "nexus")
