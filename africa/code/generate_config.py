# generate_config.py
import pandas as pd
import os
import yaml

serotypes = ['Dengue_1', 'Dengue_2', 'Dengue_3', 'Dengue_4']
major_lineages = []

for st in serotypes:
    lineage_info_path = f"results/Unaligned_output/Unaligned_{st}_infoTbl_with_genotype_filtered.csv"
    if os.path.exists(lineage_info_path):
        df = pd.read_csv(lineage_info_path)
        major_lineages.extend(df['Major_Lineage'].unique().tolist())
    else:
        print(f"Warning: {lineage_info_path} does not exist. Skipping.")

major_lineages = list(set(major_lineages))

config = {
    'serotypes': serotypes,
    'major_lineages': major_lineages
}

with open('config.yaml', 'w') as f:
    yaml.dump(config, f)
