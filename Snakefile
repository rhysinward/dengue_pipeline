# Serotypes
serotype = ["Dengue_1", "Dengue_2", "Dengue_3", "Dengue_4"]

# Sequence sampling sizes
subsampling_counts = {
    "Dengue_1": {
        "local": 100,
        "background": 2000,
        "background_method": "even"
    },
    "Dengue_2": {
        "local": "all",
        "background": 2000,
        "background_method": "even"
    },
    "Dengue_3": {
        "local": "all",
        "background": 2000,
        "background_method": "even"
    },
    "Dengue_4": {
        "local": "all",
        "background": 1500,
        "background_method": "even"
    }
}


rule all:
    input:
        expand("results/imports_{serotype}.csv", serotype=serotype)

# Rule for creating necessary directories
rule create_directories:
    output:
        directory("logs"),
        directory("data"),
        directory("results")
    shell:
        "mkdir -p logs data results"

# Step 1: Acquisition of Genomic Data and Metadata from GenBank
rule acquire_data:
    output:
        zip="data/genbank_data.zip",
        metadata="data/metadata.tsv",
        fasta="data/ncbi_dataset/data/genomic.fna"
    params:
        date = "01/01/2010",
    log:
        "logs/acquire_data.log"
    message:
        "Acquiring Genomic Data and Metadata from NCBI database"
    conda:
        "config/ncbi_datasets_env.yaml"
    shell:
        """
        datasets download virus genome taxon "Dengue Virus" --filename {output.zip} --released-after {params.date} &> {log}
        datasets summary virus genome taxon "Dengue Virus" --released-after {params.date} --as-json-lines | dataformat tsv virus-genome > {output.metadata} 2>> {log}
        
        # Unzip the downloaded data
        unzip -o -d data {output.zip}
        """

# Step 2: Clean metadata and FASTA
rule process_genbank_data:
    input:
        script="code/clean_metadata_and_fasta_general.R",
        metadata="data/metadata.tsv",
        fasta="data/ncbi_dataset/data/genomic.fna",
        extra_metadata="data/extra/metadata_extra.tsv",
        extra_fasta="data/extra/genomic_extra.fna"
    output:
        fasta_files ="results/Unaligned.fasta",
        info_tables_txt ="results/infoTbl.txt",
        info_tables_csv ="results/infoTbl.csv"
    params:
        start_date = "2010-01-01",
        end_date = "2024-02-28",
        host = "Homo sapiens"

    log:
        "logs/process_genomic_data.log"
    message:
        "Processing and cleaning data downloaded from NCBI"
    shell:
        """
        Rscript {input.script} \
            --metadata {input.metadata} \
            --fasta {input.fasta} \
            --extra_metadata {input.extra_metadata} \
            --extra_fasta {input.extra_fasta} \
            --start_date {params.start_date} \
            --end_date {params.end_date} \
            --host "{params.host}" \
            --outfile_fasta {output.fasta_files} \
            --outfile_csv {output.info_tables_csv} \
            --outfile_tsv {output.info_tables_txt} \
            > {log} 2>&1
        """

# Step 3: Split into serotype, add serotypes to sequence name and generate sequence specific metadata
rule process_dengue_data:
    input:
        script="code/split_dengue.R",
        metadata="results/infoTbl.csv",
        fasta="results/Unaligned.fasta"
    output:
        fasta = "results/Unaligned_output/Unaligned_{serotype}.fasta",
        csv = "results/Unaligned_output/Unaligned_{serotype}_infoTbl.csv",
        tsv = "results/Unaligned_output/Unaligned_{serotype}_infoTbl.txt"
    params:
        out_prefix="_results/s04_process_dengue_data/Unaligned_",
        serotype = "{serotype}"
    log:
        "logs/process_data_{serotype}.log"
    message:
        "Split into {wildcards.serotype}, add serotypes to sequence name and generate sequence specific metadata"
    shell:
        """
        mkdir -p results/Unaligned_output/
        Rscript {input.script} \
            --metadata {input.metadata} \
            --fasta {input.fasta} \
            --serotype {params.serotype} \
            --outfile {params.out_prefix} \
            > {log} 2>&1
        """

# Step 4: Sequence alignment
rule sequence_alignment:
    input:
        sequences="results/Unaligned_output/Unaligned_{serotype}.fasta",
        reference="reference_genomes/reference_{serotype}.fasta",
        genemap="genemap/genemap_{serotype}.gff"
    output:
        fasta="results/Aligned_{serotype}/nextalign.aligned.fasta"
    log:
        "logs/sequence_alignment_{serotype}.log"
    conda:
        "config/alignment.yaml"
    message:
        "Running sequence alignment for Dengue serotypes"
    shell:
        """
        mkdir -p results/Aligned_{wildcards.serotype}/data
        nextalign run \
            {input.sequences} \
            --reference={input.reference} \
            --genemap={input.genemap} \
            --output-fasta={output.fasta} > {log} 2>&1
        """


# Step 5: Segregating E gene and Whole Genomes and performing quaility control
rule split_genome_and_QC:
    input:
        script="code/Seperate_EG_and_WG.R",
        fasta = "results/Aligned_{serotype}/nextalign.aligned.fasta"
    output:
        E_gene_dir = "results/{serotype}_EG.fasta",
        WG_gene_dir = "results/{serotype}_WG.fasta"
    params:
        wg_threshold = 0.29,
        eg_threshold = 0.29
    log:
        "logs/process_dengue_data.log"
    message:
        "Segregating E gene and whole genomes from aligned Dengue virus sequences and performing quality control."
    shell:
        """
        Rscript {input.script} \
            --fasta {input.fasta} \
            --WG_threshold {params.wg_threshold} \
            --EG_threshold {params.eg_threshold} \
            --outfile_fasta_eg {output.E_gene_dir} \
            --outfile_fasta_wg {output.WG_gene_dir} \
            > {log} 2>&1
        """

# Step 6: Subsampling DENVs
rule subsample_denv:
    input:
        script = "code/subsampler.R",
        metadata_file = "results/Unaligned_output/Unaligned_{serotype}_infoTbl.csv",
        fasta_file = "results/{serotype}_EG.fasta",
        location_local = "data/number_of_sequences.csv"
    output:
        subsample_fasta = "results/subsampled_{serotype}.fasta",
        subsample_csv = "results/subsampled_{serotype}_infoTbl.csv",
        subsample_txt = "results/subsampled_{serotype}_infoTbl.tsv"
    params:
        number_sequences_local = lambda wc: subsampling_counts[wc.serotype]["local"],
        number_sequences_background = lambda wc: subsampling_counts[wc.serotype]["background"],
        sampling_method = lambda wc: subsampling_counts[wc.serotype]["background_method"],
        time_interval = "Year",
        serotype = "{serotype}",
        out_prefix = "results/subsampled_{serotype}",
        out_dir = "_results/subsampling_plots/"
    log:
        "logs/subsample_{serotype}.log"
    message:
        "Subsampling {wildcards.serotype} virus E gene sequences based on specified criteria."
    shell:
        """
        Rscript {input.script} \
            --metadata {input.metadata_file} \
            --fasta {input.fasta_file} \
            --location_local {input.location_local} \
            --time_interval {params.time_interval} \
            --number_sequences_local {params.number_sequences_local} \
            --number_sequences_background {params.number_sequences_background} \
            --sampling_method {params.sampling_method} \
            --serotype {params.serotype} \
            --out_prefix {params.out_prefix} \
            --output_dir {params.out_dir} \
            > {log} 2>&1
        """

# Step 7: Correct metadata and fasta files into the correct format for iqtree and treetime  
rule reformatting:
    input:
        script = "code/reformatting_iqtree_treetime.R",
        fasta_file = "results/subsampled_{serotype}.fasta",
        metadata_file = "results/subsampled_{serotype}_infoTbl.csv",
    output:
        cleaned_fasta = "results/subsampled_{serotype}_cleaned.fasta",
        cleaned_metadata = "results/subsampled_{serotype}_cleaned_infoTbl.csv"    
    log:
        "logs/reformatting_{serotype}.log"    
    message:
        "Correct metadata and fasta files into the correct format for iqtree and treetime"
    shell:
        """
        Rscript {input.script} \
            --metadata {input.metadata_file} \
            --fasta {input.fasta_file} \
            --output_dir_fasta {output.cleaned_fasta} \
            --output_dir_csv {output.cleaned_metadata} \
            > {log} 2>&1
        """

# Step 8: Treebuilding
rule treebuilding:
    input:
        aln = "results/subsampled_{serotype}_cleaned.fasta" 
    output:
        tree = "results/subsampled_{serotype}_cleaned.fasta.treefile" 
    conda:
        "config/iqtree.yaml"
    params:
        model = "GTR+F+I"
    log:
        "logs/iqtree_{serotype}.log"
    message:
        "Inferring maximum likelihood phylogenetic trees for {wildcards.serotype} using IQ-TREE."
    shell:
        """
        iqtree2 -s {input.aln} -m {params.model} -redo > {log} 2>&1
        """

# Step 9: Building Time-Calibrated Trees
rule treetime:
    input:
        aln = "results/subsampled_{serotype}_cleaned.fasta",
        ml_tree = "results/subsampled_{serotype}_cleaned.fasta.treefile",
        metadata = "results/subsampled_{serotype}_cleaned_infoTbl.csv"
    output:
        tree = "results/timetrees/timetree_{serotype}.tree",
        branch_lengths = "results/timetrees/{serotype}_branch_lengths.json"
    conda:
        "config/nextstrain_all.yaml"
    params:
        clock_filter = 3
    log:
        "logs/treetime_{serotype}.log"
    message:
        "Inferring time-calibrated trees for each Dengue virus serotype."
    shell:
        """
        augur refine \
            --tree {input.ml_tree} \
            --alignment {input.aln} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.branch_lengths} \
            --timetree \
            --coalescent opt \
            --date-confidence \
            --clock-filter-iqd {params.clock_filter} \
            --root best \
            > {log} 2>&1
        """


# Step 10: Infer "ancestral" mutations across the tree
rule mutations:
    input:
        aln = "results/subsampled_{serotype}_cleaned.fasta",
        time_tree = "results/timetrees/timetree_{serotype}.tree"
    output:
        mutations = "results/subsampled_{serotype}_cleaned_mutations.json"
    conda:
        "config/nextstrain_all.yaml"
    log:
        "logs/mutations_{serotype}.log"
    message:
        "Infer ancestral mutations across the tree"
    shell:
        """
        augur ancestral \
            --tree {input.time_tree} \
            --alignment {input.aln} \
            --inference joint \
            > {log} 2>&1
        """

# Step 11: Translate sequences
rule translation:
    input:
        mutations = "results/subsampled_{serotype}_cleaned_mutations.json",
        time_tree = "results/timetrees/timetree_{serotype}.tree",
        ref_genomes = "reference_genomes/reference_{serotype}.gb"
    output:
        amino = "results/aa_muts_{serotype}.json"
    conda:
        "config/nextstrain_all.yaml"
    log:
        "logs/translations_{serotype}.log"
    params:
        genes = "E"
    message:
        "Translate sequences"
    shell:
        """
        augur translate \
            --tree {input.time_tree} \
            --ancestral-sequences {input.mutations} \
            --reference-sequence {input.ref_genomes} \
            --genes {params.genes} \
            --output {output.amino} \
            > {log} 2>&1
        """

# Step 12: Discrete trait reconstruction
rule mugration:
    input:
        time_tree = "results/timetrees/timetree_{serotype}.tree",
        metadata = "results/subsampled_{serotype}_cleaned_infoTbl.csv"
    output:
        traits = "results/traits_{serotype}.json"
    conda:
        "config/nextstrain_all.yaml"
    log:
        "logs/mugration_{serotype}.log"
    message:
        "Use Discrete Trait Reconstruction for country and state"
    shell:
        """
        augur traits \
            --tree {input.time_tree} \
            --metadata {input.metadata} \
            --columns Country State \
            --confidence \
            --output {output.traits} \
            > {log} 2>&1
        """

# Step 13: Export for visualisation in Auspice

rule export:
    """Exporting data files for auspice"""
    input:
        tree = "results/timetrees/timetree_{serotype}.tree",
        metadata = "results/subsampled_{serotype}_cleaned_infoTbl.csv",
        branch_lengths = "results/timetrees/{serotype}_branch_lengths.json",
        traits = "results/traits_{serotype}.json",
        aa_muts = "results/aa_muts_{serotype}.json",
        nt_muts = "results/subsampled_{serotype}_cleaned_mutations.json",
        auspice_config = "config/auspice_config.json",
    output:
        auspice_json = "auspice/dengue_{serotype}.json"

    conda:
        "config/nextstrain_all.yaml"
    log:
        "logs/export_{serotype}.log"
    message:
        "Export for visualisation in Auspice"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json} \
            > {log} 2>&1
        """

# Step 14: Extract annotated tree from nextstrain JSON format 

rule extract_phylogenetic_tree:
    input:
        script = "code/extract_tree_from_json.py",
        json_file = "auspice/dengue_{serotype}.json"
    output:
        nexus_file = "results/dengue_{serotype}_timetree.nexus"
    conda:
        "config/python_env.yaml"
    log:
        "logs/extract_tree_{serotype}.log"
    message:
        "Extract annotated tree from nextstrain JSON format "
    shell:
        """
        python {input.script} \
            {input.json_file} \
            {output.nexus_file} \
            > {log} 2>&1
        """

# Step 15: Extract information from tree 

rule extract_information_from_phylogenetic_tree:
    input:
        script = "code/DENV_tree_breakdown.py",
        tree_file = "results/dengue_{serotype}_timetree.nexus"
    output:
        tsv_file = "results/dengue_{serotype}_timetree_extracted.tsv"
    conda:
        "config/python_env.yaml"
    log:
        "logs/extract_information_from_tree_{serotype}.log"
    message:
        "Extract information from annotated tree"
    shell:
        """
        python {input.script} \
            {input.tree_file} \
            {output.tsv_file} \
            > {log} 2>&1
        """

# Step 16: Quantify number of exports and imports from desired country

rule plot_exports_and_imports:
    input:
        script = "code/plot_exports_and_imports.R",
        metadata = "results/dengue_{serotype}_timetree_extracted.tsv"
    output:
        imports = "results/imports_{serotype}.csv",
        exports = "results/exports_{serotype}.csv"
    params:
        country = "Vietnam",
        serotype = "{serotype}"
    log:
        "logs/exports_and_imports_{serotype}.log"
    message:
        "Quantify number of exports and imports from desired country"
    shell:
        """
        Rscript {input.script} \
            --metadata {input.metadata} \
            --output_dir results/  \
            --output_dir_export {output.exports} \
            --output_dir_import {output.imports} \
            --country {params.country} \
            --serotype {params.serotype} \
            > {log} 2>&1
        """

# Step 17: Create updated metadata and fasta file based of prunning from treetime

# Step 18: BEASTGEN 

# Step 19: BEAST 



