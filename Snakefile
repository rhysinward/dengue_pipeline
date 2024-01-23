# Serotypes
serotype = ['Dengue_1', 'Dengue_2', 'Dengue_3', 'Dengue_4']

rule all:
    input:
       expand("auspice/dengue_{serotype}.json", serotype=serotype)

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
    log:
        "logs/acquire_data.log"
    message:
        "Acquiring Genomic Data and Metadata from NCBI database"
    conda:
        "config/ncbi_datasets_env.yaml"
    shell:
        """
        datasets download virus genome taxon "Dengue Virus" --filename {output.zip} &> {log}
        datasets summary virus genome taxon "Dengue Virus" --released-after 2000-01-01 --as-json-lines | dataformat tsv virus-genome > {output.metadata} 2>> {log}
        
        # Unzip the downloaded data
        unzip -d data {output.zip}
        """

# Step 2: Clean metadata and generate FASTA
rule process_dengue_data:
    input:
        script="code/Clean_metadata_and_fasta.R",
        metadata="data/metadata.tsv",
        fasta="data/ncbi_dataset/data/genomic.fna"
    output:
        fasta_files = expand("results/Unaligned_{serotype}.fasta", serotype=serotype),
        info_tables_txt = expand("results/{serotype}_infoTbl.txt", serotype=serotype),
        info_tables_csv = expand("results/{serotype}_infoTbl.csv", serotype=serotype)
    params:
        start_date = "2000-01-01",
        end_date = "2023-12-24"
    log:
        "logs/process_dengue_data.log"
    message:
        "Processing and cleaning Dengue data downloaded from NCBI"
    shell:
        """
        rscript {input.script} --metadata {input.metadata} --fasta {input.fasta} --start-date {params.start_date} --end-date {params.end_date}
        """

# Step 3: Sequence alignment
rule sequence_alignment:
    input:
        fasta_files = expand("results/Unaligned_{serotype}.fasta", serotype=serotype)
    output:
        align_dir = expand("results/Aligned_{serotype}/nextalign.aligned.fasta", serotype=serotype)
    log:
        "logs/sequence_alignment.log"
    conda:
        "config/alignment.yaml"
    message:
        "Running sequence alignment for Dengue serotypes"
    shell:
        """
        bash bash/align_sequences.sh
        """

# Step 4: Segregating E gene and Whole Genomes and performing quaility control
rule split_genome_and_QC:
    input:
        fasta_files = expand("results/Aligned_{serotype}/nextalign.aligned.fasta", serotype=serotype)
    output:
        E_gene_dir = expand("results/{serotype}_EG.fasta", serotype=serotype),
        WG_gene_dir = expand("results/{serotype}_WG.fasta", serotype=serotype)
    params:
        wg_threshold = 0.29,
        eg_threshold = 0.05
    log:
        "logs/process_dengue_data.log"
    message:
        "Segregating E gene and whole genomes from aligned Dengue virus sequences and performing quality control."
    shell:
        """
        Rscript Code/Seperate_EG_and_WG.R --WG_threshold {params.wg_threshold} --EG_threshold {params.eg_threshold}
        """

# Step 5a: Subsampler DENV1
rule subsample_denv1:
    input:
        fasta_files = "results/Dengue_1_EG.fasta",
        metadata_files = "results/Dengue_1_infoTbl.csv"
    output:
        subsample_fasta = "results/subsampled_Dengue_1.fasta",
        subsample_csv = "results/subsampled_Dengue_1_infoTbl.csv",
        subsample_txt = "results/subsampled_Dengue_1_infoTbl.tsv"
    params:
        country = "Colombia",
        number_sequences = 500,
        prop_rd = 0.8,
        prop_or = 0.2
    log:
        "logs/subsample.log"
    message:
        "Subsampling Dengue 1 virus E gene sequences based on specified criteria."
    shell:
        """
        rscript Code/subsampler.R --metadata '{input.metadata_files}' --fasta '{input.fasta_files}' --country {params.country} --number_sequences {params.number_sequences} --prop_rd {params.prop_rd} --prop_or {params.prop_or}
        """

# Step 5b: Subsampler DENV2
rule subsample_denv2:
    input:
        fasta_files = "results/Dengue_2_EG.fasta",
        metadata_files = "results/Dengue_2_infoTbl.csv"
    output:
        subsample_fasta = "results/subsampled_Dengue_2.fasta",
        subsample_csv = "results/subsampled_Dengue_2_infoTbl.csv",
        subsample_txt = "results/subsampled_Dengue_2_infoTbl.tsv"
    params:
        country = "Colombia",
        number_sequences = 500,
        prop_rd = 0.8,
        prop_or = 0.2
    log:
        "logs/subsample.log"
    message:
        "Subsampling Dengue 2 virus E gene sequences based on specified criteria."
    shell:
        """
        rscript Code/subsampler.R --metadata '{input.metadata_files}' --fasta '{input.fasta_files}' --country {params.country} --number_sequences {params.number_sequences} --prop_rd {params.prop_rd} --prop_or {params.prop_or}
        """

# Step 5c: Subsampler DENV3
rule subsample_denv3:
    input:
        fasta_files = "results/Dengue_3_EG.fasta",
        metadata_files = "results/Dengue_3_infoTbl.csv"
    output:
        subsample_fasta = "results/subsampled_Dengue_3.fasta",
        subsample_csv = "results/subsampled_Dengue_3_infoTbl.csv",
        subsample_txt = "results/subsampled_Dengue_3_infoTbl.tsv"
    params:
        country = "Colombia",
        number_sequences = 500,
        prop_rd = 0.8,
        prop_or = 0.2
    log:
        "logs/subsample.log"
    message:
        "Subsampling Dengue 3 virus E gene sequences based on specified criteria."
    shell:
        """
        rscript Code/subsampler.R --metadata '{input.metadata_files}' --fasta '{input.fasta_files}' --country {params.country} --number_sequences {params.number_sequences} --prop_rd {params.prop_rd} --prop_or {params.prop_or}
        """

# Step 5D: Subsampler DENV4
rule subsample_denv4:
    input:
        fasta_files = "results/Dengue_4_EG.fasta",
        metadata_files = "results/Dengue_4_infoTbl.csv"
    output:
        subsample_fasta = "results/subsampled_Dengue_4.fasta",
        subsample_csv = "results/subsampled_Dengue_4_infoTbl.csv",
        subsample_txt = "results/subsampled_Dengue_4_infoTbl.tsv"
    params:
        country = "Colombia",
        number_sequences = 500,
        prop_rd = 0.8,
        prop_or = 0.2
    log:
        "logs/subsample.log"
    message:
        "Subsampling Dengue 4 virus E gene sequences based on specified criteria."
    shell:
        """
        rscript Code/subsampler.R --metadata '{input.metadata_files}' --fasta '{input.fasta_files}' --country {params.country} --number_sequences {params.number_sequences} --prop_rd {params.prop_rd} --prop_or {params.prop_or}
        """

# Step 6: Treebuilding
rule iqtree:
    input:
        aln = "results/subsampled_{serotype}.fasta"
    output:
        tree = "results/subsampled_{serotype}.fasta.treefile"
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
        iqtree2 -s {input.aln} -m {params.model} -redo
        """


# Step 7: Building Time-Calibrated Trees
rule treetime:
    input:
        aln = "results/subsampled_{serotype}.fasta",
        ml_tree = "results/subsampled_{serotype}.fasta.treefile",
        metadata = "results/subsampled_{serotype}_infoTbl.csv"
    output:
        tree = "results/timetrees/timetree_{serotype}.nwk",
        branch_lengths = "results/timetrees/{serotype}_branch_lengths.json"
    conda:
        "config/nextstrain_all.yaml"
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
            --date-inference marginal \
            --clock-filter-iqd 3 \
            --root best
        """
# Step 7: Infer "ancestral" mutations across the tree
rule mutations:
    input:
        aln = "results/subsampled_{serotype}.fasta",
        time_tree = "results/timetrees/timetree_{serotype}.nwk"
    output:
        mutations = "results/subsampled_{serotype}_mutations.json"
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
            --inference joint
        """

# Step 8: Translate sequences
rule translation:
    input:
        mutations = "results/subsampled_{serotype}_mutations.json",
        time_tree = "results/timetrees/timetree_{serotype}.nwk",
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
            --output {output.amino} 
        """

# Step 9: Discrete trait reconstruction
rule mugration:
    input:
        time_tree = "results/timetrees/timetree_{serotype}.nwk",
        metadata = "results/subsampled_{serotype}_infoTbl.csv"
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
            --output {output.traits} 
        """

# Step 10: Export for visualisation in Auspice

rule export:
    """Exporting data files for auspice"""
    input:
        tree = "results/timetrees/timetree_{serotype}.nwk",
        metadata = "results/subsampled_{serotype}_infoTbl.csv",
        branch_lengths = "results/timetrees/{serotype}_branch_lengths.json",
        traits = "results/traits_{serotype}.json",
        aa_muts = "results/aa_muts_{serotype}.json",
        nt_muts = "results/subsampled_{serotype}_mutations.json",
        auspice_config = "config/auspice_config.json",
    output:
        auspice_json = "auspice/dengue_{serotype}.json",
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
            --output {output.auspice_json}
        """
