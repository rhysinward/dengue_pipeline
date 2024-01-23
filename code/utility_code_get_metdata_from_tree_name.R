# Load required packages
required_packages <- c("dplyr", "lubridate", "tidyr", "ape", "seqinr")
suppressMessages(
  lapply(required_packages, function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, repos = "http://cran.us.r-project.org")
    }
    library(pkg, character.only = TRUE)
  })
)

process_dengue_data <- function(serotype) {
  tree_file_path <- sprintf("/Users/rhysinward/Documents/Dengue_anaysis/results/treetime_round2/Dengue_%d_cleaned.fasta.treefile", serotype)
  tree <- read.tree(tree_file_path)
    taxa_data <- data.frame(tree$tip.label)
  taxa_split_data <- data.frame(do.call('rbind', strsplit(as.character(taxa_data$tree.tip.label), '|', fixed = TRUE)))
  taxa_split_data$name <- taxa_data$tree.tip.label
  colnames(taxa_split_data) <- c('GenBank_ID', "Country", "State", "City", "Serotype", "date","Decimal_Date","name")
    output_file_path <- sprintf("/Users/rhysinward/Documents/Dengue_anaysis/results/treetime_round2/Dengue_%d_cleaned_metadata.csv", serotype)
  write.table(taxa_split_data, file = output_file_path, row.names = FALSE, sep = ",", quote = FALSE)
}

# Process data for each Dengue serotype
for (i in 1:4) {
  process_dengue_data(i)
}

