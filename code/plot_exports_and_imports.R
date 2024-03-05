################################################################################
################### DENV phylogenetic summaries ##########################
################################################################################

########################### Bernardo Gutierrez - adapated by Rhys Inward #################################

# List of required packages
required_packages <- c("optparse", "tidyverse", "janitor", "lubridate",
                       "stringr", "patchwork","NatParksPalettes","RColorBrewer",
                       "ape","phylobase","scales")

options(repos = c(CRAN = "http://cran.us.r-project.org"))

# Function to check and install packages
suppressMessages(
  for (package in required_packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, repos = "http://cran.us.r-project.org")
    }
    library(package, character.only = TRUE)
  }
)

# Define and parse command-line options
opt_parser <- OptionParser(
  option_list = list(
    make_option(c("-m", "--metadata"), type="character", help="Input csv file containing cleaned metadata, with sequence identifiers matching those in the input fasta file."),
    make_option(c("-o", "--outfile"), type="character", default="subsampled", help="Base name for output files."),
    make_option(c("-c", "--country"), type="character", help="Country to plot viral movements from."),
    make_option(c("-s", "--serotype"), type="character", help="Serotype to plot viral movements from."),
    make_option(c("-d", "--output_dir"), type="character", default="subsampled", help="Output directory"),
    make_option(c("-a", "--output_dir_export"), type="character", default="subsampled", help="Output directory for export metadata"),
    make_option(c("-b", "--output_dir_import"), type="character", default="subsampled", help="Output directory for import metadata")
    
    
  )
)

opt = parse_args(opt_parser)


############################## Data wrangling ##################################
# Deconstructed ph# Deconstructed ph# Deconstructed phylogenetic trees as data frames
# Generated through custom Python scripts (by Joseph Tsui) from phylogenetic
# pipeline (by Rhys Inward).
if (!is.null(opt$metadata)) {
  metadata.df <- read.csv(opt$metadata,sep ="\t")
} else {
  cat("Input metadata file. Exiting now...")
  quit()
}


####################### Cross-border viral movements ###########################
## Generate data frame counting viral movements between neighbours ####

# Remove branches with two nodes in the same country
metadata.df_intl <- metadata.df |> filter(head_country != tail_country)

# Plot numbers of international imports into specific countries
country <- opt$country
serotype <- opt$serotype

#add serotype 

metadata.df_intl$serotype <- serotype

plot_data_2 <- metadata.df_intl |>
  filter(tail_country == country) |>
  group_by(serotype, head_country) |>
  summarise(count = n())

plot_data_2 |> 
  group_by(serotype) |>
  mutate(total = sum(count))

write.csv(plot_data_2,
          file = paste0(opt$output_dir_import),
          row.names = FALSE,quote=FALSE)

total_imports <- ggplot(plot_data_2, aes(x = head_country, y = count, fill = count)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  labs(y = "No. of inferred imports", x = "Source country", fill = "Number of Imports",
       title = paste0("DENV imports to ", country)) +
  theme_minimal() +
  scale_fill_gradient(low = "lightblue", high = "darkblue")

ggsave(filename = paste0(opt$output_dir, "total_imports_", serotype, ".pdf"), plot = total_imports)

#imports over time 

plot_data_over_time <- metadata.df_intl |>
  filter(tail_country == country) |>
  mutate(Year = year(tail_date)) %>%  
  group_by(serotype,Year) |>
  summarise(count = n())

imports_over_time <- ggplot(plot_data_over_time, aes(x = Year, y = count, fill = count)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  labs(y = "No. of inferred imports", x = "Date", fill = "Number of Imports",
       title = paste0("DENV imports to ", country)) +
  theme_minimal() +
  scale_fill_gradient(low = "lightblue", high = "darkblue") 

ggsave(filename = paste0(opt$output_dir, "total_imports_over_time", serotype, ".pdf"), plot = imports_over_time)

# Plot numbers of international export from specific countries

plot_data_3 <- metadata.df_intl |>
  filter(head_country == country) |>
  group_by(serotype, tail_country) |>
  summarise(count = n())

write.csv(plot_data_3,
          file = paste0(opt$output_dir_export),
          row.names = FALSE,quote=FALSE)


total_exports <- ggplot(plot_data_3, aes(x = tail_country, y = count, fill = count)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  labs(y = "No. of inferred exports", x = "Source country", fill = "Number of export",
       title = paste0("DENV export from ", country)) +
  theme_minimal() +
  scale_fill_gradient(low = "lightblue", high = "darkblue")

ggsave(filename = paste0(opt$output_dir,"total_export_", serotype, ".pdf"), plot = total_exports)

#export over time 

plot_data_over_time_exports <- metadata.df_intl |>
  filter(head_country == country) |>
  mutate(Year = year(tail_date)) %>%  
  group_by(serotype,Year) |>
  summarise(count = n())

exports_over_time <- ggplot(plot_data_over_time_exports, aes(x = Year, y = count, fill = count)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  labs(y = "No. of inferred export", x = "Date", fill = "Number of export",
       title = paste0("DENV export from ", country)) +
  theme_minimal() +
  scale_fill_gradient(low = "lightblue", high = "darkblue") 

ggsave(filename = paste0(opt$output_dir, "total_exports_over_time", serotype, ".pdf"), plot = exports_over_time)


