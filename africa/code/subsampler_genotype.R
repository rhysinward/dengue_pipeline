## load packages
required_packages <- c("optparse", "dplyr","lubridate","tidyr",
                       "readr","ape","seqinr","countrycode",
                       "ggplot2","purrr","zoo","rlang","wrswoR")

# Function to load required packages
load_packages <- function(packages) {
  # Iterate over each package
  for (pkg in packages) {
    # Try to load the package quietly
    if (!suppressPackageStartupMessages(require(pkg, character.only = TRUE))) {
      # If the package is not installed, stop and inform the user
      stop(paste0(
        "The required package '", pkg, "' is not installed. ",
        "Please install it using install.packages('", pkg, "') and try again."
      ), call. = FALSE)
    }
  }
}

load_packages(required_packages)

# Define and parse command-line options
opt_parser <- OptionParser(
  option_list = list(
    make_option(c("-m", "--metadata"), type="character", help="Input metadata file."),
    make_option(c("-f", "--fasta"), type="character", help="Input FASTA file."),
    make_option(c("-c", "--location_local"), type="character", help="CSV for location granularity."),
    make_option(c("-x", "--location_background"), type="character", help="CSV for background locations."),
    make_option(c("-t", "--time_interval"), type="character", default='Year', help="Sampling interval (Year, Month, Week)."),
    make_option(c("-n", "--number_sequences_local"), type="character", help="Number of sequences from local."),
    make_option(c("-e", "--number_sequences_background"), type="numeric", help="Number of sequences from background."),
    make_option(c("-w", "--sampling_method"), type="character", default='Even', help="Sampling method: Even or Proportional."),
    make_option(c("-p", "--outfile_plots"), type="character", default="subsampled", help="Base name for output plots."),
    make_option(c("-i", "--outfile_fasta_csv"), type="character", default="subsampled", help="Base name for output FASTA and CSV files.")
  )
)

opt = parse_args(opt_parser)

##########################################################
## main
##########################################################

# SET SEED FOR CONSISTENT SAMPLING REPRODUCTION
set.seed(123)

##########################################################
# Step 1: read in data
##########################################################
## read in input metadata file
if (!is.null(opt$metadata)) {
  metadata.df <- read.csv(opt$metadata)
} else {
  cat("Input metadata file. Exiting now...")
  quit()
}

## read in input fasta file
if (!is.null(opt$fasta)) {
  seqs <- read.fasta(opt$fasta)
} else {
  cat("Input fasta file. Exiting now...")
  quit()
}

if (!is.null(opt$location_background)) {
  proportional_variable <- read.csv(opt$location_background)
}

##########################################################
# Step 2: processing of metadata
##########################################################

#do metadata and fasta file match? If not remove any in the metadata that do not match fasta file 

# Extract identifiers from FASTA file
seq_name <- as.data.frame(as.matrix(attributes(seqs)$names))

# Filter metadata to only include rows with identifiers that match the FASTA file
length <- nrow(metadata.df)

metadata.df <- metadata.df[metadata.df$Sequence_name %in% seq_name$V1, ]
metadata.df <- as.data.frame(metadata.df)
#add year, month, week to metadata

metadata.df <- metadata.df %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
  mutate(Year = year(ymd(Date))) %>%
  mutate(Month = floor_date(ymd(metadata.df$Date), "month")) %>%
  mutate(Week = floor_date(ymd(metadata.df$Date), "week"))

# Check if any rows were removed from the metadata
if (nrow(metadata.df) < length) {
  cat("Some identifiers in the metadata did not match the FASTA file. They have been removed.\n")
}

#assign region (SA, NA, Africa, Asia, Europe, Oceania) 
metadata.df$Country <- ifelse(metadata.df$Country == "Micronesia", "Micronesia (Federated States of)", metadata.df$Country)
metadata.df$Country <- ifelse(metadata.df$Country == "Borneo", "Indonesia", metadata.df$Country)
metadata.df$Country <- ifelse(metadata.df$Country == "South_Korea", "Korea (Republic of)", metadata.df$Country)
metadata.df$Country <- ifelse(metadata.df$Country == "Republic_of_the_Congo", "Congo (Republic of)", metadata.df$Country)

#Pacific ocean is not a country so remove it

remove_countries <- c("Pacific_Ocean")

metadata.df <- metadata.df %>%
  filter(!Country %in% remove_countries)

#add continent to metadata

metadata.df$Continent <- countrycode(metadata.df$Country, origin = "country.name", destination = "continent")

#add Saint_Martin to Americas 

metadata.df$Continent <- ifelse(metadata.df$Country == "Saint_Martin", "Americas", metadata.df$Continent)

# add csv containing the column (which can be country, continent etc) and the value of interest (e.g. Australia, Asia etc)

if (!is.null(opt$location_local)) {
  location_of_interest <- read.csv(opt$location_local)
} else {
  cat("Input metadata file. Exiting now...")
  quit()
}

# Get the column name and value of interest
column_of_interest <- colnames(location_of_interest)[1]
values_of_interest <- location_of_interest[, 1]

#add time interval

time_interval <- opt$time_interval

#add serotype

Serotype <- unique(metadata.df$Genotype)

# Filter metadata.df where the dynamic column equals the value of interest
metadata_countries_of_interest <- metadata.df %>%
  filter(!!sym(column_of_interest) %in% values_of_interest)

##############################################################################
# Step 3: Select sample from Location of Interest using weighted random sampling
##############################################################################

#add sampling weights 1/sample for each letter at each time

metadata_countries_of_interest <- metadata_countries_of_interest %>%
  group_by(!!sym(time_interval)) %>%
  mutate(weights = 1/n())

#plot histogram with ggplot this will help us understand if the sample has worked

metadata_countries_of_interest_plot <- metadata_countries_of_interest %>%
  group_by(!!sym(time_interval)) %>%
  summarise(count = n()) %>%
  ungroup()

plot_1 <- ggplot(metadata_countries_of_interest_plot, aes(x = !!sym(time_interval), y = count)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Year", y = "Number of Sequences", title = paste0('Number of Sequences per Year in ', paste0(values_of_interest, collapse = ";"))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set1")

#save as pdf

ggsave(filename = paste0(opt$outfile_plots,"metadata_countries_of_interest_plot_", Serotype, ".pdf"), plot = plot_1)

#conduct random weighted sampling

if (!is.null(opt$number_sequences_local)) {
  if (opt$number_sequences_local == "all") {
    number_sequences_local <- nrow(metadata_countries_of_interest)
  } else {
    number_sequences_local <- opt[["number_sequences_local"]]
  }
} else {
  number_sequences_local <- nrow(metadata_countries_of_interest)
}

df_sample <- sample_int_crank(n = nrow(metadata_countries_of_interest),number_sequences_local,
                              metadata_countries_of_interest$weights)

df_sample_location_of_interest <- metadata_countries_of_interest[df_sample,]

df_sample_location_of_interest_plot <- df_sample_location_of_interest %>%
  group_by(!!sym(time_interval)) %>%
  summarise(count = n()) %>%
  ungroup()

#plot data with good color scheme

plot_2 <- ggplot(df_sample_location_of_interest_plot, aes(x = !!sym(time_interval), y = count, fill)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Year", y = "Number of Sequences", title = paste0('Number of Sequences per Year in ', paste0(values_of_interest, collapse = ";"))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set1")

#save as pdf

ggsave(filename = paste0(opt$outfile_plots,"df_sample_location_of_interest_plot_", Serotype,".pdf"), plot = plot_2)

##############################################################################
# Step 4: Select even sample from background using weighted random sampling
##############################################################################

if (opt$sampling_method == 'Even') {
  
  if (exists("proportional_variable")) {
    # 'proportional_variable' exists, proceed with filtering
    background_metadata_even <- metadata.df %>%
      filter(!(!!sym(column_of_interest) %in% values_of_interest)) %>%
      filter(!!sym(column_of_interest) %in% unique(proportional_variable$Location))
  } else {
    # 'proportional_variable' does not exist, skip the second filter
    background_metadata_even <- metadata.df %>%
      filter(!(!!sym(column_of_interest) %in% values_of_interest))
  }
  
  #add sampling weights 1/sample for each letter at each time
  
  background_metadata_even <- background_metadata_even %>%
    group_by(!!sym(time_interval),!!sym(column_of_interest)) %>%
    mutate(weights = 1/n())
  
  #plot histogram with ggplot
  
  background_metadata_plot_even <- background_metadata_even %>%
    group_by(!!sym(time_interval),!!sym(column_of_interest)) %>%
    summarise(count = n()) %>%
    ungroup()
  
  #plot data with good color scheme
  
  plot_3 <- ggplot(background_metadata_plot_even, aes(x = !!sym(time_interval), y = count, fill = !!sym(column_of_interest))) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "Year", y = "Number of Sequences", title = "Number of Sequences per Year, Background Dataset") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 5),  # Reduce text size
          legend.key.size = unit(0.2, "cm"),  # Reduce key size
          legend.spacing.y = unit(0.1, "cm"),  # Reduce spacing between legend keys
          legend.margin = margin(2, 2, 2, 2)  # Reduce margin around the legend
    ) +
    guides(fill = guide_legend(ncol = 3, byrow = FALSE))  # Organize legend items in rows
  
  #save as pdf
  
  ggsave(filename = paste0(opt$outfile_plots,"background_metadata_plot_even_", Serotype,".pdf"), plot = plot_3)
  
  
  #conduct random weighted sampling
  
  # -------------------------------
  # Handling number_sequences_background
  # -------------------------------
  
  # Ensure the option is provided
  if (!is.null(opt$number_sequences_background)) {
    
    # Attempt to convert to numeric
    number_sequences_background <- suppressWarnings(as.numeric(opt$number_sequences_background))
    
    # Check if conversion was successful
    if (is.na(number_sequences_background)) {
      stop("Error: 'number_sequences_background' must be a numeric value representing either a fixed number or a percentage (as a decimal).")
    }
    
    # Determine if the value is a percentage or a fixed number
    if (number_sequences_background > 0 && number_sequences_background < 1) {
      # Interpret as a percentage of the total background sequences
      number_sequences_background_calculated <- floor(number_sequences_background * nrow(background_metadata_even))
      
      # Ensure at least one sequence is selected
      if (number_sequences_background_calculated < 1) {
        number_sequences_background_calculated <- 1
        warning("Calculated number of background sequences is less than 1. Setting to 1.")
      }
      
      message(sprintf("Sampling %.2f%% of background sequences (%d sequences).",
                      number_sequences_background * 100,
                      number_sequences_background_calculated))
      
    } else if (number_sequences_background >= 1) {
      # Interpret as a fixed number of sequences
      number_sequences_background_calculated <- floor(number_sequences_background)
      
      # Ensure the requested number does not exceed available sequences
      if (number_sequences_background_calculated > nrow(background_metadata_even)) {
        warning(sprintf("Requested number of background sequences (%d) exceeds available sequences (%d). Using all available sequences.",
                        number_sequences_background_calculated,
                        nrow(background_metadata_even)))
        number_sequences_background_calculated <- nrow(background_metadata_even)
      }
      
      message(sprintf("Sampling a fixed number of background sequences: %d sequences.", 
                      number_sequences_background_calculated))
      
    } else {
      stop("Error: 'number_sequences_background' must be a positive number.")
    }
    
  } else {
    # Default behavior: sample 10% of background sequences
    number_sequences_background_calculated <- floor(0.1 * nrow(background_metadata_even))
    
    # Ensure at least one sequence is selected
    if (number_sequences_background_calculated < 1) {
      number_sequences_background_calculated <- 1
      warning("Default sampling resulted in less than 1 background sequence. Setting to 1.")
    }
    
    message(sprintf("Default sampling: 10%% of background sequences (%d sequences).",
                    number_sequences_background_calculated))
  }
  
  # -------------------------------
  # Handling number_sequences_background
  # -------------------------------
  
  # Ensure the option is provided
  if (!is.null(opt$number_sequences_background)) {
    
    # Attempt to convert to numeric
    number_sequences_background <- suppressWarnings(as.numeric(opt$number_sequences_background))
    
    # Check if conversion was successful
    if (is.na(number_sequences_background)) {
      stop("Error: 'number_sequences_background' must be a numeric value representing either a fixed number or a percentage (as a decimal).")
    }
    
    # Determine if the value is a percentage or a fixed number
    if (number_sequences_background > 0 && number_sequences_background < 1) {
      # Interpret as a percentage of the total background sequences
      number_sequences_background_calculated <- floor(number_sequences_background * nrow(background_metadata_even))
      
      # Ensure at least one sequence is selected
      if (number_sequences_background_calculated < 1) {
        number_sequences_background_calculated <- 1
        warning("Calculated number of background sequences is less than 1. Setting to 1.")
      }
      
      message(sprintf("Sampling %.2f%% of background sequences (%d sequences).",
                      number_sequences_background * 100,
                      number_sequences_background_calculated))
      
    } else if (number_sequences_background >= 1) {
      # Interpret as a fixed number of sequences
      number_sequences_background_calculated <- floor(number_sequences_background)
      
      # Ensure the requested number does not exceed available sequences
      if (number_sequences_background_calculated > nrow(background_metadata_even)) {
        warning(sprintf("Requested number of background sequences (%d) exceeds available sequences (%d). Using all available sequences.",
                        number_sequences_background_calculated,
                        nrow(background_metadata_even)))
        number_sequences_background_calculated <- nrow(background_metadata_even)
      }
      
      message(sprintf("Sampling a fixed number of background sequences: %d sequences.", 
                      number_sequences_background_calculated))
      
    } else {
      stop("Error: 'number_sequences_background' must be a positive number.")
    }
    
  } else {
    # Default behavior: sample 10% of background sequences
    number_sequences_background_calculated <- floor(0.1 * nrow(background_metadata_even))
    
    # Ensure at least one sequence is selected
    if (number_sequences_background_calculated < 1) {
      number_sequences_background_calculated <- 1
      warning("Default sampling resulted in less than 1 background sequence. Setting to 1.")
    }
    
    message(sprintf("Default sampling: 10%% of background sequences (%d sequences).",
                    number_sequences_background_calculated))
  }
  
  # -------------------------------
  # Sampling Background Sequences
  # -------------------------------
  
  # Perform the sampling using the calculated number of sequences
  df_sample <- sample_int_crank(
    nrow(background_metadata_even),             # Total number of background sequences
    number_sequences_background_calculated,     # Number of sequences to sample
    background_metadata_even$weights            # Weights for sampling
  )
  
  # Subset the background metadata to obtain the sampled sequences
  background_sample_even <- background_metadata_even[df_sample, ]
  
  # -------------------------------
  # Preparing Data for Plotting
  # -------------------------------
  
  background_sample_plot_even <- background_sample_even %>%
    group_by(!!sym(time_interval), !!sym(column_of_interest)) %>%
    summarise(count = n(), .groups = 'drop')  # Added .groups = 'drop' for clarity
  
  #plot data with good color scheme
  
  plot_4 <- ggplot(background_sample_plot_even, aes(x = Year, y = count,fill = !!sym(column_of_interest))) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "Year", y = "Number of Sequences", title = "Number of Sequences per Year, Sampled Background") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 5),  # Reduce text size
          legend.key.size = unit(0.2, "cm"),  # Reduce key size
          legend.spacing.y = unit(0.1, "cm"),  # Reduce spacing between legend keys
          legend.margin = margin(2, 2, 2, 2)  # Reduce margin around the legend
    ) +
    guides(fill = guide_legend(ncol = 3, byrow = FALSE))  # Organize legend items in rows
  
  #save as pdf
  
  ggsave(filename = paste0(opt$outfile_plots,"background_sample_plot_even_", Serotype,".pdf"), plot = plot_4)
  
} else {
  
  ##############################################################################
  #Step 5: Proportional sampling of sequences from countries of interest and background
  ##############################################################################
  background_metadata_proportional <- metadata.df %>%
    filter(!(!!sym(column_of_interest) %in% values_of_interest))  %>%
    filter(!!sym(column_of_interest) %in% unique(proportional_variable$Location))
  
  #add sampling weights at each time
  join_conditions <- setNames(c("Date", "Location"), c(time_interval, column_of_interest))
  
  background_metadata_proportional_with_variable <- left_join(background_metadata_proportional,proportional_variable,
                                                              by = join_conditions)
  
  # List column names of proportional_variable
  
  background_metadata_proportional_with_variable <- background_metadata_proportional_with_variable %>%
    group_by(!!sym(time_interval),!!sym(column_of_interest)) %>%
    mutate(weights = Variable)
  
  #plot histogram with ggplot
  
  background_metadata_proportional_with_variable_plot <- background_metadata_proportional_with_variable %>%
    group_by(!!sym(time_interval),!!sym(column_of_interest)) %>%
    summarise(count = n()) %>%
    ungroup()
  
  #plot data with good color scheme
  
  plot_5 <- ggplot(background_metadata_proportional_with_variable_plot, aes(x = !!sym(time_interval), y = count, fill = !!sym(column_of_interest))) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "Year", y = "Number of Sequences", title = "Number of Sequences per Year") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 5),  # Reduce text size
          legend.key.size = unit(0.2, "cm"),  # Reduce key size
          legend.spacing.y = unit(0.1, "cm"),  # Reduce spacing between legend keys
          legend.margin = margin(2, 2, 2, 2)  # Reduce margin around the legend
    ) +
    guides(fill = guide_legend(ncol = 3, byrow = FALSE))  # Organize legend items in rows
  
  #save as pdf
  
  ggsave(filename = paste0(opt$outfile_plots,"background_metadata_plot_proportional_", Serotype, ".pdf"), plot = plot_5)
  
  df_sample <- sample_int_crank(n = nrow(background_metadata_proportional_with_variable),opt[["number_sequences_background"]],
                                background_metadata_proportional_with_variable$weights)
  
  background_sample_proportional <- background_metadata_proportional_with_variable[df_sample,]
  
  background_sample_plot_proportional <- background_sample_proportional %>%
    group_by(!!sym(time_interval),!!sym(column_of_interest)) %>%
    summarise(count = n()) %>%
    ungroup()
  
  #plot data with good color scheme
  
  plot_6 <- ggplot(background_sample_plot_proportional, aes(x = !!sym(time_interval), y = count,fill = !!sym(column_of_interest))) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "Year", y = "Number of Sequences", title = "Number of Sequences per Year") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 5),  # Reduce text size
          legend.key.size = unit(0.2, "cm"),  # Reduce key size
          legend.spacing.y = unit(0.1, "cm"),  # Reduce spacing between legend keys
          legend.margin = margin(2, 2, 2, 2)  # Reduce margin around the legend
    ) +
    guides(fill = guide_legend(ncol = 3, byrow = FALSE))  # Organize legend items in rows
  
  #save as pdf
  
  ggsave(filename = paste0(opt$outfile_plots,"background_sample_plot_proportional_", Serotype, ".pdf"), plot = plot_6)
  
}

##############################################################################
#Step 6: Merge the local and background datasets
##############################################################################

if (opt$sampling_method == 'Even') {
  final_sample <- rbind(df_sample_location_of_interest, background_sample_even)
} else {
  final_sample <- rbind(df_sample_location_of_interest, background_sample_proportional)
}

##############################################################################
#Step 7: Save the final sample
##############################################################################


seq_name <- as.data.frame(as.matrix(attributes(seqs)$names))

keep <- seq_name %>%
  filter(V1 %in% final_sample$Sequence_name)

taxa_split <- data.frame(do.call('rbind',strsplit(as.character(keep$V1),'|',fixed = TRUE)))
taxa_split$name <- keep$V1
species.to.keep <- taxa_split$name

vec.names<-unlist(lapply(strsplit(names(seqs), ";"), function(x)x[length(x)]))

vec.tokeep <-which(vec.names %in%  species.to.keep)

# Write the sequences to a new FASTA file

write.fasta(sequences=seqs[vec.tokeep], names=names(seqs)[vec.tokeep],
            file.out=paste0(opt$outfile_fasta_csv, '.fasta'))

write.csv(final_sample,
          file = paste0(opt$outfile_fasta_csv, "_infoTbl.csv"),
          row.names = FALSE, quote=FALSE)
