#!/usr/bin/env Rscript

# Function to check and install packages
check_and_install_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
    }
  }
}

# Check and install necessary packages
check_and_install_packages(c("argparse", "dplyr", "ggplot2"))

# Load the packages
library(argparse)
library(dplyr)
library(ggplot2)

# Argument parsing
parser <- ArgumentParser()
parser$add_argument("-i", "--infile", type="character", help="Input file (output of fertree)")
parser$add_argument("-s", "--source", type="character", help="Name of source location (case-sensitive), for identifying local transmission-lineages")
parser$add_argument("-m", "--max_ntaxa", type="numeric", help="Maximum lineage size to bin lineages by (all lineages of size greater than max_ntaxa are aggregated)")
parser$add_argument("-t", "--title", type="character", default="Lineage size distribution", help="Title of output figure")
parser$add_argument("-x", "--width", type="numeric", default=7.5, help="Width of output figure")
parser$add_argument("-y", "--height", type="numeric", default=4.5, help="Height of output figure")
parser$add_argument("-p", "--pdf", action="store_true", default=FALSE, help="Save plot as PDF, otherwise PNG (default)")
parser$add_argument("-o", "--outfile", type="character", help="Output figure")
args <- parser$parse_args()

# Read in output from fertree transmission-lineages
fertree_out.df <- read.csv(args$infile, sep='\t')
# Take only lineages 
fertree_out.df <- fertree_out.df %>%
  filter(source == args$source & ntaxa > 0)

# Calculate total number of genomes
total_seqs <- sum(fertree_out.df$ntaxa)

# Specify max_ntaxa
max_ntaxa <- min(args$max_ntaxa, max(fertree_out.df$ntaxa))

# Bin by lineage size and calculate cumulative number of genomes across bins
fertree_out.binned.df <- fertree_out.df %>%
  mutate(capped_ntaxa=ifelse(ntaxa > max_ntaxa, max_ntaxa+1, ntaxa)) %>%
  group_by(capped_ntaxa) %>%
  summarise(count=n(), total_seqs_count=sum(ntaxa)) %>%
  mutate(proportion=total_seqs_count/total_seqs)
# Make sure all consecutive lineage sizes are there
all_consecutive_ntaxa <- seq(min(fertree_out.binned.df$capped_ntaxa), max(fertree_out.binned.df$capped_ntaxa))
empty.binned.df <- data.frame(
  capped_ntaxa=all_consecutive_ntaxa[!all_consecutive_ntaxa %in% fertree_out.binned.df$capped_ntaxa],
  count=0, total_seqs_count=0, proportion=0
)
fertree_out.binned.df <- rbind(fertree_out.binned.df, empty.binned.df) %>%
  arrange(capped_ntaxa) %>%
  mutate(cumulative_proportion=cumsum(proportion))

# Draw plot
vert_offset.p <- 0.05
sf <- 1/max(fertree_out.binned.df$count)
plot <- ggplot() +
  geom_line(dat=fertree_out.binned.df, aes(x=capped_ntaxa, y=cumulative_proportion/sf), color='#0C4C5F', linewidth=1, alpha=0.9) +
  geom_bar(dat=fertree_out.binned.df, aes(x=capped_ntaxa, y=count), stat='identity',
           fill=c(rep('#3E4855', nrow(fertree_out.binned.df)-1), 'lightgrey'),
           color=c(rep('black', nrow(fertree_out.binned.df)-1), 'darkgrey'), alpha=0.9, width=0.8) +
  scale_y_continuous(
    expand=c(0, 0),
    limits=c(0, max(fertree_out.binned.df$count) * (1+vert_offset.p)),
    sec.axis = sec_axis(trans=~.*sf, name='Cumulative prop. of total genomes')
  ) +
  labs(title=sprintf('%s (n=%d)', args$title, total_seqs), x='Lineage size', y='Count') +
  coord_cartesian(clip='off') +
  theme_bw()

# Save plot to file
ggsave(args$outfile, plot, device=ifelse(args$pdf, 'pdf', 'png'), width=args$width, height=args$height)
