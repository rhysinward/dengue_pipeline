#!/usr/bin/env Rscript

# -----------------------------
# plot_viral_movements.R
# -----------------------------

# **1. Load Necessary Packages**

# Install and load required packages
# Note: For reproducibility and better dependency management, it's recommended to handle package installations via Conda environments in Snakemake.

suppressPackageStartupMessages({
  library(optparse)        # For command-line argument parsing
  library(ggplot2)         # For plotting
  library(dplyr)           # For data manipulation
  library(grid)            # For graphical functions like unit()
  library(countrycode)     # For country-to-continent conversion
  library(rnaturalearth)   # For natural earth map data
  library(rnaturalearthdata)
  library(sf)              # For handling spatial data
  library(readr)           # For reading CSV files
  library(scales)          # For scaling functions in plots
})

# **2. Define Command-Line Arguments**

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to input CSV file containing annotated tree events", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Path to output PNG file for the viral movements map", metavar = "character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# **3. Validate Command-Line Arguments**

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Error: Input file must be specified (--input).", call. = FALSE)
}

if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Error: Output file must be specified (--output).", call. = FALSE)
}

# **4. Define Decimal Date Conversion Function**

decimal2Date <- function(decimal_date) {
  year <- floor(decimal_date)
  remainder <- decimal_date - year
  # Determine if the year is a leap year
  is_leap <- (year %% 4 == 0 & year %% 100 != 0) | (year %% 400 == 0)
  days_in_year <- ifelse(is_leap, 366, 365)
  # Convert decimal to days
  date <- as.Date(paste0(year, "-01-01")) + round(remainder * days_in_year)
  return(date)
}

# **5. Load and Process Data**

# Read the CSV data and replace underscores with spaces in 'Origin' and 'Destination'
DENV_1IIIA_introductions <- read_csv(opt$input) %>%
  mutate(
    Origin = gsub("_", " ", Origin),
    Destination = gsub("_", " ", Destination)
  )

# **6. Obtain World Country Centroids and Continent Information**

# Retrieve country data with centroids and continent information
world_centroids <- ne_countries(scale = "medium", returnclass = "sf") %>%
  # Select relevant columns
  select(name = name_long, continent, geometry) %>%
  # Compute centroids
  st_centroid() %>%
  # Extract latitude and longitude from centroids
  mutate(
    latitude = st_coordinates(geometry)[,2],
    longitude = st_coordinates(geometry)[,1]
  ) %>%
  # Remove geometry as it's no longer needed
  st_set_geometry(NULL)

# **7. Standardize Country Names in the Dataset**

# Function to manually correct unmatched country names
manual_country_corrections <- function(name_vector){
  # Define manual mappings: names as they appear in the dataset to standardized names
  corrections <- c(
    "Dominican Rep." = "Dominican Republic",
    "USA" = "United States",
    "UK" = "United Kingdom",
    "Côte d'Ivoire" = "Ivory Coast",
    "DRC" = "Democratic Republic of the Congo",
    "S. Korea" = "South Korea",
    "N. Korea" = "North Korea",
    "Laos" = "Lao People's Democratic Republic",
    "Russia" = "Russian Federation",
    "Viet Nam" = "Vietnam",
    "Czech Republic" = "Czechia",
    "Bolivia" = "Bolivia (Plurinational State of)",
    "Moldova" = "Republic of Moldova"
  )
  
  # Apply corrections
  corrected_names <- sapply(name_vector, function(x){
    if(x %in% names(corrections)){
      return(corrections[[x]])
    } else {
      return(x)
    }
  }, USE.NAMES = FALSE)
  
  return(corrected_names)
}

# Apply manual corrections to 'Origin' and 'Destination' columns
DENV_1IIIA_introductions <- DENV_1IIIA_introductions %>%
  mutate(
    Origin_Standardized = manual_country_corrections(Origin),
    Destination_Standardized = manual_country_corrections(Destination)
  )

# **8. Define Mapping Functions**

country2continent <- function(country_name, data = world_centroids){
  continent <- data$continent[match(country_name, data$name)]
  return(as.character(continent))
}

country2lat <- function(country_name, data = world_centroids){
  lat <- data$latitude[match(country_name, data$name)]
  return(as.numeric(lat))
}

country2long <- function(country_name, data = world_centroids){
  lon <- data$longitude[match(country_name, data$name)]
  return(as.numeric(lon))
}

# **9. Map Countries to Continents, Latitudes, and Longitudes**

DENV_1IIIA_introductions <- DENV_1IIIA_introductions %>%
  mutate(
    # Convert decimal date to standard date
    date = decimal2Date(EventTime),
    
    # Map standardized country names to continents
    Origin_Continent = sapply(Origin_Standardized, country2continent),
    Destination_Continent = sapply(Destination_Standardized, country2continent),
    
    # Map standardized country names to latitudes and longitudes
    origin_lat = sapply(Origin_Standardized, country2lat),
    origin_long = sapply(Origin_Standardized, country2long),
    destination_lat = sapply(Destination_Standardized, country2lat),
    destination_long = sapply(Destination_Standardized, country2long)
  )

# **10. Handle Missing or Unmatched Country Names**

# Identify rows with any NA values in crucial mapping fields
missing_coords <- DENV_1IIIA_introductions %>%
  filter(
    is.na(origin_lat) | is.na(origin_long) |
      is.na(destination_lat) | is.na(destination_long)
  )

# Print unmatched countries for inspection and remove them
if(nrow(missing_coords) > 0){
  message("Warning: The following rows have unmatched country names and will be removed:")
  print(missing_coords %>% select(Origin, Destination))
  
  # Remove rows with any NA coordinates
  DENV_1IIIA_introductions <- DENV_1IIIA_introductions %>%
    filter(
      !is.na(origin_lat) & !is.na(origin_long) &
        !is.na(destination_lat) & !is.na(destination_long)
    )
}

# **11. Remove Identical Origin-Destination Pairs**

# Identify and remove rows where origin and destination coordinates are identical
identical_pairs <- DENV_1IIIA_introductions %>%
  filter(
    origin_lat == destination_lat & origin_long == destination_long
  )

if(nrow(identical_pairs) > 0){
  message("Info: Removing rows with identical origin and destination coordinates.")
  DENV_1IIIA_introductions <- DENV_1IIIA_introductions %>%
    filter(!(origin_lat == destination_lat & origin_long == destination_long))
}

# **12. Aggregate Data for Movements**

# Aggregate data to get total movements between origin-destination pairs
DENV_1IIIA_aggregated <- DENV_1IIIA_introductions %>%
  group_by(origin_lat, origin_long, destination_lat, destination_long) %>%
  summarize(total_movements = n(), .groups = 'drop')

# Join aggregated data back to the main dataset for plotting
DENV_1IIIA_plot_data <- DENV_1IIIA_introductions %>%
  left_join(DENV_1IIIA_aggregated, 
            by = c("origin_lat", "origin_long", "destination_lat", "destination_long"))

# **13. Prepare Map Data for South America**

# Define South American countries for filtering
south_american_countries <- c(
  "Argentina", "Bolivia", "Brazil", "Chile", "Colombia", 
  "Ecuador", "Guyana", "Paraguay", "Peru", "Suriname", 
  "Uruguay", "Venezuela", "French Guiana"
)

# Retrieve South American map data using rnaturalearth
map <- ne_countries(scale = "medium", returnclass = "sf") 

# **14. Final Data Verification**

# Verify no NA values remain in plotting data
if(any(is.na(DENV_1IIIA_plot_data$origin_long) |
       is.na(DENV_1IIIA_plot_data$origin_lat) |
       is.na(DENV_1IIIA_plot_data$destination_long) |
       is.na(DENV_1IIIA_plot_data$destination_lat))){
  stop("Error: NA values detected in plotting coordinates after cleaning.")
}

# **15. Plotting the Map with Viral Movements**

viral_movements_map <- ggplot() +
  # Add South America map layer
  geom_sf(data = map, fill = 'thistle4', color = "white", alpha = 0.7) +
  
  # Add viral movement curves
  geom_curve(
    data = DENV_1IIIA_plot_data,
    aes(
      x = origin_long,
      y = origin_lat,
      xend = destination_long,
      yend = destination_lat,
      color = EventTime,             # Color based on EventTime (converted to date)
      size = total_movements    # Line width based on total movements
    ),
    alpha = 0.8,
    curvature = 0.2,
    arrow = arrow(length = unit(0.2, "cm"), type = "closed")
  ) +
  
  # Define color scale for dates
  scale_color_gradientn(
    colors = c("#003049", "deepskyblue4", "skyblue3", "darkseagreen3", "lightyellow2"),
    name = 'Date of Viral Movement',
    labels = date_format("%Y-%m-%d"),
    breaks = pretty_breaks(n = 5)
  ) +
  
  # Define size scale for total movements
  scale_size(range = c(0.5, 3), name = "Number of Movements") +
  
  # Customize legends
  guides(
    color = guide_colourbar(
      barwidth = 20, 
      barheight = 0.5,
      title.position = 'top', 
      title.hjust = 0.5
    ),
    size = guide_legend(
      title.position = 'top',
      title.hjust = 0.5
    )
  ) +
  
  # Apply theme settings
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = 'top',
    legend.direction = 'horizontal',
    plot.title = element_text(family = "Helvetica", size = 16, face = "bold", margin = margin(b = 20)),
    legend.box = 'vertical',
    legend.title.align = 0.5
  ) 

# **16. Display the Map**

print(viral_movements_map)

# **17. Save the Map**

ggsave(filename = opt$output, plot = viral_movements_map, width = 12, height = 8, dpi = 300)

# -----------------------------
# End of Script
# -----------------------------
