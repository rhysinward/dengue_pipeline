#write tsv containing latlongs for all my countries in the DTA

# Define the data frame with countries and their latitude and longitude
countries <- data.frame(
    Location = c("Cambodia", "Central_Vietnam", "Fujian", "Guangdong", 
              "Hainan", "Henan", "India", "Indonesia", 
              "Jiangxi", "Laos","Malaysia", "Myanmar", "Philippines",
              "Northern_Vietnam", "Singapore", "Southern_Vietnam","Sichuan",
              "Taiwan", "Thailand", "Yunnan", "Zhejiang","Shanghai","Shandong"),
  Latitude = c(11.5564, 16.4637, 26.0998, 23.1317, 20.0200, 
               34.7657, 28.6139, -6.1944, 28.6742, 
               17.9757, 3.1319, 16.8409, 14.5995,21.0278,
               1.3521, 10.8231,30.6509,23.6978,13.7563,
               25.0453,30.2655,31.2304,36.6683),
  Longitude = c(104.9282, 107.5909, 119.2966, 113.2663, 110.3486, 
                113.7532, 77.2090, 106.8229, 115.9100, 
                102.6331, 101.6841, 96.1735, 120.9842,105.8342,
                103.8198, 106.6297,104.0757,120.9605,100.5018,
                102.7097,120.1536,121.4737,117.0208)
)

# Save the data frame to a TSV file
write.table(countries, file = "countries_sea_lat_long.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
