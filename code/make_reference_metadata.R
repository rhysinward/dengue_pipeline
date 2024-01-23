#dengue 1
df <- data.frame(
  name = "EF457905.1_Dengue_virus_type1_isolate_P72-1244_complete_genome",
  date = "1972-06-15",
  GenBank_ID = "EF457905.1", 
  Country = "Malaysia", 
  State = NA,
  City = NA, 
  Serotype = NA, 
  Decimal_Date = NA,
  stringsAsFactors = FALSE 
)

#export to csv

write.csv(df, "sylvatic_dengue1.csv", row.names = FALSE)

#dengue 2

df <- data.frame(
  name = "EU003591.1_Dengue_virus_type_2_isolate_IBH11234_polyprotein_gene_complete_cds",
  date = "1966-08-18",
  GenBank_ID = "EU003591.1", 
  Country = "Nigeria", 
  State = NA,
  City = NA, 
  Serotype = NA, 
  Decimal_Date = NA,
  stringsAsFactors = FALSE 
)

#export to csv

write.csv(df, "aligned_root_denv2.csv", row.names = FALSE)

#dengue 3

df <- data.frame(
  name = "KU050695.1_Dengue_virus_3_complete_genome",
  date = "1956-06-15",
  GenBank_ID = "KU050695.1", 
  Country = "Philippines", 
  State = NA,
  City = NA, 
  Serotype = NA, 
  Decimal_Date = NA,
  stringsAsFactors = FALSE 
)

#export to csv

write.csv(df, "aligned_root_dengue_3.csv", row.names = FALSE)

#dengue 4

df <- data.frame(
  name = "JF262780.1_Dengue_virus_4_isolate_P73-1120,_complete_genome",
  date = "1973-06-15",
  GenBank_ID = "JF262780.1", 
  Country = "Malaysia", 
  State = NA,
  City = NA, 
  Serotype = NA, 
  Decimal_Date = NA,
  stringsAsFactors = FALSE 
)

#export to csv

write.csv(df, "sylvatic_dengue4.csv", row.names = FALSE)

