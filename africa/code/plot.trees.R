lapply(c("ggplot2", "ggnewscale", "ggtree", "ggtreeExtra", "treeio", "tidytree", "dplyr", "phangorn", 
         "tidyverse", "ape", "seqinr", "dispRity", "rstudioapi", "grDevices", "ggstar"), library, character.only = TRUE)

#Dengue 1

#Load in MCC

tree <- read.nexus("/Users/rhysinward/Documents/grapevne_modules/Beast_Vietnam/visulisation/dta_dengue_1.tre")
genotype <- read.csv("/Users/rhysinward/Documents/grapevne_modules/dengue_pipeline_vietnam/results/nextclade_dta_dengue_1.csv",sep =";")
DTA <- read.csv("/Users/rhysinward/Documents/grapevne_modules/dengue_pipeline_vietnam/results/dta_metadata_dengue_1.txt", sep = "\t")

#Get metadata from tip name

metadata_df <- as.data.frame(tree[["tip.label"]])

genotype <- dplyr :: select(genotype,c(2,3))

metadata_df <- filter(genotype, seqName %in% metadata_df$`tree[["tip.label"]]`)

metadata_df <- left_join(metadata_df, DTA, by = c("seqName" = "traits"))

#merge china into 1 for plotting

chinese_provinces <- c("Guangdong", "Sichuan", "Yunnan", "Hainan", "Henan", "Jiangxi", "Fujian", "Zhejiang")
metadata_df$location <- ifelse(metadata_df$location %in% chinese_provinces, "China", metadata_df$location)


for (i in 1:nrow(metadata_df)){
  if (metadata_df$location[i] == 'China'){
    metadata_df$location[i] <- 'Border'
  } else if (metadata_df$location[i] == 'Laos'){
    metadata_df$location[i] <- 'Border'
  } else if (metadata_df$location[i] == 'Cambodia'){
    metadata_df$location[i] <- 'Border'
  } else if (metadata_df$location[i] == 'Central_Vietnam'){
    metadata_df$location[i] <- 'Central Vietnam'
  } else if (metadata_df$location[i] == 'Northern_Vietnam'){
    metadata_df$location[i] <- 'Northern Vietnam'
  } else if (metadata_df$location[i] == 'Southern_Vietnam'){
    metadata_df$location[i] <- 'Southern Vietnam'
  } else {
    metadata_df$location[i] <- 'Other Southern and SEA'
  }
}

#plot tree
p <- ggtree(tree,size=1, mrsd="2023-08-01") %<+% metadata_df

  p1 <- p +
  geom_tippoint(aes(color=location), size=3, alpha=.75) +
  scale_color_manual(values = c(
    '#C2AFF0',  '#C97064',
    '#FFD700', '#FF6347',  '#7FFFD4'
  )) +  theme(legend.position="right") +
  geom_text2(aes(subset = !isTip, label=label), # subset=!isTip
             size = 10,
             color = "black",
             hjust = 1, 
             vjust = -1.5
  ) + theme_tree2()

  p2 <-p1 +
    new_scale_fill() +
    geom_fruit(
      geom=geom_star,
      mapping=aes(fill=clade),
      starshape=26,
      color=NA,
      size=3,
      starstroke=0,
      offset=.02,
    ) +
    scale_fill_manual(values = c(
      '#C2AFF0',  '#C97064',
      '#FFD700'),
      name="Genotype",
      guide=guide_legend(keywidth=0.3, keyheight=0.3, ncol=2, order=2)
    )
  
  p2


p1 <- p1 + geom_text2(aes(subset = !isTip, label=label))


p <- p1 + geom_hilight(mapping = aes(subset = (genotype == "DENV1/I")), fill = "red", alpha = 0.5)

#Dengue 2

#Load in MCC

tree <- read.nexus("/Users/rhysinward/Documents/grapevne_modules/Beast_Vietnam/visulisation/denv2_pruned.tre")
genotype <- read.csv("/Users/rhysinward/Documents/grapevne_modules/dengue_pipeline_vietnam/results/nextclade_dta_dengue_2.csv",sep =";")
DTA <- read.csv("/Users/rhysinward/Documents/grapevne_modules/dengue_pipeline_vietnam/results/dta_metadata_dengue_2.txt", sep = "\t")

#Get metadata from tip name

tree[["tip.label"]] <- gsub("'", "",tree[["tip.label"]])


metadata_df <- as.data.frame(tree[["tip.label"]])

#remove ' from tree name 

metadata_df$`tree[["tip.label"]]` <- gsub("'", "", metadata_df$`tree[["tip.label"]]`)

genotype <- dplyr :: select(genotype,c(2,3))

metadata_df <- filter(genotype, seqName %in% metadata_df$`tree[["tip.label"]]`)

metadata_df <- left_join(metadata_df, DTA, by = c("seqName" = "traits"))

#merge china into 1 for plotting

chinese_provinces <- c("Guangdong", "Sichuan", "Yunnan", "Hainan", "Henan", "Jiangxi", "Fujian", "Zhejiang")
metadata_df$location <- ifelse(metadata_df$location %in% chinese_provinces, "China", metadata_df$location)


for (i in 1:nrow(metadata_df)){
  if (metadata_df$location[i] == 'China'){
    metadata_df$location[i] <- 'Border'
  } else if (metadata_df$location[i] == 'Laos'){
    metadata_df$location[i] <- 'Border'
  } else if (metadata_df$location[i] == 'Cambodia'){
    metadata_df$location[i] <- 'Border'
  } else if (metadata_df$location[i] == 'Central_Vietnam'){
    metadata_df$location[i] <- 'Central Vietnam'
  } else if (metadata_df$location[i] == 'Northern_Vietnam'){
    metadata_df$location[i] <- 'Northern Vietnam'
  } else if (metadata_df$location[i] == 'Southern_Vietnam'){
    metadata_df$location[i] <- 'Southern Vietnam'
  } else {
    metadata_df$location[i] <- 'Other Southern and SEA'
  }
}

#plot tree
p <- ggtree(tree,size=1, mrsd="2023-08-03") %<+% metadata_df

p1 <- p +
  geom_tippoint(aes(color=location), size=3, alpha=.75) +
  scale_color_manual(values = c(
    '#C2AFF0',  '#C97064',
    '#FFD700', '#FF6347',  '#7FFFD4'
  )) +  theme(legend.position="right") +
  geom_text2(aes(subset = !isTip, label=label), # subset=!isTip
             size = 10,
             color = "black",
             hjust = 1, 
             vjust = -1.5
  ) + theme_tree2()

p2 <-p1 +
  new_scale_fill() +
  geom_fruit(
    geom=geom_star,
    mapping=aes(fill=clade),
    starshape=26,
    color=NA,
    size=3,
    starstroke=0,
    offset=.02,
  ) +
  scale_fill_manual(values = c(
    '#C2AFF0',  '#C97064',
    '#FFD700','#7FFFD4'),
    name="Genotype",
    guide=guide_legend(keywidth=0.3, keyheight=0.3, ncol=2, order=2)
  )

p2
