# FIGURE 3
```{r}
# Import libraries 
library(tidyverse)
library(ggthemes)
library(sf)
library(scales)
library(patchwork)
library(dplyr)
library(ggplot2)
library(mclust)
library(sjmisc)
```
## FIGURE 3a
```{r}
# Import LINCodes and create strings of relevant groupings
lincodes <- read.csv("FINAL_DATA/short_LINcodes_cgST.csv", header=TRUE) %>% na.omit()
lincodes$gx <- paste(lincodes$L1, lincodes$L2, sep = "_")
lincodes$gy <- paste(lincodes$L1, lincodes$L2,lincodes$L3, sep = "_")
lincodes$gz <- paste(lincodes$L1, lincodes$L2,lincodes$L3,lincodes$L4, sep = "_")

# Import phylogeny
tree3 <- midpoint.root(read.newick("19k_tree.nwk"))

# Import list of samples
sample_list <- read.table("shortlist2.txt")

# Subset LINcode list to relevant ids
tx22 <- lincodes[lincodes$id %in% sample_list$V1, ]

# Find shared node for entire lineage (e.g. Linegae 15 shown here)
liny <- tx22 %>% subset(gy=="0_21_0") %>% select(id) 
zinny <- as.character(unlist(tail(liny$id,200)))
findMRCA(tree3, tips=zinny, type=c("node"))

# Plot
ggtree(tree3, layout="circular", size=0.1, color="grey30")  %<+% subset(lincodes) +  
  scale_color_manual(values = cols, na.value = "grey50") +
  geom_tippoint(aes(color=gy)) +
  theme(legend.position = "none") + geom_treescale(width = 0.001) +
  geom_cladelabel(node=3851, color='black', label = "PS", offset = 0.00015, barsize = 0.45)  +
  geom_cladelabel(node=3816, color='black', label = "OG", offset = 0.00015, barsize = 0.45)  +
  geom_cladelabel(node=2518, color='black', label = "0_58_0", offset = 0.00015, barsize = 0.45,hjust = 1) + geom_hilight(node =2518 ,alpha = .5, fill="#f44336") +
  geom_cladelabel(node=2927, color='black', label = "0_15_4", offset = 0.00015, barsize = 0.45,hjust = 1) + geom_hilight(node =2927 ,alpha = .5, fill="#e81e63") +
  geom_cladelabel(node=2871, color='black', label = "0_15_0", offset = 0.00015, barsize = 0.45,hjust = 1) + geom_hilight(node =2871 ,alpha = .5, fill="#9c27b0") +
  geom_cladelabel(node=2647, color='black', label = "0_0_0", offset = 0.00015, barsize = 0.45) + geom_hilight(node =2647,alpha = .5, fill="#673ab7") +
  geom_cladelabel(node=3439, color='black', label = "0_20_1", offset = 0.00015, barsize = 0.45) + geom_hilight(node =3439,alpha = .5, fill="#3f51b5") +
  geom_cladelabel(node=1968, color='black', label = '0_11_3', offset = 0.00015, barsize = 0.45) + geom_hilight(node = 1968,alpha = .5, fill="#2196f3") +
  geom_cladelabel(node=2496, color='black', label = '0_55_1', offset = 0.00015, barsize = 0.45) + geom_hilight(node = 2496,alpha = .5, fill="#03a9f4") +
  geom_cladelabel(node=2029, color='black', label = '0_11_8', offset = 0.00015, barsize = 0.45) + geom_hilight(node = 2029,alpha = .5, fill="#00bcd4") +
  geom_cladelabel(node=2193, color='black', label = "0_75_0", offset = 0.00015, barsize = 0.45) + geom_hilight(node =2193 ,alpha = .5, fill="#009688") +
  geom_cladelabel(node=2783, color='black', label = "0_57_0", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 2783,alpha = .5, fill="#4caf50") +
  geom_cladelabel(node=3752, color='black', label = "0_145_2", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 3752,alpha = .5, fill="#8bc34a") +
  geom_cladelabel(node=3079, color='black', label = "0_130_1", offset = 0.00015, barsize = 0.45) + geom_hilight(node =3079 ,alpha = .5, fill="#cddc39") +
  geom_cladelabel(node=3019, color='black', label = "0_10_0", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 3019,alpha = .5, fill="#ffeb3b") +
  geom_cladelabel(node=3619, color='black', label = "0_169_0", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 3619,alpha = .5, fill="#ffc107") +
  geom_cladelabel(node=2430, color='black', label = "0_232_1", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 2430,alpha = .5, fill="#ff9800") +
  geom_cladelabel(node=2147, color='black', label = "0_74_0", offset = 0.00015, barsize = 0.45) + geom_hilight(node =2147 ,alpha = .5, fill="#ff5722") +
  geom_cladelabel(node=3524, color='black', label = "0_21_0", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 3524,alpha = .5, fill="#795548") +
  geom_cladelabel(node=2386, color='black', label = "0_237_0", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 2386 ,alpha = .5, fill="#9e9e9e") +
  geom_cladelabel(node=3201, color='black', label = "0_265_0", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 3201,alpha = .5, fill="#E76BF3") +
  geom_cladelabel(node=3259, color='black', label = "0_310_0", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 3259,alpha = .5, fill="#607d8b") 
````
## Figure 3b
```{r}
# Import tree (from Figure 2A)
tree3 <- midpoint.root(read.newick("19k_tree.nwk"))

# Calculate pairwise cophentic distance 
coph1 <- (cophenetic(tree3) %>% reshape2::melt())

# Get allelic mismatch stats from MSTclust
mst_1950 <- read.csv("mst_cluster.19k.input.out.2.d.csv", header=TRUE, check.names = FALSE) %>% reshape2::melt(id.vars = c("ID")) %>% na.omit()

# Get LINcodes 
lincodes <- read.csv("FINAL_DATA/short_LINcodes_cgST.csv", header=TRUE) %>% na.omit()

# Merge dataframes
mrg1 <- merge(coph1, mst_1950, by.y=c("ID","variable"), by.x=c("Var1","Var2"), all.x=TRUE)
mrg2 <- merge(mrg1, mst_1950, by.y=c("variable","ID"), by.x=c("Var1","Var2"), all.x=TRUE) %>% subset(!is.na(value) | !is.na(value.y))

# Copy across columns to get consistency for pairwise comparisons
mrg2$mst3 <- rowSums(mrg2[,c("value", "value.y")], na.rm=TRUE)

# Merge LINcodes for each pair
mrg5 <- merge(lincodes, mrg2, by.y=c("Var1"),by.x=c("id"), all.y=TRUE)
mrg6 <- merge(lincodes, mrg5, by.y=c("Var2"),by.x=c("id"), all.y=TRUE)

# Create strings for each level of LINcode

mrg6$SUPERLINEAGE.x <- paste(mrg6$L1.x, mrg6$L2.x, sep = "_")
mrg6$SUPERLINEAGE.y <- paste(mrg6$L1.y, mrg6$L2.y, sep = "_")
mrg6$LINEAGE.x <- paste(mrg6$L1.x, mrg6$L2.x,mrg6$L3.x, sep = "_")
mrg6$LINEAGE.y <- paste(mrg6$L1.y, mrg6$L2.y,mrg6$L3.y, sep = "_")
mrg6$SUBLINEAGE.x <- paste(mrg6$L1.x, mrg6$L2.x, mrg6$L3.x, mrg6$L4.x, sep = "_")
mrg6$SUBLINEAGE.y <- paste(mrg6$L1.y, mrg6$L2.y,mrg6$L3.y,mrg6$L4.y, sep = "_")
mrg6$CLONALGROUP.x <- paste(mrg6$L1.x, mrg6$L2.x,mrg6$L3.x,mrg6$L4.x, mrg6$L5.x, sep = "_")
mrg6$CLONALGROUP.y <- paste(mrg6$L1.y, mrg6$L2.y,mrg6$L3.y,mrg6$L4.y, mrg6$L5.y, sep = "_")
mrg6$CLONALSUBGROUP.x <- paste(mrg6$L1.x, mrg6$L2.x,mrg6$L3.x,mrg6$L4.x, mrg6$L5.x,mrg6$L6.y, sep = "_")
mrg6$CLONALSUBGROUP.y <- paste(mrg6$L1.y, mrg6$L2.y,mrg6$L3.y,mrg6$L4.y, mrg6$L5.y,mrg6$L6.y, sep = "_")
mrg6$sp.x <- ifelse(grepl(pattern = "ps_", mrg6$id), "SSP","NULL")
mrg6$sp.y <- ifelse(grepl(pattern = "ps_", mrg6$id.y), "SSP","NULL")

# Replace the false NA_* strings
mrg6b <- mrg6  %>% sample_n(150000) 
mrg9 <- mrg6b %>%  rowwise() %>%  mutate(across(everything(), ~ifelse(str_contains(., "NA_"), NA, .)))

# Define relationships between pairs
mrg9$test <- ifelse(mrg9$SUPERLINEAGE.x!=mrg9$SUPERLINEAGE.y,"DIF SUPERLINEAGE",
                    ifelse(mrg9$sp.x!=mrg9$sp.y,"DIF SPECIES",
                           ifelse(mrg9$LINEAGE.x!=mrg9$LINEAGE.y, "DIF LINEAGE",
                                  ifelse(mrg9$SUBLINEAGE.x!=mrg9$SUBLINEAGE.y, "DIF SUBLINEAGE",
                                         ifelse(mrg9$CLONALGROUP.x!=mrg9$CLONALGROUP.y, "DIF CG",
                                                ifelse(mrg9$CLONALSUBGROUP.y==mrg9$CLONALSUBGROUP.y, "DIF CSG",
                                                       ifelse(mrg9$SUPERLINEAGE.x==mrg9$SUPERLINEAGE.y,"SAME SUPERLINEAGE",
                                                              ifelse(mrg9$LINEAGE.x==mrg9$LINEAGE.y,"SAME LINEAGE",
                                                                     ifelse(mrg9$SUBLINEAGE.x==mrg9$SUBLINEAGE.y, "SAME SUBLINEAGE",
                                                                            ifelse(mrg9$CLONALGROUP.x==mrg9$CLONALGROUP.y, "SAME CG",
                                                                                   ifelse(mrg9$CLONALSUBGROUP.x==mrg9$CLONALSUBGROUP.y, "SAME CSG","UNKNOWN","ELSE")))))))))))

# Subset and rename groupings
subset1 <- subset(mrg9, sp.x!=sp.y)
subset2 <- subset(mrg9, sp.x==sp.y & sp.x!="SSP") %>%  subset(!is.na(L2.x) & !is.na(L2.y))
subset1$test <- "DIF Sp."
allset4 <- rbind(subset1, subset2)

# Randomly subset for plotting
#mrg4 <- allset4  %>% sample_n(1000) 

# Plot
ad_coph_plot <- ggplot(data=subset(allset4), aes(y=mst3*100, x=log10(value.x))) + 
  geom_point(aes(color=test), shape=16, alpha=0.5, stroke = 0, size=0.75) + 
  scale_x_continuous(expand=c(0,0), limits=c(-6,-1)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,100)) + theme_bw() +
  scale_fill_manual(values=c("#332288","#88CCEE","#bfaf41","#AA4499","#44AA99","#CC6677","grey50")) +
  scale_color_manual(values=c("#332288","#88CCEE","#bfaf41","#AA4499","#44AA99","#CC6677","grey50")) +
  xlab("Nucleotide divergence (cophenetic distance)") + ylab("Allelic mismatches (% of loci)") +
  theme( panel.grid = element_blank(), plot.title = element_text(face="bold", color="black"), 
         axis.line = element_line(colour = "black"), panel.border = element_blank(),
         legend.title = element_blank(), legend.position="none",
         axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=10),
         axis.title.x = element_text(face="bold", color="black", size=10)) + coord_flip()


```
## Figure 3C
```{r}
# Import metadata
linmeta <- read.csv("lineage_meta.csv", header=TRUE, na.strings=c(""," ","NA"))

# Set LINcode strings
linmeta$gy <- paste(linmeta$L1, linmeta$L2,linmeta$L3, sep="_")

# Get counts by country and LINcode
country_count <- linmeta %>% subset(!is.na(country)) %>% group_by(gy) %>% select(gy,country) %>% unique() %>% summarise(countries=n())

# Get counts by serotype and LINcode
sero_count <-linmeta %>% subset(!is.na(serotype)) %>% subset(serotype!="inconclusive" & serotype!="genetic variant") %>% group_by(gy,serotype) %>% 
  #filter(n() > 2) %>% 
  ungroup %>% group_by(gy) %>% select(gy,serotype) %>% unique() %>% summarise(serotype2=n())

# Get counts by Lineage
all_count <- linmeta %>% select(id,gy) %>% group_by(gy) %>% summarise(count=n())

# Merge dataframes
testmergecg_LIN_country_serotype <- merge(country_count,sero_count, by.x=c("gy"), by.y=c("gy"), all=TRUE)
testmergecg_LIN_country_serotype_all <- merge(all_count,testmergecg_LIN_country_serotype, by.x=c("gy"), by.y=c("gy"), all=TRUE) %>% na.omit()  %>% subset(gy!="NA_NA_NA")

# Plot
bubl_div <- ggplot(data=(testmergecg_LIN_country_serotype_all)) + 
  geom_point(aes(y=countries,x=count, size=(serotype2)),shape=21,stroke=0,color="black",fill="grey50", alpha=0.25) + 
  theme_bw() + xlab("Number of isolates") + ylab("Number of countries") +
  scale_x_continuous(limits=c(0,1500), expand=c(0,0)) +
  scale_size_area(max_size = 7, limits=c(1,25), breaks=c(1,5,10,15,20,25)) +
  scale_y_continuous(limits=c(0,50), expand=c(0,0))  + 
  theme( panel.grid = element_blank(), legend.position = c(0.87, 0.45),
         axis.line = element_line(colour = "black"), panel.border = element_blank(),
         legend.title = element_blank(),axis.text.x=element_text(color="black", size=6),
         axis.text.y=element_text(color="black", size=6),axis.title.y = element_text(face="bold", color="black", size=10),
         axis.title.x = element_text(face="bold", color="black", size=10))
````
## Figure 3D
```{r}
# Read in LIN codes
k5_lin <- read.table("cgMLST_groupings/lincodes_all.list", header=TRUE, sep=",", fill=TRUE) %>% na.omit()

# Create LIN codes for clonal groups
k5_lin$gx <- paste(k5_lin$L1, k5_lin$L2,k5_lin$L3,k5_lin$L4,k5_lin$L5, sep = "_")

# Create LIN codes for lineages
k5_lin$gy <- paste(k5_lin$L1, k5_lin$L2,k5_lin$L3, sep = "_")

# Get counts of number isolates for 10 most prevalent lineages in the PGL
lincounts <- k5_lin %>% select(gy,id) %>% group_by(gy) %>% summarise(count=n()) %>% arrange(-count) %>% head(10)

# Select relevant columns
lin2x <- k5_lin %>% select(gy,gx) %>% na.omit() 

# Create empty dataframe to add values to
L_data2 <- data.frame(size = numeric(), counts = numeric(),ID = character())

# Calculate rarefaction
zz <- 0
repeat { # Repeat calculations for desired number of replicates 
  zz = zz+1  # Replicate count
  for(i in 1:10) { # Run for 10 lineages
    for(j in 1:as.numeric(lincounts[i,2])) { # Run over the total number of isolates for each lineage
      x <- subset(lin2x, gy==as.character(lincounts[i,1])) %>% sample_n(j) %>% unique() # For each lineage, calculate the number of unique clonal groups
      xy <- data.frame(size = j, counts=nrow(x), as.character(lincounts[i,1])) # Add label to dataframe
      L_data2[nrow(L_data2) + 1, ] <- xy  # Add data to empty dataframe
    }}
  if (zz == 3){ # zz = number of replicates
    break # Finish when desired number of replicates has been reached
  }
}

# Plot data
lin_rare <- ggplot(data=subset(L_data2)) + 
  geom_smooth(aes(x=size, y=counts, color=ID), size=0.7)+ 
  scale_x_continuous(expand=c(0,0), limits=c(0,1500)) + 
  scale_y_continuous(expand=c(0,0), limits=c(0,500)) +
  theme_bw() + 
  theme( panel.grid = element_blank(), plot.title = element_text(face="bold", color="black"), axis.line = element_line(colour = "black"), panel.border = element_blank(),
         legend.title = element_blank(),axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=10),axis.title.x = element_text(face="bold", color="black", size=10),
         legend.position = "none") 
```
## Figure 3E
```{r}
# Import sample data
metadata_short <- read.csv("FINAL_DATA/temp_summary.csv", header=TRUE, check.names = FALSE)
metadata_mandrake <- read.csv("def_knn_5k.embedding_hdbscan_clusters.csv", header=TRUE, check.names = FALSE)
metadata_merged_b <- merge(metadata_short, metadata_mandrake, by.x=c("id"), by.y=c("id"), all=TRUE)

# Create LINcode strings
metadata_merged_b$LL <- paste(metadata_merged_b$L1,  sep = "_")
metadata_merged_b$LK <- paste(metadata_merged_b$L1,metadata_merged_b$L2,  sep = "_")
metadata_merged_b$LJ <- paste(metadata_merged_b$L1,metadata_merged_b$L2,metadata_merged_b$L3,  sep = "_")
metadata_merged_b$LI <- paste(metadata_merged_b$L1,metadata_merged_b$L2,metadata_merged_b$L3,metadata_merged_b$L4,  sep = "_")
metadata_merged_b$LH <- paste(metadata_merged_b$L1,metadata_merged_b$L2,metadata_merged_b$L3,metadata_merged_b$L4,metadata_merged_b$L5,  sep = "_")
metadata_merged_b$LG <- paste(metadata_merged_b$L1,metadata_merged_b$L2,metadata_merged_b$L3,metadata_merged_b$L4,metadata_merged_b$L5,metadata_merged_b$L6,  sep = "_")
metadata_merged_b$LF <- paste(metadata_merged_b$L1,metadata_merged_b$L2,metadata_merged_b$L3,metadata_merged_b$L4,metadata_merged_b$L5,metadata_merged_b$L6,metadata_merged_b$L7,  sep = "_")
metadata_merged_b$LE <- paste(metadata_merged_b$L1,metadata_merged_b$L2,metadata_merged_b$L3,metadata_merged_b$L4,metadata_merged_b$L5,metadata_merged_b$L6,metadata_merged_b$L7,metadata_merged_b$L8,  sep = "_")
metadata_merged_b$LD <- paste(metadata_merged_b$L1,metadata_merged_b$L2,metadata_merged_b$L3,metadata_merged_b$L4,metadata_merged_b$L5,metadata_merged_b$L6,metadata_merged_b$L7,metadata_merged_b$L8,metadata_merged_b$L9,  sep = "_")
metadata_merged_b$LC <- paste(metadata_merged_b$L1,metadata_merged_b$L2,metadata_merged_b$L3,metadata_merged_b$L4,metadata_merged_b$L5,metadata_merged_b$L6,metadata_merged_b$L7,metadata_merged_b$L8,metadata_merged_b$L9,metadata_merged_b$L10,   sep = "_")
metadata_merged_b$LB <- paste(metadata_merged_b$L1,metadata_merged_b$L2,metadata_merged_b$L3,metadata_merged_b$L4,metadata_merged_b$L5,metadata_merged_b$L6,metadata_merged_b$L7,metadata_merged_b$L8,metadata_merged_b$L9,metadata_merged_b$L10,metadata_merged_b$L11,   sep = "_")

# Subset, removing rows with missing values
metadata_merged_b_CC <- metadata_merged_b %>% subset(!is.na(`Clonal complex`)) %>% subset(!is.na(LK)) 
metadata_merged_b_GPSC <- metadata_merged_b %>% subset(!is.na(GPSC)) %>% subset(!is.na(LK)) 
metadata_merged_b_ST <- metadata_merged_b %>% subset(!is.na(ST)) %>% subset(!is.na(LK))  %>% subset(!is.na(LK)) 
metadata_merged_b_rST <- metadata_merged_b %>% subset(!is.na(rST)) %>% subset(!is.na(LK))  %>% subset(!is.na(LK)) 
metadata_merged_b_mandrake <- metadata_merged_b %>% subset(!is.na(mandrake_cluster)) %>% subset(!is.na(LK)) 

# Calculate adjusted Rand Index (examples below, repeat for each pairwise comparison)
adjustedRandIndex(metadata_merged_b_CC$`Clonal complex`, metadata_merged_b_CC$LF)
adjustedRandIndex(metadata_merged_b_GPSC$GPSC, metadata_merged_b_GPSC$LK)
adjustedRandIndex(metadata_merged_b_ST$ST, metadata_merged_b_ST$LK)
adjustedRandIndex(metadata_merged_b_rST$rST, metadata_merged_b_rST$LK)
adjustedRandIndex(metadata_merged_b_mandrake$mandrake_cluster, metadata_merged_b_mandrake$LK)
#rand <- read.csv("rand.csv", header=TRUE)

# Plot
rand_plot <- ggplot(data=subset(rand, Level>1)) + 
  geom_point(aes(x=Level, y=value, color=group, shape=group), size=2.5) +
  geom_line(aes(x=Level, y=value, color=group)) +
  theme_bw() + xlab("") + ylab("Adjusted Rand index") +
  scale_x_continuous(expand=c(0,0), limits=c(1,11), breaks=c(2,3,4,5,6,7,8,9,10,11)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,1)) +
  scale_shape_manual(values=c(15,16,17,18,19)) +
  scale_color_manual(values=c("grey20","black","grey40","grey80","grey60")) +
  scale_fill_manual(values=c("grey20","black","grey40","grey80","grey60")) +
  theme( panel.grid.major =  element_blank(), plot.title = element_text(face="bold", color="black"),
         axis.line = element_line(colour = "black"),
         panel.border = element_blank(),panel.grid.minor =element_blank(),legend.title = element_blank(),
         axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=8),
         axis.title.x = element_text(face="bold", color="black", size=8), legend.position = "none") 
```
