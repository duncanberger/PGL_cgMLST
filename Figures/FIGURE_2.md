

# FIGURE 2
```{r}
# Import libraries 
library(tidyverse)
library(ggthemes)
library(sf)
library(scales)
library(patchwork)
library(dplyr)
library(ggplot2)
library(ggridges)
````
## Figure 2A
```{r}
# Remove pairwise relationships between Streptococcus pseudopneumoniae
m5000_ALL_sub <- mstclust_5000_melted_b  %>% subset(!is.na(species.x) && !is.na(species.y)) %>% select(variable, ID, value) %>% sample_n(200000) 

# Plot
borders <- ggplot() + theme_bw() +
  geom_density(data=subset(m5000_ALL_sub), aes(x=(value*100)), bw=0.00001, color="black", fill="#5F84A2", alpha=1) + 
  scale_y_continuous(expand=c(0,0), limits=c(0,2400), breaks=c(0,300,600,900,1200,1500,1800,2100,2400)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
  geom_vline(xintercept = (1180/1222)*100, color="grey60", linetype="dashed") + 
  geom_vline(xintercept = (750/1222)*100, color="grey60", linetype="dashed") +
  geom_vline(xintercept = (540/1222)*100, color="grey60", linetype="dashed") +
  geom_vline(xintercept = (160/1222)*100, color="grey60", linetype="dashed") +
  geom_vline(xintercept = (25/1222)*100, color="grey60", linetype="dashed") +
  xlab("Pairwise allelic mismatches (%)") + ylab("Frequency") +
  theme( panel.grid = element_blank(), plot.title = element_text(face="bold", color="black"), axis.line = element_line(colour = "black"), panel.border = element_blank(),
         legend.title = element_blank(),legend.position="none",axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=10),axis.title.x = element_text(face="bold", color="black", size=10)) 
```
## FIGURE 2B
```{r}
# Read MSTclust output
mstclust_5000 <- read.csv("mstclust.5000.d", header=TRUE, check.names = FALSE)

# Convert dataframe to columns of pairwise relationships
mstclust_5000_melted <- reshape2::melt(mstclust_5000, id.vars = c("ID")) %>% na.omit()

# Read in metadata (short version), merge tables
metadata_short <- read.csv("short_metadata.csv", header=TRUE, check.names = FALSE)
metadata_mandrake <- read.csv("def_knn_5k.embedding_hdbscan_clusters.csv", header=TRUE, check.names = FALSE)
metadata_merged <- merge(metadata_short, metadata_mandrake, by.x=c("ID"), by.y=c("id"))

# Add metadata from both isolates for each row
mstclust_5000_melted_a <- merge(mstclust_5000_melted, metadata_merged, by.x=c("ID"), by.y=c("ID"), all=TRUE)
mstclust_5000_melted_b <- merge(mstclust_5000_melted_a, metadata_merged, by.x=c("variable"), by.y=c("ID"), all=TRUE) %>% subset(!is.na(value))

#write.table(mstclust_5000_melted_b,"mstclust_5000_melted_b.csv", sep=",", row.names=FALSE, quote=FALSE)
#mstclust_5000_melted_b <- read.csv("mstclust_5000_melted_b.csv", header=TRUE)

# Create a new dataframe for each subset of samples with matching clustering metric
m5000_GPSC <- subset(mstclust_5000_melted_b, GPSC.x==GPSC.y & !is.na(GPSC.x)) %>% filter(!grepl(';', GPSC.x)) %>% select(variable, ID, value) 
m5000_CC <- subset(mstclust_5000_melted_b, CC.x==CC.y & !is.na(CC.x)) %>% select(variable, ID, value) 
m5000_ST <- subset(mstclust_5000_melted_b, ST.x==ST.y & !is.na(ST.x)) %>% select(variable, ID, value) 
m5000_SERO <- subset(mstclust_5000_melted_b, serotype.x==serotype.y & !is.na(serotype.x) & serotype.x!="inconclusive" & serotype.x!="nontypable" & serotype.x!="genetic variant") %>% select(variable, ID, value) 
m5000_MANDRAKE <- subset(mstclust_5000_melted_b, mandrake_cluster.x==mandrake_cluster.y & !is.na(mandrake_cluster.x)) %>% select(variable, ID, value) 
m5000_rST <- subset(mstclust_5000_melted_b, rST.x==rST.y & !is.na(rST.x)) %>% select(variable, ID, value) 
m5000_ALL <- mstclust_5000_melted_b  %>% subset(species.x==species.y & species.x=="S_pneumoniae") %>% select(variable, ID, value) 

# Create a column labelling the relevant metric
m5000_ALL$TAG <- "1_ALL"
m5000_SERO$TAG <- "4_sero"
m5000_CC$TAG <- "2_CC"
m5000_GPSC$TAG <- "3_GPS"
m5000_ST$TAG <- "6_ST"
m5000_rST$TAG <- "7_rST"
m5000_MANDRAKE$TAG <- "5_MAND"

# Merge dataframes
m5000_MERGED <- (rbind(m5000_ALL, m5000_SERO, m5000_CC, m5000_GPSC, m5000_ST, m5000_rST,m5000_MANDRAKE))
#write.table(m5000_MERGED,"m5000_MERGED.csv", sep=",", row.names=FALSE, quote=FALSE)

# Subset for plotting
m5000_MERGED_test <- m5000_MERGED %>% sample_n(2000000)

# Plot
densities <- ggplot(m5000_MERGED_test, aes(x = value*100, y = TAG, height = stat(density), fill=TAG)) +
  geom_density_ridges(stat = "density",bw=0.0001, alpha=0.35,scale =2.5,color="black") + theme_bw() +
  scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
  scale_y_discrete(expand=c(0,0)) +
  geom_vline(xintercept = (1180/1222)*100, color="grey60", linetype="dashed") + 
  geom_vline(xintercept = (750/1222)*100, color="grey60", linetype="dashed") +
  geom_vline(xintercept = (540/1222)*100, color="grey60", linetype="dashed") +
  geom_vline(xintercept = (160/1222)*100, color="grey60", linetype="dashed") +
  geom_vline(xintercept = (25/1222)*100, color="grey60", linetype="dashed") +
  xlab("Pairwise allelic mismatches (%)") + ylab("") +
  scale_fill_manual(values=c("#023858","#194B6A","#305E7D","#487190","#5F84A2","#7797B5","#8EAAC8","#A6BDDB")) +
  theme( panel.grid.major =  element_blank(), plot.title = element_text(face="bold", color="black"),
         axis.line = element_line(colour = "black"),
         panel.border = element_blank(),panel.grid.minor =element_blank(),legend.title = element_blank(),
         legend.position="none",axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=10),
         axis.title.x = element_text(face="bold", color="black", size=10)) 
```
## Figure 2C
### Chaguza et al. (2020)
```{r}
# Import and merge metadata from multiple sources
meta1 <- read.csv("chaguza_meta1.csv", header=TRUE)
meta2 <- read.csv("chaguza_meta2.csv", header=TRUE)
meta_merge <- merge(meta1,meta2, by.x=c("ENA_run_accession"), by.y=c("Lane.accession"), all=TRUE) %>% subset(!is.na(id)) %>% subset(!is.na(Sample.ID)) 

# Read in table of pairwise allelic differences and convert to row format
chaguza_difs <- read.csv("mst_cluster.chaguza.input.out.2.d", header=TRUE, check.names = FALSE)
chaguza_difs_melted <- reshape2::melt(chaguza_difs, id.vars = c("ID")) %>% na.omit()

# Merge allelic differences and metadata
chaguza_difs_melted_a <- merge(chaguza_difs_melted, meta_merge, by.x=c("ID"), by.y=c("id"), all=TRUE)
chaguza_difs_melted_b <- merge(chaguza_difs_melted_a, meta_merge, by.x=c("variable"), by.y=c("id"))

# Label comparisons within and between subjects
chaguza_difs_melted_b$out1 <- ifelse(chaguza_difs_melted_b$Subject.ID.x==chaguza_difs_melted_b$Subject.ID.y, "Same subject",
                                     ifelse(chaguza_difs_melted_b$Subject.ID.x!=chaguza_difs_melted_b$Subject.ID.y,"Different subject","ERROR"))

# Subset columns for merging
chaguza_difs_melted_c <- chaguza_difs_melted_b %>% select(variable, ID, value, out1)

# Add label
chaguza_difs_melted_c$TAG <- "1CHAGUZA"
```
### Lees et al. (2017)
```{r}
# Import and merge metadata 
meta1_lees <- read.csv("lees_meta.csv", header=TRUE)

# Read in table of pairwise allelic differences and convert to row format
lees_difs <- read.csv("mst_cluster.lees.input.out.2.d", header=TRUE, check.names = FALSE)
lees_difs_melted <- reshape2::melt(lees_difs, id.vars = c("ID")) %>% na.omit()

# Merge allelic differences and metadata
lees_difs_melted_a <- merge(lees_difs_melted, meta1_lees, by.x=c("ID"), by.y=c("id"), all=TRUE)
lees_difs_melted_b <- merge(lees_difs_melted_a, meta1_lees, by.x=c("variable"), by.y=c("id"))

# Label comparisons within and between subjects
lees_difs_melted_b$out1 <- ifelse(lees_difs_melted_b$isolate.x==lees_difs_melted_b$isolate.y, "Same subject",
                                  ifelse(lees_difs_melted_b$isolate.x!=lees_difs_melted_b$isolate.y,"Different subject","BLAH"))

# Subset columns for merging
lees_difs_melted_c <- lees_difs_melted_b %>% select(variable, ID, value, out1) %>% subset(!is.na(out1))

# Add label
lees_difs_melted_c$TAG <- "2LEES"
```
### Merge datasets
```{r}
# Merge
merged_d <- rbind(chaguza_difs_melted_c, lees_difs_melted_c)

# Plot
ch_b <- ggplot(merged_d, aes(x = value*100, y = TAG, height = stat(density), fill=out1)) +
  geom_density_ridges(stat = "density",bw=0.05, alpha=0.35,scale =2,color=NA) + theme_bw() +
  scale_x_continuous(expand=c(0,0), limits=c(0,5)) +
  scale_y_discrete(expand=c(0,0)) +
  geom_vline(xintercept = (1180/1222)*100, color="grey60", linetype="dashed") + 
  geom_vline(xintercept = (750/1222)*100, color="grey60", linetype="dashed") +
  geom_vline(xintercept = (540/1222)*100, color="grey60", linetype="dashed") +
  geom_vline(xintercept = (160/1222)*100, color="grey60", linetype="dashed") +
  geom_vline(xintercept = (25/1222)*100, color="grey60", linetype="dashed") +
  geom_vline(xintercept = (15/1222)*100, color="grey60", linetype="dashed") +
  geom_vline(xintercept = (8/1222)*100, color="grey60", linetype="dashed") +
  geom_vline(xintercept = (4/1222)*100, color="grey60", linetype="dashed") +
  geom_vline(xintercept = (2/1222)*100, color="grey60", linetype="dashed") +
  geom_vline(xintercept = (1/1222)*100, color="grey60", linetype="dashed") +
  xlab("Pairwise allelic mismatches (%)") + ylab("") +
  scale_fill_manual(values=c("red","blue")) +
  theme( panel.grid.major =  element_blank(), plot.title = element_text(face="bold", color="black"),
         axis.line = element_line(colour = "black"),
         panel.border = element_blank(),panel.grid.minor =element_blank(),legend.title = element_blank(), 
         legend.position="none",axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=10),
         axis.title.x = element_text(face="bold", color="black", size=10)) 
  
ch_a <- ggplot(merged_d, aes(x = value*100, y = TAG, height = stat(density), fill=out1)) +
  geom_density_ridges(stat = "density",bw=0.25, alpha=0.35,scale = 2.3,color="black") + theme_bw() +
  scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
  scale_y_discrete(expand=c(0,0)) +
  geom_vline(xintercept = (1180/1222)*100, color="grey60", linetype="dashed") + 
  geom_vline(xintercept = (750/1222)*100, color="grey60", linetype="dashed") +
  geom_vline(xintercept = (540/1222)*100, color="grey60", linetype="dashed") +
  geom_vline(xintercept = (160/1222)*100, color="grey60", linetype="dashed") +
  geom_vline(xintercept = (25/1222)*100, color="grey60", linetype="dashed") +
  scale_fill_manual(values=c("red","blue")) +
  theme( panel.grid.major =  element_blank(), plot.title = element_text(face="bold", color="black"),
         axis.line = element_line(colour = "black"),
         panel.border = element_blank(),panel.grid.minor =element_blank(),legend.title = element_blank(),
         legend.position="none",axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=10),
         axis.title.x = element_text(face="bold", color="black", size=10)) 
```
