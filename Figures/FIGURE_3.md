# FIGURE 3
## Figure 3A: 
```{r}
# Import libraries
library(ggplot2)
library(dplyr)
library(ggsankey)
library(ggtree)
library(ape)
library(phytools)
library(phangorn)
library(mclust)

# Import table with basic example of LINcodes
sub_lin2 <- read.csv("subset_lin.csv", header=TRUE, check.names = FALSE)

# Create LINcode strings for each level
sub_lin2$LL <- paste(sub_lin2$L1,  sep = "_")
sub_lin2$LK <- paste(sub_lin2$L1,sub_lin2$L2,  sep = "_")
sub_lin2$LJ <- paste(sub_lin2$L1,sub_lin2$L2,sub_lin2$L3,  sep = "_")
sub_lin2$LI <- paste(sub_lin2$L1,sub_lin2$L2,sub_lin2$L3,sub_lin2$L4,  sep = "_")
sub_lin2$LH <- paste(sub_lin2$L1,sub_lin2$L2,sub_lin2$L3,sub_lin2$L4,sub_lin2$L5,  sep = "_")
sub_lin2$LG <- paste(sub_lin2$L1,sub_lin2$L2,sub_lin2$L3,sub_lin2$L4,sub_lin2$L5,sub_lin2$L6,  sep = "_")
sub_lin2$LF <- paste(sub_lin2$L1,sub_lin2$L2,sub_lin2$L3,sub_lin2$L4,sub_lin2$L5,sub_lin2$L6,sub_lin2$L7,  sep = "_")
sub_lin2$LE <- paste(sub_lin2$L1,sub_lin2$L2,sub_lin2$L3,sub_lin2$L4,sub_lin2$L5,sub_lin2$L6,sub_lin2$L7,sub_lin2$L8,  sep = "_")
sub_lin2$LD <- paste(sub_lin2$L1,sub_lin2$L2,sub_lin2$L3,sub_lin2$L4,sub_lin2$L5,sub_lin2$L6,sub_lin2$L7,sub_lin2$L8,sub_lin2$L9,  sep = "_")
sub_lin2$LC <- paste(sub_lin2$L1,sub_lin2$L2,sub_lin2$L3,sub_lin2$L4,sub_lin2$L5,sub_lin2$L6,sub_lin2$L7,sub_lin2$L8,sub_lin2$L9,sub_lin2$L10,   sep = "_")
sub_lin2$LB <- paste(sub_lin2$L1,sub_lin2$L2,sub_lin2$L3,sub_lin2$L4,sub_lin2$L5,sub_lin2$L6,sub_lin2$L7,sub_lin2$L8,sub_lin2$L9,sub_lin2$L10,sub_lin2$L11,   sep = "_")

# Convert to input for ggsankey
df2_q <- sub_lin2 %>% make_long(LL,LK,LJ,LI,LH,LG,LF,LE,LD,LC,LB)

# Plot
ggplot(df2_q, aes(x = x,next_x = next_x, node = (node),next_node = next_node, label = node)) + xlab("") + 
  geom_sankey(flow.alpha = 0.65, space = 6, node.color = 1, node.fill="black", width=0.01, smooth = 4) +
  scale_fill_viridis_d() + 
  scale_fill_manual(values=c("#3288bd","#5e4fa2","#AA4499", "#d53e4f", "#d5a33e" ,"#117733" ,"#d53e4f", "#AA4499","lightblue","orange")) +
  theme_sankey() + theme(legend.position = "none")
```
## Figure 3B: 
```{r}
# Import tree (from Figure 2A)
tree3 <- midpoint.root(read.newick("18k_tree.nwk"))

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
mrg6$LINEAGE.x <- paste(mrg6$L1.x, mrg6$L2.x, sep = "_")
mrg6$LINEAGE.y <- paste(mrg6$L1.y, mrg6$L2.y, sep = "_")
mrg6$SUBLINEAGE.x <- paste(mrg6$L1.x, mrg6$L2.x,mrg6$L3.x, sep = "_")
mrg6$SUBLINEAGE.y <- paste(mrg6$L1.y, mrg6$L2.y,mrg6$L3.y, sep = "_")
mrg6$CLONALGROUP.x <- paste(mrg6$L1.x, mrg6$L2.x,mrg6$L3.x,mrg6$L4.x, sep = "_")
mrg6$CLONALGROUP.y <- paste(mrg6$L1.y, mrg6$L2.y,mrg6$L3.y,mrg6$L4.y, sep = "_")
mrg6$STRAIN.x <- paste(mrg6$L1.x, mrg6$L2.x,mrg6$L3.x,mrg6$L4.x,mrg6$L5.x, sep = "_")
mrg6$STRAIN.y <- paste(mrg6$L1.y, mrg6$L2.y,mrg6$L3.y,mrg6$L4.y,mrg6$L5.y, sep = "_")
mrg6$SUBSTRAIN.x <- paste(mrg6$L1.x, mrg6$L2.x,mrg6$L3.x,mrg6$L4.x,mrg6$L5.x,mrg6$L6.x, sep = "_")
mrg6$SUBSTRAIN.y <- paste(mrg6$L1.y, mrg6$L2.y,mrg6$L3.y,mrg6$L4.y,mrg6$L5.y,mrg6$L6.y, sep = "_")
mrg6$sp.x <- ifelse(grepl(pattern = "ps_", mrg6$id), "SSP","NULL")
mrg6$sp.y <- ifelse(grepl(pattern = "ps_", mrg6$id.y), "SSP","NULL")

# Replace the false NA_* strings
mrg6[mrg6 == "NA"] <- NA
mrg6[mrg6 == "NA_NA"] <- NA
mrg6[mrg6 == "NA_NA_NA"] <- NA
mrg6[mrg6 == "NA_NA_NA_NA"] <- NA
mrg6[mrg6 == "NA_NA_NA_NA_NA"] <- NA
mrg6[mrg6 == "NA_NA_NA_NA_NA_NA"] <- NA

# Define relationships between pairs
mrg6$test <- ifelse(mrg6$LINEAGE.x!=mrg6$LINEAGE.y,"DIF LINEAGE",
                    ifelse(mrg6$sp.x!=mrg6$sp.y,"DIF S",
                           ifelse(mrg6$SUBLINEAGE.x!=mrg6$SUBLINEAGE.y, "DIF SUBLINEAGE",
                                  ifelse(mrg6$CLONALGROUP.x!=mrg6$CLONALGROUP.y, "DIF CG",
                                         ifelse(mrg6$STRAIN.x!=mrg6$STRAIN.y, "DIF STRAIN",
                                                ifelse(mrg6$LINEAGE.x==mrg6$LINEAGE.y,"SAME LINEAGE",
                                                       ifelse(mrg6$SUBLINEAGE.x==mrg6$SUBLINEAGE.y, "SAME SUBLINEAGE",
                                                              ifelse(mrg6$CLONALGROUP.x==mrg6$CLONALGROUP.y, "SAME CG",
                                                                     ifelse(mrg6$STRAIN.x==mrg6$STRAIN.y, "SAME STRAIN",
                                                                            ifelse(mrg6$STRAIN.x==mrg6$STRAIN.y, "SAME SUBSTRAIN",
                                                                                   ifelse(is.na(mrg6$LINEAGE.x) | is.na(mrg6$LINEAGE.y), "UNKNOWN","ELSE")))))))))))

# Subset and rename groupings
subset1 <- subset(mrg6, sp.x!=sp.y)
subset2 <- subset(mrg6, sp.x==sp.y & sp.x!="SSP") %>%  subset(!is.na(L2.x) & !is.na(L2.y))
subset1$test <- "DIF Sp."
allset4 <- rbind(subset1, subset2)

# Randomly subset for plotting
mrg4 <- allset4  %>% sample_n(40000) 

# Plot
ad_coph_plot <- ggplot(data=subset(mrg4), aes(y=mst3*100, x=value.x)) + 
  geom_point(aes(color=test), shape=16, alpha=0.5, stroke = 0, size=0.75) + 
  scale_x_continuous(expand=c(0,0), limits=c(0,0.0125)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,100)) + theme_bw() +
  scale_fill_manual(values=c("#44AA99","#CC6677","#AA4499","#332288","#bfaf41","#88CCEE","grey50")) +
  scale_color_manual(values=c("#44AA99","#CC6677","#AA4499","#332288","#bfaf41","#88CCEE","grey50")) +
  xlab("Nucleotide divergence (cophenetic distance)") + ylab("Allelic mismatches (% of loci)") +
  theme( panel.grid = element_blank(), plot.title = element_text(face="bold", color="black"), 
         axis.line = element_line(colour = "black"), panel.border = element_blank(),
         legend.title = element_blank(), legend.position="none",
         axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=10),
         axis.title.x = element_text(face="bold", color="black", size=10)) 
```
## Figure 3C:
```{r}
# 
# Import sample data
temp_gps <- read.table("test/gpos.short.txt", header=TRUE)
temp_cc <- read.table("test/cc.st.list", header=TRUE)
temp_mand <- read.csv("mandrake/mandrake.embedding_hdbscan_clusters.csv", header=TRUE)
temp_rst <- read.table("test/rst.list", header=TRUE)
lin <- read.csv("FINAL_DATA/short_LINcodes_cgST.csv", header=TRUE) %>% na.omit() %>% select(-cgST)

# Merge dataframes
temp_merge1a <- merge(temp_gps, temp_cc, by.x=c("ID"), by.y=c("ID"), all = TRUE)
temp_merge1b <- merge(temp_merge1a, temp_mand, by.x=c("ID"), by.y=c("id"))
temp_merge1c <- merge(temp_merge1b, temp_rst, by.x=c("ID"), by.y=c("id"))
temp_merge2 <- merge(temp_merge1c,lin, by.x=c("ID"), by.y=c("id"), all.y = TRUE ) %>% na.omit()

# Create LINcode strings
temp_merge2$LL <- paste(temp_merge2$L1,  sep = "_")
temp_merge2$LK <- paste(temp_merge2$L1,temp_merge2$L2,  sep = "_")
temp_merge2$LJ <- paste(temp_merge2$L1,temp_merge2$L2,temp_merge2$L3,  sep = "_")
temp_merge2$LI <- paste(temp_merge2$L1,temp_merge2$L2,temp_merge2$L3,temp_merge2$L4,  sep = "_")
temp_merge2$LH <- paste(temp_merge2$L1,temp_merge2$L2,temp_merge2$L3,temp_merge2$L4,temp_merge2$L5,  sep = "_")
temp_merge2$LG <- paste(temp_merge2$L1,temp_merge2$L2,temp_merge2$L3,temp_merge2$L4,temp_merge2$L5,temp_merge2$L6,  sep = "_")
temp_merge2$LF <- paste(temp_merge2$L1,temp_merge2$L2,temp_merge2$L3,temp_merge2$L4,temp_merge2$L5,temp_merge2$L6,temp_merge2$L7,  sep = "_")
temp_merge2$LE <- paste(temp_merge2$L1,temp_merge2$L2,temp_merge2$L3,temp_merge2$L4,temp_merge2$L5,temp_merge2$L6,temp_merge2$L7,temp_merge2$L8,  sep = "_")
temp_merge2$LD <- paste(temp_merge2$L1,temp_merge2$L2,temp_merge2$L3,temp_merge2$L4,temp_merge2$L5,temp_merge2$L6,temp_merge2$L7,temp_merge2$L8,temp_merge2$L9,  sep = "_")
temp_merge2$LC <- paste(temp_merge2$L1,temp_merge2$L2,temp_merge2$L3,temp_merge2$L4,temp_merge2$L5,temp_merge2$L6,temp_merge2$L7,temp_merge2$L8,temp_merge2$L9,temp_merge2$L10,   sep = "_")
temp_merge2$LB <- paste(temp_merge2$L1,temp_merge2$L2,temp_merge2$L3,temp_merge2$L4,temp_merge2$L5,temp_merge2$L6,temp_merge2$L7,temp_merge2$L8,temp_merge2$L9,temp_merge2$L10,temp_merge2$L11,   sep = "_")

# Calculate adjusted Rand Index (examples below, repeat for each pairwise comparison)
adjustedRandIndex(temp_merge2$CC, temp_merge2$LK)
adjustedRandIndex(temp_merge2$GPS, temp_merge2$LK)
adjustedRandIndex(temp_merge2$rST, temp_merge2$LK)
adjustedRandIndex(temp_merge2$ST, temp_merge2$LK)

rand <- read.csv("rand.csv", header=TRUE)

# Plot
rand_plot <- ggplot(data=rand) + 
  geom_point(aes(x=Level, y=value, color=group)) +
  geom_line(aes(x=Level, y=value, color=group)) +
  theme_bw() + xlab("") + ylab("Adjusted Rand index") +
  scale_x_continuous(expand=c(0,0), limits=c(0,11), breaks=c(0,1,2,3,4,5,6,7,8,9,10,11)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,1)) +
  scale_fill_manual(values=c("red","#117733","#332288","orange","lightblue")) +
  scale_color_manual(values=c("red","#117733","#332288","orange","lightblue")) +
  theme( panel.grid.major =  element_blank(), plot.title = element_text(face="bold", color="black"),
         axis.line = element_line(colour = "black"),
         panel.border = element_blank(),panel.grid.minor =element_blank(),legend.title = element_blank(),
         axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=8),
         axis.title.x = element_text(face="bold", color="black", size=8), legend.position = "none") 
```
## Figure 3D:
```{r}
# Import tree
nk27 <- midpoint(read.newick("10e.nwk"))

# Import metadata, add 
linmeta <- read.csv("lineage_meta.csv", header=TRUE, na.strings=c(""," ","NA"))
lincodes$gx <- paste(lincodes$L1, lincodes$L2, sep="_")
lincodes$gy <- paste(lincodes$L1,lincodes$L2,lincodes$L3, sep = "_")
lincodes$gz <- paste(lincodes$L1, lincodes$L2,lincodes$L3,lincodes$L4, sep="_")

# Create a pallette of colors in a random order
n <- 52
cols <- hue_pal(h = c(0, 360) + 15, 
                c = 100, l = 65, 
                h.start = 0, direction = 1)(n)[order(sample(1:n, n))]

# Plot lineage tree
px1 <- ggtree(nk27, layout="rectangular", size=0.1, ladderize = TRUE) %<+% subset(lincodes) +  
  scale_color_manual(values = cols, na.value = "grey50")  + 
  geom_tippoint(aes(color=as.factor(gy))) + geom_treescale() +  theme(legend.position="none")

# Plot specific clade of lineage tree
viewClade(px1, 204) + 
  scale_color_manual(values = cols, na.value = "grey50") +
  geom_tippoint(aes(color=gz)) + geom_treescale() + geom_tiplab(aes(label=gz), hjust=-0.75, size=3.5)
```
