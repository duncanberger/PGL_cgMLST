# FIGURE 4
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
```
## Figure 4a
```{r}
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
## Figure 4b
```{r}
# Import tree
nk27 <- midpoint(read.newick("10e.nwk"))

# Import metadata, add LIN codes
linmeta <- read.csv("lineage_meta.csv", header=TRUE, na.strings=c(""," ","NA"))
lincodes$gx <- paste(lincodes$L1, lincodes$L2, sep="_")
lincodes$gy <- paste(lincodes$L1,lincodes$L2,lincodes$L3, sep = "_")
lincodes$gz <- paste(lincodes$L1, lincodes$L2,lincodes$L3,lincodes$L4, sep="_")
lincodes$gz2 <- paste(lincodes$L1, lincodes$L2,lincodes$L3,lincodes$L4,lincodes$L5, sep="_")

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
px2 <- viewClade(px1, 204) + 
  scale_color_manual(values = cols, na.value = "grey50") +
  geom_tippoint(aes(color=gz)) + geom_treescale() + geom_tiplab(aes(label=gz2), hjust=-0.75, size=2) +
  geom_tippoint(aes(color=as.factor(gz)))

# Plot specific subset of clade
px3 <- viewClade(px1, 210) + 
  scale_color_manual(values = cols, na.value = "grey50") +
  geom_tippoint(aes(color=gz)) + geom_treescale() + geom_tiplab(aes(label=gz2), hjust=-0.75, size=2) +
  geom_tippoint(aes(color=as.factor(gz))) 
```
