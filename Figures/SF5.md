# Supplementary Figure 5
```{r}
library(ggtree)
library(phytools)
library(ggplot2)
library(ape)

# Import table of per loci allele missingness
missing_list <- read.csv("missingness.csv", header=TRUE)

# Import phylogeny (from Figure 2a)
tree3 <- midpoint.root(read.newick("18k_tree.nwk"))

# Plot simple phylogeny
treezz <- ggtree(tree3, layout="rectangular", size=0.1, color="grey30")

# Plot merged phylogeny and missingness data
facet_plot(treezz, panel="x", 
           data=missing_list, 
           geom_segment, mapping=aes(x=0,xend=MISSING, y=y, yend=y), 
             color='grey80', lwd=.1) + theme_tree2()
```
