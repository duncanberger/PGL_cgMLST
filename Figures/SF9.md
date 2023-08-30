# Supplementary Figure 9
```{r}
# Load libraries
library(ggtree)

# Import lineage metadata
linmeta <- read.csv("lineage_meta.csv", header=TRUE, na.strings=c(""," ","NA"))

# Plot 20 most prevalent serotypes on trees
ggtree(tree3, layout="circular", size=0.1, color="grey30")  %<+% subset(linmeta, serotype=="19F" | serotype=="19A" | 
                                                                          serotype=="23F" | serotype=="6A" | serotype=="1" | serotype=="14" | 
                                                                          serotype=="15BC" | serotype=="6B" | serotype=="3" | serotype=="35B" ) +  
  geom_tippoint(aes(color=serotype)) +
  geom_cladelabel(node=2518, color='black', label = "0_58_0", offset = 0.00015, barsize = 0.45, hjust = 1.5) +
  geom_cladelabel(node=2927, color='black', label = "0_15_4", offset = 0.00015, barsize = 0.45,hjust = 1.5)  +
  geom_cladelabel(node=2871, color='black', label = "0_15_0", offset = 0.00015, barsize = 0.45,hjust = 1.5) +
  geom_cladelabel(node=2647, color='black', label = "0_0_0", offset = 0.00015, barsize = 0.45,hjust = 1.5) + 
  geom_cladelabel(node=2193, color='black', label = "0_75_0", offset = 0.00015, barsize = 0.45,hjust = 1.5) +
  geom_cladelabel(node=3439, color='black', label = "0_20_1", offset = 0.00015, barsize = 0.45,hjust = 1) +
  geom_cladelabel(node=1968, color='black', label = '0_11_3', offset = 0.00015, barsize = 0.45,hjust = 1.5) + 
  geom_cladelabel(node=2496, color='black', label = '0_55_1', offset = 0.00015, barsize = 0.45,hjust = 1.5) +
  geom_cladelabel(node=2029, color='black', label = '0_11_8', offset = 0.00015, barsize = 0.45,hjust = 1.5)  +
  geom_cladelabel(node=2193, color='black', label = "0_75_0", offset = 0.00015, barsize = 0.45,hjust = 1.5) +
  geom_cladelabel(node=2783, color='black', label = "0_57_0", offset = 0.00015, barsize = 0.45,hjust = 1.5)  +
  geom_cladelabel(node=3752, color='black', label = "0_145_2", offset = 0.00015, barsize = 0.45, hjust=-1) +
  geom_cladelabel(node=3077, color='black', label = "0_130_1", offset = 0.00015, barsize = 0.45) + 
  geom_cladelabel(node=3019, color='black', label = "0_10_0", offset = 0.00015, barsize = 0.45, hjust=1) +
  geom_cladelabel(node=3619, color='black', label = "0_169_0", offset = 0.00015, barsize = 0.45, hjust=1) +
  geom_cladelabel(node=2430, color='black', label = "0_232_1", offset = 0.00015, barsize = 0.45, hjust=1) +
  geom_cladelabel(node=2147, color='black', label = "0_74_0", offset = 0.00015, barsize = 0.45, hjust=1) +
  geom_cladelabel(node=2386, color='black', label = "0_237_0", offset = 0.00015, barsize = 0.45, hjust=1) +
  geom_cladelabel(node=3524, color='black', label = "0_21_0", offset = 0.00015, barsize = 0.45) +
  geom_cladelabel(node=3201, color='black', label = "0_265_0", offset = 0.00015, barsize = 0.45, hjust=1) +
  geom_cladelabel(node=3259, color='black', label = "0_310_0", offset = 0.00015, barsize = 0.45, hjust=1)

ggtree(tree3, layout="circular", size=0.1, color="grey30")  %<+% subset(linmeta, serotype=="15A" | serotype=="11A" | 
                                                                                  serotype=="12F" | serotype=="5" | serotype=="16F" | serotype=="23A" | 
                                                                                  serotype=="9V" | serotype=="34" | serotype=="6C" | serotype=="18C" ) +  
  geom_tippoint(aes(color=serotype)) +
  geom_cladelabel(node=2518, color='black', label = "0_58_0", offset = 0.00015, barsize = 0.45, hjust = 1.5) +
  geom_cladelabel(node=2927, color='black', label = "0_15_4", offset = 0.00015, barsize = 0.45,hjust = 1.5)  +
  geom_cladelabel(node=2871, color='black', label = "0_15_0", offset = 0.00015, barsize = 0.45,hjust = 1.5) +
  geom_cladelabel(node=2647, color='black', label = "0_0_0", offset = 0.00015, barsize = 0.45,hjust = 1.5) + 
  geom_cladelabel(node=2193, color='black', label = "0_75_0", offset = 0.00015, barsize = 0.45,hjust = 1.5) +
  geom_cladelabel(node=3439, color='black', label = "0_20_1", offset = 0.00015, barsize = 0.45,hjust = 1) +
  geom_cladelabel(node=1968, color='black', label = '0_11_3', offset = 0.00015, barsize = 0.45,hjust = 1.5) + 
  geom_cladelabel(node=2496, color='black', label = '0_55_1', offset = 0.00015, barsize = 0.45,hjust = 1.5) +
  geom_cladelabel(node=2029, color='black', label = '0_11_8', offset = 0.00015, barsize = 0.45,hjust = 1.5)  +
  geom_cladelabel(node=2193, color='black', label = "0_75_0", offset = 0.00015, barsize = 0.45,hjust = 1.5) +
  geom_cladelabel(node=2783, color='black', label = "0_57_0", offset = 0.00015, barsize = 0.45,hjust = 1.5)  +
  geom_cladelabel(node=3752, color='black', label = "0_145_2", offset = 0.00015, barsize = 0.45, hjust=-1) +
  geom_cladelabel(node=3077, color='black', label = "0_130_1", offset = 0.00015, barsize = 0.45) + 
  geom_cladelabel(node=3019, color='black', label = "0_10_0", offset = 0.00015, barsize = 0.45, hjust=1) +
  geom_cladelabel(node=3619, color='black', label = "0_169_0", offset = 0.00015, barsize = 0.45, hjust=1) +
  geom_cladelabel(node=2430, color='black', label = "0_232_1", offset = 0.00015, barsize = 0.45, hjust=1) +
  geom_cladelabel(node=2147, color='black', label = "0_74_0", offset = 0.00015, barsize = 0.45, hjust=1) +
  geom_cladelabel(node=2386, color='black', label = "0_237_0", offset = 0.00015, barsize = 0.45, hjust=1) +
  geom_cladelabel(node=3524, color='black', label = "0_21_0", offset = 0.00015, barsize = 0.45) +
  geom_cladelabel(node=3201, color='black', label = "0_265_0", offset = 0.00015, barsize = 0.45, hjust=1) +
  geom_cladelabel(node=3259, color='black', label = "0_310_0", offset = 0.00015, barsize = 0.45, hjust=1)
```
