Supplementary Figure 8
```{r}
# Load libraries
library(ggplot2)
library(dplyr)

# Import table of per loci allele missingness
missing_list <- read.csv("missingness.csv", header=TRUE)

# Import table of per isolate clonal complex assignment
cc_list <- read.table("test/cc.list", header = FALSE)

# Merge tables
missing_merge <- merge(cc_list,missing_list, by.x=c("V1"), by.y=c("id"))

# Get mean missing counts per clonal complex
missing_merge_counts <- missing_merge %>% group_by(V2) %>% select(V2, MISSING) %>% summarise(count=n(), mean=mean(MISSING), sstdev=sd(MISSING))

# Plot
ggplot(data=subset(missing_merge_counts, count>=10)) + 
  geom_pointrange(aes(x=mean, y=reorder(V2, as.numeric(V2)), xmin=mean-sstdev, xmax=mean+sstdev)) + theme_bw() +
  xlab("Mean number of missing alleles") +
  ylab("Clonal complex") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        plot.title = element_text(face="bold", color="black"),
        legend.title = element_blank(), legend.position = "none",
        axis.text.x=element_text(color="black", size=7),
        axis.text.y=element_text(color="black", size=3),
        axis.title.y = element_text(face="bold", color="black", size=8),
        axis.title.x = element_text(face="bold", color="black", size=8))
```
