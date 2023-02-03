# Supplementary Figure 6
```{r}
# Load libraries
library(ggpubr)
library(ggplot2)

# Read in data
loci_des <- read.csv("loci_description.csv", header=TRUE)

# Plot allele length distribution
supl1x_A <- ggplot(data=loci_des) + 
  geom_histogram(aes(x=length),bins=80, size=0.3, color="black", fill="#3182bd") +
  xlab("Median allelle length (bp)") + theme_bw() + 
  scale_y_continuous(expand=c(0,0), limits=c(0,100)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,6000)) +
  ylab("Frequency") +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face="bold", color="black"),
        legend.title = element_blank(),
        axis.text.x=element_text(color="black", size=6),
        axis.text.y=element_text(color="black", size=6),
        axis.title.y = element_text(face="bold", color="black", size=10),
        axis.title.x = element_text(face="bold", color="black", size=10))

# Plot unqiue alleles per loci count distribution
supl1x_B <- ggplot(data=loci_des) + 
  geom_histogram(aes(x=count),bins=80, size=0.3, color="black", fill="#3182bd") +
  xlab("Number of unique alleles") + theme_bw() + 
  scale_x_continuous(expand=c(0,0), limits=c(0,1200)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,100)) +
  ylab("Frequency") +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face="bold", color="black"),
        legend.title = element_blank(),
        axis.text.x=element_text(color="black", size=6),
        axis.text.y=element_text(color="black", size=6),
        axis.title.y = element_text(face="bold", color="black", size=10),
        axis.title.x = element_text(face="bold", color="black", size=10))

# Plot correlation between allele length and number of unique alleles per loci
supl1x_C <- ggscatter(loci_des, x = "length", y = "count",add = "reg.line",shape = 19,color="#3182bd", size=0.75,
                      add.params = list(color = "black", fill = "lightgray"),
                      conf.int = TRUE) +
  stat_cor(method = "pearson", label.x = 3000, label.y = 100) +
  ylab("Number of unique alleles") +
  xlab("Loci length (bp)") + theme_bw() + 
  scale_y_continuous(expand=c(0,0), limits=c(0,1200)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,6000)) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face="bold", color="black"),
        legend.title = element_blank(),
        axis.text.x=element_text(color="black", size=6),
        axis.text.y=element_text(color="black", size=6),
        axis.title.y = element_text(face="bold", color="black", size=10),
        axis.title.x = element_text(face="bold", color="black", size=10))
```
