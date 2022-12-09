# Loci characteristics split by significant pairwise homoplasty index values
```{r}
# Load libraries
library(ggplot2)

# Read in loci characteristics and phi scores
loci_des <- read.csv("loci_description.csv", header=TRUE)
phi <- read.csv("phi.csv", header=TRUE, na.strings = "--")

# Filter for significant phi scores (after Bonferroni correction)
phi_sign <- subset(phi, phi<4.091653e-05)  
phi_sign$x <- "A"

# Merge dataframes
phi2 <- merge(phi_sign, loci_des, by.x=c("loci"), by.y=c("loci"), all.y=TRUE)

# Plot
ggplot(data=phi2) + geom_density(aes(x=count, fill=x, color=x), alpha=0.6) +
  xlab("Number of unique alleles") +
  ylab("Density") + theme_bw() + 
  scale_y_continuous(expand=c(0,0), limits=c(0,0.005)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,3000)) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face="bold", color="black"),
        legend.title = element_blank(),
        axis.text.x=element_text(color="black", size=6),
        axis.text.y=element_text(color="black", size=6),
        axis.title.y = element_text(face="bold", color="black", size=10),
        axis.title.x = element_text(face="bold", color="black", size=10))
        
ggplot(data=phi2) + geom_density(aes(x=length, fill=x, color=x), alpha=0.6) +
  xlab("Allele length (bp)") +
  ylab("Density") + theme_bw() + 
  scale_y_continuous(expand=c(0,0), limits=c(0,0.0015)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,6000)) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face="bold", color="black"),
        legend.title = element_blank(),
        axis.text.x=element_text(color="black", size=6),
        axis.text.y=element_text(color="black", size=6),
        axis.title.y = element_text(face="bold", color="black", size=10),
        axis.title.x = element_text(face="bold", color="black", size=10))
        
# Merge plots
supl3x_A|supl3x_B
```
