# Supplementary Figure 3
```{r}
Import libraries
library(ggplot2)
library(dplyr)

# Import table of loci lengths
loci_des <- read.csv("loci_description.csv", header=TRUE)

# Import table of phi scores
phi <- read.csv("phi.csv", header=TRUE, na.strings = "--")

# Import table of coordination for each loci/reference assembly combination
reps <- read.table("filt.blast.hits", header=FALSE) %>% select(V1,V2,V3,V4)

# Subset loci with significant phi scores
phi_sign <- subset(phi, phi<4.091653e-05)  
phi_sign$x <- "A"

# Merge dataframes
phi2 <- merge(phi_sign, loci_des, by.x=c("loci"), by.y=c("loci"), all.y=TRUE)

# Import full PGL cgMLST allele profiles
profiles_list <- read.csv("FINAL_DATA/PGL_cgMLST_profiles.csv", header=TRUE, na.strings=c(""," ","NA")) %>% select(-cgST,-LINcode)

# Make the description of missing allele uniform
profiles_list[profiles_list==0] <- NA

# Count missing values
profiles_list$count_na <- rowSums(is.na(profiles_list))

#profiles_list2 <- profiles_list %>% select(-id, -count_na)
#counts_distinct <- as.data.frame(sapply(profiles_list2, n_distinct))

# Get allele frequencies for all samples in the PGL
profiles_list_melt <- profiles_list  %>% reshape2::melt((id.vars = c("id")))
profiles_list_melt_loci_count <- profiles_list_melt %>% select(variable, value) %>% unique() %>% group_by(variable) %>% summarize(count=n())
profiles_list_melt_loci_count_b <- profiles_list_melt_loci_count %>% subset(variable!="count_na")

# Merge allele count and phi data
freq_rep <- merge(reps, profiles_list_melt_loci_count_b, by.x=c("V4"), by.y=c("variable"),all.x=TRUE)
freq_rep_phi <- merge(freq_rep, phi, by.x=c("V4"), by.y=c("loci"),all.x=TRUE)

# Plot allele frequencies across the genome
all_recomb <- ggplot() + 
  geom_point(data=subset(freq_rep_phi, !is.na(phi) & phi<4.091653e-05), aes(x=V2/1000000,y=count), color="red", size=0.75) +
  geom_point(data=subset(freq_rep_phi, is.na(phi) | phi>=4.091653e-05), aes(x=V2/1000000,y=count), color="grey50", size=0.75) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3000)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,2.230000)) +
  facet_grid(rows = vars(V1), scales = "free_x") + 
  theme_bw() + xlab("Position (Mb)") + ylab("Allele frequency") +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face="bold", color="black"),
        legend.title = element_blank(),
        axis.text.x=element_text(color="black", size=6),
        axis.text.y=element_text(color="black", size=6),
        axis.title.y = element_text(face="bold", color="black", size=8),
        axis.title.x = element_text(face="bold", color="black", size=8))

# Plot unique alleles per locus
recomb_A1 <- ggplot(data=phi2) + geom_density(aes(x=count, fill=x, color=x), alpha=0.6) +
  xlab("Number of unique alleles") +
  ylab("Density") + theme_bw() + 
  scale_color_manual(values=c("red","red"),na.value="grey50") +
  scale_fill_manual(values=c("red","red"),na.value="grey50") +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.005)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,3000)) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face="bold", color="black"),
        legend.title = element_blank(), legend.position = "none",
        axis.text.x=element_text(color="black", size=6),
        axis.text.y=element_text(color="black", size=6),
        axis.title.y = element_text(face="bold", color="black", size=8),
        axis.title.x = element_text(face="bold", color="black", size=8))

# Plot allele length per locus
recomb_A2 <- ggplot(data=phi2) + geom_density(aes(x=length, fill=x, color=x), alpha=0.6) +
  xlab("Allele length (bp)") +
  ylab("Density") + theme_bw() + 
  scale_color_manual(values=c("red","red"),na.value="grey50") +
  scale_fill_manual(values=c("red","red"),na.value="grey50") +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.0015)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,6000)) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face="bold", color="black"),
        legend.title = element_blank(), legend.position = "none",
        axis.text.x=element_text(color="black", size=6),
        axis.text.y=element_text(color="black", size=6),
        axis.title.y = element_text(face="bold", color="black", size=8),
        axis.title.x = element_text(face="bold", color="black", size=8))
