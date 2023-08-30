# FIGURE 1
```{r}
# Import libraries
library(tidyverse)
library(ggthemes)
library(sf)
library(scales)
library(patchwork)
library(dplyr)
library(ggplot2)
```
## FIGURE 1b-e: Assembly statistics
```{r}
# Plot N50
A1 <- ggplot(data=seq_stats) +
  geom_histogram(aes(x=((N50/1000))),bins=65, color="black", fill="grey50") +
  theme_bw()+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits=c(10,3500)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3000)) +
  annotation_logticks(sides = 'b',outside = TRUE, short=unit(0.05,'cm'),mid=unit(0.1,'cm'),long=unit(0.15,'cm'))  +
  coord_cartesian(clip = 'off') +
  ylab("") +
  labs(x=expression(bold("N"["50 "]+"(kb)"))) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face="bold", color="black"),
    legend.title = element_blank(),
    legend.position = c(0.75, 0.35),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    axis.text.x=element_text(color="black", size=6),
    axis.text.y=element_text(color="black", size=6),
    axis.title.y = element_text(face="bold", color="black", size=8),
    axis.title.x = element_text(face="bold", color="black", size=8)) 

# Plot L90
A2 <- ggplot(data=seq_stats) +
  geom_histogram(aes(x=L90),bins=60, color="black", fill="grey50") +
  theme_bw()+
  scale_y_continuous(expand=c(0,0), limits=c(0,2500)) +  
  scale_x_continuous(limits=c(0,120), expand=c(0,0)) +
  coord_cartesian(clip = 'off') +
  ylab("") +  labs(x=expression(bold("L"["90"]))) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face="bold", color="black"),
    legend.title = element_blank(),
    legend.position = c(0.75, 0.35),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    axis.text.x=element_text(color="black", size=6),
    axis.text.y=element_text(color="black", size=6),
    axis.title.y = element_text(face="bold", color="black", size=8),
    axis.title.x = element_text(face="bold", color="black", size=8)) 

# Plot assembly sizes 

A3 <- ggplot(data=seq_stats) +
  geom_histogram(aes(x=Total.length/1000000),bins=60, color="black", fill="grey50") +
  theme_bw()+
  scale_y_continuous(expand=c(0,0), limits=c(0,2000)) +  
  scale_x_continuous(limits=c(1.9, 2.3), breaks=c(1.9,2.0,2.1,2.2,2.3)) +
  coord_cartesian(clip = 'off') +
  ylab("") + xlab("Assembly size (Mb)") +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face="bold", color="black"),
    legend.title = element_blank(),
    legend.position = c(0.75, 0.35),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    axis.text.x=element_text(color="black", size=6),
    axis.text.y=element_text(color="black", size=6),
    axis.title.y = element_text(face="bold", color="black", size=8),
    axis.title.x = element_text(face="bold", color="black", size=8)) 

# Plot BUSCO scores
A4 <- ggplot(data=busco) +
  geom_histogram(aes(x=V2),binwidth=0.25, color="black", fill="grey50") +
  theme_bw()+
  scale_y_continuous(expand=c(0,0), limits=c(0,2000)) +  
  scale_x_continuous(limits=c(90,101), breaks=c(90,92.5,95,97.5,100)) +
  coord_cartesian(clip = 'off') +
  ylab("") + xlab("Complete (%)") +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face="bold", color="black"),
    legend.title = element_blank(),
    legend.position = c(0.75, 0.35),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    axis.text.x=element_text(color="black", size=6),
    axis.text.y=element_text(color="black", size=6),
    axis.title.y = element_text(face="bold", color="black", size=8),
    axis.title.x = element_text(face="bold", color="black", size=8)) + coord_cartesian(xlim=c(90,100))

# Merge plots

A1/A2/A3/A4
```
## FIGURE 1f: WORLD MAP
```{r}
## Import PubMLST data  
metadata_PGL <- read.csv("PGL_FREEZE_21072023 (1).csv" , na.strings=c(" ","NA",""), comment.char = "", header=TRUE)

## Import list of vaccine/non-vaccine serotypes
vac_sero <- read.table("vaccine_serotypes.csv", header=TRUE, sep=",")

# Import table of sequencing statistics 
seq_stats <- read.csv("sequence_bins_stats.csv", sep=",")
busco <- read.table("csc.busco.txt")

# Import clonal complex/serotypes
cc_st <- read.table("cc.st.list", header=TRUE)

## Import world map data
world <- map_data("world")

## Change country name format to match map data, then summarize the results for plotting
country_counts <- metadata_PGL %>%
  mutate( across(.cols = everything(),~str_replace( ., 'UK \\[England\\]', "UK"))) %>%
  mutate( across(.cols = everything(),~str_replace( ., "UK \\[Scotland\\]", "UK"))) %>%
  mutate( across(.cols = everything(),~str_replace( ., "The Gambia", "Gambia"))) %>%
  mutate( across(.cols = everything(),~str_replace( ., "The Netherlands", "Netherlands"))) %>%
  mutate( across(.cols = everything(),~str_replace( ., "Trinidad and Tobago", "Trinidad"))) %>%
  mutate( across(.cols = everything(),~str_replace( ., "China \\[Hong Kong\\]", "China"))) %>%
  mutate( across(.cols = everything(),~str_replace( ., "Congo \\[DRC\\]", "Democratic Republic of the Congo"))) %>% 
  select(country) %>% dplyr::group_by(country) %>% dplyr::summarize(count = n())

# Convert map data to sf format set reference system
world.sf <- sf::st_as_sf(world, coords = c("long", "lat"), crs = 4326 ) %>% 
  group_by(group) %>% 
  summarize(do_union = FALSE) %>%
  st_cast("POLYGON") %>% 
  ungroup()

# Get list of countries by rowid

countries = map_data("world")  %>% distinct(region) %>% rowid_to_column()

# Merge map data and metadata
merge_country <- merge(countries, country_counts, by.x=("region"), by.y=("country"), all=TRUE) %>% subset(region!="Unknown")
world_x1 <- world %>% select(group,region)
world_x2 <- merge(world_x1, merge_country, by.x=("region"), by.y=("region"))
world_x3 <- merge(world_x2, world.sf, by.x=("group"), by.y=("group"))

# Group sample counts by country into categories
world_x3$gcount <- ifelse(world_x3$count<1000,"A",
                          ifelse(world_x3$count>=1000 & world_x3$count<2000,"B",
                                 ifelse(world_x3$count>=2000 & world_x3$count<3000,"C",
                                        ifelse(world_x3$count>=3000 & world_x3$count<4000,"D",
                                               ifelse(world_x3$count>=4000 & world_x3$count<5000,"E",NA)))))

# Plot
world_plotb <- world_x3 %>% unique() %>% ggplot() +
  geom_sf(colour = "black", aes(fill = gcount, geometry=geometry)) + 
  coord_sf(ylim = c(-50, 90), datum = NA) +
  theme(panel.background = element_rect(fill = 'white')) + theme_void() +
  scale_fill_manual(values=c("#c7a2cd","#ad76b5","#8d5196","#64396a","#3a213e"), na.value="white") +
  theme_void() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.3),
        plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position = "none")
```

## FIGURE 1g: Counts by year
```{r}
# Summarize counts by year
metadata_PGL_year <- metadata_PGL %>% 
  dplyr::group_by(year) %>% dplyr::summarize(count = n()) %>% select(year, count) %>% na.omit()

# Plot
PGL_year <- ggplot(data=metadata_PGL_year) +
  geom_col(aes(x=year,y=(count/26411)*100), color="black", fill="grey50") + theme_bw()  +
  scale_y_continuous(expand=c(0,0), limits=c(0,25)) +
  scale_x_continuous(expand=c(0,0), limits=c(1910,2020)) +
  xlab("Year collected") + ylab("Genomes (%)") +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face="bold", color="black"),
    legend.title = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    axis.text.x=element_text(color="black", size=6),
    axis.text.y=element_text(color="black", size=6),
    axis.title.y = element_text(face="bold", color="black", size=8),
    axis.title.x = element_text(face="bold", color="black", size=8))
```

## FIGURE 1h: 
```{r}
# Get counts by serotype
metadata_PGL_serotype <- metadata_PGL %>% dplyr::mutate(across('serotype', str_replace, "6E\\(6Bii\\)", 'genetic variant'))  %>%
  dplyr::group_by(serotype) %>% 
  dplyr::summarize(count = n()) %>% 
  select(serotype, count) %>% na.omit() %>% 
  subset(serotype!="inconclusive"& serotype!="nontypable" ) 

# Select columns
vac_sero2 <- vac_sero %>% select(Serotype_renamed) %>% unique()
vac_sero2$VS <- "Vaccine serotype"

# Merge vaccine serotypes and metadata
metadata_PGL_serotype_vs <- merge(metadata_PGL_serotype, vac_sero2, by.x=("serotype"), by.y=("Serotype_renamed"), all=TRUE)
metadata_PGL_serotype_vs$VS[is.na(metadata_PGL_serotype_vs$VS)] <- "Non-vaccine serotype"

# Plot
PGL_sero <- ggplot(data=metadata_PGL_serotype_vs %>% arrange(-count) %>% subset(serotype!="genetic variant") %>% head(24)) +
  geom_col(aes(x=reorder(serotype, -count),y=(count/25290)*100, fill=VS), color="black") + theme_bw()  +
  scale_y_continuous(expand=c(0,0), limits=c(0,10)) +
  xlab("") + ylab("Genomes (%)") +
  scale_fill_manual(na.value="white", values=c("grey85","grey35")) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face="bold", color="black"),
    legend.title = element_blank(),
    legend.position = c(0.75, 0.75),
    axis.text.x=element_text(color="black", size=6),
    axis.text.y=element_text(color="black", size=6),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    axis.title.y = element_text(face="bold", color="black", size=8),
    axis.title.x = element_text(face="bold", color="black", size=8)) +
  guides(shape = guide_legend(override.aes = list(size = 0.3))) 
```

## FIGURE 1i: 
```{r}
# Select sequence types
tf <- metadata_PGL %>% select(ST..MLST.) %>% na.omit() 

# Create an empty data frame 
data_ST_u <- data.frame(size = numeric(),            
                   counts = numeric())

# Subsample for each sample size and calculate unique sequence types
for(i in 1:30894) {
  x <- sample_n(tf, i) %>% unique()
  xy <- data.frame(size = i, counts=nrow(x))
  data_ST_u[nrow(data_ST_u) + 1, ] <- xy 
}

# Select GPSCs
tree_meta <- read.csv("../PGL/cgMLST/tree_meta1.csv")
tree_meta$nGPS <- ifelse(grepl(";", tree_meta$GPS), NA, tree_meta$GPS)
tf2 <- tree_meta %>% select(nGPS) %>% na.omit() 

# Create an empty data frame 
data_GPS_u <- data.frame(size = numeric(),            
                    counts = numeric())

# Subsample for each sample size and calculate unique GPSCs
for(i in 1:27989) {
  x <- sample_n(tf2, i) %>% unique()
  xy <- data.frame(size = i, counts=nrow(x))
  data_GPS_u[nrow(data_GPS_u) + 1, ] <- xy 
}

# Select ribosomal sequence types
tf3 <- metadata_PGL %>% select(rST..Ribosomal.MLST.) %>% na.omit() #%>% group_by(ST..MLST.) %>% dplyr::summarize(count = n())

# Create an empty data frame 
data_rst_u <- data.frame(size = numeric(),            
                    counts = numeric())

# Subsample for each sample size and calculate ribosomal sequence type
for(i in 1:30648) {
  x <- sample_n(tf3, i) %>% unique()
  xy <- data.frame(size = i, counts=nrow(x))
  data_rst_u[nrow(data_rst_u) + 1, ] <- xy 
}


# Select clonal complexes
cc_list <- read.table("test/cc.list", header = FALSE)
tf5 <- cc_list %>% select(V2) %>% na.omit() 

# Create an empty data frame 
data_CC_u <- data.frame(size = numeric(),            
                    counts = numeric())

# Subsample for each sample size and calculate ribosomal sequence type
for(i in 1:30892) {
  x <- sample_n(tf5, i) %>% unique()
  xy <- data.frame(size = i, counts=nrow(x))
  data_CC_u[nrow(data_CC_u) + 1, ] <- xy 
}

#data_ST_u <- read.csv("data_ST_u.csv", sep=" ")
#data_rst_u <- read.csv("data_rst_u.csv", sep=" ")
#data_GPS_u <- read.csv("data_GPS_u.csv", sep=" ")
#data_CC_u <- read.table("data_CC_u.csv", header=TRUE)

# Plot
rarefaction <- ggplot() + 
  geom_smooth(data=subset(data_ST_u, size!="A"), aes(x=as.numeric(size), y=as.numeric(counts)), color="grey60") + 
  geom_smooth(data=subset(data_GPS_u, size!="A"), aes(x=as.numeric(size), y=as.numeric(counts)), color="grey20") +
  geom_smooth(data=subset(data_rst_u, size!="A"), aes(x=as.numeric(size), y=as.numeric(counts)), color="grey80") + 
  geom_smooth(data=subset(data_CC_u, size!="A"), aes(x=as.numeric(size), y=as.numeric(counts)), color="grey40") + 
  theme_bw() +
  xlab("Sample size") + ylab("Sequence types") +
  scale_x_continuous(expand=c(0,0), limits=c(0,32000)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,9000), breaks=c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000)) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face="bold", color="black"),
    legend.title = element_blank(),
    legend.position = c(0.75, 0.85),
    axis.text.x=element_text(color="black", size=6),
    axis.text.y=element_text(color="black", size=6),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    axis.title.y = element_text(face="bold", color="black", size=8),
    axis.title.x = element_text(face="bold", color="black", size=8))
```

## FIGURE 1j: 
```{r}
# Count number of samples per clonal complex
tf4x <- cc_st %>% select(CC,ST) %>% na.omit() %>% group_by(CC) %>% summarise(count_CC = n())

# Count number sequence types per clonal complex
tf4y <- cc_st %>% select(CC,ST) %>% na.omit() %>% group_by(CC,ST) %>% unique() %>% group_by(CC) %>% summarise(count_CC = n())

# Merge two dataframes
tf4z <- merge(tf4x,tf4y, by=c("CC"))

# Plot
tx <- ggplot(data=tf4z) + 
  geom_point(aes(x=count_CC.x, y=count_CC.y), color="black",fill="grey10", alpha=0.3, shape=21) +
  theme_bw() +
  scale_y_continuous(expand=c(0,0), limits=c(0,160)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,1600)) +
  xlab("No. of isolates") + ylab("Unique sequence types") + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face="bold", color="black"),
    legend.title = element_blank(),
    legend.position = c(0.75, 0.85),
    axis.text.x=element_text(color="black", size=6),
    axis.text.y=element_text(color="black", size=6),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    axis.title.y = element_text(face="bold", color="black", size=8),
    axis.title.x = element_text(face="bold", color="black", size=8))
```
