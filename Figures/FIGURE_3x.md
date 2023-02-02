
# FIGURE 2
## FIGURE 2A: 
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

# Create random list of unique colors
n <- 160# number of colors
cols <- hue_pal(h = c(0, 360) + 15, 
                c = 100, l = 65, 
                h.start = 0, direction = 1)(n)[order(sample(1:n, n))] # color palette in random order

# Import LINCodes and create strings of relevant groupings
lincodes <- read.csv("FINAL_DATA/short_LINcodes_cgST.csv", header=TRUE) %>% na.omit()
lincodes$gx <- paste(lincodes$L1, lincodes$L2, sep = "_")
lincodes$gy <- paste(lincodes$L1, lincodes$L2,lincodes$L3, sep = "_")

# Import phylogeny
tree3 <- midpoint.root(read.newick("18k_tree.nwk"))

# Import list of samples
sample_list <- read.table("shortlist.txt")

# Subset LINcode list to relevant ids
tx22 <- lincodes[lincodes$id %in% sample_list$V1, ]

# Subset relevant lineages to be highlighted
rownames(lincodes) <- lincodes[,1]
linx <-  lincodes %>% select(gx) %>% 
  subset(gx=="0_15" | gx=="0_58" |gx=="0_0" |gx=="0_11"  | gx=="0_75" |gx=="0_21"| 
           gx=="0_57" | gx=="0_145" | gx=="0_130" | gx=="0_10" | gx=="0_169" |
           gx=="0_232" | gx=="0_74" | gx=="0_154" | gx=="0_237" | gx=="0_265" |
           gx=="0_23" | gx=="0_310" | gx=="0_20" | gx=="0_346" | gx=="0_70" | gx=="0_55" )

# Find shared node for entire lineage (e.g. Linegae 15 shown here)
liny <- tx22 %>% subset(gx=="0_15") %>% select(id) 
zinny <- as.character(unlist(tail(liny$id,200)))
findMRCA(tree3, tips=zinny, type=c("node"))

# Plot
ggtree(tree3, layout="circular", size=0.1, color="grey30")  %<+% subset(lincodes) +  
  scale_color_manual(values = cols, na.value = NA) +
  theme(legend.position = "none") +
  geom_cladelabel(node=3851, color='black', label = "PS", offset = 0.00015, barsize = 0.45)  +
  geom_cladelabel(node=3816, color='black', label = "OG", offset = 0.00015, barsize = 0.45)  +
  geom_cladelabel(node=2848, color='black', label = "0_15", offset = 0.00015, barsize = 0.45,hjust = 1) + geom_hilight(node =2848 ,alpha = .5, fill="#D89000") +
  geom_cladelabel(node=2518, color='black', label = "0_58", offset = 0.00015, barsize = 0.45) + geom_hilight(node =2518,alpha = .5, fill="#FF62BC") +
  geom_cladelabel(node=2647, color='black', label = '0_0', offset = 0.00015, barsize = 0.45) + geom_hilight(node = 2647,alpha = .5, fill="#00BF7D") +
  geom_cladelabel(node=1965, color='black', label = '0_11', offset = 0.00015, barsize = 0.45) + geom_hilight(node = 1965,alpha = .5, fill="#FA62DB") +
  geom_cladelabel(node=2193, color='black', label = "0_75", offset = 0.00015, barsize = 0.45) + geom_hilight(node =2193 ,alpha = .5, fill="#C09B00") +
  geom_cladelabel(node=3490, color='black', label = "0_21", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 3490,alpha = .5, fill="#FF6A98") +
  geom_cladelabel(node=2783, color='black', label = "0_57", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 2783,alpha = .5, fill="#00BB4E") +
  geom_cladelabel(node=3752, color='black', label = "0_145", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 3752,alpha = .5, fill="#00B0F6") +
  geom_cladelabel(node=3059, color='black', label = "0_130", offset = 0.00015, barsize = 0.45) + geom_hilight(node =3059 ,alpha = .5, fill="#A3A500") +
  geom_cladelabel(node=2995, color='black', label = "0_10", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 2995,alpha = .5, fill="#35A2FF") +
  geom_cladelabel(node=3619, color='black', label = "0_169", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 3619,alpha = .5, fill="#9590FF") +
  geom_cladelabel(node=2427, color='black', label = "0_232", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 2427,alpha = .5, fill="#C77CFF") +
  geom_cladelabel(node=2147, color='black', label = "0_74", offset = 0.00015, barsize = 0.45) + geom_hilight(node =2147 ,alpha = .5, fill="#EA8331") +
  geom_cladelabel(node=3150, color='black', label = "0_154", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 3150,alpha = .5, fill="#7CAE00") +
  geom_cladelabel(node=2386, color='black', label = "0_237", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 2386 ,alpha = .5, fill="#00C1A3") +
  geom_cladelabel(node=3201, color='black', label = "0_265", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 3201,alpha = .5, fill="#E76BF3") +
  geom_cladelabel(node=3259, color='black', label = "0_310", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 3259,alpha = .5, fill="#39B600") +
  geom_cladelabel(node=3433, color='black', label = "0_20", offset = 0.00015, barsize = 0.45) + geom_hilight(node = 3433,alpha = .5, fill="#00BFC4") 

````
## FIGURE 2B: 
```{r}
tx3 <- read.table("melted_ax_5k6_sub.csv", sep=",", header=TRUE)
fro <- (subset(tx3, match_GPS=="MATCH" & !is.na(GPS.x)) %>% filter(!grepl(';', GPS.x)))
tsvv <- read.table("txvv.csv", sep=",", header=TRUE)
tx3 <- melted_ax_5k6sub

nps <- read.table("ML_5k_test/x.tsv", header=TRUE, na.strings=c(" ","NA",""), comment.char = "", sep="\t")
nps_melted <- reshape2::melt(nps, id.vars = c("ID")) %>% na.omit()
nps_melted$C <- ifelse(grepl("PS", nps_melted$ID), "yes", "no")
nps_melted$D <- ifelse(grepl("PS", nps_melted$variable), "yes", "no")
nps_melted_diff <- subset(nps_melted, C!=D)
nps_melted_diff$comps <- ""
nps_melted$comps <- ""

all_tx3 <- tx3 %>% subset(is.na(match_ST) & is.na(match_GPS) & is.na(match_rST) & is.na(MANDRAKE)) %>% select(variable, ID, comps, value)
all_tsvv <- tsvv %>% select(variable, ID, comps, value)
all_fro <- fro %>% select(variable, ID, comps, value)
all_ST <- tx3 %>% subset(match_ST=="MATCH") %>% select(variable, ID, comps, value) 
all_rST <- tx3 %>% subset(match_rST=="MATCH") %>% select(variable, ID, comps, value) 
all_MAND <- tx3 %>% subset(MANDRAKE=="MATCH") %>% select(variable, ID, comps, value) 
nps_test <- nps_melted %>% subset(C=="yes" | D=="yes") %>% subset(C=="no" | D=="no")  %>% select(variable, ID, comps, value) 

all_tx3$TAG <- "1_ALL"
all_tsvv$TAG <- "2_CC"
all_fro$TAG <- "3_GPS"
all_ST$TAG <- "4_ST"
all_rST$TAG <- "5_rST"
all_MAND$TAG <- "35_MAND"
nps_test$TAG <- "0_NPS"
#all_tx3_fb2$TAG <- "x_fb"
merged_prelim_dens2 <- (rbind(all_tx3, all_tsvv, all_fro, all_ST, all_rST, all_MAND))

densities <- ggplot(merged_prelim_dens2, aes(x = value*100, y = TAG, height = stat(density), fill=TAG)) +
  geom_density_ridges(stat = "density",bw=0.00001, alpha=0.35,scale =3.5,color=NA) + theme_bw() +
  scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
  scale_y_discrete(expand=c(0,0)) +
  xlab("Pairwise allelic mismatches (%)") + ylab("") +
  scale_fill_manual(values=c("#023858","#194B6A","#305E7D","#487190","#5F84A2","#7797B5","#8EAAC8","#A6BDDB")) +
  theme( panel.grid.major =  element_blank(), plot.title = element_text(face="bold", color="black"),
         axis.line = element_line(colour = "black"),
         panel.border = element_blank(),panel.grid.minor =element_blank(),legend.title = element_blank(),
         legend.position="none",axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=10),
         axis.title.x = element_text(face="bold", color="black", size=10)) 
```
## Figure 2C:
```{r}
tx3x <- tx3 %>% sample_n(205000) %>% select(variable, ID, comps, value) 
nps_test$TAG <- "0_NPS"
all_tx3b$TAG <- "1_ALL"
all_tsvv$TAG <- "1x_ALL"
all_tx3d$TAG <- "WHO"
tx3x$TAG <- "other"
merged_tx3x <- (rbind(tx3x,  nps_test))

borders <- ggplot() + theme_bw() +
  geom_density(data=subset(merged_tx3x), aes(x=(value*100)), bw=0.00001, color="#5F84A2", fill="#5F84A2", alpha=1) + 
  scale_y_continuous(expand=c(0,0)) +
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
## Figure 2D: 
```{r}
# Import metadata
linmeta <- read.csv("lineage_meta.csv", header=TRUE, na.strings=c(""," ","NA"))

# Set LINcode strings
linmeta$gx <- paste(linmeta$L1, linmeta$L2, sep="_")

# Get counts by country and LINcode
country_count <- linmeta %>% subset(!is.na(country)) %>% group_by(gx) %>% select(gx,country) %>% unique() %>% summarise(countries=n())

# Get counts by serotype and LINcode
sero_count <-linmeta %>% subset(!is.na(serotype)) %>% subset(serotype!="inconclusive") %>% group_by(gx,serotype) %>% 
  #filter(n() > 2) %>% 
  ungroup %>% group_by(gx) %>% select(gx,serotype) %>% unique() %>% summarise(serotype2=n())

# Get counts by Lineage
all_count <- linmeta %>% select(id,gx) %>% group_by(gx) %>% summarise(count=n())

# Merge dataframes
testmergecg_LIN_country_serotype <- merge(country_count,sero_count, by.x=c("gx"), by.y=c("gx"), all=TRUE)
testmergecg_LIN_country_serotype_all <- merge(all_count,testmergecg_LIN_country_serotype, by.x=c("gx"), by.y=c("gx"), all=TRUE) %>% na.omit()  %>% subset(gx!="NA_NA")

# Plot
ggplot(data=(testmergecg_LIN_country_serotype_all)) + 
  geom_point(aes(y=countries,x=count, size=(serotype2)),shape=21,stroke=0,color="black",fill="grey50", alpha=0.5) + 
  theme_bw() + xlab("Number of isolates") + ylab("Number of countries") +
  scale_fill_gradientn(breaks = c(0,1,2,3,4,5,6),colors=c("#fed976","#fd8d3c","#b10026")) +
  scale_x_continuous(limits=c(0,2500), expand=c(0,0)) +
  scale_size_area(max_size = 6, limits=c(1,30)) +
  scale_y_continuous(limits=c(0,60), expand=c(0,0))  + 
  theme( panel.grid = element_blank(), legend.position = c(0.87, 0.45),
         axis.line = element_line(colour = "black"), panel.border = element_blank(),
         legend.title = element_blank(),axis.text.x=element_text(color="black", size=6),
         axis.text.y=element_text(color="black", size=6),axis.title.y = element_text(face="bold", color="black", size=10),
         axis.title.x = element_text(face="bold", color="black", size=10))
````
