# Supplementary Figure 3
```{r}
# Load libraries
library(circlize)

# Import data
reps <- read.table("filt.blast.hits", header=FALSE)
eggnog <- read.table("eggnog_annotations_short.csv", sep=",", header=TRUE)
egg_rep <- merge(reps, eggnog, by.x=c("V4"), by.y=c("query"),all.x=TRUE)

# Subset data for each reference genome
chr <- read.csv("ref.lengths.csv", header=FALSE) %>% subset(V1=="FM211187.1")
gc <- read.table("FM211187.1.ATCC700669.fasta.gc.bed", header=FALSE) %>% subset(V1=="FM211187.1")
reps2 <- subset(egg_rep, V1=="FM211187.1") %>% select(V2,V3,X) %>% subset(X=="Information storage and processing") 
reps3 <- subset(egg_rep, V1=="FM211187.1") %>% select(V2,V3,X) %>% subset(X=="Cellular processes and signaling") 
reps4 <- subset(egg_rep, V1=="FM211187.1") %>% select(V2,V3,X) %>% subset(X=="Metabolism") 
reps5 <- subset(egg_rep, V1=="FM211187.1") %>% select(V2,V3,X) %>% subset(X=="Poorly characterized") 

# Plot figure
svg('FM211187.svg', width = 5, height=4)
circos.clear()
circos.par("track.height" = 0.15, start.degree = 90, gap.degree=c(20))
circos.genomicInitialize(data = chr, major.by=200000, sector.names = NA)
circos.track(ylim = c(0,1))
circos.genomicRect(region = reps2, ytop = 1.1, ybottom = 0.80, border=NA, fill=c("#872154"), col=c("#872154"))
circos.genomicRect(region = reps3, ytop = 0.80, ybottom = 0.50, border=NA, fill=c("#42a899"), col=c("#42a899"))
circos.genomicRect(region = reps4, ytop = 0.50, ybottom = 0.2, border=NA, fill=c("#87ccee"), col=c("#87ccee"))
circos.genomicRect(region = reps5, ytop = 0.2, ybottom = -0.1, border=NA, fill=c("grey50"), col=c("grey50"))
circos.track(ylim = c(0.2,  0.6), factors = gc$V1,  x=gc$V2, y=gc$V5)
circos.trackLines(x=gc$V2, y=gc$V5, sectors = gc$V1, col=c("black"))
circos.yaxis(labels.cex=0.35, side = "left", tick = T, sector.index ="FM211187.1",at = c(0.2,0.3,0.4,0.5,0.6) )
text(0, 0, "ATCC 700669", cex = 1, font=2)
dev.off()

# Subset data for each reference genome
reps2 <- subset(egg_rep, V1=="AKVY01000001.1") %>% select(V2,V3,X) %>% subset(X=="Information storage and processing") 
reps3 <- subset(egg_rep, V1=="AKVY01000001.1") %>% select(V2,V3,X) %>% subset(X=="Cellular processes and signaling") 
reps4 <- subset(egg_rep, V1=="AKVY01000001.1") %>% select(V2,V3,X) %>% subset(X=="Metabolism") 
reps5 <- subset(egg_rep, V1=="AKVY01000001.1") %>% select(V2,V3,X) %>% subset(X=="Poorly characterized") 
chr <- read.csv("ref.lengths.csv", header=FALSE) %>% subset(V1=="AKVY01000001.1")
gc <- read.table("all.gc.bed", header=FALSE) %>% subset(V1=="AKVY01000001.1")
circos.clear()

# Plot figure
svg('AKVY01000001.svg', width = 5, height=4)
circos.par("track.height" = 0.15, start.degree = 90, gap.degree=c(20))
circos.genomicInitialize(data = chr, major.by=200000, sector.names = NA)
circos.track(ylim = c(0,1))
circos.genomicRect(region = reps2, ytop = 1.1, ybottom = 0.80, border=NA, fill=c("#872154"), col=c("#872154"))
circos.genomicRect(region = reps3, ytop = 0.80, ybottom = 0.50, border=NA, fill=c("#42a899"), col=c("#42a899"))
circos.genomicRect(region = reps4, ytop = 0.50, ybottom = 0.2, border=NA, fill=c("#87ccee"), col=c("#87ccee"))
circos.genomicRect(region = reps5, ytop = 0.2, ybottom = -0.1, border=NA, fill=c("grey50"), col=c("grey50"))
circos.track(ylim = c(0.2,  0.6), factors = gc$V1,  x=gc$V2, y=gc$V5)
circos.trackLines(x=gc$V2, y=gc$V5, sectors = gc$V1, col=c("black"))
circos.yaxis(labels.cex=0.35, side = "left", tick = T, sector.index ="AKVY01000001.1",at = c(0.2,0.3,0.4,0.5,0.6) )
text(0, 0, "TIGR4", cex = 1, font=2)
dev.off()

# Subset data for each reference genome
svg('CP035238.svg', width = 5, height=4)
reps2 <- subset(egg_rep, V1=="CP035238.1") %>% select(V2,V3,X) %>% subset(X=="Information storage and processing") 
reps3 <- subset(egg_rep, V1=="CP035238.1") %>% select(V2,V3,X) %>% subset(X=="Cellular processes and signaling") 
reps4 <- subset(egg_rep, V1=="CP035238.1") %>% select(V2,V3,X) %>% subset(X=="Metabolism") 
reps5 <- subset(egg_rep, V1=="CP035238.1") %>% select(V2,V3,X) %>% subset(X=="Poorly characterized") 
chr <- read.csv("ref.lengths.csv", header=FALSE) %>% subset(V1=="CP035238.1")
gc <- read.table("all.gc.bed", header=FALSE) %>% subset(V1=="CP035238.1")
circos.clear()

# Plot figure
circos.par("track.height" = 0.15, start.degree = 90, gap.degree=c(20))
circos.genomicInitialize(data = chr, major.by=200000, sector.names = NA)
circos.track(ylim = c(0,1))
circos.genomicRect(region = reps2, ytop = 1.1, ybottom = 0.80, border=NA, fill=c("#872154"), col=c("#872154"))
circos.genomicRect(region = reps3, ytop = 0.80, ybottom = 0.50, border=NA, fill=c("#42a899"), col=c("#42a899"))
circos.genomicRect(region = reps4, ytop = 0.50, ybottom = 0.2, border=NA, fill=c("#87ccee"), col=c("#87ccee"))
circos.genomicRect(region = reps5, ytop = 0.2, ybottom = -0.1, border=NA, fill=c("grey50"), col=c("grey50"))
circos.track(ylim = c(0.2,  0.6), factors = gc$V1,  x=gc$V2, y=gc$V5)
circos.trackLines(x=gc$V2, y=gc$V5, sectors = gc$V1, col=c("black"))
circos.yaxis(labels.cex=0.35, side = "left", tick = T, sector.index ="CP035238.1",at = c(0.2,0.3,0.4,0.5,0.6) )
text(0, 0, "D39", cex = 1, font=2)
dev.off()
```
