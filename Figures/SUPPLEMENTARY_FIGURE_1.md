library(ggpubr)
loci_des <- read.csv("loci_description.csv", header=TRUE)


supl1x_A <- ggscatter(loci_des, x = "length", y = "count",add = "reg.line",shape = 19,color="#3182bd", size=0.75,
                      add.params = list(color = "black", fill = "lightgray"),
                      conf.int = TRUE) +
  stat_cor(method = "pearson", label.x = 3000, label.y = 100) +
  ylab("Number of unique loci") +
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

supl1x_B <- ggplot(data=loci_des) + 
  geom_histogram(aes(x=length),bins=80, size=0.3, color="black", fill="#3182bd") +
  xlab("Loci length (bp)") + theme_bw() + 
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

supl1x_C <- ggplot(data=loci_des) + 
  geom_histogram(aes(x=count),bins=80, size=0.3, color="black", fill="#3182bd") +
  xlab("Number of unique loci") + theme_bw() + 
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


(supl1x_A|supl1x_B)/(supl1x_C|void)

void <- ggplot() + theme_void()


## CC cluster size vs cluster diversity
######### Circlize

reps <- read.table("filt.blast.hits", header=FALSE)
eggnog <- read.table("eggnog_annotations_short.csv", sep=",", header=TRUE)
egg_rep <- merge(reps, eggnog, by.x=c("V4"), by.y=c("query"),all.x=TRUE)

svg('test1.svg', width = 5, height=4)
chr <- read.csv("ref.lengths.csv", header=FALSE) %>% subset(V1=="FM211187.1")
gc <- read.table("FM211187.1.ATCC700669.fasta.gc.bed", header=FALSE) %>% subset(V1=="FM211187.1")
reps2 <- subset(egg_rep, V1=="FM211187.1") %>% select(V2,V3,X) %>% subset(X=="Information storage and processing") 
reps3 <- subset(egg_rep, V1=="FM211187.1") %>% select(V2,V3,X) %>% subset(X=="Cellular processes and signaling") 
reps4 <- subset(egg_rep, V1=="FM211187.1") %>% select(V2,V3,X) %>% subset(X=="Metabolism") 
reps5 <- subset(egg_rep, V1=="FM211187.1") %>% select(V2,V3,X) %>% subset(X=="Poorly characterized") 
svg('test1.svg', width = 5, height=4)
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

svg('test2.svg', width = 5, height=4)
reps2 <- subset(egg_rep, V1=="AKVY01000001.1") %>% select(V2,V3,X) %>% subset(X=="Information storage and processing") 
reps3 <- subset(egg_rep, V1=="AKVY01000001.1") %>% select(V2,V3,X) %>% subset(X=="Cellular processes and signaling") 
reps4 <- subset(egg_rep, V1=="AKVY01000001.1") %>% select(V2,V3,X) %>% subset(X=="Metabolism") 
reps5 <- subset(egg_rep, V1=="AKVY01000001.1") %>% select(V2,V3,X) %>% subset(X=="Poorly characterized") 
chr <- read.csv("ref.lengths.csv", header=FALSE) %>% subset(V1=="AKVY01000001.1")
gc <- read.table("all.gc.bed", header=FALSE) %>% subset(V1=="AKVY01000001.1")
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
circos.yaxis(labels.cex=0.35, side = "left", tick = T, sector.index ="AKVY01000001.1",at = c(0.2,0.3,0.4,0.5,0.6) )
text(0, 0, "TIGR4", cex = 1, font=2)
dev.off()

svg('test3.svg', width = 5, height=4)
reps2 <- subset(egg_rep, V1=="CP035238.1") %>% select(V2,V3,X) %>% subset(X=="Information storage and processing") 
reps3 <- subset(egg_rep, V1=="CP035238.1") %>% select(V2,V3,X) %>% subset(X=="Cellular processes and signaling") 
reps4 <- subset(egg_rep, V1=="CP035238.1") %>% select(V2,V3,X) %>% subset(X=="Metabolism") 
reps5 <- subset(egg_rep, V1=="CP035238.1") %>% select(V2,V3,X) %>% subset(X=="Poorly characterized") 
chr <- read.csv("ref.lengths.csv", header=FALSE) %>% subset(V1=="CP035238.1")
gc <- read.table("all.gc.bed", header=FALSE) %>% subset(V1=="CP035238.1")
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
circos.yaxis(labels.cex=0.35, side = "left", tick = T, sector.index ="CP035238.1",at = c(0.2,0.3,0.4,0.5,0.6) )
text(0, 0, "D39", cex = 1, font=2)
dev.off()



