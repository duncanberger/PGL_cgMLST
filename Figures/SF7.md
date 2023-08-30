# Supplementary Figure 7
```{r}
# Load libraries 
library(ggplot2)

# Import HierCC results
hcc <- read.table("hcc.txt", header=FALSE, sep=",")

# Plot
p_hcc <- ggplot(data=hcc) + geom_line(aes(x=V1, y=V2), size=0.75, color="darkred") + theme_bw() +
  xlab("Pairwise allelic mismatches") + ylab("Silhouette score") +
  scale_fill_manual(values=c("#009688"), na.value = "grey75") +
  scale_x_continuous(expand=c(0,0), limits=c(0,1222)) +
  geom_vline(xintercept = 1180, color="grey60", linetype="dashed") + 
  geom_vline(xintercept = 750, color="grey60", linetype="dashed") + 
  geom_vline(xintercept = 540, color="grey60", linetype="dashed") +
  geom_vline(xintercept = 160, color="grey60", linetype="dashed") +
  geom_vline(xintercept = 25, color="grey60", linetype="dashed") +
  geom_vline(xintercept = 15, color="grey60", linetype="dashed") +
  geom_vline(xintercept = 8, color="grey60", linetype="dashed") +
  geom_vline(xintercept = 4, color="grey60", linetype="dashed") + 
  geom_vline(xintercept = 2, color="grey60", linetype="dashed") +
  geom_vline(xintercept = 1, color="grey60", linetype="dashed") +
  theme( panel.grid = element_blank(), plot.title = element_text(face="bold", color="black"),
         legend.title = element_blank(),legend.position="none",axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=10),axis.title.x = element_text(face="bold", color="black", size=10)) 
```
