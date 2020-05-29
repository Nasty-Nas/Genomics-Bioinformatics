library(DESeq2)
library(FactoMineR)
library("factoextra")

final_2 <- read.csv("data/final_2_CMP.CSV", row.names = 1)
coldata <- read.csv("data/samples.CSV")
final_2_mat <- data.matrix(final_2, rownames.force = NA)
vsd_2<- vst(round(final_2_mat))
res.pca_2 <- PCA(t(vsd_2), graph = FALSE)

fviz_pca_ind(res.pca_2, geom="point", alpha=I(0)) +
  geom_point(aes(shape=coldata$species, colour=coldata$organ)) + #fill
  scale_shape_manual(values = c("Human"=15, "Mouse"=16,"Opossum"=17, "Rabbit"=18, "Chicken"=8, "Rhesus"=22, "Rat"=21)) +
  scale_colour_manual(breaks = c("Brain","Cerebellum", "Heart", "Kidney", "Liver", "Ovary", "Testis"),
                     values=c("#0494c6","#00cfff", "#e10000", "#db9600", "#009c00", "#e6009d", "#ff5900")) +
  scale_y_continuous(breaks=c(-50,0, 50), limits=c(-55, 55)) + scale_x_continuous(breaks=c(-40,0, 40)) +
  labs(color="Organ", shape="Species", x="PC1 (17.6% variance explained)", y="PC2 (11.4% variance explained)",
       title="Global PCA")


#### Fig 1.c ####
library(ggplot2)
library(reshape2)

corr_brain <- read.csv("data/correlation_brain.csv", row.names = 1)
corr_brain$time = c(1:22)

d <- melt(corr_brain, id.vars="time")

corr_brain_tf <- read.csv("data/correlation_brain_tf.csv", row.names = 1)
corr_brain_tf$time = c(1:22)
d_tf <- melt(corr_brain_tf, id.vars="time")


ggplot(d, aes(time,value, col=variable)) + ylim(0.7,1) + 
  geom_point() + 
  geom_smooth(se = FALSE) +
  geom_point(data = d_tf, shape=17) + 
  geom_smooth(data= d_tf, se = FALSE, linetype = "dashed") +
  scale_color_manual(breaks = c("Cerebellum", "Heart", "Kidney", "Liver", "Ovary", "Testis"),
                     values=c("#00cfff", "#e10000", "#db9600", "#009c00", "#e6009d", "#ff5900")) +
  labs(title ="Human developmental expression similarity to the brain",
       x ="Developmental stages", y="Transcriptome similarity (Spearman's correlation)") +
  coord_fixed(ratio=100) +
  scale_x_continuous(breaks=c(1,5,9,14,19,22), labels=c("4w","9w","13w","new","ya","sen")) +
  labs(color="Organ")

#### Differential expression ####
library(limma)
library(edgeR)
hum_vs_rhesus <- read.csv("data/hum_vs_rhesus.csv", row.names = 1) #Gene ID orthologues imported form python notebook
hum_vs_rhesus = subset(hum_vs_rhesus, select = -c(Rhesus_ID) )

library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)
library(gplots)   # contains the heatmap.2 package
library(tidyr)
library(stringi)

library(limma) 
library(edgeR) 

# RNAseq data:

human <- fread("data/human/Human.CPM.txt")
rhesus <- fread("data/macaque/Rhesus.CPM.txt")

hum_vs_rhesus <- fread("data/hum_vs_rhesus.csv")

hum_vs_rhesus$Brain.newborn <- NULL
hum_vs_rhesus$Brain.P0 <- NULL

rhesus <- merge(rhesus, hum_vs_rhesus, by.x = "V1",by.y = "Rhesus_ID",sort = TRUE)
human <- merge(human, hum_vs_rhesus, by.x = "V1",by.y = "Gene stable ID",sort = TRUE)

colnames(human) <- paste("human", colnames(human), sep = "_")
colnames(rhesus) <- paste("rhesus", colnames(rhesus), sep = "_")

all_data <- merge(human,rhesus, by.x = 'human_V1',by.y ='rhesus_Gene stable ID',sort = TRUE)
all_data <- all_data %>% select(human_V1, contains(c("newborn", "P0"))) #for DE at birth only

human_samples <- colnames(human)[2:296]
rhesus_samples <- colnames(rhesus)[2:165]

brain_data <- all_data[,(grepl( "Brain" , names( all_data ) ))|(grepl( "human_V1" , names( all_data ) )),with=FALSE]

rnaseq_mat <- as.matrix(brain_data[,grepl( "Brain" , names( brain_data ) ),with=FALSE]) 

row.names(rnaseq_mat) <- brain_data$human_V1

colnames(rnaseq_mat) 

sample_types <- data.table(colnames(rnaseq_mat))

sample_types[sample_types$V1 %in% human_samples ,"type"] <- "human"
sample_types[sample_types$V1 %in% rhesus_samples ,"type"] <- "rhesus"


x <- DGEList(counts=rnaseq_mat,samples=sample_types, group=sample_types$type )
class(x)

dim(x)

samplenames <- colnames(x)

group <- x$samples$group 

#Removing genes that are lowly expressed
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

#Creating a design matrix and contrasts
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design

contr.matrix <- makeContrasts(
  HvsR = human-rhesus,
  levels = colnames(design))
contr.matrix

#Removing heteroscedascity from count data
v <- voom(x, design, plot=TRUE) 
v

#Fitting linear models for comparisons of interest
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit) 


#Examining the number of DE genes
summary(decideTests(efit))
# HvsR
# Down   3105
# NotSig 2405
# Up     3693

#HvsR at birth
#Down   1199
#NotSig 5626
#Up     1566

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
write.fit(tfit, dt, file="DE.txt") #export results for enrichment
# HvsR
# Down    701
# NotSig 7163
# Up     1339

#HvsR at birth
#Down    293
#NotSig 7613
#Up      485

#Examining individual DE genes from top to bottom
h.vs.r <- topTreat(tfit, coef=1, n=Inf)
head(h.vs.r)

#Useful graphical representations of differential expression results
plotMD(tfit, column=1, status=dt[,1], main="Human vs. Rhesus differentially expressed genes", xlim=c(0,13))
