library(ggsci)
library(ggplot2)
library(MutationalPatterns)

existing = get_known_signatures(incl_poss_artifacts = F, genome = "GRCh38")

#### load HSPC signature, available on https://github.com/ToolsVanBox/MutationalPatterns
hspc = read.table("HSPC.txt", header = T)

existing_hspc = cbind(existing, hspc$HSPC)
colnames(existing_hspc) = c(colnames(existing), "hspc")

#### load mutational matrices for all samples
load("allMutMats.Rdata")

#### mutation count per sample 
Muts482 = apply(allMutMats[[1]], 1, sum)
Muts483 = apply(allMutMats[[2]], 1, sum)
Muts486 = apply(allMutMats[[3]], 1, sum)
Muts489 = apply(allMutMats[[4]], 1, sum)
Muts495 = apply(allMutMats[[5]], 1, sum)
Muts497 = apply(allMutMats[[6]], 1, sum)
Muts501 = apply(allMutMats[[7]], 1, sum)
Muts504 = apply(allMutMats[[8]], 1, sum)
Muts494 = apply(allMutMats[[9]], 1, sum)
Muts498 = apply(allMutMats[[10]], 1, sum)

Muts493 = apply(allMutMats[[11]], 1, sum)
Muts515 = apply(allMutMats[[12]], 1, sum)
Muts516 = apply(allMutMats[[13]], 1, sum)
Muts517 = apply(allMutMats[[14]], 1, sum)
Muts518 = apply(allMutMats[[15]], 1, sum)
Muts519 = apply(allMutMats[[16]], 1, sum)

Muts9 = apply(allMutMats[[17]], 1, sum)
Muts10 = apply(allMutMats[[18]], 1, sum)
Muts11 = apply(allMutMats[[19]], 1, sum)
Muts12 = apply(allMutMats[[20]], 1, sum)
Muts13 = apply(allMutMats[[21]], 1, sum)
Muts17 = apply(allMutMats[[22]], 1, sum)

Muts18 = apply(allMutMats[[23]], 1, sum)
Muts19 = apply(allMutMats[[24]], 1, sum)

Muts20 = apply(allMutMats[[25]], 1, sum)
Muts21 = apply(allMutMats[[26]], 1, sum)

mutMatsSummed = cbind(Muts482,
Muts483 ,
Muts486 ,
Muts489,
Muts495 ,
Muts497,
Muts501 ,
Muts504,
 
Muts494,
Muts498,

Muts493 ,
Muts515,
Muts516,
Muts517,
Muts518,
Muts519,
Muts9,
Muts10,
Muts11,
Muts12,
Muts13,
Muts17,
Muts18,
Muts19,
Muts20,
Muts21)

### load metadata

load("sampleAnnotation.Rdata")

### counts per sample 


#####1b, SNV counts
### make mutation count for ggplot 

mutationCountTreated  = unlist(allMuts[which(annot$treated == "Yes")]) 
mutationCountUntreated = unlist(allMuts[which(annot$treated == "No")]) 
mutationCountHealthy = unlist(allMuts[17:26])

counts = c(mutationCountTreated, mutationCountUntreated, mutationCountHealthy)
label = c(rep("Treated", length(mutationCountTreated)), rep("Untreated", length(mutationCountUntreated)), rep("Healthy", length(mutationCountHealthy)))
dataToPlot = data.frame(counts, label)

level_order = factor(c("Healthy", "Untreated", "Treated"))


pdf("Figure1/MutationCount_3groups.pdf", width = 5, height = 6)
ggplot(dataToPlot, aes(x=factor(label, level = level_order), y=counts)) + ylab("SNV count") + xlab(NA)+
  geom_boxplot() + geom_jitter(aes(color = label),size = 0.6, alpha =1)+  scale_fill_npg(alpha=0.3)  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 12, angle = 0)) + scale_color_npg()

dev.off()


###### 1c, divide treated groups 

mutationCountPlatin  = unlist(allMuts[which(annot$treatmentGroup == "Carboplatin")]) 
mutationCountnonplatin = unlist(allMuts[which(annot$treatmentGroup == "EC")]) 
mutationCountABVD = unlist(allMuts[which(annot$treatmentGroup == "ABVD")]) 

counts = c(mutationCountHealthy, mutationCountUntreated, mutationCountPlatin, mutationCountnonplatin, mutationCountABVD)
label = c(rep("Healthy", length(mutationCountHealthy)),rep("Untreated",length(mutationCountUntreated)),rep("Carboplatin", length(mutationCountPlatin)), rep("EC", length(mutationCountnonplatin)), rep("ABVD", length(mutationCountABVD)))

treatmentColor = c(rep("#00BA38", length(mutationCountHealthy)),rep("darkolivegreen",length(mutationCountUntreated)),rep("#F8766D", length(mutationCountPlatin)), rep("dodgerblue", length(mutationCountnonplatin)), rep("darkorchid4", length(mutationCountABVD)))

samples = NULL
for(i in 1:length(counts)){
  samples[i] = strsplit(strsplit(names(counts)[i], "[.]")[[1]][1], "_")[[1]][1]
}


dataToPlot = data.frame(counts, samples, label)

level_order = factor(c("Healthy", "Untreated", "EC", "ABVD", "Carboplatin"))

pdf("Figure1/boxplot_perTreatmentGroup.pdf", width = 5, height = 6)

ggplot(dataToPlot, aes(x=factor(label, level = level_order), y=counts, fill = treatmentColor, alpha = 0.5)) + ylab("SNV count") + xlab(NA)+
  geom_boxplot() + geom_jitter(aes(colour =treatmentColor, alpha =1))  +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 12, angle = 90)) 
dev.off()


#### 1d, regression with between snv count and administrations

annot$meanSNVcount = rep(NA, nrow(annot))
for(i in 1:nrow(annot)){
  annot$meanSNVcount[i] = mean(allMuts[[i]])
  annot$lowerBound[i] = min(allMuts[[i]])
  annot$upperBound[i] = max(allMuts[[i]])
}

HealthyMutCounts = NULL
for(i in 1:10){
  HealthyMutCounts[i] = mean(allMuts[[i+15]])
}


treatedAnnot = annot[which(annot$treated == "Yes"),]


datacycles = data.frame(treatedAnnot$meanSNVcount, treatedAnnot$lowerBound, treatedAnnot$upperBound, treatedAnnot$cycles, treatedAnnot$treatmentGroup)
colnames(datacycles) = c("meanCount", "lower", "upper", "cycles", "treatment")
colorToPlot = rep("darkolivegreen3", nrow(datacycles))
colorToPlot[which(datacycles$treatment == "EC")] = "darkorchid2"
colorToPlot[which(datacycles$treatment == "ABVD")] = "dodgerblue"


pdf("Figure1/Regression_cycles.pdf")
ggplot(datacycles, aes(y = meanCount, x = cycles, color = colorToPlot)) + 
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper, color = colorToPlot), color = colorToPlot) +  
  stat_smooth(formula = y ~ x, geom = "smooth", method = 'lm', colour = "azure4") +  
  ylab("SNV count") + xlab("Administrations") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 16, angle = 0), legend.position = "top") 
dev.off()

#### 1e, regression with washout 

dataWashout = data.frame(treatedAnnot$meanSNVcount, treatedAnnot$lowerBound, treatedAnnot$upperBound, treatedAnnot$washout, treatedAnnot$treatmentGroup)
colnames(dataWashout) = c("meanCount", "lower", "upper", "washout", "treatment")
colorToPlot = rep("darkolivegreen3", nrow(dataWashout))
colorToPlot[which(dataWashout$treatment == "EC")] = "darkorchid2"
colorToPlot[which(dataWashout$treatment == "ABVD")] = "dodgerblue"


pdf("Figure1/Regression_washout.pdf")
ggplot(dataWashout, aes(y = meanCount, x = washout, color = colorToPlot)) + 
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper, color = colorToPlot), color = colorToPlot) +  
  stat_smooth(formula = y ~ x, geom = "smooth", method = 'lm', colour = "azure4") +  
  ylab("SNV count") + xlab("Washout (days)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 16, angle = 0), legend.position = "top") 
dev.off()

