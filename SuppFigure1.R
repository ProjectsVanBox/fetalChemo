### 1A plot indels per treatment group 
load("allIndels.Rdata")

mutationCountHealthy = unlist(allIndels[17:26])
mutationCountUntreated = unlist(allIndels[which(annot$treated == "No")]) 
mutationCountPlatin  = unlist(allIndels[which(annot$treatmentGroup == "Carboplatin")]) 
mutationCountnonplatin = unlist(allIndels[which(annot$treatmentGroup == "non-platinum")]) 
mutationCountABVD = unlist(allIndels[which(annot$treatmentGroup == "ABVD")]) 

counts = c(mutationCountHealthy, mutationCountUntreated, mutationCountPlatin, mutationCountnonplatin, mutationCountABVD)
label = c(rep("Healthy", length(mutationCountHealthy)),rep("Untreated",length(mutationCountUntreated)),rep("(EC +) TCpt", length(mutationCountPlatin)), rep("EC", length(mutationCountnonplatin)), rep("ABVD", length(mutationCountABVD)))

treatmentColor = c(rep("#00BA38", length(mutationCountHealthy)),rep("darkolivegreen",length(mutationCountUntreated)),rep("#F8766D", length(mutationCountPlatin)), rep("dodgerblue", length(mutationCountnonplatin)), rep("darkorchid4", length(mutationCountABVD)))

samples = NULL
for(i in 1:length(counts)){
  samples[i] = strsplit(strsplit(names(counts)[i], "[.]")[[1]][1], "_")[[1]][1]
}


dataToPlot = data.frame(counts, samples, label)


level_order = sampleOrderAll

pdf("boxplot_persample_indels.pdf", width = 5, height = 6)
ggplot(dataToPlot, aes(x=factor(samples, level = level_order), y=counts, fill = treatmentColor, alpha = 0.5)) + ylab("SNV count") + xlab(NA)+
  geom_boxplot() + geom_jitter(aes(colour =treatmentColor, alpha =1))  +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 12, angle = 90)) 

dev.off(
  
  
)
### 1B plot per treatment group

level_order = factor(c("Healthy", "Untreated", "EC", "ABVD", "(EC +) TCpt"))

pdf("boxplot_pertreatment_indels.pdf", width = 5, height = 6)

ggplot(dataToPlot, aes(x=factor(label, level = level_order), y=counts, fill = treatmentColor, alpha = 0.5)) + ylab("SNV count") + xlab(NA)+
  geom_boxplot() + geom_jitter(aes(colour =treatmentColor, alpha =1))  +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 12, angle = 90)) 
dev.off()

### 1C karyotype plots are produced using FREEC