### 2A plot SNVs per sample group 
load("allMuts.Rdata")

mutationCountHealthy = unlist(allMuts[17:26])
mutationCountUntreated = unlist(allMuts[which(annot$treated == "No")]) 
mutationCountPlatin  = unlist(allMuts[which(annot$treatmentGroup == "Carboplatin")]) 
mutationCountnonplatin = unlist(allMuts[which(annot$treatmentGroup == "non-platinum")]) 
mutationCountABVD = unlist(allMuts[which(annot$treatmentGroup == "ABVD")]) 

counts = c(mutationCountHealthy, mutationCountUntreated, mutationCountPlatin, mutationCountnonplatin, mutationCountABVD)
label = c(rep("Healthy", length(mutationCountHealthy)),rep("Untreated",length(mutationCountUntreated)),rep("(EC +) TCpt", length(mutationCountPlatin)), rep("EC", length(mutationCountnonplatin)), rep("ABVD", length(mutationCountABVD)))

treatmentColor = c(rep("#00BA38", length(mutationCountHealthy)),rep("darkolivegreen",length(mutationCountUntreated)),rep("#F8766D", length(mutationCountPlatin)), rep("dodgerblue", length(mutationCountnonplatin)), rep("darkorchid4", length(mutationCountABVD)))

samples = NULL
for(i in 1:length(counts)){
  samples[i] = strsplit(strsplit(names(counts)[i], "[.]")[[1]][1], "_")[[1]][1]
}


dataToPlot = data.frame(counts, samples, label)


level_order = sampleOrderAll

pdf("boxplot_persample_snvs.pdf", width = 5, height = 6)
ggplot(dataToPlot, aes(x=factor(samples, level = level_order), y=counts, fill = treatmentColor, alpha = 0.5)) + ylab("SNV count") + xlab(NA)+
  geom_boxplot() + geom_jitter(aes(colour =treatmentColor, alpha =1))  +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 12, angle = 90)) 

dev.off()


### 2b cosine simalarity between groups 

cosSim = cos_sim_matrix(forRefit, forRefit)

pdf("Cosine_simalarity.pdf")
plot_cosine_heatmap(cosSim)
dev.off()

### de novo extraction of signatures 
## mutMatExisting - healthy reference data 


allMutMats_forRefit = cbind(mutMatsSummed, mutMatExisting)
set.seed(42)
allMutMats_forRefit = allMutMats_forRefit + 0.00001
estimate <- nmf(allMutMats_forRefit, rank = 2:10, method = "brunet", 
                nrun = 10, seed = 42, .opt = "v-p")

nmf_res <- extract_signatures(allMutMats_forRefit, rank = 5, nrun = 50, single_core = TRUE)

### 2c fit to extracted signatures 

fitToExtracted = fit_to_signatures(forRefit, nmf_res$signatures)

pdf("Signature_contribution_extracted.pdf")
plot_contribution(fitToExtracted$contribution, mode = "relative")
dev.off()

### 2d plot cosine similarity with existing signatures 
cosSimSignatures = cos_sim_matrix(nmf_res$signatures, existing_hspc)
getHighSimilarity = which(cosSimSignatures > 0.8, arr.ind = T)

selectedSigs = colnames(existing_hspc)[getHighSimilarity[,2]]

pdf("CosineSimilarity_existingSignatures.pdf")
plot_cosine_heatmap(cosSimSignatures[,selectedSigs])
dev.off()

###2e bootstrapped refit 
bootstrappedRefit = fit_to_signatures_bootstrapped(forRefit, existing_hspc[,c(1,5,37,46,61)])

pdf("bootstrappedRefit_selectedSigs.pdf")
plot_bootstrapped_contribution(bootstrappedRefit, mode = "relative", plot_type = "dotplot")
dev.off()










