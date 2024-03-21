library(MutationalPatterns)

existing = get_known_signatures(incl_poss_artifacts = F, genome = "GRCh38")

#### load HSPC signature, available on https://github.com/ToolsVanBox/MutationalPatterns
hspc = read.table("HSPC.txt", header = T)

existing_hspc = cbind(existing, hspc$HSPC)
colnames(existing_hspc) = c(colnames(existing), "hspc")


#### Figure 2a, plot 96 channel plots 
nonplatinumMuts = mutMatsSummed[,which(annot$treatmentGroup == "EC")]
platinumMuts =mutMatsSummed[,which(annot$treatmentGroup == "Carboplatin")]
ABVDMuts = mutMatsSummed[,which(annot$treatmentGroup == "ABVD")]
untreatedMuts =  mutMatsSummed[,which(annot$treated == "No")]
healthyMuts =  mutMatsSummed[,18:26]

forRefit = cbind(healthyMuts, untreatedMuts,  nonplatinumMuts,platinumMuts, ABVDMuts )
colnames(forRefit) = c("Healthy","Unexposed", "Non-platinum", "Carboplatin", "ABVD")

pdf("96_channel_treatments.pdf", width =12)
plot_96_profile(forRefit, ymax = 0.1)
dev.off()



### Figure 2b, mutational signatures 


getSig_treatment = fit_to_signatures(forRefit, existing_hspc[,c(1,5,37,46,61)])

pdf("Signature_refit.pdf", width =12)
plot_contribution(getSig_treatment$contribution)
dev.off()


### Figure 2c
load("allDBS.Rdata")


dbs_platinum =allDBS[,which(annot$treatmentGroup == "Carboplatin")]
dbs_untreated =  allDBS[,which(annot$treated == "No")]

pdf("DBS_occurences.pdf")
plot_main_dbs_contexts(dbs_platinum)
dev.off()



### Figure 2d strand bias

load("MutMats_stranded.Rdata")

### make mutational matrix

Muts_stranded482 = apply(allMutMats_stranded[[1]], 1, sum)
Muts_stranded483 = apply(allMutMats_stranded[[2]], 1, sum)
Muts_stranded486 = apply(allMutMats_stranded[[3]], 1, sum)
Muts_stranded489 = apply(allMutMats_stranded[[4]], 1, sum)
Muts_stranded495 = apply(allMutMats_stranded[[5]], 1, sum)
Muts_stranded497 = apply(allMutMats_stranded[[6]], 1, sum)
Muts_stranded501 = apply(allMutMats_stranded[[7]], 1, sum)
Muts_stranded504 = apply(allMutMats_stranded[[8]], 1, sum)
Muts_stranded494 = apply(allMutMats_stranded[[9]], 1, sum)
Muts_stranded498 = apply(allMutMats_stranded[[10]], 1, sum)

Muts_stranded493 = apply(allMutMats_stranded[[11]], 1, sum)
Muts_stranded515 = apply(allMutMats_stranded[[12]], 1, sum)
Muts_stranded516 = apply(allMutMats_stranded[[13]], 1, sum)
Muts_stranded517 = apply(allMutMats_stranded[[14]], 1, sum)
Muts_stranded518 = apply(allMutMats_stranded[[15]], 1, sum)
Muts_stranded519 = apply(allMutMats_stranded[[16]], 1, sum)

Muts_stranded9 = apply(allMutMats_stranded[[17]], 1, sum)
Muts_stranded10 = apply(allMutMats_stranded[[18]], 1, sum)
Muts_stranded11 = apply(allMutMats_stranded[[19]], 1, sum)
Muts_stranded12 = apply(allMutMats_stranded[[20]], 1, sum)
Muts_stranded13 = apply(allMutMats_stranded[[21]], 1, sum)
Muts_stranded17 = apply(allMutMats_stranded[[22]], 1, sum)

Muts_stranded18 = apply(allMutMats_stranded[[23]], 1, sum)
Muts_stranded19 = apply(allMutMats_stranded[[24]], 1, sum)

Muts_stranded20 = apply(allMutMats_stranded[[25]], 1, sum)
Muts_stranded21 = apply(allMutMats_stranded[[26]], 1, sum)

mutMatsSummed_stranded = cbind(Muts_stranded482,
Muts_stranded483 ,
Muts_stranded486 ,
Muts_stranded489,
Muts_stranded495 ,
Muts_stranded497,
Muts_stranded501 ,
Muts_stranded504,
 
Muts_stranded494,
Muts_stranded498,

Muts_stranded493 ,
Muts_stranded515,
Muts_stranded516,
Muts_stranded517,
Muts_stranded518,
Muts_stranded519,
Muts_stranded9,
Muts_stranded10,
Muts_stranded11,
Muts_stranded12,
Muts_stranded13,
Muts_stranded17,
Muts_stranded18,
Muts_stranded19,
Muts_stranded20,
Muts_stranded21)

#### test for strand bias

platinumMuts_stranded =mutMatsSummed_stranded[,which(annot$treatmentGroup == "Carboplatin")]
untreatedMuts_stranded =  mutMatsSummed_stranded[,which(annot$treated == "No")]


strandBias = strand_occurences(cbind(platinumMuts_stranded, untreatedMuts_stranded))

strand_bias_test(strandBias)

pdf("StrandBias_occurence.pdf")
plot_strand(strandBias)
dev.off()

