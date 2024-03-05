library(phyloseq)
library(NetCoMi)
library(genefilter)
library(tidyverse)
library(GGally)
library(igraph)
library(network)
library(qgraph)

#My Data
physeq.tree <- readRDS("Data/Objetos/physeq.tree.rds")
physeq <- tax_glom(physeq.tree, "Species")

#Keep taxa with relative abundance greater than 0.05% in at least 1 sample
# k = 1; A = 0.0005
flist <- filterfun(kOverA(1, 0.05))
# Filter taxa
physeq.fil <- filter_taxa(physeq, flist)
physeq.species <- prune_taxa(physeq.fil, physeq)

#Construction of the two networks
#Cluster 1 
physeqC1 <- subset_samples(physeq.species, Cluster=="1")

#Filter to eliminate zeros
flistC <- filterfun(kOverA(1, 0))
physeq.filC1 <- filter_taxa(physeqC1, flistC)
physeqC1_fil <- prune_taxa(physeq.filC1, physeqC1)

#Cluster 2 
physeqC2 <- subset_samples(physeq.species, Cluster=="2")
#Filter to eliminate zeros
flistC <- filterfun(kOverA(1, 0))
physeq.filC2<- filter_taxa(physeqC2, flistC)
physeqC2_fil <- prune_taxa(physeq.filC2, physeqC2)


#Bug's network Cluster 1
netC1_fil <- netConstruct(data = physeqC1_fil,data2=NULL, dataType = "counts", group = NULL,
                       matchDesign = NULL, measure = "spieceasi", measurePar = NULL,
                       jointPrepro = NULL, filtTax = "none", filtTaxPar = NULL,
                       filtSamp = "none", filtSampPar = NULL,  zeroMethod = "none",
                       zeroPar = NULL, normMethod = "none", normPar = NULL, sparsMethod = "t-test",
                       thresh = 0.3, alpha = 0.05, adjust = "adaptBH", trueNullMethod = "convest",
                       lfdrThresh = 0.2, nboot = 1000L,  cores = 8,  logFile = "log.txt",
                       softThreshType = "signed",  softThreshPower = NULL, softThreshCut = 0.8,
                       kNeighbor = 3L, knnMutual = FALSE,  dissFunc = "signed",  dissFuncPar = NULL,
                       simFunc = NULL, simFuncPar = NULL,  scaleDiss = TRUE, weighted = TRUE,
                       sampleSize = NULL,  verbose = 2, seed = NULL)

props_clusC1.fil <- netAnalyze(netC1_fil, 
                           centrLCC = FALSE,
                           avDissIgnoreInf = TRUE,
                           sPathNorm = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = c("degree", "between"),
                           hubQuant = 0.90,
                           lnormFit = TRUE,
                           normDeg = FALSE,
                           normBetw = FALSE,
                           normClose = FALSE,
                           normEigen = FALSE)

#Plotting nicely
taxtab <- physeqC1_fil@tax_table@.Data
phyla <- as.factor(gsub("p__", "", taxtab[, "Phylum"]))
colors <- MetBrewer::met.brewer("Gauguin",n=10)

plot(props_clusC1.fil, 
     sameLayout = F, 
     nodeColor = "feature", 
     featVecCol = phyla, 
     colorVec =  colors,
     borderCol = "gray40",
     nodeSize = "clr", 
     cexNodes = 2,
     nodeSizeSpread = 4, 
     edgeTranspLow = 30, 
     edgeTranspHigh = 10,
     mar = c(3,3,3,3), 
     repulsion = 1, 
     rmSingles = "inboth",
     nodeFilter = "clustMin", 
     nodeFilterPar = 3, 
     nodeTransp = 10, 
     hubTransp = 10, 
     highlightHubs = T,
     cexLabels=0)


summary(props_clusC1.fil)

#Bug's network Cluster 2
netC2_fil <- netConstruct(data = physeqC2_fil,data2=NULL, dataType = "counts", group = NULL,
                       matchDesign = NULL, measure = "spieceasi", measurePar = NULL,
                       jointPrepro = NULL, filtTax = "none", filtTaxPar = NULL,
                       filtSamp = "none", filtSampPar = NULL,  zeroMethod = "none",
                       zeroPar = NULL, normMethod = "none", normPar = NULL, sparsMethod = "t-test",
                       thresh = 0.3, alpha = 0.05, adjust = "adaptBH", trueNullMethod = "convest",
                       lfdrThresh = 0.2, nboot = 1000L,  cores = 8,  logFile = "log.txt",
                       softThreshType = "signed",  softThreshPower = NULL, softThreshCut = 0.8,
                       kNeighbor = 3L, knnMutual = FALSE,  dissFunc = "signed",  dissFuncPar = NULL,
                       simFunc = NULL, simFuncPar = NULL,  scaleDiss = TRUE, weighted = TRUE,
                       sampleSize = NULL,  verbose = 2, seed = NULL)

props_clusC2.fil  <- netAnalyze(netC2_fil, 
                           centrLCC = FALSE,
                           avDissIgnoreInf = TRUE,
                           sPathNorm = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = c("degree", "between"),
                           hubQuant = 0.9,
                           lnormFit = TRUE,
                           normDeg = FALSE,
                           normBetw = FALSE,
                           normClose = FALSE,
                           normEigen = FALSE)

#Plotting nicely
taxtab <- physeqC2_fil@tax_table@.Data
phyla <- as.factor(gsub("p__", "", taxtab[, "Phylum"]))
colors <- MetBrewer::met.brewer("Gauguin",n=10)

plot(props_clusC2.fil, 
     sameLayout = F, 
     nodeColor = "feature", 
     featVecCol = phyla, 
     colorVec =  colors,
     borderCol = "gray40",
     nodeSize = "clr", 
     cexNodes = 2,
     nodeSizeSpread = 4, 
     edgeTranspLow = 30, 
     edgeTranspHigh = 10,
     mar = c(3,3,3,3), 
     repulsion = 1, 
     rmSingles = "inboth",
     nodeFilter = "clustMin", 
     nodeFilterPar = 3, 
     nodeTransp = 10, 
     hubTransp = 10, 
     highlightHubs = T,
     cexLabels=0)

summary(props_clusC2.fil)


### Comparisons

#Network comparison including the two clusters
net_CeD <- netConstruct(data = physeqC1_fil, data2=physeqC2_fil, dataType = "counts", group = NULL,
                          matchDesign = NULL, measure = "spieceasi", measurePar = NULL,
                          jointPrepro = NULL, filtTax = "none", filtTaxPar = NULL,
                          filtSamp = "none", filtSampPar = NULL,  zeroMethod = "none",
                          zeroPar = NULL, normMethod = "none", normPar = NULL, sparsMethod = "t-test",
                          thresh = 0.3, alpha = 0.05, adjust = "adaptBH", trueNullMethod = "convest",
                          lfdrThresh = 0.2, nboot = 1000L,  cores = 8,  logFile = "log.txt",
                          softThreshType = "signed",  softThreshPower = NULL, softThreshCut = 0.8,
                          kNeighbor = 3L, knnMutual = FALSE,  dissFunc = "signed",  dissFuncPar = NULL,
                          simFunc = NULL, simFuncPar = NULL,  scaleDiss = TRUE, weighted = TRUE,
                          sampleSize = NULL,  verbose = 2, seed = NULL)

props_clusCeD <- netAnalyze(net_CeD, 
                               centrLCC = FALSE,
                               avDissIgnoreInf = TRUE,
                               sPathNorm = FALSE,
                               clustMethod = "cluster_fast_greedy",
                               hubPar = c("degree", "between"),
                               hubQuant = 0.9,
                               lnormFit = TRUE,
                               normDeg = FALSE,
                               normBetw = FALSE,
                               normClose = FALSE,
                               normEigen = FALSE)
summary(props_clusCeD)

#Plotting
taxtab <- physeqC1_fil@tax_table@.Data
phyla <- as.factor(gsub("p__", "", taxtab[, "Phylum"]))
colors <- MetBrewer::met.brewer("Gauguin",n=12)

plot(props_clusCeD, 
     sameLayout = F, 
     nodeColor = "feature", 
     featVecCol = phyla, 
     colorVec =  colors,
     borderCol = "gray40",
     nodeSize = "clr", 
     cexNodes = 2,
     nodeSizeSpread = 4, 
     edgeTranspLow = 30, 
     edgeTranspHigh = 10,
     mar = c(3,3,3,3), 
     repulsion = 1, 
     rmSingles = "inboth",
     nodeFilter = "clustMin", 
     nodeFilterPar = 3, 
     nodeTransp = 10, 
     hubTransp = 10, 
     highlightHubs = T,
     cexLabels=0)


plot(props_clusCeD, 
     nodeColor = "cluster", 
     nodeSize = "eigenvector",
     title1 = "Network on OTU level with SPRING associations", 
     showTitle = TRUE,
     cexTitle = 2.3)

summary(props_clusCeD)


#Comparisson with Perm Test *Takes time*

comp_clus_perm <- netCompare(props_clusCeD, permTest = TRUE, nPerm = 100, cores=9L,
                             storeAssoPerm = TRUE, verbose=T,
                             fileStoreAssoPerm = "assoPerm_comp",
                             storeCountsPerm = FALSE, 
                             seed = 123456)

summary(comp_clus_perm)

summary(comp_clus_perm, 
        groupNames = c("Cluster 1", "Cluster 2"),
        showCentr = c("degree", "between", "closeness"), 
        numbNodes = 10)

#Differential networks 
net_clus_pears <- netConstruct(data = physeqC1_fil, 
                               data2 = physeqC2_fil, 
                               measure = "pearson", 
                               filtTax = "highestFreq",
                               filtTaxPar = list(highestFreq = 70),
                               zeroMethod = "pseudo",
                               normMethod = "clr",
                               sparsMethod = "none", 
                               thresh = 0.2,
                               verbose = 3)

diff_clus <- diffnet(net_clus_pears,
                     diffMethod = "fisherTest", 
                     adjust = "lfdr")

summary(diff_clus)

plot(diff_clus,
     cexLegend = 0.7,cexLabels=3)


#Differentially associated species

props_clus_pears <- netAnalyze(net_clus_pears, 
                               clustMethod = "cluster_fast_greedy",
                               weightDeg = TRUE,
                               normDeg = FALSE)

# Identify the differentially associated species
diffmat_sums <- rowSums(diff_clus$diffAdjustMat)
diff_asso_names <- names(diffmat_sums[diffmat_sums > 0])

plot(props_clus_pears, 
     nodeFilter = "names",
     nodeFilterPar = diff_asso_names,
     nodeColor = "gray",
     highlightHubs = FALSE,
     sameLayout = TRUE, 
     layoutGroup = "union",
     rmSingles = FALSE, 
     nodeSize = "mclr",
     edgeTranspHigh = 10,
     edgeTranspLow = 10, 
     shortenLabels="none",
     labelScale = F,
     cexNodes = 1, 
     cexLabels = 1,
     cexTitle = 1,
     mar = c(5, 5, 5, 5),
     posCol = "#b2df8a", 
     negCol = "#B2ABD2",
     groupNames = c("Cluster 1", "Cluster 2"),
     hubBorderCol  = "gray40",
     nodeTransp = 20)

plot(1,0.6)
legend("bottom", title = "estimated correlation:", legend = c("+","-"), 
       col = c("#b2df8a","#B2ABD2"), cex = 4, lty = 1, lwd = 4, 
       bty = "n", horiz = TRUE)

