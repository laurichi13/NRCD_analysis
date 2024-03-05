#Loading libraries
library(tidyverse)
library(phyloseq)
library(curatedMetagenomicData)
library(ape)
library(microbiome)
library(factoextra)

#Abundance table
Abundance <- read_csv("Data/Metagenomica/merged_abundance_table.txt")
Abundance <- column_to_rownames(Abundance, var="clade_name")
Abundance$NCBI_tax_id <- NULL

#Converte to phyloseq object
#Metadata after clustering to add cluster as a variable
DATA_plus <- read_csv("Data/Objetos/data_plus.csv")
DATAmfa.cluster <- read_csv("Data/Objetos/DATAmfa.cluster.csv")
DATA <- DATAmfa.cluster %>% select(ID,Cluster)
names(DATA) <- c("Sample_ID", "Cluster")
DATA <- DATA %>% left_join(DATA_plus)
DATA$Cluster <- factor(DATA$Cluster, levels = c("1", "2"), labels=c("1","2"))
row.names(DATA) <- DATA$Sample_ID
sample <- sample_data(DATA)
sample_names(sample) <- DATA$Sample_ID

#Load script Function metaphlantophyloseq
phyloseqin= metaphlanToPhyloseq(Abundance, metadat = sample)
phyloseqin

#Tree 
random_tree=rtree(ntaxa(phyloseqin), rooted = T, tip.label = taxa_names(phyloseqin))
#Objeto phyloseq con arbol
physeq.tree <- merge_phyloseq(phyloseqin, random_tree)

#Transfom counts for diversity analysis
#You might need to transform to absolute counts by multiplying each taxa abundance by one million - counts per million:
#because of how MetaPhlAn works it is not possible to have read counts. You can however have "pseudo" read counts by multiplying relative abundances by a constant and rounding to the closest integer.
NRCD <- transform_sample_counts(physeq.tree, function(x) 1E6 * x)
plot_richness(NRCD)
rich <- richness(NRCD)
plot_rich <- plot_richness(NRCD, color = "Cluster", x = "Cluster", measures = c("Chao1", "Simpson", "Shannon")) + geom_boxplot(aes(fill = Cluster), alpha=.7) + scale_color_manual(values = c("#B2ABD2", "#b2df8a")) + scale_fill_manual(values = c("#B2ABD2", "#b2df8a")) + theme_bw()


#Compositional analísis transform to CLR 
pseq.compositional <- transform(physeq.tree, "compositional")

#Unifrac
ordinated_taxa_unifrac = ordinate(physeq.tree, method="MDS", distance="unifrac",weighted=TRUE)
plot_ordination(physeq.tree, ordinated_taxa_unifrac, color="Cluster", title = "Unifrac MDS Analysis") + theme_bw()

#Análisis abundancia#### 
#Convert to relative abundance
#colores
library(randomcoloR)

physeq.fil = filter_taxa(pseq.compositional, function(x) mean(x) > 1e-5, TRUE)
Sp <- tax_glom(pseq.compositional, taxrank = "Species") 
tb_sp <- psmelt(Sp) 

#Plot species (Not included in the mansucript)
paleta308 <- distinctColorPalette(k = 308, altCol = FALSE, runTsne = FALSE)
PlotSp <- ggplot(tb_sp, aes(Sample, Abundance ,fill=OTU)) + geom_col(position="fill") + scale_fill_manual(values=paleta308)
PlotSp <- PlotSp + theme_minimal() + xlab("Sample") + ylab("Relative Abundance") +
  theme(axis.text.x = element_text(color = "black", size = 14, angle = 45, hjust = .5, vjust = .5, face = "plain")) +
  theme(axis.text.y = element_text(color = "black", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain")) + 
  theme(legend.title=element_text(color= "black", size=12), legend.text=element_text(size=11, face="italic"))
#PlotSp <- PlotSp +  ggtitle ("15 most abundant Species") + theme(plot.title = element_text(family="Arial", size=rel(2), vjust=2,face="plain", color="black", lineheight=1.5))  
PlotSp + theme(legend.position = "none") + coord_flip()

#Data to MetOrigins analysis####
pseq.compositional <- transform(physeq.tree, "compositional")
physeq.fil = filter_taxa(pseq.compositional, function(x) mean(x) > 1e-5, TRUE)
Sp <- tax_glom(pseq.compositional, taxrank = "Species") 
tb_sp <- psmelt(Sp)

gat <- tb_sp %>% select(OTU,Abundance, Sample_ID)
gat <- gat %>% spread(Sample_ID, Abundance)

tabla <- tb_sp %>% select(OTU, Kingdom,Phylum,Class,Order,Family,Genus,Species) 
tabla <- tabla %>% left_join(gat) %>% select(-OTU)

write_csv(tabla, "Data/Objetos/metOrigins.csv")
