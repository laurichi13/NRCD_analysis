#Heatmap correlation 
library(ComplexHeatmap)
library(tidyverse)
library(microbiome)

#My data
meta <- read_csv("Corr/meta1.csv")
sintomas <- read_csv("Corr/sintomas_tocor.csv")

#Pearson 
correlation.table <- associate(meta, sintomas, method = "pearson", mode = "table", p.adj.threshold = 0.05, p.adj.method ="BH")

my_mat <- correlation.table %>% select(-p.adj) %>% spread(X1,Correlation) %>% 
          column_to_rownames(var="X2") %>% as.matrix()

my_padj_mat <- correlation.table %>% select(-Correlation) %>% spread(X1,p.adj) %>% 
               column_to_rownames(var="X2") %>% as.matrix()

Heatmap(my_mat, row_order = sort(rownames(my_mat)), 
        column_order = sort(colnames(my_mat)))

Heatmap(my_mat, name = "mat", cluster_rows = FALSE) # turn off row clustering

#Adding P-values
Heatmap(my_mat, cell_fun = function(j, i, x, y, w, h, fill) {
  if(my_padj_mat[i, j] < 0.01) {
    grid.text("**", x, y)
  } else if(my_padj_mat[i, j] < 0.05) {
    grid.text("*", x, y)
  }
})


library(cluster)
Heatmap(my_mat, 
        row_order = sort(rownames(my_mat)),
        column_km = 2,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(my_padj_mat[i, j] < 0.01) {
            grid.text("**", x, y)
          } else if(my_padj_mat[i, j] < 0.05) {
            grid.text("*", x, y)
          }
        })


#Adding anotations
note <- read_csv("Corr/Data2/name_byCLus.csv")
anot <- as.factor(note$Enriched)
column_ha <- HeatmapAnnotation(Cluster = anot, col = list(Cluster = c("1" = "#CA8A98", "2" = "#4C3B33")))

library(MetBrewer)

Heatmap(my_mat, name="Pearson corr.", col = met.brewer("Archambault", 15),
        row_order = sort(rownames(my_mat)),
        top_annotation = column_ha,
        show_parent_dend_line = FALSE,
        column_km = 2,
        column_title = NULL,
        column_names_gp = grid::gpar(fontsize = 10),
        row_names_gp = grid::gpar(fontsize = 10),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(my_padj_mat[i, j] < 0.01) {
            grid.text("**", x, y)
          } else if(my_padj_mat[i, j] < 0.05) {
            grid.text("*", x, y)
          }
        })


