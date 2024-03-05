library(tidyverse)
library(FactoMineR)
library(factoextra)


DATA_plus <- read_csv("Data/Objetos/data_plus.csv")

#Multiple Factor Analisys 
res.mfa <- MFA(DATA_plus, 
               group = c(10,11,1,3,3,1,1,3,1,1,10,1,15,1,1,1,1,5,7,8,1), 
               type = c("s","s","s","s","s","s","s","s","n", "s", "n", "s","n","s","s","s","s","s","s","s","s"),
               name.group = c("Haemathology","Biochemestry","Age","Anthropometrics","Fat-muscle","Fat-visc.","Metabolism","Blood pressure","CDAT","GIP.feaces","Symptoms1", "Symptoms1a","Symptoms","Symptoms2a", "Symptoms2b", "CeD", "GSRS", "Inflammation","Cytokines","Mucosal integrity","IFAP"),
               num.group.sup = NULL,
               graph = TRUE)

#Eigenvalues_Variances
eig.val <- get_eigenvalue(res.mfa)
head(eig.val)

fviz_screeplot(res.mfa)

group <- get_mfa_var(res.mfa, "group")
group

fviz_mfa_var(res.mfa, "group")

#According to scree plot the first 3 dimensions are more relevant 

# Contribution to the first dimension
fviz_contrib(res.mfa, "group", axes = 1)

# Contribution to the second dimension
fviz_contrib(res.mfa, "group", axes = 2)

# Contribution to the third dimension
fviz_contrib(res.mfa, "group", axes = 3)


#Quantitave variables
quanti.var <- get_mfa_var(res.mfa, "quanti.var")
quanti.var 

paleta = randomcoloR::distinctColorPalette(k=18)

fviz_mfa_var(res.mfa, "quanti.var", palette = "paleta", 
             col.var.sup = "violet", repel = TRUE,
             geom = c("point", "text"), legend = "bottom")

# Contributions to dimension 1
fviz_contrib(res.mfa, choice = "quanti.var", axes = 1, top = 20,
             palette = "paleta")

# Contributions to dimension 2
fviz_contrib(res.mfa, choice = "quanti.var", axes = 2, top = 20,
             palette = "paleta")

# Contributions to dimension 3
fviz_contrib(res.mfa, choice = "quanti.var", axes = 3, top = 20,
             palette = "paleta")


fviz_mfa_var(res.mfa, "quanti.var", col.var = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             col.var.sup = "violet", repel = TRUE,
             geom = c("point", "text"))

fviz_cos2(res.mfa, choice = "quanti.var", axes = 1)

#Individuos
ind <- get_mfa_ind(res.mfa)
ind

fviz_mfa_ind(res.mfa, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)
