###############################################################################
# 
#  codes to reproduce analyses and figures of the Dalongeville et al. article accepted for publication in Journal of Applied Ecology: 
#  "Benchmarking eleven biodiversity indicators based on environmental DNA surveys: more diverse functional traits and evolutionary lineages inside marine reserves" 
#  
#   Code author : Alicia Dalongeville
#   Date : August 2022
###############################################################################
### Load the libraries
library(dplyr)
library(ggpubr)
library(FactoMineR)
library(factoextra)
library(missMDA)
library(corrplot)

### Load the environmental data
env <- read.csv("Data/environment.csv", header=T)

##################################################################################################
### PCA with environmental variables
##################################################################################################
# select variables to be included in the PCA
env_pca <- env %>%
  dplyr::select(Site, dist_land, dist_port, min_bathy,
                habitat_div,
                chloroYear, chloroWeek,
                tempYear, tempWeek)

## Impute the missing values
res.impute <- imputePCA(env_pca[,-1], ncp=3) 
## The output can be used as an input of the PCA function of the FactoMineR package 
  ##to perform the PCA on the incomplete data  
res.pca <- FactoMineR::PCA(res.impute$completeObs, scale.unit = TRUE, ncp = 5, graph = FALSE)

## Get eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val

# Draw screeplot
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

#################################################################
## Variable contributions to axes
var <- get_pca_var(res.pca)
# Coordinates of variables
head(var$coord)
# Cos2: quality of representation on the factore map
head(var$cos2)
# Contributions to the  dimensions
head(var$contrib)

# Corrplot of the cos2 of variables
corrplot(var$cos2, is.corr=FALSE)

# Total cos2 of variables on Dim.1 to Dim.3
fviz_cos2(res.pca, choice = "var", axes = 1:3)

# Contributions of variables to PC1 : chlorophyll + temperature + bathy
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2 : distland + ports + temp week
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)


# Plot the variables with Color by contribution to the first 2 axes
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

### Dimension description
res.desc <- dimdesc(res.pca, axes = c(1:3), proba = 0.05)
# Description of dimension 1
res.desc$Dim.1
# Description of dimension 2
res.desc$Dim.2
# Description of dimension 3
res.desc$Dim.3

#####################################
## Get the PCA results 
ind <- get_pca_ind(res.pca)
s_coord <- as.data.frame(ind$coord[,1:4])        # Coordinates for the first 3 dimensions
s_coord$Samples <- env$Sample

write.csv(s_coord, file="02_PCA/PCA_sites_coord.csv", row.names=F)

## biplot with color by protection status
p1 <- fviz_pca_biplot(res.pca, 
                col.ind = env$protection, palette = c("cyan4","darkblue"), 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Protection") 
p2 <- fviz_pca_biplot(res.pca, axes = c(3,4),
                      col.ind = env$protection, palette = c("cyan4","darkblue"), 
                      addEllipses = TRUE, label = "var",
                      col.var = "black", repel = TRUE,
                      legend.title = "Protection") 

## biplot with color by site
p3 <- fviz_pca_ind(res.pca, 
                pointshape = 21,
                pointsize = 2,
                fill.ind = env$Site,
                col.ind = "black", palette = hcl.colors(n=9, "Dynamic"), 
                addEllipses = TRUE, label = "none",
                repel = TRUE, 
                legend.title = "Sites") 
p4 <- fviz_pca_ind(res.pca, axes=c(3,4),
                      pointshape = 21,
                      pointsize = 2,
                      fill.ind = env$Site,
                      col.ind = "black", palette = hcl.colors(n=9, "Dynamic"), 
                      addEllipses = TRUE, label = "none",
                      col.var = "black", repel = T, 
                      legend.title = "Sites") 
library(gridExtra)
fig <- grid.arrange(p1,p2,p3,p4, ncol=2)
### Save the figure
ggsave(filename = "02_PCA/PCA_biplot.png", 
       plot = fig,
       width = 10, 
       height = 10,
       dpi = 300)
