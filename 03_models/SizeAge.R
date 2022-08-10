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
library(rsq)
library(fastDummies)
library(tidyr)
library(purrr)
library(scales)
library(margins)
library(ggplot2)
library(gridExtra)
library(viridis)
library(spdep)

### Load the data
# indicators
indicators <- read.csv("01_indicators/indicators_updated.csv", header=T, na.strings=c("Inf","NA"))

# Environmental PCA data
env <- read.csv("02_PCA/PCA_sites_coord.csv", header=T) 

# Protection and site
metadata <- read.csv("Data/environment.csv", header=T) %>%
  dplyr::select(Sample, Site, protection, habitat_principal, latitude_start_DD, longitude_start_DD) 

# Reserve data
reserve <- read.table("Data/reserve.csv", sep=",", header=T)

## Combine all the data to be used in the models in one dataframe
mydata <- indicators %>% 
  inner_join(env, by="Samples") %>% # create a new dataframe combining indicators & environment
  # only rows that are in both dataframes are kept
  inner_join(metadata, by=c("Samples" = "Sample")) %>% # Add the protection
  filter(R >= 10)

## Add the size & age of the reserve to the dataframe 'mydata'
mydata <- mydata %>%
  left_join(reserve, by="Site") %>%
  dplyr::select(-c("Reserve", "Year", "N_reserve", "N_out")) %>%
  # Select only data that are within reserves
  #filter(protection == "reserve")

### Create a data set of response variables (indicators)
Y <- mydata %>%
  # log-transform variables
  #mutate(across(c(fric,PD_perc, DP_B_ratio,Vulner), log10)) %>%
  # express vqriables as proportion of species richness
  mutate(across(c(FD,LFI,Crypto, Chondri,  Commercial, High_commerc, IUCN_noweigths), ~ .x/R)) %>%
  dplyr::select(R, FD, High_commerc, PD_perc,LFI, Crypto, DP_B_ratio,
                Chondri,  Commercial, Vulner, IUCN_noweigths) 

### Create a data set of explanatory variables (environment + protection)
X <- mydata %>%
  dplyr::select(Dim.1, Dim.2, Dim.3, Dim.4, Age, Size_ha, protection) %>%
  # transform all character columns to factors
  mutate_if(sapply(., is.character), as.factor)

### Create a vector of site names
sites <- mydata %>%
  dplyr::pull(Site)

### Create a vector of sample names
samples <- mydata %>%
  dplyr::pull(Samples)

### Matrix of spatial coordinates
coords<- mydata %>%
  select(longitude_start_DD, latitude_start_DD) %>%
  as.matrix(.)

### Create a vector of indicator names
ind_names <- c("R",  "FD", "High_commerc", 
               "PD", "LFI", "Crypto", 
               "DeBRa", "Chondri",  "Commercial", 
               "Vulner", "RedList")

# Create a vector of distributions : Poisson for RedList and Chondri, 
# and Gaussian for all the others
distri <- c("gaussian", "gaussian", "gaussian",
            "gaussian", "gaussian", "gaussian", 
            "gaussian", "quasipoisson", "gaussian",
            "gaussian", "quasipoisson")

##################
## Fit models

#####################################################################################
## Fit the GLMs
#####################################################################################
# create a list to store the models
glm_results <- list()

# Create a vector to store the R-squared
r2 <- rep(NA, ncol(Y), names=ind_names)

# Create a matrix to store the coefficient and pvalues
coef <- matrix(NA, ncol(X), (2*ncol(Y)),
               dimnames= list(colnames(X),
                              c(paste(ind_names,"coef", sep="_" ), paste(ind_names,"pval", sep="_" ))))

style.autocov <- c("B","B","B","B","B","W", "W", "B","B","B","B")

# Make a loop to do the glm for each indicator
for (i in 1:ncol(Y)) {
  data_i <- cbind.data.frame(Y[,i], X) # create a dataframe with the data for indicator i
  m_i <- glm(data_i[,1] ~ Age + Size_ha + Dim.1 + Dim.2 +Dim.3 + Dim.4 , 
             data=data_i , family = distri[i], na.action = na.omit )
  
  # calcul de l'autocov
  acinv2a = autocov_dist(residuals(m_i),coords, nbs = 100, 
                         type="inverse", style=style.autocov[i], zero.policy = T, longlat=T)
  acinv2a[is.infinite(acinv2a) | is.na(acinv2a)] <- 0
  # re-fit model with autocov
  data_i=cbind.data.frame(data_i,acinv2a)
  m_i <- glm(data_i[,1] ~ Age + Size_ha  + Dim.1 + Dim.2 +Dim.3 + Dim.4 + acinv2a , 
             data=data_i , family = distri[i], na.action = na.omit )
  
  # store the model in the list
  glm_results[[i]] <- m_i
  
  # Calculate r2 and store it in the result vector
  r2[i] <- rsq(m_i, adj=T)
  
  # Store coefficient in the result matrix
  coef[,i] <- summary(m_i)$coefficients[2:7, 1]
  # Store ANOVA pval in the results
  coef[,(i+ncol(Y))] <- car::Anova(m_i,test.statistic="F")$"Pr(>F)"[1:6]
}


#######################################
## Plots
######################################
# list to store plots
p<-list()
p2 <- list()


for (i in 1:ncol(Y)) {
  data_i <- cbind.data.frame(Y[,i], X) # create a dataframe with the data for indicator i
  m_i <- glm(data_i[,1] ~ Age + Size_ha + Dim.1 + Dim.2 +Dim.3 + Dim.4 , 
             data=data_i , family = distri[i], na.action = na.omit )
  
  # calcul de l'autocov
  acinv2a = autocov_dist(residuals(m_i),coords, nbs = 100, 
                         type="inverse", style=style.autocov[i], zero.policy = T, longlat=T)
  acinv2a[is.infinite(acinv2a) | is.na(acinv2a)] <- 0
  # re-fit model with autocov
  data_i=cbind.data.frame(data_i,acinv2a)
  x <- glm(data_i[,1] ~ Age + Size_ha  + Dim.1 + Dim.2 +Dim.3 + Dim.4 + acinv2a , 
             data=data_i , family = distri[i], na.action = na.omit )
  
  # Calculate marginal effects at the mean (i.e., partial effects at the average of all covariates)
  m <- margins(model = x, data = data_i, type = "response")
  ame_age <- round(c(summary(m)[2,2],summary(m)[2,5]), 3)
  ame_size <- round(c(summary(m)[7,2],summary(m)[7,5]),3)
  
  # Gather data to plot
  to_plot <- cbind.data.frame(m$fitted, X$Age, X$Size_ha)
  colnames(to_plot) <- c("Y", "Age", "Size")
  
  ## Draw the plot for Age
  p[[i]]<-ggplot(data=to_plot, aes(x=Age, y=Y)) +
    geom_point() +
    geom_smooth(method='lm') +
    ylab(paste("Predicted ",ind_names[i], sep=" ")) +
    xlab("Age of reserve (years)") + 
    annotate("text",  x=Inf, y = Inf, 
             label = paste0("AME=",ame_age[1],"; pval=",ame_age[2]), 
             vjust=1.05, hjust=1.05, size=6) +
    theme_bw() +
    theme(legend.position="none") +
    theme(axis.text.y=element_text(colour="black",size=13)) +
    theme(axis.text.x=element_text(colour="black",size=16)) +
    theme(axis.title=element_text(colour="black",size=16)) 
  
  ## Draw the plot for size
  p2[[i]]<-ggplot(data=to_plot, aes(x=Size, y=Y)) +
    geom_point() +
    geom_smooth(method='lm') +
    ylab(paste("Predicted ",ind_names[i], sep=" ")) +
    xlab("Size of no-take zone (ha)") + 
    annotate("text",  x=Inf, y = Inf, 
             label = paste0("AME=",ame_size[1],"; pval=",ame_size[2]), 
             vjust=1, hjust=1, size=6) +
    theme_bw() +
    theme(legend.position="none") +
    theme(axis.text.y=element_text(colour="black",size=13)) +
    theme(axis.text.x=element_text(colour="black",size=16)) +
    theme(axis.title=element_text(colour="black",size=16)) 
  

}

png("03_models/Plot_indicators_AGE.png", width = 1000, height = 1000)
do.call(grid.arrange,p)
dev.off()

png("03_models/Plot_indicators_SIZE.png", width = 1000, height = 1000)
do.call(grid.arrange,p2)
dev.off()

#######################################################################
## Check the effect of size and age on the residuals of main models
#######################################################################
## Fit the main GLMs
#######################################################################
# create a list to store the models
glm_results <- list()

style.autocov <- c("B","B","B","B","B","W", "W", "B","B","B","B")

# Make a loop to do the glm for each indicator
for (i in 1:ncol(Y)) {
  data_i <- cbind.data.frame(Y[,i], X) # create a dataframe with the data for indicator i
  m_i <- glm(data_i[,1] ~ protection + Dim.1 + Dim.2 +Dim.3 + Dim.4 , 
             data=data_i , family = distri[i], na.action = na.omit )
  
  # calcul de l'autocov
  acinv2a = autocov_dist(residuals(m_i),coords, nbs = 100, 
                         type="inverse", style=style.autocov[i], zero.policy = T, longlat=T)
  acinv2a[is.infinite(acinv2a) | is.na(acinv2a)] <- 0
  # re-fit model with autocov
  data_i=cbind.data.frame(data_i,acinv2a)
  m_i <- glm(data_i[,1] ~ protection + Dim.1 + Dim.2 +Dim.3 + Dim.4 + acinv2a , 
             data=data_i , family = distri[i], na.action = na.omit )
  
  # store the model in the list
  glm_results[[i]] <- m_i
  
}


