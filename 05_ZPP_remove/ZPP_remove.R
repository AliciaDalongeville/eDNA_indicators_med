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
  dplyr::select(Sample, Site, protection,ZPP_remove, habitat_principal, latitude_start_DD, longitude_start_DD) 

## Combine all the data to be used in the models in one dataframe
mydata <- indicators %>% 
  inner_join(env, by="Samples") %>% # create a new dataframe combining indicators & environment
  # only rows that are in both dataframes are kept
  inner_join(metadata, by=c("Samples" = "Sample")) %>% # Add the protection
  filter(R >= 10 & ZPP_remove=="No")


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
  dplyr::select(protection, Dim.1, Dim.2, Dim.3, Dim.4) %>%
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
  
  # Calculate r2 and store it in the result vector
  r2[i] <- rsq(m_i, adj=T)
  
  # Store coefficient in the result matrix
  coef[,i] <- summary(m_i)$coefficients[2:6, 1]
  # Store ANOVA pval in the results
  coef[,(i+ncol(Y))] <- Anova(m_i,test.statistic="F")$"Pr(>F)"[1:5]
}


## For which indicator is protection significant ?
pval <- t(coef[,12:22])
sign <- pval[which(pval[,1] < 0.1),]

####################################################################
#check Spatial autocorrelation 
lstw <- nb2listw((knn2nb(knearneigh(coords, k=1, longlat = T, 
                                    use_kd_tree=F))))

# Compute Moran's I using residuals of model and also raw data
png("05_ZPP_remove/QQplot_residuals_ZPP.png", width = 800, height = 900)
layout(
  matrix(1:12, ncol=3, byrow=TRUE),  # plot 12 graphs in 3 columns
  widths=c(1,1,1), # widths of each column
  heights=c(1,1,1,1)) # height of each row

moran_p <- vector()
for (i in 1:11) {
  m_i <- glm_results[[i]]
  moran_p[i] <-moran.test(m_i$residuals, lstw, 
                          alternative="two.sided")$p.value 
  
  # QQ plot of residuals
  m.stdres = rstandard(m_i)
  qqnorm(m.stdres, 
         ylab="Standardized Residuals", 
         xlab="Normal Scores", 
         main=paste0(ind_names[i], " - Moran.I test p=", round(moran_p[i], 3)),
         cex=1.2, cex.axis=1.3,cex.lab=1.3, cex.main=1.5)
  qqline(m.stdres)
}

dev.off()

moran_res <- cbind.data.frame(ind_names,moran_p)
moran_res[which(moran_res$moran_p < 0.05),] # No SAS
mean(moran_p) 

#####################################################################################
## Fit the marginal models
#####################################################################################
## Rescale the response varibles between 0 and 1
Y <- Y %>%
  mutate(across(everything(), ~ rescale(.x, to=c(0,1))))


## prepare the result table
ame_all <- data.frame(Y=rep(NA,(6*ncol(Y))),
                      X=rep(NA,(6*ncol(Y))),
                      AME=rep(NA,(6*ncol(Y))),
                      SE=rep(NA,(6*ncol(Y))),
                      z=rep(NA,(6*ncol(Y))),
                      p=rep(NA,(6*ncol(Y))),
                      lower=rep(NA,(6*ncol(Y))),
                      upper=rep(NA,(6*ncol(Y))))
e=1 # index to fill result matrix

for (i in 1:ncol(Y)) {
  # Gather data
  data_i <- cbind.data.frame(Y[,i], X)
  
  # compute model
  m_i <- glm(data_i[,1] ~ protection + Dim.1 + Dim.2 +Dim.3 + Dim.4 , 
             data=data_i , family = distri[i], na.action = na.omit )
  # calcul de l'autocov
  acinv2a = autocov_dist(residuals(m_i),coords, nbs = 100, 
                         type="inverse", style=style.autocov[i], zero.policy = T, longlat=T)
  acinv2a[is.infinite(acinv2a) | is.na(acinv2a)] <- 0
  # re-fit model with autocov
  data_i=cbind.data.frame(data_i,acinv2a)
  x <- glm(data_i[,1] ~ protection + Dim.1 + Dim.2 +Dim.3 + Dim.4 + acinv2a , 
           data=data_i , family = distri[i], na.action = na.omit )
  
  # Calculate marginal effects at the mean (i.e., partial effects at the average of all covariates)
  m <- margins(model = x, data = data_i, type = "response")
  
  # Print the summary
  print(ind_names[i])
  print(summary(m))
  
  ## Fill result table
  ame_all[e:(e+5),"Y"] <- rep(ind_names[i],6)
  ame_all[e:(e+5),2:8] <- summary(m)
  
  rm(x)
  rm(m)
  rm(data_i)
  e=e+6
}

#######################################################
## FIGURE 2 - Plot the predicted values inside and outside reserve 
# for each indicators
#######################################################
## vector of p-values
signif <- ame_all %>%
  filter(X=="protectionreserve") %>%
  mutate(p = round(p,3)) %>%
  pull(p)
print(signif)

# list to store plots
p <-list()

# vector of effect strength according to the scale of Muff et al. 2021
evidence <- c("Little", "Strong", "No", "No", "Strong",  "No",
              "Moderate", "No", "No", "No", "No")

names <- c("Species richness", gsub("_perc", "",ind_names[-1]))

## Draw plot for the significant indicators
l=1
for (i in c(1,2,4,7,3,5,6,8:11)) {
  # Gather data
  data_i <- cbind.data.frame(Y[,i], X)
  # compute model
  m_i <- glm(data_i[,1] ~ protection + Dim.1 + Dim.2 +Dim.3 + Dim.4 , 
             data=data_i , family = distri[i], na.action = na.omit )
  # calcul de l'autocov
  acinv2a = autocov_dist(residuals(m_i),coords, nbs = 100, 
                         type="inverse", style=style.autocov[i], zero.policy = T, longlat=T)
  acinv2a[is.infinite(acinv2a) | is.na(acinv2a)] <- 0
  # re-fit model with autocov
  data_i=cbind.data.frame(data_i,acinv2a)
  x <- glm(data_i[,1] ~ protection + Dim.1 + Dim.2 +Dim.3 + Dim.4 + acinv2a , 
           data=data_i , family = distri[i], na.action = na.omit )
  
  # Calculate marginal effects at the mean (i.e., partial effects at the average of all covariates)
  m <- margins(model = x, data = data_i, type = "response")
  
  # Gather data to plot
  to_plot <- cbind.data.frame(m$fitted, X$protection, m$se.fitted)
  colnames(to_plot) <- c("Y", "Protection", "SE")
  
  ## Draw the plot
  p[[l]]<-ggplot(data=to_plot, aes(x=Protection, y=Y, color=Protection)) +
    geom_boxplot(notch=F) +
    scale_color_manual(values=c("cyan4", "darkblue")) +
    ylab(paste("Predicted ",names[i], sep=" ")) +
    ggtitle(paste("ME test, p = ",signif[i], " - ", evidence[i], " effect", sep="")) +
    theme_bw() +
    theme(plot.title = element_text(size=15, 
                                    colour = "black",
                                    face = "bold")) +
    theme(legend.position="none") +
    theme(axis.text.y=element_text(colour="black",size=13)) +
    theme(axis.text.x=element_text(colour="black",size=16)) +
    theme(axis.title.x=element_blank()) +
    theme(axis.title.y=element_text(colour="black",size=16)) 
  
  l=l+1
}

## Save plot
png("05_ZPP_remove/Boxplot_predicted_ZPP_remove.png", 
    width = 1000, height = 1000) 
do.call(grid.arrange,p)
dev.off()

#######################################################
## FIGURE 3 - Plot the AME of reserve protection on a bar chart
#######################################################
## Gather data to plot
to_plot <- ame_all %>%
  filter(X=="protectionreserve") %>%
  mutate(Indicator=gsub("_perc", "",Y)) %>%
  dplyr::select(Indicator,AME,SE,p) %>%
  arrange(desc(AME)) %>%
  mutate(signif = (p<0.05)) 
to_plot$Indicator[10] <- "Species Richness"

# Change the order of factors to plot in ascending order
to_plot$Indicator <- factor(to_plot$Indicator,                                    # Factor levels in decreasing order
                            levels = to_plot$Indicator[order(to_plot$AME)])

my_breaks <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 1)
breaks_lab <- c("0.0001", "0.001", "0.01", "0.05", "0.1", "1")

# Basic barplot
p<-ggplot(data=to_plot, aes(x=Indicator, y=AME, width=0.8, fill=p)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=AME-SE, ymax=AME+SE), width=.2) +
  ylim(-0.33,0.35) +
  ylab("Partial effect of reserve protection (AME)") +
  coord_flip() +
  scale_fill_viridis(limits = c(0,1),breaks = my_breaks, labels = breaks_lab, 
                     trans = scales::pseudo_log_trans(sigma = 0.001)) +
  labs(fill="ME test \np-value") +
  theme_bw() +
  theme(legend.position="right") +
  theme(axis.text=element_text(colour="black",size=17))+
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x=element_text(colour="black",size=18)) +
  theme(legend.title=element_text(colour="black",size=18)) +
  theme(legend.text=element_text(colour="black",size=18)) +
  theme(legend.key=element_rect(size=23)) +
  guides(fill = guide_colourbar(barwidth = 1.2, barheight = 20))

p

### Save the plot
ggsave("05_ZPP_remove/Protection_effects_ZPP_remove.png",
       width=12, height=10,
       dpi = 300)

