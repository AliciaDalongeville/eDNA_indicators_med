###############################################################################
# 
#  codes to reproduce analyses and figures of the Dalongeville et al. article accepted for publication in Journal of Applied Ecology: 
#  "Benchmarking eleven biodiversity indicators based on environmental DNA surveys: more diverse functional traits and evolutionary lineages inside marine reserves" 
#  
#   Code author : Alicia Dalongeville
#   Date : August 2022
###############################################################################
## load packages
library(dplyr)
library(tidyr)
library(fastDummies)
library(ape)
library(fishtree)
library(picante)
library(geiger)
loadNamespace("rfishbase")

## Load the eDNA data (matrix species per sample)
adne <- read.csv("Data/eDNA.csv", header=T, row.names = 1 )

# list of species
species <- rownames(adne)
# list of samples
samples <- colnames(adne)

## Load the species traits matrix
traits <- read.csv("Data/functional_data.csv", header=T)


#######################################################################################################
## Calculate the indicators for each sample
#######################################################################################################
## Create the result matrix
indicators <- matrix(NA,ncol(adne), 11,
                     dimnames=list(colnames(adne),
                                  c("R",  "FD", "LFI", "Crypto", "DP_B_ratio",  "RedList", "Chondri", "Commercial", "High_commerc", "PD", "Vulner")))

##########################################################
## 1 - Species Richness R
##########################################################
indicators[,1] <- apply(adne, 2, sum)

##########################################################
## 2 - Functional diversity FD
##########################################################
for (i in 1:nrow(indicators)) { # for each sample
  # list species present in the sample
  s_i <- adne %>%
    filter(adne[,i] == 1) %>%
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)
  
  # Get the functional groups of these species
  fd_i <- as.factor(traits[which(traits$Species %in% s_i), "GF"])
  
  # Number of unique functional groups
  indicators[i,2] <- nlevels(fd_i)
}

##########################################################
## 3 - Large Reef Fish Index - LFI
## 5 - Ratio Demerso-pelagic / benthic
## 7 - Chondrichtyen species
## 8 - Commercial species
## 9 - Highly commercial species 
##########################################################
for (i in 1:nrow(indicators)) { # for each sample
  # list species present in the sample
  s_i <- adne %>%
    filter(adne[,i] == 1) %>%
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)

  # Calculate indicators
  indicators[i,"LFI"] <- sum(traits[which(traits$Species %in% s_i), "LRFI"]) 
  indicators[i,"DP_B_ratio"] <- sum(traits[which(traits$Species %in% s_i), "DP"]) / (sum(traits[which(traits$Species %in% s_i), "B"])+1)
  indicators[i,"Chondri"] <- sum(traits[which(traits$Species %in% s_i), "SHarK"])
  indicators[i,"Commercial"] <- sum(traits[which(traits$Species %in% s_i), "all_commercial_level"])
  indicators[i,"High_commerc"] <- sum(traits[which(traits$Species %in% s_i), "highly_commercial_only"])
}

###########################################################
## 4 Cryptobenthic (definition Brandl et al. 2018)
                       # Brandl SJ, Goatley CHR, Bellwood DR, Tornabene L. 2018 The hidden half: ecology and evolution of cryptobenthic fishes on coral reefs. Biol. Rev. 93, . (doi:10.1111/brv.12423))
###########################################################
crypto_families = c("Tripterygiidae", "Grammatidae", "Creediidae", "Aploactinidae", "Gobiidae", 
                    "Chaenopsidae", "Gobiesocidae", "Labrisomidae", "Pseudochromidae", "Bythitidae", 
                    "Plesiopidae", "Blenniidae", "Apogonidae", "Callionymidae", "Opistognathidae", "Syngnathidae")

traits <- traits %>%
  mutate(crypto_Brandl = if_else(Family %in% crypto_families, 1,0))

for (i in 1:nrow(indicators)) { # for each sample
  # list species present in the sample
  s_i <- adne %>%
    filter(adne[,i] == 1) %>%
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)
  
  # Calculate indicator
  indicators[i,"Crypto"] <- sum(traits[which(traits$Species %in% s_i), "crypto_Brandl"]) 
}

##########################################################
## 11 - Vulnerability
##########################################################
for (i in 1:nrow(indicators)) { # for each sample
  # list species present in the sample
  s_i <- adne %>%
    filter(adne[,i] == 1) %>%
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)
  
  # calculate indicators
 indicators[i,"Vulner"] <- mean(traits[which(traits$Species %in% s_i), "Vulnerability"], na.rm=T)
}

##########################################################
## 6 - Red List IUCN
##########################################################
## Count all species listed VU, EN or CR on the Red List of Threatened Species

### Make a dummy variable for IUCN categories
traits <- dummy_cols(traits, select_columns = 'IUCN_Red_List_Category')

## Calculate indicator
for (i in 1:nrow(indicators)) { # for each sample
  # list species present in the sample
  s_i <- adne %>%
    filter(adne[,i] == 1) %>%
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)
  
  # Number of species per category
  VU <- sum(traits[which(traits$Species %in% s_i), "IUCN_Red_List_Category_VU"], na.rm=T)
  EN <- sum(traits[which(traits$Species %in% s_i), "IUCN_Red_List_Category_EN"], na.rm=T)
  CR <- sum(traits[which(traits$Species %in% s_i), "IUCN_Red_List_Category_CR"], na.rm=T)
  
  # Calculate indicator
    indicators[i,"RedList"] <- VU + EN + CR
}

##########################################################
## 10 - Phylogenetic Diversity - PD
##########################################################
# Retrieve the phylogeny of only native reef species across all three oceans.
phy <- fishtree_phylogeny(species = species)

plot(phy, show.tip.label = FALSE)
 tiplabels(tip = which(phy$tip.label %in% species),
             pch=19, cex=2)

rownames(adne) <- gsub(" ", "_", rownames(adne), fixed = TRUE)

## check that phylogeny and data have matching names
#nc <- geiger::name.check(phy, adne) # 24 species not in tree (chondrychtyans + synonyms)

# Manually check synonyms and find the species
species[species == "Mullus barbatus"] <- "Mullus barbatus barbatus"
rownames(adne)[rownames(adne) == "Mullus_barbatus"] <- "Mullus_barbatus_barbatus"

species[species == "Diplodus sargus"] <- "Diplodus sargus sargus"
rownames(adne)[rownames(adne) == "Diplodus_sargus"] <- "Diplodus_sargus_sargus"

species[species == "Diplodus cervinus"] <- "Diplodus cervinus cervinus"
rownames(adne)[rownames(adne) == "Diplodus_cervinus"] <- "Diplodus_cervinus_cervinus"

species[species == "Chelon auratus"] <- "Liza aurata"
rownames(adne)[rownames(adne) == "Chelon_auratus"] <- "Liza_aurata"

species[species == "Chelon ramada"] <- "Liza ramada"
rownames(adne)[rownames(adne) == "Chelon_ramada"] <- "Liza_ramada"

# Retrieve the missing phylogeny 
phy <- fishtree_phylogeny(species = species, type="phylogram")
nc <- geiger::name.check(phy, adne) # 19 species not in the tree

# Remove from the data the species that are not in the tree
adne2 <- adne[which(rownames(adne) %in% nc$data_not_tree == F),]

# Transpose the ADNe matrix 
adne2 <- t(adne2)

# prune the tree
prunedTree <- prune.sample(adne2,phy)

# Calculate Faith's PD
pd.result <- pd(adne2, prunedTree, include.root=T)

# Add PD to indicator dataframe
indicators[,"PD"] <- pd.result$PD/pd.result$SR  #percentage of species richness with only species that are in the tree

# Save the indicator table
write.csv(indicators, file="01_indicators/indicators_updated.csv")

