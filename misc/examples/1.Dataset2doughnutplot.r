# (C) Copyright 2016 Sur Herrera Paredes
#
#    This file is part of pbi.
#
#    pbi is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    pbi is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with pbi.  If not, see <http://www.gnu.org/licenses/>.

# This is a self contained demonstration of how to use the monte_carlo functions.
# It uses a dataset contained in the AMOR package (https://github.com/surh/AMOR).
# It is redundant in the sense that some functions have been copied from other parts
# of the repository into here

############ LOAD FUNCTIONS AND AMOR PACKAGE #############
library(AMOR)

#' Fit four ZINB models
#' 
#' Fits the 4 zinb models, and computes the best fit
#' for each taxon based on AIC
#' 
#' First it removes any sample that is not in the specified experiments. Requires
#' a Experiment column in the Map
#' 
#' Then it rmeoves OTUs under the measurable threshold (default 25in5)
#' 
#' It reorders the Genotype variable. It also creates Fracgen and LDepth.
#' It requires Genotype, Fraction and Depth columns to be present.
#' 
#' # Fits four ZINB models and computes the best for each taxa according to AIC
#' 
#' @param Dat a Dataset object
#' @param f1 right hand side formula for the model
#' @param group character vector of experiments to include
#' @param mutant_order Order if variables in Genotype column of the Map
#' @param min_reads_otu Minimum count for an otu in a sample for that sample to
#' be included in the measurable threshold
#' @param min_samples_otu Minimum samples that pass the min_reads_otu threshold for
#' a given otu to be counted as measurable.
#' 
fit_4_models <- function(Dat, f1, group, mutant_order, min_reads_otu = 25,
                         min_samples_otu = 5){
  
  # remove samples in different experiments
  Dat.group <- subset(Dat, Experiment %in% group, clean = TRUE, drop = TRUE)
  
  # Remove otus under the  measurable threshold
  # to_remove <- !findGoodOTUs(OTUtab = Dat.group$Tab,min_reads_otu = min_reads_otu,
  #                            min_samples_otu = min_samples_otu)
  # to_remove <- row.names(Dat.group$Tab)[ to_remove ]
  # Dat.group <- remove_taxons(Dat = Dat.group, to_remove)
  # Dat.group <- clean(Dat.group)
  Dat.group <- measurable_taxa(Dat,min_reads_otu = min_reads_otu,
                               min_samples_otu = min_samples_otu,
                               method = "absolute",table = TRUE, clean = TRUE)
  
  
  # Provess mappinfg file
  Dat.group$Map$Genotype <- factor(Dat.group$Map$Genotype, levels = mutant_order)
  Dat.group$Map$Fracgen <- factor(Dat.group$Map$Fraction, levels = c("EC","R","Soil"))
  Dat.group$Map$Fracgen <- droplevels(Dat.group$Map$Fracgen:Dat.group$Map$Genotype)
  Dat.group$Map$LDepth <- log10(Dat.group$Map$Depth)
  
  # Fit 4 models
  m1 <- matrix_glm(x = Dat.group, formula = f1, family = poisson(link = "log"),
                   verbose = TRUE)
  m2 <- matrix_glmNB(x = Dat.group, formula = f1,
                     verbose = TRUE)
  m3 <- matrix_zeroinfl(x = Dat.group, formula = f1, family = "poisson",
                        link = "log", verbose = TRUE)
  m4 <- matrix_zeroinfl(x = Dat.group, formula = f1, family = "negbin",
                        link = "log",verbose = TRUE)
  
  aic <- data.frame(glm = m1$AIC, glm.NB = m2$AIC, zip = m3$AIC, zinb = m4$AIC)
  best.index <- apply(aic,1,which.min)
  to_remove <- sapply(best.index,length)
  to_remove <- which(to_remove  == 0)
  best.index[ to_remove ] <- 5
  best.index <- unlist(best.index)
  
  aic$best <-  colnames(aic)[ best.index ]
  
  return(list(m1 = m1, m2 = m2, m3 = m3, m4 = m4, aic = aic))
}


##################3

# Create Dataset object from AMOR data
data(Rhizo)
data(Rhizo.map)
data(Rhizo.tax)

# IMPORTANT!!! We need to change the name of some variables to make them consistent with the
# functions in his Dataset.
Rhizo.map$Fraction <- Rhizo.map$fraction
Rhizo.map$Depth <- colSums(Rhizo)
Rhizo.map$Experiment <- 'A'
Rhizo.map$Genotype <- Rhizo.map$accession
Rhizo.map$fraction <- Rhizo.map$accession <- NULL

Dat <- create_dataset(Tab = Rhizo, Map = Rhizo.map, Tax = Rhizo.tax)
rm(Rhizo, Rhizo.map, Rhizo.tax)
Dat


# First step is to fit linear models.
# IMPORTANT!!: The following function was written with the specific dataset of this project in mind.
# It makes strong assumptions about the metadata available.


