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

#' Plot fracgen heatmap
#' 
#' First selects the variables that match the egular expression and then
#' finds the significant ones using sig_thres. Requires columns with the
#' appropriate names in Full, meaning Taxon, Variable and q.value
#' 
#' The it processes the dataset object. Rescales the abundance to percent,
#' requires a Depth variable in the Map, and removes any taxa that are never
#' significant. Removes Soil and R fraction samples (requires a Fraction
#' variable with R and Soil values in it), and pools samples by genotype
#' (requires a Genotype variable in the Map)
#' 
#' 
#' @param Full output from combine_models
#' @param Dat Dataset object
#' @param group character vector indicating experiments to keep.
#' @param variable regular expression used to extract variable of interest
#' @param sig_thres significance threshold
#' @param pool.FUN function to use when pooling observations for the heatmap
#' @param phyl.level level at which phylum is found in Tax
#' @param mut_order
plot_fragen_heatmap <- function(Full,Dat,group, variable = "FracgenEC:",
                                sig_thres = 0.05, pool.FUN = mean, phyl.level = 3,
                                mut_order){
  
  # Get significant
  Gen <- Full[ grep(variable,Full$Variable), ]
  fracgen_effect <- as.character(unique(Gen$Taxon[ Gen$q.value < sig_thres ]))
  Gen <- Gen[ Gen$Taxon %in% fracgen_effect, ]
  Gen <- droplevels(Gen)
  Gen$Taxon <- factor(Gen$Taxon)
  Gen$significant <- Gen$q.value < sig_thres
  
  
  # Process Dat
  pool.FUN <- match.fun(pool.FUN)
  # Get samples in group and rescale
  to_remove <- row.names(Dat$Map)[ !(Dat$Map$Experiment %in% group) ]
  Dat.group <- remove_samples(Dat, samples = to_remove)
  #   Dat.group <- create_dataset(Tab = t(100 * t(Dat.group$Tab) / Dat.group$Map$Depth),
  #                               Map = Dat.group$Map,
  #                               Tax = data.frame(ID = row.names(Dat.group$Tab), 
  #                                                Taxonomy = row.names(Dat.group$Tab),
  #                                                row.names = row.names(Dat.group$Tab)))
  
  Dat.group <- create_dataset(Tab = t(100 * t(Dat.group$Tab) / Dat.group$Map$Depth),
                              Map = Dat.group$Map,
                              Tax = Dat.group$Tax)
  
  
  Dat.group$Map$Genotype <- factor(Dat.group$Map$Genotype, levels = mut_order)
  # Remove non significant taxa
  Dat.group <- remove_taxons(Dat = Dat.group,
                             taxons = row.names(Dat.group$Tab)[ !(row.names(Dat.group$Tab) %in% levels(Gen$Taxon)) ])
  # Remove Soil and R samples
  Dat.group <- remove_samples(Dat = Dat.group,
                              samples = colnames(Dat.group$Tab)[ Dat.group$Map$Fraction %in% c("R","Soil") ])
  # Pool by genotype
  Dat.group <- pool_samples(Dat = Dat.group, groups = "Genotype", FUN = mean)
  
  # Create plotting object
  dat <- as.data.frame(Dat.group$Tab)
  dat$OTU <- row.names(dat)
  dat <- melt(dat,id.vars=c("OTU"),value.name="Abundance",variable.name="Genotype")
  dat$Abundance <- log2(dat$Abundance + 1)
  dat$Enrichment <- 0
  dat$Enrichment[ paste(dat$OTU,variable,dat$Genotype,sep = "") %in%
                    paste(Gen$Taxon,Gen$Variable,sep="")[ Gen$significant & Gen$Estimate > 0 ] ] <- 1
  dat$Enrichment[ paste(dat$OTU,variable,dat$Genotype,sep = "") %in%
                    paste(Gen$Taxon,Gen$Variable,sep="")[ Gen$significant & Gen$Estimate < 0 ] ] <- -1
  dat$Phylum <- get_tax_level(Dat.group$Tax[as.character(dat$OTU),],level = phyl.level)
  dat$Type <- "Mutant"
  dat$Type[ dat$Genotype == "Col-0" ] <- "WT"
  dat <- dat[order(dat$Phylum,dat$Enrichment),]
  dat$OTU <- factor(dat$OTU,levels = unique(dat$OTU))
  dat$Type <- factor(dat$Type,levels=c("WT","Mutant"))
  dat$Enrichment <- factor(dat$Enrichment,levels=c("-1","1","0"))
  dat <- dat[order(dat$Enrichment),]
  
  p1 <- ggplot(dat,aes(x=Genotype,y=OTU)) +
    #   geom_tile(aes(fill=Abundance,col=Enrichment),size=2) +
    geom_tile(aes(fill=Abundance)) +
    geom_tile(aes(col=Enrichment),size=2,alpha=1,fill = NA) +
    facet_grid(Phylum ~ Type,scales="free",space="free") +
    scale_fill_continuous(low="white",high="black") +
    scale_colour_manual(values = c("blue","red",NA)) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(face="bold",colour="black",
                                     angle = 90, size = 20),
          panel.background = element_blank(),
          strip.text.y = element_text(angle = 360))
  
  return(p1)
}


#' Doughnut plot
#' 
#' @param dat The output data attribute from plot_fracgen_heatmap
#' @param top_phyla Phyla to show explicitly. The rest will be lumped into other
#' @param mut_order Order of genotypes to show
#' @param text_size Size of text for numbers inside doughnuts
#' @param phyla_colors Vector of colors to use
doughnut_plot <- function(dat,top_phyla, mut_order, text_size = 10, phyla_colors = phyla_colors){
  # mut_order <- c('Col','Ler','mygen','Soil')
  # top_phyla <- c("Root;Root;k__Bacteria;p__Actinobacteria",
  #               "Root;Root;k__Bacteria;p__Bacteroidetes",
  #               "Root;Root;k__Bacteria;p__Proteobacteria")
  
  count <- matrix(ftable(Enrichment~Genotype,dat)[-1,-3], ncol = 2)
  row.names(count) <- mut_order[c(-1,-length(mut_order))]
  colnames(count) <- c("Depletion", "Enrichment")
  count <- melt(count,varnames = c("Genotype","Type"),value.name = "Total")
  count$Phylum <- NA
  count$ymin <- NA
  count$ymax <- NA
  
  enrich <- dat[ dat$Enrichment == 1, ]
  enrich <- calculate_enrich(enrich = enrich,top_phyla = top_phyla,
                             mut_order = mut_order)
  enrich$Type <- "Enrichment"
  deplete <- dat[ dat$Enrichment == -1, ]
  deplete <- calculate_enrich(enrich = deplete,top_phyla = top_phyla,
                              mut_order = mut_order)
  deplete$Type <- "Depletion"
  enrich <- rbind(enrich,deplete)
  enrich$Type <- factor(enrich$Type, levels = c("Enrichment","Depletion"))
  
  enrich$Total <- NA
  enrich <- rbind(enrich,count)
  
  p1 <- ggplot(enrich, aes(xmax = 4, xmin = 3)) +
    facet_grid(Type ~ Genotype) +
    geom_rect(data = subset(enrich, !is.na(ymin)), aes(ymax = ymax, ymin = ymin, fill = Phylum), color = NA) +
    xlim(c(0,4)) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = phyla_colors) +
    geom_text(data = subset(enrich, is.na(ymin)),
              aes(label = Total), x = 0, y = 0, size = text_size,
              fontface = "bold") +
    theme_blackbox +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
  #p1
  return(p1)
}

#' Take a set of enrichments or depletions and make
#' The object for geom_rect
calculate_enrich <- function(enrich, top_phyla = top_phyla, mut_order){
  #  mut_order <- pepr_order
  
  enrich$Enrichment <- 1
  enrich$Phylum[ !(enrich$Phylum %in% top_phyla) ] <- "Low Abundance"
  enrich$Phylum <- factor(enrich$Phylum, levels = top_phyla )
  enrich$Genotype <- factor(enrich$Genotype, levels = mut_order[c(-1, -length(mut_order))])
  enrich <- acast(data = enrich, formula = Phylum ~ Genotype, fun.aggregate = sum, value.var = "Enrichment", drop = FALSE)
  
  enrich <- t(t(enrich) / colSums(enrich))
  enrich[ is.na(enrich) ] <- 0
  enrich.ymax <- apply(enrich,2,cumsum)
  enrich.ymin <- apply(enrich,2,function(vec) c(0,cumsum(vec)[-length(vec)]))
  row.names(enrich.ymin) <- row.names(enrich)
  
  enrich.ymin <- melt(enrich.ymin,varnames = c("Phylum","Genotype"),value.name = "ymin")
  enrich.ymax <- melt(enrich.ymax,varnames = c("Phylum","Genotype"),value.name = "ymax")
  enrich <- enrich.ymin
  enrich$ymax <- enrich.ymax$ymax
  
  return(enrich)
}

#' Test similarity of enrichments and depletions
#' 
#' @param dat output data from plot_fracgen_heatmap
#' @param wt reference genotype
#' @param dist.method distance method from dist function
#' @param nperm Number of permutations
monte_carlo_test <- function(dat, wt = "Col-0",dist.method = "manhattan",
                             nperm = 1000){
  dat$Enrichment <- as.numeric(as.character(dat$Enrichment))
  dat <- subset(dat, Genotype != wt)
  dat$Genotype <- droplevels(dat$Genotype)
  mat <- acast(data = dat, formula = OTU ~ Genotype,
               value.var = "Enrichment")
  
  all <- mat
  #all <- all[ , -which(colnames(all) == wt) ]
  
  all.d <- dist(t(all),method=dist.method)
  all.d <- as.matrix(all.d)
  D.reps <- array(dim=c(nrow(all.d),nrow(all.d),nperm),
                  dimnames = list(colnames(all),
                                  colnames(all),
                                  1:nperm))
  for (i in 1:nperm){
    #i <- 1
    all.rep <- apply(all,2,sample)
    all.rep.d <- dist(t(all.rep),method=dist.method)
    all.rep.d <- as.matrix(all.rep.d)
    D.reps[,,i] <- all.d < all.rep.d
  }
  pval <- as.matrix(as.dist(1 - apply(D.reps,MARGIN=c(1,2),sum) / nperm))
  
}

##################

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
Rhizo.map$Fracgen <- interaction(Rhizo.map$Genotype,
                                 factor(Rhizo.map$Fraction, levels = c('E','R','Soil')),
                                 drop = TRUE)

Dat <- create_dataset(Tab = Rhizo, Map = Rhizo.map, Tax = Rhizo.tax)
rm(Rhizo, Rhizo.map, Rhizo.tax)
Dat


# First step is to fit linear models.
# IMPORTANT!!: The following function was written with the specific dataset of this project in mind.
# It makes strong assumptions about the metadata available.
# m1 <- fit_4_models(Dat = Dat, f1 = ~ soil + Genotype + Fraction + plate,
#                    group = c('A'), 
#                    mutant_order = levels(Dat$Map$Genotype),
#                    min_reads_otu = 25,
#                    min_samples_otu = 5)
# m1 <- fit_4_models(Dat = Dat, f1 = ~ soil + Genotype + Fraction + plate,
#                    group = c('A'), 
#                    mutant_order = levels(Dat$Map$Genotype),
#                    min_reads_otu = 25,
#                    min_samples_otu = 5)

# This is not an appropriate model, and it will throw a bunch of errors, but it serves as
# demonstration

# levels(Dat$Map$Fracgen)
# # Dat$Map$Fracgen <- relevel(Dat$Map$Fracgen, ref = 'E.Col')
# levels(Dat$Map$Fracgen)

m1 <- matrix_glm(x = Dat,formula = ~ soil + Fraction + Genotype + Fracgen + plate)
Full <- summary(m1)$coefficients
Full <- droplevels(subset(Full, !is.na(p.value)))
Full$q.value <- p.adjust(Full$p.value, method = 'fdr')
head(Full)

# This function, though it generates a plot, can also be used to reformat data
# NOTE WE ARE USING AN ABSURD Q-VALUE THRESHOLD (0.5) FOR DEMONSTRATION PURPOSES ONLY
p1 <- plot_fragen_heatmap(Full = Full, Dat = Dat, group = c("A"),
                          variable = "Genotype", sig_thres = 0.5,
                          pool.FUN = mean, phyl.level = 4,
                          mut_order = levels(Dat$Map$Genotype))
dat <- p1$data
head(dat)

# The way the model was specified, there is only one genotype to compare (Ler).
# I'm going create a fake genotype and assign some enrichments and depletions
# because we need at least two groups for comparisons
set.seed(193)
t <- subset(dat, Genotype == 'Col')
t$Enrichment[ sample(nrow(t),size = 5) ] <- '-1'
t$Enrichment[ sample(nrow(t),size = 5) ] <- '1'
t$Genotype <- 'mygen'
t
dat <- rbind(dat,t)

# Plot doughnut
p1 <- doughnut_plot(dat = dat,
                    top_phyla = c("Root;Root;k__Bacteria;p__Actinobacteria",
                                  "Root;Root;k__Bacteria;p__Bacteroidetes",
                                  "Root;Root;k__Bacteria;p__Proteobacteria",
                                  "Low Abundance"),
                    mut_order = c('Col','Ler','mygen','Soil'), phyla_colors = c('red', 'blue', 'orange','black'))
p1

set.seed(7151)
pval <- monte_carlo_test(dat = dat,wt = 'Col', nperm = 5000)
pval
