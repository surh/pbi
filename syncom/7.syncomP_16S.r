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

library(AMOR)

# Read data
#Dat <- read.am("~/rhizogenomics/experiments/2016/2016-03-17.syncomP/ref_tab.txt",format = "qiime") 

Dat <- read.am("reftab/otu_table.txt",format = "qiime") 
Dat.otus <- read.am("otutab/otu_table.txt", format = "qiime", taxonomy = "taxonomy")
Map <- read.table("~/rhizogenomics/data/phosphate/colonization_synthetic/metadata.txt",
                  sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# Sort taxa
Dat$Tab <- Dat$Tab[ names(sort(rowSums(Dat$Tab),decreasing = TRUE)), ]
Dat <- create_dataset(Dat$Tab)
Dat.otus$Tab <- Dat.otus$Tab[ names(sort(rowSums(Dat.otus$Tab),decreasing = TRUE)), ]
Dat.otus <- create_dataset(Tab = Dat.otus$Tab,
                           Tax = Dat.otus$Tax[row.names(Dat.otus$Tab),])
all(colnames(Dat$Tab) %in% colnames(Dat.otus$Tab))
all(colnames(Dat.otus$Tab) %in% colnames(Dat$Tab))

# Sort samples
Dat.otus$Tab <- Dat.otus$Tab[,colnames(Dat$Tab)]
row.names(Map) <- Map$ID
Map <- Map[ colnames(Dat$Tab), ]

# Taxonomy
Dat$Tax <- data.frame(ID = row.names(Dat$Tab), Taxonomy = "Isolate",
                      stringsAsFactors = FALSE)
row.names(Dat$Tax) <- Dat$Tax$ID
Dat$Tax$Taxonomy[ Dat$Tax$ID == "Atplastid" ] <- "Plastid"
Dat$Tax$Taxonomy[ Dat$Tax$ID == "Atmitochondria" ] <- "Mitochondria"
Dat$Tax$Taxonomy[ Dat$Tax$ID == "47Yellow" ] <- "Contaminant"
Dat.otus$Tax$Taxonomy <- paste("OTU; ",as.character(Dat.otus$Tax$Taxonomy), sep = "")


# Master Dataset
Dat <- create_dataset(Tab = rbind(Dat$Tab,Dat.otus$Tab), Map = Map,
                      Tax = rbind(Dat$Tax, Dat.otus$Tax))

Dat$Map$Depth <- colSums(Dat$Tab)
Dat$Map$Fraction <- factor(Dat$Map$Fraction,
                                levels = c("Inoculum","Agar","Root", "Blank"))
Dat$Map$Pi[ Dat$Map$Pi == "LowP" ] <- "50uM"
Dat$Map$Pi[ Dat$Map$Pi == "FullP" ] <- "625uM"
Dat <- clean(Dat)
rm(Dat.otus)

########## Plots full  ############
# Plot depth
p1 <- plotgg_var(Dat,var.name = "Depth", x = "Fraction", col = "Bacteria") +
  scale_y_log10()
p1

p1 <- phylogram(collapse_by_taxonomy(Dat,level = 1),facet = ~ Fraction + Bacteria) +
  theme(strip.text = element_text(angle = 90))
ggsave("phylogram_full_type.svg",p1,width = 10, height = 5)

# Plot usable reads
Dat <- remove_taxons(Dat, taxons = c("Atmitochondria","Atplastid"))
Dat <- remove_taxons(Dat, taxons = Dat$Tax$ID[ Dat$Tax$Taxonomy == "OTU; Root; k__Bacteria"])
Dat$Map$Usable <- colSums(Dat$Tab)
Dat <- clean(Dat)

p1 <- plotgg_var(Dat,var.name = "Usable", x = "Fraction", col = "Bacteria") +
  scale_y_log10()
p1
ggsave("usable_reads.svg",p1,width = 5, height = 4)

p1 <- phylogram(collapse_by_taxonomy(Dat,level = 1),facet = ~ Fraction + Bacteria) +
  theme(strip.text = element_text(angle = 90))
ggsave("phylogram_usable_type.svg",p1,width = 10, height = 5)

# Heatmap blank norm
# Normalize by blanks
Dat.freq <- Dat
Dat.freq$Tab <- t(t(Dat.freq$Tab) - round(rowMeans(subset(Dat.freq,Fraction == "Blank")$Tab)))
Dat.freq$Tab[ Dat.freq$Tab < 0 ] <- 0
Dat.freq$Map$Usable.norm <- colSums(Dat.freq$Tab)

Dat.freq$Map$Fraction <- factor(Dat.freq$Map$Fraction,
                                levels = c("Inoculum","Agar","Root"))
Dat.freq$Map$Genotype <- factor(Dat.freq$Map$Genotype,
                                levels = c("Col-0","phf1","phr1/phl1","Agar","Inoculum"))
Dat.perc <- normalize(Dat.freq, norm = "Usable.norm")

write.qiime(Dat.freq,file = "syncomP_counts.txt")
write.table(Dat.freq$Map,file = "syncomP_metadata", col.names = NA,
            sep = "\t", quote = F)

p1 <- plotgg_var(subset(Dat.freq,Fraction %in% c("Inoculum","Agar","Root")),
                 var.name = "Usable.norm", x = "Fraction", col = "Bacteria") +
  scale_y_log10()
#p1
ggsave("syncomP_usable.norm.svg",p1,width = 5,height = 4)

p1 <- phylogram(collapse_by_taxonomy(Dat.freq,level = 1),
                facet = ~ Fraction + Bacteria) +
  theme(strip.text = element_text(angle = 90))
p1

# Sort
Dat.perc <- create_dataset(Tab = Dat.perc$Tab[ names(sort(rowMeans(Dat.perc$Tab),decreasing = TRUE)), ],
                           Map = Dat.perc$Map,
                           Tax = Dat.perc$Tax[ names(sort(rowMeans(Dat.perc$Tab),decreasing = TRUE)), ])

# Identify colonizers
col_thres <- mean(subset(Dat.perc, Fraction == "Root")$Tab["Ecoli",])
nsam <- ncol(subset(Dat.perc,Fraction == "Root")$Tab)
count <- rowSums(subset(Dat.perc,Fraction == "Root")$Tab > col_thres)
colonizers <- 1 - pbinom(q = count - 1, size = nsam, prob = 0.5)
#colonizers <- p.adjust(colonizers,method = "fdr")
colonizers <- colonizers[ colonizers < 0.05 ]
colonizers <- names(colonizers)
colonizers

p1 <- heatgg(subset(Dat.perc, Bacteria == "+Bacteria", drop = TRUE))
p1$data$Colonizer <- "Non-Colonizer"
p1$data$Colonizer[ p1$data$Taxon %in% colonizers ] <- "Colonizer"
p1 <- p1 + facet_grid(Genotype + posHarvest ~ Colonizer, scales = "free",space = "free") +
  theme(axis.text.x = element_text(angle = 90, color = "black", face = "bold"),
        strip.text.y = element_text(angle = 0, face = "bold"),
        strip.text.x = element_text(angle = 0, face = "bold"),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.y = element_blank())
#p1
ggsave("syncomP_heatmap_relative_abundances_all.svg",p1,width = 30, height = 8)

p1 <- heatgg(subset(x = remove_taxons(Dat.perc,
                                      taxons = Dat.perc$Tax$ID[ Dat.perc$Tax$Taxonomy != "Isolate" ]),
                    Bacteria == "+Bacteria", drop = TRUE))
p1$data$Colonizer <- "Non-Colonizer"
p1$data$Colonizer[ p1$data$Taxon %in% colonizers ] <- "Colonizer"
p1 <- p1 + facet_grid(Genotype + posHarvest ~ Colonizer, scales = "free",space = "free") +
  ylab(label = "Sample") +
  xlab(label = "Isolate") +
  theme(axis.text.x = element_text(angle = 90, color = "black", face = "bold",
                                   size = 12),
        strip.text.y = element_text(angle = 0, face = "bold"),
        strip.text.x = element_text(angle = 0, face = "bold", size = 12),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.y = element_blank())
#p1
ggsave("syncomP_heatmap_relative_abundances_isolates.svg",p1,
       width = 10, height = 8)


################# Model
# Select data for modes
Dat.freq <- subset(Dat.freq, Fraction != "Blank")
Dat.freq <- subset(Dat.freq, Bacteria == "+Bacteria")
Dat.freq <- subset(Dat.freq, posHarvest != "Overnight")
Dat.freq <- subset(Dat.freq, !(Fraction == "Agar" & Day == 0),
                   clean = TRUE,drop = TRUE)
ftable(Genotype + Experiment ~ Fraction,Dat.freq$Map)

# Subset percentage table for plotting
Dat.perc <- subset(Dat.perc, ID %in% colnames(Dat.freq$Tab),
                   clean = TRUE, drop = TRUE)
p1 <- heatgg(Dat.perc)
p1$data$Colonizer <- "Non-Colonizer"
p1$data$Colonizer[ p1$data$Taxon %in% colonizers ] <- "Colonizer"
p1 <- p1 + facet_grid(Genotype + Pi~ Colonizer, scales = "free",space = "free") +
  ylab(label = "Sample") +
  xlab(label = "Isolate") +
  theme(axis.text.x = element_text(angle = 90, color = "black", face = "bold",
                                   size = 12),
        strip.text.y = element_text(angle = 0, face = "bold"),
        strip.text.x = element_text(angle = 0, face = "bold", size = 12),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.y = element_blank())
#p1
ggsave("syncomP_heatmap_relative_abundances_all_model.svg",
       p1,width = 30, height = 8)

p1 <- heatgg(remove_taxons(Dat.perc,
                           taxons = Dat.perc$Tax$ID[ Dat.perc$Tax$Taxonomy != "Isolate" ]))
p1$data$Colonizer <- "Non-Colonizer"
p1$data$Colonizer[ p1$data$Taxon %in% colonizers ] <- "Colonizer"
p1 <- p1 + facet_grid(Genotype + Pi~ Colonizer, scales = "free",space = "free") +
  ylab(label = "Sample") +
  xlab(label = "Isolate") +
  theme(axis.text.x = element_text(angle = 90, color = "black", face = "bold",
                                   size = 12),
        strip.text.y = element_text(angle = 0, face = "bold"),
        strip.text.x = element_text(angle = 0, face = "bold", size = 12),
        strip.background = element_rect(fill = NA, color = "black"),
        axis.text.y = element_blank())
#p1
ggsave("syncomP_heatmap_relative_abundances_isolates_model.svg",
       p1,width = 10, height = 8)


# Add 1 for scale plotting 
Dat.perc$Tab <- Dat.perc$Tab + 1

# Transform depths for offset
# Not needed for edgeR
# Dat.freq$Map$Usable <- log(Dat.freq$Map$Usable)
# Dat.freq$Map$Usable.norm <- log(Dat.freq$Map$Usable.norm)

# Define formulas
# Not needed for edgeR
# f1 <- formula(~ Fraction + Genotype + Pi + Experiment + offset(Usable))
# f2 <- formula(~ Fraction + Genotype + Pi + Experiment + offset(Usable.norm))

library(edgeR)
# Prepare design matrix
f1 <- formula(~ Fraction + Genotype + Pi + Experiment)
design1 <- model.matrix(f1,data = Dat.freq$Map)
design1 <- design1[ , !(colnames(design1) %in% c("GenotypeAgar","GenotypeInoculum","PiInoculum"))]

# Try different normalization methods and estimate overdispersion parameters
dge1 <- DGEList(counts = Dat.freq$Tab,norm.factors = Dat.freq$Map$Usable.norm)
dge1 <- estimateDisp(dge1,design = design1)

# dge2 <- DGEList(counts = Dat.freq$Tab,norm.factors = Dat.freq$Map$Usable)
# dge2 <- estimateDisp(dge2,design = design1)

# Fit models with Quasi-likelihood
dge1.fit <- glmQLFit(dge1,design = design1)
# dge2.fit <- glmQLFit(dge2,design = design1)

# Test results using blank normalized offset
write.table(topTags(glmQLFTest(dge1.fit,coef = 4),n = nrow(dge1$counts)),
            file = "phf1.col.edgeR.txt",sep = "\t", quote = FALSE, col.names = NA)
write.table(topTags(glmQLFTest(dge1.fit,coef = 5),n = nrow(dge1$counts)),
            file = "phr1phl1.col.edgeR.txt",sep = "\t", quote = FALSE, col.names = NA)
write.table(topTags(glmQLFTest(dge1.fit,contrast = c(0,0,0,1,-1,0,0)),n = nrow(dge1$counts)),
            file = "phf1.phr1phl1.edgeR.txt",sep = "\t", quote = FALSE, col.names = NA)
write.table(topTags(glmQLFTest(dge1.fit,contrast = c(0,1,-1,0,0,0,0)),n = nrow(dge1$counts)),
            file = "agar.root.edgeR.txt",sep = "\t", quote = FALSE, col.names = NA)
write.table(topTags(glmQLFTest(dge1.fit,coef = 6),n = nrow(dge1$counts)),
            file = "full.low.edgeR.txt",sep = "\t", quote = FALSE, col.names = NA)

# Test results using usable offset
# topTags(glmQLFTest(dge2.fit,coef = 4),n = 33)
# topTags(glmQLFTest(dge2.fit,coef = 5),n = 33)
# topTags(glmQLFTest(dge2.fit,coef = 6),n = 33)

# Plotting all isolates
for(strain in Dat.perc$Tax$ID[ Dat.perc$Tax$Taxonomy == "Isolate" ]){
  #strain <- "376"
  p1 <- plotgg_taxon(subset(Dat.perc,Fraction == "Root"),
                     taxon = strain, x = "Genotype", col = "Pi") +
    scale_y_log10(breaks = seq(from = 1,
                               to = max(subset(Dat.perc,Fraction == "Root")$Tab[strain,]),
                               length.out = 5)) +
    ylab("% Abundance") +
    #scale_y_log10(breaks = c(1,5,10,15,20,30,40,50,60,70)) +
    ggtitle(strain)
  #p1
  filename <- paste("syncomP_",strain,'.svg', sep = "")
  ggsave(filename,p1,width = 5,height = 4)
}

strain <- "2"
p1 <- plotgg_taxon(subset(Dat.perc, Fraction != "Inoculum"),
                   taxon = strain, x = "Genotype", col = "Pi") +
  scale_y_log10(breaks = seq(from = 1,
                             to = max(subset(Dat.perc,Fraction == "Root")$Tab[strain,]),
                             length.out = 5)) +
  ylab("% Abundance") +
  #scale_y_log10(breaks = c(1,5,10,15,20,30,40,50,60,70)) +
  ggtitle(strain)
#p1
filename <- paste("syncomP_",strain,'_fraction.svg', sep = "")
ggsave(filename,p1,width = 5,height = 4)

############## Ordination
source("~/rhizogenomics/src/trunk/phosphate_code/functions.r")
# Rarefy for estimating richness and PCA/PCoA
Dat.rar <- Dat.freq
Dat.rar <- remove_taxons(Dat.rar,Dat.rar$Tax$ID[ Dat.rar$Tax$Taxonomy != "Isolate" ])
Dat.rar$Map$rarcounts <- colSums(Dat.rar$Tab)
Dat.rar <- subset(Dat.rar, rarcounts > 1500, drop = TRUE, clean = TRUE)
set.seed(35241)
Dat.rar <- rarefaction(Dat.rar,sample = 1500)

# Calculate, plot and test richness
Dat.rar$Map$Richness <- colSums(Dat.rar$Tab > 0)
Dat.rar$Map$Shannon <- vegan::diversity(t(Dat.rar$Tab))
p1 <- plotgg_var(Dat.rar, var.name = "Richness", x = "Fraction",col = "Genotype") +
  scale_colour_manual(values = c(gen_colors[1], gen_colors[4], gen_colors[7],
                                 "grey15","magenta")) +
  scale_fill_manual(values = c(gen_colors[1], gen_colors[4], gen_colors[7],
                               "grey15","magenta"))
p1
ggsave("syncomP_richness_model.svg",p1,width = 5,height = 5)


Dat.rar$Map$new <- Dat.rar$Map$Genotype
Dat.rar$Map$new[ Dat.rar$Map$new == "phf1" ] <- "phr1/phl1"  

summary(aov(Richness ~ new + Experiment, data = subset(Dat.rar$Map,Fraction == "Root")))
TukeyHSD(aov(Richness ~ new + Experiment, data = subset(Dat.rar$Map,Fraction == "Root")))

summary(aov(Shannon ~ new + Experiment, data = subset(Dat.rar$Map,Fraction == "Root")))
TukeyHSD(aov(Shannon ~ new + Experiment, data = subset(Dat.rar$Map,Fraction == "Root")))


  
summary(aov(Richness ~ Genotype + Experiment, data = subset(Dat.rar$Map,Fraction == "Root")))
TukeyHSD(aov(Richness ~ Genotype + Experiment, data = subset(Dat.rar$Map,Fraction == "Root")))

summary(aov(Richness ~ Genotype*Pi, data = subset(Dat.rar$Map,Fraction == "Root")))


summary(aov(Richness ~ Genotype, data = subset(Dat.rar$Map,Fraction == "Root")))
TukeyHSD(aov(Richness ~ Genotype, data = subset(Dat.rar$Map,Fraction == "Root")))


p1 <- plotgg_var(Dat.rar, var.name = "Shannon", x = "Fraction",col = "Genotype") +
  scale_colour_manual(values = c(gen_colors[1], gen_colors[4], gen_colors[7],
                                 "grey15","magenta")) +
  scale_fill_manual(values = c(gen_colors[1], gen_colors[4], gen_colors[7],
                               "grey15","magenta"))
p1
ggsave("syncomP_shannon_model.svg",p1,width = 5,height = 5)
summary(aov(Shannon ~ Genotype + Experiment, data = subset(Dat.rar$Map,Fraction == "Root")))
TukeyHSD(aov(Shannon ~ Genotype + Experiment, data = subset(Dat.rar$Map,Fraction == "Root")))

summary(aov(Shannon ~ Genotype, data = subset(Dat.rar$Map,Fraction == "Root")))
TukeyHSD(aov(Shannon ~ Genotype, data = subset(Dat.rar$Map,Fraction == "Root")))

# PCA
p1 <- plotgg(PCA(Dat.rar, cor = TRUE),col="Fraction",point_size = 3)
(100*(PCA(Dat.rar, cor = TRUE)$sdev^2) / PCA(Dat.rar, cor = TRUE)$totvar)[1:5]
p1 <- p1 + scale_color_manual(values = c("magenta","grey15","chartreuse4")) +
  xlab("PC1 (28.52%)") + ylab("PC2 (11.51%)")
p1
ggsave("syncomP_cor_model_rar1500.svg",p1, width = 5,height = 5)

#PCoA
p1 <- plotgg(PCO(Dat.rar,dim = 2,distfun = function(x) vegan::vegdist(x,method = "cao")),
             col = "Fraction",point_size = 3)
eig <- PCO(Dat.rar,dim = 2,distfun = function(x) vegan::vegdist(x,method = "cao"))$eig
(100*eig[eig > 0] / sum(eig[eig>0]))[1:10]
p1 <- p1 + scale_color_manual(values = c("magenta","grey15","chartreuse4")) +
  xlab("PCo1 (31.22%)") + ylab("PC2 (13.44%)")
p1
ggsave("syncomP_pco_model_rar1500_cao.svg",p1, width = 5,height = 5)

# CAP analysis
Dat.cap <- subset(Dat.rar, Fraction == "Root",drop = TRUE, clean = TRUE)
cap.gen <- vegan::capscale(vegan::vegdist(t(Dat.cap$Tab),method = "cao") ~ Genotype + Condition(Pi) + Condition(Fraction) + Condition(Experiment) + Condition(Usable.norm),
                           data = Dat.cap$Map)
p1 <- ggplot(as.data.frame(vegan::scores(cap.gen)$sites),
             aes(x=CAP1,y=CAP2)) +
  geom_point(aes(col = Dat.cap$Map$Genotype),size = 3) +
  scale_color_manual(values = c(gen_colors[c(1,4,7)])) +
  ggtitle("Constrained variance = 5.7%") +
  xlab("CAP1 (67.42%)") + ylab("CAP2 (32.58%)") +
  theme_blackbox
p1

ggsave("syncomP_cap_gen_rar1500_cao.svg",p1,width = 7,height = 5)


pal <- colorRampPalette(colors = rev(c("#00441b","#99d8c9")))
Dat.cap <- subset(Dat.rar, Fraction != "Inoculum",drop = TRUE, clean = TRUE)
cap.pi <- vegan::capscale(vegan::vegdist(t(Dat.cap$Tab),method="cao") ~ Pi + Condition(Genotype) + Condition(Fraction) + Condition(Experiment) + Condition(Usable.norm),
                           data = Dat.cap$Map)
p1 <- ggplot(as.data.frame(vegan::scores(cap.pi)$sites),
             aes(x=CAP1,y=MDS1)) +
  geom_point(aes(col = Dat.cap$Map$Pi),size = 3) +
  scale_color_manual(values = pal(5)[c(3,5)]) +
  ggtitle("Constrained variance = 5.5%") +
  xlab("CAP1 (100%)") + ylab("MDS1 (19.06%)") +
  theme_blackbox
p1
ggsave("syncomP_cap_pi_rar1500_cao.svg",p1,width = 6,height = 5)


Dat.cap <- subset(Dat.rar, Fraction != "Inoculum",drop = TRUE, clean = TRUE)
cap.frac <- vegan::capscale(vegan::vegdist(t(Dat.cap$Tab),method = "cao") ~ Fraction + Condition(Pi) + Condition(Experiment) + Condition(Usable.norm),
                          data = Dat.cap$Map)
p1 <- ggplot(as.data.frame(vegan::scores(cap.frac)$sites),
             aes(x=CAP1,y=MDS1)) +
  geom_point(aes(col = Dat.cap$Map$Frac),size = 3) +
  scale_color_manual(values = c("grey15","chartreuse4")) +
  ggtitle("Constrained variance = 9.1%") +
  xlab("CAP1 (100%)") + ylab("MDS1 (17.92%)") +
  theme_blackbox
p1

ggsave("syncomP_cap_frac_rar1500_cao.svg",p1,width = 6,height = 5)
