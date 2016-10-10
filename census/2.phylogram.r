library(AMOR)
datadir <- "data/"

top_phyla <- c("Root; k__Bacteria; p__Acidobacteria",
               "Root; k__Bacteria; p__Actinobacteria",
               "Root; k__Bacteria; p__Bacteroidetes",
               "Root; k__Bacteria; p__Chloroflexi",
               "Root; k__Bacteria; p__Firmicutes",
               "Root; k__Bacteria; p__Proteobacteria",
               "Root; k__Bacteria; p__Verrucomicrobia",
               "Low Abundance",
               "unclassified")
phyla_names <- c("Acidobacteria","Actinobacteria","Bacteroidetes",
                 "Chloroflexi","Firmicutes","Proteobacteria",
                 "Verrucomicrobia","Low Abundance","Unknown")
phyla_colors <- c("#41F0AC","#0000C0","#FFFF00","#008000",
                  "#B856D7","#FF8000","#8D8DFF","#000000",
                  "#FFFFFF")


# Read and process
filename <- paste(datadir,"/Dat.filter.processed.rdat",sep = "")
load(filename)

# Get phylum dist anr rename
Dat.phyl <- collapse_by_taxonomy(Dat = Dat,level = 3)
Tax <- data.frame(ID = row.names(Dat.phyl$Tab), stringsAsFactors = FALSE)
Tax$Taxonomy <- Tax$ID
Tax$Taxonomy[ !(Tax$Taxonomy %in% top_phyla) ] <- "Low Abundance"
names(phyla_names) <- top_phyla
Tax$Taxonomy <- phyla_names[Tax$Taxonomy]
Tax$Taxonomy <- factor(Tax$Taxonomy,levels = phyla_names)
row.names(Tax) <- Tax$ID
Dat.phyl <- create_dataset(Tab = Dat.phyl$Tab, Map = Dat$Map, Tax = Tax)
Dat.phyl <- collapse_by_taxonomy(Dat = Dat.phyl, level = 1)

p1 <- phylogram(Dat = Dat.phyl, facet = . ~ Genotype,ntaxa = 9,nrow.legend = 9)
p1$data$Taxon <- factor(p1$data$Taxon,levels = phyla_names)
p1 <- p1 + scale_fill_manual(values = phyla_colors)
p1 <- p1 + theme(strip.text = element_text(angle = 90),
                 axis.text.x = element_blank())
#p1
ggsave("phylogram_allsamples.png",width = 10,height = 5)
ggsave("phylogram_allsamples.svg",width = 10,height = 5)

Dat.phyl.pool <- pool_samples(Dat = Dat.phyl,groups = "Genotype",FUN = sum)
Dat.phyl.pool$Map <- data.frame(Genotype = factor(colnames(Dat.phyl.pool$Tab),
                                                  levels = c("Col-0","pht1;1",
                                                             "pht1;1/pht1;4","phf1",
                                                             "nla","pho2","pho1:2",
                                                             "35","phr1","spx1/spx2",
                                                             "Soil")))
row.names(Dat.phyl.pool$Map) <- colnames(Dat.phyl.pool$Tab)
p1 <- phylogram(Dat = Dat.phyl.pool, ntaxa = 11, nrow.legend = 12)
p1$data$Taxon <- factor(p1$data$Taxon,levels = phyla_names)
p1 <- p1 + scale_fill_manual(values = phyla_colors)
p1$data$Sample <- factor(p1$data$Sample,levels = levels(Dat.phyl.pool$Map$Genotype))
p1 <- p1 + facet_grid(facets = . ~ Sample,scales = "free")
p1 <- p1 + theme(strip.text = element_text(angle = 90),
                 axis.text.x = element_blank())
#p1
ggsave("phylogram_pooled.png",width = 10,height = 5)
ggsave("phylogram_pooled.svg",width = 10,height = 5)
