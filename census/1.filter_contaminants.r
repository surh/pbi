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
outdir <- "data/"

otutab_file <- "~/rhizogenomics/data/phosphate/census/data/otu_table.txt"

# Read data
Dat <- read.am(file = otutab_file,taxonomy = "taxonomy", format = "qiime")

contam_otus <- c(grep(pattern = "chloroplast", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "mitochondri", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "oomycete", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 which(Dat$Tax$Taxonomy == "Root; k__Bacteria"),
                 which(Dat$Tax$Taxonomy == "Root"))
contam_otus <- row.names(Dat$Tax)[contam_otus]

# Filter Dataset
Dat.filter <- remove_taxons(Dat = Dat, taxons = contam_otus)
Dat.filter <- clean(Dat = Dat.filter,verbose = TRUE)
contam_otus <- Dat$Tax[ contam_otus, ]

# Write resutls of filtering
dir.create(outdir)

# Log of removed files
filename <- paste(outdir,"/tax_removed.txt",sep = "")
write.table(contam_otus,file = filename, col.names = NA, quote = FALSE, sep = "\t")

# Write filtered table
filename <- paste(outdir,"/otu_table_filter.txt",sep = "")
write.qiime(Tab = Dat.filter$Tab, file = filename)
filename <- paste(outdir,"/tax_filtered.txt",sep = "")
write.table(Dat.filter$Tax,file = filename)

# Clean
rm(Dat,Dat.filter,contam_otus,filename,otutab_file)
gc()

################## Create dataset ################
filename <- paste(outdir,"/otu_table_filter.txt",sep = "")
Tab <- read.am(filename,format = "qiime",taxonomy = FALSE)$Tab
map_file <- "~/rhizogenomics/data/phosphate/census/data/map.txt"
Map <- read.table(map_file,sep="\t",row.names = 1,header = TRUE)
filename <- paste(outdir,"/tax_filtered.txt",sep = "")
Tax <- read.table(file = filename)

# convert into dataset
Map <- Map[ (row.names(Map) %in% colnames(Tab)), ]
Tab <- Tab[ ,colnames(Tab) %in% row.names(Map)]
Map <- Map[colnames(Tab),]
Dat <- create_dataset(Tab = Tab,Map = Map, Tax = Tax)
Dat$Map$Depth <- colSums(Dat$Tab)

# Process
# We remove genotypes that won't be used and blanks and unknowns
Dat <- remove_samples(Dat = Dat,
                      samples = row.names(Map)[ Dat$Map$Genotype %in% c("unk","blank","pho1:2","35")],
                      droplevels = TRUE)
Dat <- clean(Dat)
Dat$Map$Genotype <- factor(Dat$Map$Genotype,
                           levels = c("Col-0","pht1;1","pht1;1/pht1;4","phf1",
                                      "nla","pho2","phr1","spx1/spx2","Soil"))
Dat$Map$Fraction <- relevel(Dat$Map$Fraction,ref = "Soil")
Dat <- remove_samples(Dat,row.names(Dat$Map)[ Dat$Map$Depth < 1000 ])
Dat <- clean(Dat)
ftable(Genotype ~ Experiment, Dat$Map)
filename <- paste(outdir,"/Dat.filter.processed.rdat",sep = "")
save(Dat,file = filename)

### rarefy for UniFrac  ##
set.seed(2437)
Dat.rar <- rarefaction(x = Dat,sample = 1000)
Dat.rar <- clean(Dat.rar)
filename <- paste(outdir,"/otu_table_filter_rar1000.txt",sep = "")
write.qiime(Tab = Dat.rar$Tab,file = filename)
filename <- paste(outdir,"/Dat.filter.processed.rar1000.rdat",sep = "")
save(Dat.rar,file = filename)

