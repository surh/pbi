date()
library(AMOR)
library(AMORpscl)
library(gtable)
source("~/rhizogenomics/src/trunk/immune_code/PRR_Propep/functions.r")
source("~/rhizogenomics/src/trunk/phosphate_code//functions.r")

gtable_select <- function (x, ...) {
  matches <- c(...)
  x$layout <- x$layout[matches, , drop = FALSE]
  x$grobs <- x$grobs[matches]
  x
}

gtable_stack <- function(g1, g2){
  g1$grobs <- c(g1$grobs, g2$grobs)
  g1$layout <- transform(g1$layout, z= z-max(z), name="g2")
  g1$layout <- rbind(g1$layout, g2$layout)
  g1
}

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

###### Define global variables
#datadir <- "~/rhizogenomics/experiments/2016/2016-02-20.phosphateMFfam/data/"
datadir <- "data/"
###### Load data
filename <- paste(datadir,"/Dat.filter.processed.rdat",sep = "")
load(filename)

###### Prepare for model
Dat$Map$Fraction <- factor(Dat$Map$Fraction,levels = c("Soil","EC") )
#levels(Dat$Map$Genotype)[1] <- "Col-0"

# Uncomment for family level model
#Dat <- collapse_by_taxonomy(Dat,level = 6,FUN = sum)

Dat <- remove_taxons(Dat = Dat,
                     taxons = row.names(Dat$Tab)[ !findGoodOTUs(OTUtab = Dat$Tab,
                                                                min_reads_otu = 25,
                                                                min_samples_otu = 5,
                                                                method = "absolute") ])
Dat <- clean(Dat,verbose = TRUE)
Dat$Map$LDepth <- log(Dat$Map$Depth)

###### Run models
f1 <- formula(~ Fraction + Genotype + Experiment + offset(LDepth))
models <- fit_4_models(f1 = f1, Dat = Dat)
Full <- summary_4models(models)

Full[ Full$q.value < 1 & Full$Variable == "GenotypeSoil", ]

c("OTU_516","OTU_540")

# m1 <- glm(Dat$Tab["OTU_516",] ~ Fraction + Genotype + Experiment + offset(LDepth),
#           data = Dat$Map,
#           family = poisson(link = "log"))
# m2 <- glm.nb(Dat$Tab["OTU_516",] ~ Fraction + Genotype + Experiment + offset(LDepth),
#              data = Dat$Map)
# 
# summary(m1)

# f1 <- formula(~ Genotype + Experiment + offset(LDepth))
# models2 <- fit_4_models(f1 = f1, Dat = Dat)
# Full2 <- summary_4models(models2)

table(Full$model)
table(Full$Variable[ Full$q.value < 0.05 ])
# table(Full2$Variable[ Full2$q.value < 0.05 ])

filename <- paste("phosphate_combined.summary.txt",sep = "")
write.table(Full,filename, sep = "\t", quote = FALSE, row.names = FALSE)

###### Model plots
filename <- paste("phosphate_combined.summary.txt",sep = "")
Full <- read.table(filename, sep = "\t", header = TRUE)

p1 <- plot_fragen_heatmap(Full = Full, Dat = Dat, group = c("GC2","GC1"),
                          variable = "Genotype", sig_thres = 0.05,
                          pool.FUN = mean, phyl.level = 3,
                          mut_order = levels(Dat$Map$Genotype))
p1
dat <- p1$data
dat$Phylum[ !(dat$Phylum %in% top_phyla) ] <- "Low Abundance"
names(phyla_names) <- top_phyla
dat$Phylum <- phyla_names[ dat$Phylum ]
dat$Phylum <- factor(dat$Phylum,levels = phyla_names)
levels(dat$Enrichment) <- c("Decrease","Increase","none")
dat$Abundance <- (2^dat$Abundance) - 1


p1 <- ggplot(dat,aes(x=Genotype,y=OTU)) +
  #   geom_tile(aes(fill=Abundance,col=Enrichment),size=2) +
  geom_tile(aes(fill=Abundance)) +
  geom_tile(aes(col=Enrichment),size=0.2,alpha=1,fill = NA) +
  facet_grid(Phylum ~ Type,scales="free",space="free") +
  scale_fill_continuous(low="white",high="black",
                        trans = "sqrt",na.value = "white") +
  scale_colour_manual(values = c("blue","red",NA)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(face="bold",colour="black",
                                   angle = 90, size = 10),
        panel.background = element_blank(),
        strip.text.y = element_text(color = NA))
p1

dummy <- ggplot(dat, aes(x = Genotype, y = OTU)) +
  facet_grid(Phylum ~ Type, scales = "free",space = "free") + 
  geom_rect(aes(fill=Phylum), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  scale_fill_manual(values = phyla_colors) +
  theme_minimal() +
  theme(strip.text.y = element_text(color = NA),
        axis.text.y = element_blank())
dummy

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(dummy)
panels <- grepl(pattern="panel", g2$layout$name) & g2$layout$l == max(g2$layout$l[ g2$layout$name == "panel"])
strips <- grepl(pattern="strip-right", g2$layout$name)
g2$layout$r[panels] <- g2$layout$r[panels] + 1
g2$layout$l[panels] <- g2$layout$l[panels] + 1
new_strips <- gtable_select(g2, panels | strips)

# grid.newpage()
# grid.draw(new_strips)

new_plot <- gtable_stack(g1, new_strips)
svg("otu_fracgen_heatmap.svg",width = 7, height = 10)
grid.newpage()
grid.draw(new_plot)
dev.off()

new_plot <- gtable_stack(g1, new_strips)
png("otu_fracgen_heatmap.png",width = 1050  , height = 1500)
grid.newpage()
grid.draw(new_plot)
dev.off()

tab <- dcast(data = dat, formula = OTU ~ Genotype, value.var = "Enrichment")
filename <- "gen_heatmap_phosphate.txt"
write.table(tab, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)

# Doughnut plot
dat <- dat
levels(dat$Enrichment) <- c("-1","1","0")
p1 <- doughnut_plot(dat = dat,top_phyla = phyla_names,
                    mut_order = levels(Dat$Map$Genotype))
#p1
filename <- paste("otu_doughnut_phosphate.png",sep = "")
ggsave(filename,p1,width = 10, height = 6)
filename <- paste("otu_doughnut_phosphate.svg",sep = "")
ggsave(filename,p1,width = 10, height = 6)


# Monte carlo testing
set.seed(7151)
pval <- monte_carlo_test(dat = dat,nperm = 5000)
filename <- paste("otu_overlap_pval_phosphate.txt",sep = "")
write.table(pval,file=filename,sep="\t",
            quote=FALSE,col.names = NA)

p1 <- ggplot(melt(pval,value.name = "p.value"),aes(x = Var1, y = Var2)) + 
  geom_tile(aes(fill = p.value < 0.05)) +
  geom_text(aes(label = round(p.value,3)), size = 2) +
  theme_blackbox +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90))
#p1
filename <- paste("otu_overlap_pval_phosphate.png",sep = "")
ggsave(filename,p1, width = 5,height = 4)
filename <- paste("otu_overlap_pval_phosphate.svg",sep = "")
ggsave(filename,p1, width = 5,height = 4)


p1 <- ggplot(subset(Full[grep("Genotype",Full$Variable),],q.value < 0.05 & Variable != "GenotypeSoil"),
             aes(x = Variable, y = Estimate,color = Estimate > 0, fill = Estimate > 0)) +
  geom_boxplot(fill = NA, size = 2,
               position = position_dodge(width = 0.9),outlier.colour = NA) +
  geom_point(size = 4, shape = 21, col = "black",
             position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1)) +
  theme_blackbox +
  theme(axis.text.x = element_text(angle = 90) )
#p1
filename <- paste("otu_significant_gen_effects_phosphate.png",sep = "")
ggsave(filename,p1, width = 5,height = 4)
filename <- paste("otu_significant_gen_effects_phosphate.svg",sep = "")
ggsave(filename,p1, width = 5,height = 4)


p1 <- ggplot(subset(Full[grep("Genotype",Full$Variable),],q.value < 2 & Variable != "GenotypeSoil"),
             aes(x = Variable, y = Estimate,color = Estimate > 0, fill = Estimate > 0)) +
  geom_boxplot(fill = NA, size = 2,
               position = position_dodge(width = 0.9),outlier.colour = NA) +
  geom_point(size = 4, shape = 21, col = "black",
             position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1)) +
#  ylim(c(-3.5,3.5)) +
  theme_blackbox +
  theme(axis.text.x = element_text(angle = 90) )
p1

p1 <- ggplot(droplevels(subset(Full[grep("Genotype",Full$Variable),],q.value < 2 & Variable != "GenotypeSoil")),
             aes(x = Estimate,color = Estimate > 0, fill = Estimate > 0)) +
  facet_grid(Variable ~ .) +
  geom_density(alpha = 0.5) +
  
#  xlim(c(-3.5,3.5)) +
  theme_blackbox +
  theme(strip.text.y = element_text(angle = 360))
p1


################# For family #################
# For fam
filename <- paste("~/rhizogenomics/experiments/2016/2016-02-20.phosphateMFfam/phosphate_combined.summary.txt",sep = "")
Full <- read.table(filename, sep = "\t", header = TRUE)
Full$Variable <- gsub(pattern = "pht1:1",
                      replacement = "pht1;1",
                      as.character(Full$Variable))
Full$Variable <- gsub(pattern = "pht1:4",
                      replacement = "pht1;4",
                      as.character(Full$Variable))
Full$Variable <- gsub(pattern = "SPX1/SPX2",
                      replacement = "spx1/spx2",
                      as.character(Full$Variable))

Dat <- collapse_by_taxonomy(Dat,level = 6,FUN = sum)
Dat <- create_dataset(Dat$Tab, Dat$Map,
                      data.frame(ID = row.names(Dat$Tab),
                                 Taxonomy = row.names(Dat$Tab),
                                 row.names = row.names(Dat$Tab)))

p1 <- plot_fragen_heatmap(Full = Full, Dat = Dat, group = c("GC2","GC1"),
                          variable = "Genotype", sig_thres = 0.05,
                          pool.FUN = mean, phyl.level = 3,
                          mut_order = levels(Dat$Map$Genotype))
p1
dat <- p1$data
dat$Phylum[ !(dat$Phylum %in% top_phyla) ] <- "Low Abundance"
names(phyla_names) <- top_phyla
dat$Phylum <- phyla_names[ dat$Phylum ]
dat$Phylum <- factor(dat$Phylum,levels = phyla_names)
levels(dat$Enrichment) <- c("Decrease","Increase","none")
dat$Abundance <- (2^dat$Abundance) - 1

p1 <- ggplot(dat,aes(x=Genotype,y=OTU)) +
  #   geom_tile(aes(fill=Abundance,col=Enrichment),size=2) +
  geom_tile(aes(fill=Abundance)) +
  geom_tile(aes(col=Enrichment),size=0.5,alpha=1,fill = NA) +
  facet_grid(Phylum ~ Type,scales="free",space="free") +
  scale_fill_continuous(low="white",high="black",
                        trans = "sqrt",na.value = "white") +
  scale_colour_manual(values = c("blue","red",NA)) +
  ylab("Family") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(face="bold",colour="black",
                                   angle = 90, size = 10),
        panel.background = element_blank(),
        strip.text.y = element_text(color = NA))
p1

dummy <- ggplot(dat, aes(x = Genotype, y = OTU)) +
  facet_grid(Phylum ~ Type, scales = "free",space = "free") + 
  geom_rect(aes(fill=Phylum), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  scale_fill_manual(values = phyla_colors) +
  theme_minimal() +
  theme(strip.text.y = element_text(color = NA),
        axis.text.y = element_blank())
dummy

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(dummy)
panels <- grepl(pattern="panel", g2$layout$name) & g2$layout$l == max(g2$layout$l[ g2$layout$name == "panel"])
strips <- grepl(pattern="strip-right", g2$layout$name)
g2$layout$r[panels] <- g2$layout$r[panels] + 1
g2$layout$l[panels] <- g2$layout$l[panels] + 1
new_strips <- gtable_select(g2, panels | strips)

# grid.newpage()
# grid.draw(new_strips)

new_plot <- gtable_stack(g1, new_strips)
svg("fam_fracgen_heatmap.svg",width = 7, height = 10)
grid.newpage()
grid.draw(new_plot)
dev.off()

new_plot <- gtable_stack(g1, new_strips)
png("fam_fracgen_heatmap.png",width = 1050  , height = 1500)
grid.newpage()
grid.draw(new_plot)
dev.off()

tab <- dcast(data = dat, formula = OTU ~ Genotype, value.var = "Enrichment")
filename <- "fam_gen_heatmap_phosphate.txt"
write.table(tab, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)








# Doughnut plot
dat <- dat
levels(dat$Enrichment) <- c("-1","1","0")
p1 <- doughnut_plot(dat = dat,top_phyla = phyla_names,
                    mut_order = levels(Dat$Map$Genotype))
#p1
filename <- paste("fam_doughnut_phosphate.png",sep = "")
ggsave(filename,p1,width = 10, height = 6)
filename <- paste("fam_doughnut_phosphate.svg",sep = "")
ggsave(filename,p1,width = 10, height = 6)

# Monte carlo testing
set.seed(7151)
pval <- monte_carlo_test(dat = dat,nperm = 5000)
filename <- paste("fam_overlap_pval_phosphate.txt",sep = "")
write.table(pval,file=filename,sep="\t",
            quote=FALSE,col.names = NA)

p1 <- ggplot(melt(pval,value.name = "p.value"),aes(x = Var1, y = Var2)) + 
  geom_tile(aes(fill = p.value < 0.05)) +
  geom_text(aes(label = round(p.value,3)), size = 2) +
  theme_blackbox +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90))
#p1
filename <- paste("fam_overlap_pval_phosphate.png",sep = "")
ggsave(filename,p1, width = 5,height = 4)
filename <- paste("fam_overlap_pval_phosphate.svg",sep = "")
ggsave(filename,p1, width = 5,height = 4)

p1 <- ggplot(subset(Full[grep("Genotype",Full$Variable),],q.value < 0.05 & Variable != "GenotypeSoil"),
             aes(x = Variable, y = Estimate,color = Estimate > 0, fill = Estimate > 0)) +
  geom_boxplot(fill = NA, size = 2,
               position = position_dodge(width = 0.9),outlier.colour = NA) +
  geom_point(size = 4, shape = 21, col = "black",
             position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1)) +
  theme_blackbox +
  theme(axis.text.x = element_text(angle = 90) )
#p1
filename <- paste("fam_significant_gen_effects_phosphate.png",sep = "")
ggsave(filename,p1, width = 5,height = 4)
filename <- paste("fam_significant_gen_effects_phosphate.svg",sep = "")
ggsave(filename,p1, width = 5,height = 4)
