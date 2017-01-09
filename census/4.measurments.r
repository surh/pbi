library(AMOR)
#source("~/rhizogenomics/src/trunk/immune_code/PRR_Propep/functions.r")
source("~/rhizogenomics/src/trunk/phosphate_code//functions.r")

datadir <- "~/rhizogenomics/data/phosphate/census/data/"
outdir <- "measurments_out/"

dir.create(outdir)

############################
# Load data
filename <- paste(datadir,"/map_complete.txt",sep = "")
Map <- read.table(filename,sep="\t",quote = "", header = TRUE)

load(file = "data/Dat.filter.processed.rdat")

# Map <- read.table("~/rhizogenomics/data/phosphate/census/data/map_complete.txt",sep = "\t",quote = "", header = TRUE)
# load("~/rhizogenomics/experiments/2016/2016-03-30.phosphateMF/data/Dat.filter.processed.rdat")

# Homogenize samples
Map <- Map[row.names(Dat$Map),]
Map$Fraction <- NULL
Map$Genotype <- NULL
Map$Soil <- NULL
Map$Experiment <- NULL
Map$lowgermination <- NULL
Map$replant <- NULL
Dat$Map <- cbind(Dat$Map,Map)
levels(Dat$Map$Genotype)[1] <- "Col-0"
rm(Map)

filename <- paste(outdir,"/plant_measurments_cor.png",sep = "")
png(filename,width = 1000,height = 1000)
pairs(log(Dat$Map[,9:17]), upper.panel = panel.cor,lower.panel = panel.lm)
dev.off()

filename <- paste(outdir,"/plant_measurments_cor_wt.png",sep = "")
png(filename,width = 1000,height = 1000)
pairs(log(subset(Dat$Map, Genotype == "Col-0")[,9:17]),upper.panel = panel.cor,lower.panel = panel.lm)
dev.off()

filename <- paste(outdir,"/plant_measurments_cor_scale.png",sep = "")
png(filename,width = 1000,height = 1000)
pairs(scale(Dat$Map[,9:17]), upper.panel = panel.cor,lower.panel = panel.lm)
dev.off()

filename <- paste(outdir,"/plant_measurments_cor_scalewt.png",sep = "")
png(filename,width = 1000,height = 1000)
pairs(scale(subset(Dat$Map, Genotype == "Col-0")[,9:17]), upper.panel = panel.cor,lower.panel = panel.lm)
dev.off()

###### Plot features
Dat <- subset(Dat,Genotype != "Soil")

p1 <- plotgg_var(Dat,var.name = "Pi",x = "Genotype", col = "Genotype")
dat <- p1$data
p1 <- ggplot(dat, aes(x = Genotype, y = Pi, col = Genotype,
                      fill = Genotype)) +
  geom_boxplot(fill = NA, outlier.colour = NA, 
               position = position_dodge(width = 0.9), size = 2) +
  geom_point(aes(shape = Experiment), col = "black",
             position = position_jitterdodge(dodge.width = 0.0,
                                             jitter.width = 8),
             size = 4) +
  scale_shape_manual(values = c(21,24)) +
  scale_color_manual(values = gen_colors) +
  scale_fill_manual(values = gen_colors) +
  scale_y_log10() +
  theme_blackbox +
  theme(axis.text.x = element_text(face = "bold", angle = 90, size = 12))
p1
  
filename <- paste(outdir,"/pi.png",sep = "")
ggsave(filename,p1,width = 6,height = 6)
filename <- paste(outdir,"/pi.svg",sep = "")
# ggsave(filename,p1,width = 7,height = 6)
grDevices::svg(filename = filename,width = 6,height = 6)
print(p1)
dev.off()

m1 <- aov(log(Pi) ~ Genotype + Experiment, data = Dat$Map)
summary(m1)
m1.tuk <- TukeyHSD(m1)
m1.tuk$Genotype
gr <- multcomp::cld(multcomp::glht(m1,linfct = multcomp::mcp(Genotype = "Tukey")))
gr

p1 <- plotgg_var(Dat,var.name = "Anthoc",x = "Genotype", col = "Genotype")
p1 <- p1 + scale_color_manual(values = gen_colors) +
  scale_fill_manual(values = gen_colors) +
  scale_y_log10() +
  theme(axis.text.x = element_text(face = "bold", angle = 90, size = 12))
#p1
filename <- paste(outdir,"/anthoc.png",sep = "")
ggsave(filename,p1,width = 6,height = 6)
filename <- paste(outdir,"/anthoc.svg",sep = "")
ggsave(filename,p1,width = 6,height = 6)

m1 <- aov(log(Anthoc) ~ Genotype + Experiment, data = Dat$Map)
summary(m1)
m1.tuk <- TukeyHSD(m1)
m1.tuk$Genotype
gr <- multcomp::cld(multcomp::glht(m1,linfct = multcomp::mcp(Genotype = "Tukey")))
gr

p1 <- plotgg_var(Dat,var.name = "ACP5",x = "Genotype", col = "Genotype")
p1 <- p1 + scale_color_manual(values = gen_colors) +
  scale_fill_manual(values = gen_colors) +
  scale_y_log10() +
  theme(axis.text.x = element_text(face = "bold", angle = 90, size = 12))
#p1
filename <- paste(outdir,"/acp5.png",sep = "")
ggsave(filename,p1,width = 6,height = 6)
filename <- paste(outdir,"/acp5.svg",sep = "")
ggsave(filename,p1,width = 6,height = 6)

groups <- NULL
m1 <- aov(log(ACP5) ~ Genotype + Experiment, data = Dat$Map)
summary(m1)
m1.tuk <- TukeyHSD(m1)
m1.tuk$Genotype
gr <- multcomp::cld(multcomp::glht(m1,linfct = multcomp::mcp(Genotype = "Tukey")))
gr <- data.frame(Gene = "ACP5",Genotype = names(gr$mcletters$Letters),
                 groups = gr$mcletters$Letters)
groups <- rbind(groups,gr)

p1 <- plotgg_var(Dat,var.name = "At4",x = "Genotype", col = "Genotype")
p1 <- p1 + scale_color_manual(values = gen_colors) +
  scale_fill_manual(values = gen_colors) +
  scale_y_log10() +
  theme(axis.text.x = element_text(face = "bold", angle = 90, size = 12))
#p1
filename <- paste(outdir,"/at4.png",sep = "")
ggsave(filename,p1,width = 6,height = 6)
filename <- paste(outdir,"/at4.svg",sep = "")
ggsave(filename,p1,width = 6,height = 6)

m1 <- aov(log(At4) ~ Genotype + Experiment, data = Dat$Map)
summary(m1)
m1.tuk <- TukeyHSD(m1)
m1.tuk$Genotype
gr <- multcomp::cld(multcomp::glht(m1,linfct = multcomp::mcp(Genotype = "Tukey")))
temp <- gr$mcletters$Letters
temp <- gsub(pattern = "b",replacement = "A",x = temp)
temp <- gsub(pattern = "a",replacement = "B",x = temp)
temp <- tolower(temp)
temp
gr$mcletters$Letters
gr$mcletters$Letters <- temp
gr <- data.frame(Gene = "At4",Genotype = names(gr$mcletters$Letters),
                 groups = gr$mcletters$Letters)
groups <- rbind(groups,gr)

p1 <- plotgg_var(Dat,var.name = "PHF1",x = "Genotype", col = "Genotype")
p1 <- p1 + scale_color_manual(values = gen_colors) +
  scale_fill_manual(values = gen_colors) +
  scale_y_log10() +
  theme(axis.text.x = element_text(face = "bold", angle = 90, size = 12))
#p1
filename <- paste(outdir,"/phf1.png",sep = "")
ggsave(filename,p1,width = 6,height = 6)
filename <- paste(outdir,"/phf1.svg",sep = "")
ggsave(filename,p1,width = 6,height = 6)

m1 <- aov(log(PHF1) ~ Genotype + Experiment, data = Dat$Map)
summary(m1)
m1.tuk <- TukeyHSD(m1)
m1.tuk$Genotype
gr <- multcomp::cld(multcomp::glht(m1,linfct = multcomp::mcp(Genotype = "Tukey")))
gr <- data.frame(Gene = "PHF1",Genotype = names(gr$mcletters$Letters),
                 groups = gr$mcletters$Letters)
groups <- rbind(groups,gr)


p1 <- plotgg_var(Dat,var.name = "PHT1.1",x = "Genotype", col = "Genotype")
p1 <- p1 + scale_color_manual(values = gen_colors) +
  scale_fill_manual(values = gen_colors) +
  scale_y_log10() +
  theme(axis.text.x = element_text(face = "bold", angle = 90, size = 12))
#p1
filename <- paste(outdir,"/pht1.1.png",sep = "")
ggsave(filename,p1,width = 6,height = 6)
filename <- paste(outdir,"/pht1.1.svg",sep = "")
ggsave(filename,p1,width = 6,height = 6)

m1 <- aov(log(PHT1.1) ~ Genotype + Experiment, data = Dat$Map)
summary(m1)
m1.tuk <- TukeyHSD(m1)
m1.tuk$Genotype
gr <- multcomp::cld(multcomp::glht(m1,linfct = multcomp::mcp(Genotype = "Tukey")))
temp <- gr$mcletters$Letters
temp <- gsub(pattern = "b",replacement = "A",x = temp)
temp <- gsub(pattern = "a",replacement = "B",x = temp)
temp <- tolower(temp)
temp
gr$mcletters$Letters
gr$mcletters$Letters <- temp
gr <- data.frame(Gene = "PHT1;1",Genotype = names(gr$mcletters$Letters),
                 groups = gr$mcletters$Letters)
groups <- rbind(groups,gr)

p1 <- plotgg_var(Dat,var.name = "RNS1",x = "Genotype", col = "Genotype")
p1 <- p1 + scale_color_manual(values = gen_colors) +
  scale_fill_manual(values = gen_colors) +
  scale_y_log10() +
  theme(axis.text.x = element_text(face = "bold", angle = 90, size = 12))
#p1
filename <- paste(outdir,"/rns1.png",sep = "")
ggsave(filename,p1,width = 6,height = 6)
filename <- paste(outdir,"/rns1.svg",sep = "")
ggsave(filename,p1,width = 6,height = 6)

m1 <- aov(log(RNS1) ~ Genotype + Experiment, data = Dat$Map)
summary(m1)
m1.tuk <- TukeyHSD(m1)
m1.tuk$Genotype
gr <- multcomp::cld(multcomp::glht(m1,linfct = multcomp::mcp(Genotype = "Tukey")))
gr <- data.frame(Gene = "RNS1",Genotype = names(gr$mcletters$Letters),
                 groups = gr$mcletters$Letters)
groups <- rbind(groups,gr)

p1 <- plotgg_var(Dat,var.name = "SOD1",x = "Genotype", col = "Genotype")
p1 <- p1 + scale_color_manual(values = gen_colors) +
  scale_fill_manual(values = gen_colors) +
  scale_y_log10() +
  theme(axis.text.x = element_text(face = "bold", angle = 90, size = 12))
#p1
filename <- paste(outdir,"/sod1.png",sep = "")
ggsave(filename,p1,width = 6,height = 6)
filename <- paste(outdir,"/sod1.svg",sep = "")
ggsave(filename,p1,width = 6,height = 6)

m1 <- aov(log(SOD1) ~ Genotype + Experiment, data = Dat$Map)
summary(m1)
m1.tuk <- TukeyHSD(m1)
m1.tuk$Genotype
gr <- multcomp::cld(multcomp::glht(m1,linfct = multcomp::mcp(Genotype = "Tukey")))
gr <- data.frame(Gene = "SQD1",Genotype = names(gr$mcletters$Letters),
                 groups = gr$mcletters$Letters)
groups <- rbind(groups,gr)


p1 <- plotgg_var(Dat,var.name = "SPX1",x = "Genotype", col = "Genotype")
p1 <- p1 + scale_color_manual(values = gen_colors) +
  scale_fill_manual(values = gen_colors) +
  scale_y_log10() +
  theme(axis.text.x = element_text(face = "bold", angle = 90, size = 12))
#p1
filename <- paste(outdir,"/spx1.png",sep = "")
ggsave(filename,p1,width = 6,height = 6)
filename <- paste(outdir,"/spx1.svg",sep = "")
ggsave(filename,p1,width = 6,height = 6)

m1 <- aov(log(SPX1) ~ Genotype + Experiment, data = Dat$Map)
summary(m1)
m1.tuk <- TukeyHSD(m1)
m1.tuk$Genotype
gr <- multcomp::cld(multcomp::glht(m1,linfct = multcomp::mcp(Genotype = "Tukey")))
temp <- gr$mcletters$Letters
temp <- gsub(pattern = "b",replacement = "A",x = temp)
temp <- gsub(pattern = "a",replacement = "B",x = temp)
temp <- tolower(temp)
temp
gr$mcletters$Letters
gr$mcletters$Letters <- temp
gr <- data.frame(Gene = "SPX1",Genotype = names(gr$mcletters$Letters),
                 groups = gr$mcletters$Letters)
groups <- rbind(groups,gr)

rm(p1,m1,m1.tuk,gr,temp)

##### Plot by Genotype
dat <- melt(Dat$Map,id.vars = c("Genotype","Experiment","ID"),
            value.name = "Expression", variable.name = "Gene",
            measure.vars = c("ACP5","At4","PHF1","PHT1.1",
                             "RNS1","SOD1","SPX1"),
            factorsAsStrings = FALSE)
levels(dat$Gene)[ which(levels(dat$Gene) == "PHT1.1") ] <- "PHT1;1"
levels(dat$Gene)[ which(levels(dat$Gene) == "SOD1") ] <- "SQD1" # Naming error

p1 <- ggplot(dat,aes(x = Gene, y = Expression, color = Genotype)) +
  geom_boxplot() +
  scale_color_manual(values = gen_colors) +
  scale_y_log10() +
  theme_blackbox +
  theme(axis.text.x = element_text(angle = 90))
#p1
filename <- paste(outdir,"/expression_by_gene.png",sep = "")
ggsave(filename,p1,width = 6,height = 4)
filename <- paste(outdir,"/expression_by_gene.svg",sep = "")
ggsave(filename,p1,width = 6,height = 4)

p1 <- ggplot(dat,aes(x = Genotype, y = Expression, color = Genotype)) +
  facet_grid(Gene ~ ., scales = "free") + 
  geom_boxplot() +
  scale_color_manual(values = gen_colors) +
  scale_y_log10() +
  theme_blackbox +
  theme(axis.text.x = element_text(angle = 90))
#p1
filename <- paste(outdir,"/expression_by_gene_facet.png",sep = "")
ggsave(filename,p1,width = 4,height = 6)
filename <- paste(outdir,"/expression_by_gene_facet.svg",sep = "")
ggsave(filename,p1,width = 4,height = 6)

p1 <- ggplot(dat,aes(x = Genotype, y = Expression, color = Gene)) +
  geom_boxplot() +
  #scale_y_log10() +
  scale_y_log10(limits = c(0.0004,110)) +
  geom_label(data = groups, aes(x = Genotype, y = 100,
                               color = Gene, label = groups)) +
  theme_blackbox +
  theme(axis.text.x = element_text(angle = 90))
p1

p1 <- ggplot(dat,aes(x = Gene, y = Expression, color = Gene)) +
  facet_grid(~ Genotype,scales = "free_x", switch = "x") +
  geom_boxplot(width = 0.9) +
  geom_text(data = groups, aes(x = Gene, y = 65,
                               color = Gene, label = groups),
            size = 4.5, angle = 90) +
  scale_y_log10(limits = c(0.0004,100)) +
  ylab("Relative Expression") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        panel.background = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_blank(),
        strip.text = element_text(face = "bold", size = 14),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA,size = 2),
        panel.margin = unit(0.5,"lines"))
#p1
#min(dat$Expression,na.rm = TRUE)
filename <- paste(outdir,"/expression_by_genotype.png",sep = "")
ggsave(filename,p1,width = 15,height = 4)
filename <- paste(outdir,"/expression_by_genotype.svg",sep = "")
ggsave(filename,p1,width = 15,height = 4)

p1 <- ggplot(dat,aes(x = Gene, y = Expression, color = Gene)) +
  facet_grid(Genotype ~ ., scales = "free") + 
  geom_boxplot() +
  scale_y_log10() +
  theme_blackbox +
  theme(axis.text.x = element_text(angle = 90))
#p1
filename <- paste(outdir,"/expression_by_genotype_facet.png",sep = "")
ggsave(filename,p1,width = 4,height = 6)
filename <- paste(outdir,"/expression_by_genotype_facet.svg",sep = "")
ggsave(filename,p1,width = 4,height = 6)

# Heatmaps
dat2 <- acast(data = dat,formula = Gene ~ Genotype,value.var = "Expression",
              fun.aggregate = mean,na.rm = TRUE)
png("measurments_out/heatmap_mean_clus.png",width = 1000, height = 1000)
gplots::heatmap.2(log10(dat2), trace = "none",scale = "none",
                  col = colorRampPalette(c("white", "#67001F")),
                  margins = c(12,10),cexRow = 2, cexCol = 2)
dev.off()

dat2 <- acast(data = dat,formula = Gene ~ Genotype,value.var = "Expression",
              fun.aggregate = median,na.rm = TRUE)
png("measurments_out/heatmap_median_clus.png",width = 1000, height = 1000)
gplots::heatmap.2(log10(dat2), trace = "none",scale = "none",
                  col = colorRampPalette(c("white", "#67001F")),
                  margins = c(12,10),cexRow = 2, cexCol = 2)
dev.off()

dat2 <- aggregate(Expression ~ Gene + Genotype,
                  FUN = mean, data = dat,na.rm = TRUE)
p1 <- ggplot(dat2,aes(x = Genotype, y = Gene)) +
  geom_tile(aes(fill=Expression)) +
  scale_fill_gradientn(colours = c("white", "#67001F"), trans = "log10") +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(color = "black", angle = 90),
        axis.text.y = element_text(color = "black"),
        axis.title = element_text(face = "bold"))
#p1
filename <- paste(outdir,"/heatmap_mean.png",sep = "")
ggsave(filename,p1,width = 5,height = 5)
filename <- paste(outdir,"/heatmap_mean.svg",sep = "")
ggsave(filename,p1,width = 5,height = 5)

dat2 <- aggregate(Expression ~ Gene + Genotype,
                  FUN = median, data = dat,na.rm = TRUE)
p1 <- ggplot(dat2,aes(x = Genotype, y = Gene)) +
  geom_tile(aes(fill=Expression)) +
  scale_fill_gradientn(colours = c("white", "#67001F"), trans = "log10") +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(color = "black", angle = 90),
        axis.text.y = element_text(color = "black"),
        axis.title = element_text(face = "bold"))
#p1
filename <- paste(outdir,"/heatmap_median.png",sep = "")
ggsave(filename,p1,width = 5,height = 5)
filename <- paste(outdir,"/heatmap_median.svg",sep = "")
ggsave(filename,p1,width = 5,height = 5)


p1 <- ggplot(dat,aes(x = Genotype, y = Expression)) +
  geom_boxplot() + 
  scale_y_log10()
p1

###### Associate measurment
Dat$Map$Fraction <- factor(Dat$Map$Fraction,levels = c("Soil","EC") )
levels(Dat$Map$Genotype)[1] <- "Col-0"
Dat <- collapse_by_taxonomy(Dat,level = 6,FUN = sum)
Dat <- remove_taxons(Dat = Dat,
                     taxons = row.names(Dat$Tab)[ !findGoodOTUs(OTUtab = Dat$Tab,
                                                                min_reads_otu = 25,
                                                                min_samples_otu = 5,
                                                                method = "absolute") ])
Dat <- clean(Dat,verbose = TRUE)
Dat$Map$LDepth <- log(Dat$Map$Depth)

# Convert to percent
Dat.ra <- create_dataset(Tab = 100*t(t(Dat$Tab)/Dat$Map$Depth),Map = Dat$Map,Tax = Dat$Tax)

# Associate
res <- associate_measurment(Dat = Dat.ra,measurment = "Pi")
res <- res[order(res$pval),]
res$qval <- p.adjust(res$pval,method = "bonferroni")
head(res)

p1 <- plotgg_taxon_vs_var(Dat = Dat.ra,
                          taxon = as.character(res$Taxon[1]),
                          var = "Pi",
                          log.transform = FALSE)
p1

res <- associate_measurment(Dat = Dat.ra,measurment = "Anthoc")
res <- res[order(res$pval),]
res$qval <- p.adjust(res$pval,method = "bonferroni")
head(res)

p1 <- plotgg_taxon_vs_var(Dat = Dat.ra,
                          taxon = as.character(res$Taxon[1]),
                          var = "Anthoc",
                          log.transform = FALSE) +
  ggtitle("Bartonellaceae")
p1
ggsave("measurments_out/bartonellaceae_anthoc.png",p1,width = 6,height = 6)
ggsave("measurments_out/bartonellaceae_anthoc.svg",p1,width = 6,height = 6)

res <- associate_measurment(Dat = Dat.ra,measurment = "ACP5")
res <- res[order(res$pval),]
res$qval <- p.adjust(res$pval,method = "bonferroni")
head(res)

p1 <- plotgg_taxon_vs_var(Dat = Dat.ra,
                          taxon = as.character(res$Taxon[1]),
                          var = "ACP5",
                          log.transform = FALSE)
p1

res <- associate_measurment(Dat = Dat.ra,measurment = "At4")
res <- res[order(res$pval),]
res$qval <- p.adjust(res$pval,method = "bonferroni")
head(res)

p1 <- plotgg_taxon_vs_var(Dat = Dat.ra,
                          taxon = as.character(res$Taxon[1]),
                          var = "At4",
                          log.transform = FALSE)
p1


res <- associate_measurment(Dat = Dat.ra,measurment = "PHF1")
res <- res[order(res$pval),]
res$qval <- p.adjust(res$pval,method = "bonferroni")
head(res)

p1 <- plotgg_taxon_vs_var(Dat = Dat.ra,
                          taxon = as.character(res$Taxon[1]),
                          var = "PHF1",
                          log.transform = FALSE)
p1

res <- associate_measurment(Dat = Dat.ra,measurment = "PHT1.1")
res <- res[order(res$pval),]
res$qval <- p.adjust(res$pval,method = "bonferroni")
head(res)

p1 <- plotgg_taxon_vs_var(Dat = Dat.ra,
                          taxon = as.character(res$Taxon[1]),
                          var = "PHT1.1",
                          log.transform = FALSE)
p1

res <- associate_measurment(Dat = Dat.ra,measurment = "RNS1")
res <- res[order(res$pval),]
res$qval <- p.adjust(res$pval,method = "bonferroni")
head(res)

p1 <- plotgg_taxon_vs_var(Dat = Dat.ra,
                          taxon = as.character(res$Taxon[1]),
                          var = "RNS1",
                          log.transform = FALSE)
p1

res <- associate_measurment(Dat = Dat.ra,measurment = "SOD1")
res <- res[order(res$pval),]
res$qval <- p.adjust(res$pval,method = "bonferroni")
head(res)

p1 <- plotgg_taxon_vs_var(Dat = Dat.ra,
                          taxon = as.character(res$Taxon[1]),
                          var = "SOD1",
                          log.transform = FALSE)
p1

res <- associate_measurment(Dat = Dat.ra,measurment = "SPX1")
res <- res[order(res$pval),]
res$qval <- p.adjust(res$pval,method = "bonferroni")
head(res)

p1 <- plotgg_taxon_vs_var(Dat = Dat.ra,
                          taxon = as.character(res$Taxon[1]),
                          var = "SPX1",
                          log.transform = FALSE)
p1


