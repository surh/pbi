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
datadir <- "data/"
# gen_colors <- c("#A6CEE3",
#                 "#66BD63",
#                 scales::brewer_pal(type = "qual", pal = "Set3")(7)[2],
#                 "#E5C494",
#                 scales::brewer_pal(type = "qual", pal = "Set3")(7)[4],
#                 "#DE77AE",
#                 scales::brewer_pal(type = "qual", pal = "Set3")(7)[6:7],
#                 "#000000")


gen_colors <- c("#A6CEE3",
                "#66BD63",
                "#386cb0",
                "#E5C494",
                scales::brewer_pal(type = "qual", pal = "Set3")(7)[4],
                "#DE77AE",
                scales::brewer_pal(type = "qual", pal = "Set3")(7)[6:7],
                "#000000")


barplot(rep(1,9),col = gen_colors)

filename <- paste(datadir,"/Dat.filter.processed.rar1000.rdat",sep = "")
load(filename)

##### Bray-Curtis ######
distfun <- function(x,method) vegan::vegdist(x = x, method = "bray")
Dat.rar.pco <- PCO(Dat = Dat.rar,dim = 2,distfun = distfun)

perc_var <- (100*Dat.rar.pco$eig / sum(Dat.rar.pco$eig[ Dat.rar.pco$eig > 0 ]))[1:2]
perc_var <- round(perc_var,digits = 1)

p1 <- plotgg(Dat.rar.pco, shape = "Fraction:Experiment",
             col = "Genotype", point_size = 9)
p1 <- p1 + scale_shape_manual(values = c(15,0 , 19,1,15,0))
#p1 <- p1 + scale_shape_manual(values = c("\u25CF","\u2B58","\u25A0","\u25A1"))
p1 <- p1 + scale_color_manual(values = gen_colors, guide = guide_legend(ncol = 2))
p1 <- p1 + coord_fixed() + xlab(label = paste("PCo1 (",perc_var[1],"%)",sep="")) + 
  ylab(label = paste("PCo1 (",perc_var[2],"%)",sep=""))
p1 <- p1 + theme(axis.text = element_blank(),
                 axis.title = element_text(size = 24)) 
#p1
png("bray_pco.png",width = 1000, height = 1000)
p1
grid.force()
grid::grid.edit("geom_point.points", grep = TRUE, gp = gpar(lwd = 7))
dev.off()

svg("bray_pco.svg",width = 10,height = 10)
p1
grid.force()
grid::grid.edit("geom_point.points", grep = TRUE, gp = gpar(lwd = 7))
dev.off()



# CAP analysis
Dat.rar <- subset(Dat.rar, Fraction != "Soil", drop = TRUE, clean = TRUE)
cap <- vegan::capscale(t(Dat.rar$Tab) ~ Genotype + Condition(Experiment) + Condition(Depth),
                       data = Dat.rar$Map,dfun = distfun)
cap
cap.sum <- summary(cap)
Dat.rar$Map <- cbind(Dat.rar$Map,cap.sum$sites)
percvar <- round(100 * cap$CCA$eig / cap$CCA$tot.chi,2)

p1 <- ggplot(Dat.rar$Map, aes(x = CAP1, y = CAP2, col = Genotype)) +
  geom_point(aes(shape = Experiment), size = 3) +
  stat_ellipse(aes(fill = Genotype), geom = "polygon",
               level = 0.5, alpha = 0.3) +
  scale_color_manual(values = gen_colors) +
  scale_fill_manual(values = gen_colors) +
  xlab(paste("CAP1 (",percvar[1],"%)",sep = "")) +
  ylab(paste("CAP2 (",percvar[2],"%)",sep = "")) +
  theme(axis.text = element_text(color = "black"), 
        axis.title = element_text(face = "bold"),
        panel.background = element_rect(color = "black", size = 3, fill = NA),
        panel.grid = element_blank())
#p1
ggsave("bray_cap1_cap2.svg",p1,width = 5,height = 4)


p1 <- ggplot(Dat.rar$Map, aes(x = CAP3, y = CAP4, col = Genotype)) +
  geom_point(aes(shape = Experiment), size = 3) +
  stat_ellipse(aes(fill = Genotype), geom = "polygon",
               level = 0.5, alpha = 0.3) +
  scale_color_manual(values = gen_colors) +
  scale_fill_manual(values = gen_colors) +
  xlab(paste("CAP3 (",percvar[3],"%)",sep = "")) +
  ylab(paste("CAP4 (",percvar[4],"%)",sep = "")) +
  theme(axis.text = element_text(color = "black"), 
        axis.title = element_text(face = "bold"),
        panel.background = element_rect(color = "black", size = 3, fill = NA),
        panel.grid = element_blank())
#p1
ggsave("bray_cap3_cap4.svg",p1,width = 5,height = 4)

rm(Dat.rar.pco, distfun, p1, perc_var,filename,cap,cap.sum,percvar)

########### Weighted unifrac #########
filename <- paste(datadir,"/Dat.filter.processed.rar1000.rdat",sep = "")
load(filename)
wUnifrac <- read.table("data/weighted_unifrac_otu_table_filter_rar1000.txt")
distfun <- function(mat,method) as.dist(wUnifrac)

Dat.rar.pco <- PCO(Dat = Dat.rar,dim = 2,distfun = distfun)

perc_var <- (100*Dat.rar.pco$eig / sum(Dat.rar.pco$eig[ Dat.rar.pco$eig > 0 ]))[1:2]
perc_var <- round(perc_var,digits = 1)

p1 <- plotgg(Dat.rar.pco, shape = "Fraction:Experiment", col = "Genotype",point_size = 9)
p1 <- p1 + scale_shape_manual(values = c(19,1,17,2))
p1 <- p1 + scale_color_manual(values = gen_colors)
p1 <- p1 + coord_fixed() + xlab(label = paste("PCo1 (",perc_var[1],"%)",sep="")) + 
  ylab(label = paste("PCo1 (",perc_var[2],"%)",sep=""))
p1 <- p1 + theme(axis.text = element_blank(),
                 axis.title = element_text(size = 24)) 
  
png("wUnifrac_pco.png",width = 1000, height = 1000)
p1
grid.force()
grid.edit("geom_point.points", grep = TRUE, gp = gpar(lwd = 7))
dev.off()

# CAP analysis
Dat.rar <- subset(Dat.rar, Fraction != "Soil", drop = TRUE, clean = TRUE)
wUnifrac <- wUnifrac[ colnames(Dat.rar$Tab), colnames(Dat.rar$Tab) ]
cap <- vegan::capscale(t(Dat.rar$Tab) ~ Genotype + Condition(Experiment) + Condition(Depth),
                       data = Dat.rar$Map,dfun = distfun)
cap
cap.sum <- summary(cap)
Dat.rar$Map <- cbind(Dat.rar$Map,cap.sum$sites)
percvar <- round(100 * cap$CCA$eig / cap$CCA$tot.chi,2)

p1 <- ggplot(Dat.rar$Map, aes(x = CAP1, y = CAP2, col = Genotype)) +
  geom_point(aes(shape = Experiment), size = 3) +
  stat_ellipse(aes(fill = Genotype), geom = "polygon",
               level = 0.5, alpha = 0.3) +
  scale_color_manual(values = gen_colors) +
  scale_fill_manual(values = gen_colors) +
  xlab(paste("CAP1 (",percvar[1],"%)",sep = "")) +
  ylab(paste("CAP2 (",percvar[2],"%)",sep = "")) +
  theme(axis.text = element_text(color = "black"), 
        axis.title = element_text(face = "bold"),
        panel.background = element_rect(color = "black", size = 3, fill = NA),
        panel.grid = element_blank())
#p1
ggsave("wUnifrac_cap1_cap2.svg",p1,width = 5,height = 4)

p1 <- ggplot(Dat.rar$Map, aes(x = CAP3, y = CAP4, col = Genotype)) +
  geom_point(aes(shape = Experiment), size = 3) +
  stat_ellipse(aes(fill = Genotype), geom = "polygon",
               level = 0.5, alpha = 0.3) +
  scale_color_manual(values = gen_colors) +
  scale_fill_manual(values = gen_colors) +
  xlab(paste("CAP3 (",percvar[3],"%)",sep = "")) +
  ylab(paste("CAP4 (",percvar[4],"%)",sep = "")) +
  theme(axis.text = element_text(color = "black"), 
        axis.title = element_text(face = "bold"),
        panel.background = element_rect(color = "black", size = 3, fill = NA),
        panel.grid = element_blank())
#p1
ggsave("wUnifrac_cap3_cap4.svg",p1,width = 5,height = 4)

rm(Dat.rar.pco, distfun, p1, perc_var,wUnifrac,filename,cap,cap.sum,percvar)

################### Unweighted UniFrac ###############
filename <- paste(datadir,"/Dat.filter.processed.rar1000.rdat",sep = "")
load(filename)
uUnifrac <- read.table("data/unweighted_unifrac_otu_table_filter_rar1000.txt")
distfun <- function(mat,method) as.dist(uUnifrac)

Dat.rar.pco <- PCO(Dat = Dat.rar,dim = 2,distfun = distfun)

perc_var <- (100*Dat.rar.pco$eig / sum(Dat.rar.pco$eig[ Dat.rar.pco$eig > 0 ]))[1:2]
perc_var <- round(perc_var,digits = 1)

p1 <- plotgg(Dat.rar.pco, shape = "Fraction:Experiment", col = "Genotype",point_size = 9)
p1 <- p1 + scale_shape_manual(values = c(19,1,17,2))
p1 <- p1 + scale_color_manual(values = gen_colors)
p1 <- p1 + coord_fixed() + xlab(label = paste("PCo1 (",perc_var[1],"%)",sep="")) + 
  ylab(label = paste("PCo1 (",perc_var[2],"%)",sep=""))
p1 <- p1 + theme(axis.text = element_blank(),
                 axis.title = element_text(size = 24)) 

png("uUnifrac_pco.png",width = 1000, height = 1000)
p1
grid.force()
grid.edit("geom_point.points", grep = TRUE, gp = gpar(lwd = 7))
dev.off()

# CAP analysis
Dat.rar <- subset(Dat.rar, Fraction != "Soil", drop = TRUE, clean = TRUE)
uUnifrac <- uUnifrac[ colnames(Dat.rar$Tab), colnames(Dat.rar$Tab) ]
cap <- vegan::capscale(t(Dat.rar$Tab) ~ Genotype + Condition(Experiment) + Condition(Depth),
                       data = Dat.rar$Map,dfun = distfun)
cap
cap.sum <- summary(cap)
Dat.rar$Map <- cbind(Dat.rar$Map,cap.sum$sites)
percvar <- round(100 * cap$CCA$eig / cap$CCA$tot.chi,2)

p1 <- ggplot(Dat.rar$Map, aes(x = CAP1, y = CAP2, col = Genotype)) +
  geom_point(aes(shape = Experiment), size = 3) +
  stat_ellipse(aes(fill = Genotype), geom = "polygon",
               level = 0.5, alpha = 0.3) +
  scale_color_manual(values = gen_colors) +
  scale_fill_manual(values = gen_colors) +
  xlab(paste("CAP1 (",percvar[1],"%)",sep = "")) +
  ylab(paste("CAP2 (",percvar[2],"%)",sep = "")) +
  theme(axis.text = element_text(color = "black"), 
        axis.title = element_text(face = "bold"),
        panel.background = element_rect(color = "black", size = 3, fill = NA),
        panel.grid = element_blank())
#p1
ggsave("uUnifrac_cap1_cap2.svg",p1,width = 5,height = 4)

p1 <- ggplot(Dat.rar$Map, aes(x = CAP3, y = CAP4, col = Genotype)) +
  geom_point(aes(shape = Experiment), size = 3) +
  stat_ellipse(aes(fill = Genotype), geom = "polygon",
               level = 0.5, alpha = 0.3) +
  scale_color_manual(values = gen_colors) +
  scale_fill_manual(values = gen_colors) +
  xlab(paste("CAP3 (",percvar[3],"%)",sep = "")) +
  ylab(paste("CAP4 (",percvar[4],"%)",sep = "")) +
  theme(axis.text = element_text(color = "black"), 
        axis.title = element_text(face = "bold"),
        panel.background = element_rect(color = "black", size = 3, fill = NA),
        panel.grid = element_blank())
#p1
ggsave("uUnifrac_cap3_cap4.svg",p1,width = 5,height = 4)

rm(Dat.rar.pco, distfun, p1, perc_var,uUnifrac,filename,cap,cap.sum,percvar)

################ Between site diversity ###############
filename <- paste(datadir,"/Dat.filter.processed.rar1000.rdat",sep = "")
load(filename)

set.seed(87429)
sitediv.gen <- compare_site_diversity(Dat = Dat.rar, factor = "Genotype", nperm = 20)
sitediv.gen

p1 <- plotgg(sitediv.gen)
p1 <- p1 + scale_color_manual(values = gen_colors) +
  scale_fill_manual(values = gen_colors)
#p1
ggsave("sitediv_gen.png",p1,width = 6,height = 4)
ggsave("sitediv_gen.svg",p1,width = 6,height = 4)

p1 <- plotgg(sitediv.gen, confints = FALSE)
p1 <- p1 + scale_color_manual(values = gen_colors)
#p1
ggsave("sitediv_gen_noconfints.png",p1,width = 6,height = 4)
ggsave("sitediv_gen_noconfints.svg",p1,width = 6,height = 4)

# p1 <- plotgg(subset(sitediv.gen,group %in% c("Col","pho2")), confints = TRUE)
# p1 <- p1 + scale_color_manual(values = gen_colors) +
#   scale_fill_manual(values = gen_colors)
# p1

