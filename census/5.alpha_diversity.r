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
gen_colors <- c("#A6CEE3",
                "#66BD63",
                "#386cb0",
                "#E5C494",
                scales::brewer_pal(type = "qual", pal = "Set3")(7)[4],
                "#DE77AE",
                scales::brewer_pal(type = "qual", pal = "Set3")(7)[6:7],
                "#000000")

filename <- paste(datadir,"/Dat.filter.processed.rar1000.rdat",sep = "")
load(filename)
Dat.rar$Map$Fraction <- factor(Dat.rar$Map$Fraction,levels = c("Soil","EC") )

Dat.rar$Map$Shannon <- vegan::diversity(x = Dat.rar$Tab, index = "shannon", MARGIN = 2 )
Dat.rar$Map$Richness <- colSums(Dat.rar$Tab > 0)

p1 <- plotgg_var(Dat = Dat.rar,var.name = "Shannon", x = "Genotype", col = "Genotype")
p1 <- p1 + scale_color_manual(values = gen_colors, guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = gen_colors) +
  theme(axis.text.x = element_text(angle = 90))
#p1
ggsave("shannon.png",p1, width = 8, height = 6)
ggsave("shannon.svg",p1, width = 8, height = 6)

m1 <- aov(Shannon ~ Fraction + Genotype + Experiment, data = Dat.rar$Map)
summary(m1)
m1.tuk <- TukeyHSD(m1)
m1.tuk$Genotype

p1 <- plotgg_var(Dat = Dat.rar,var.name = "Richness", x = "Genotype", col = "Genotype")
p1 <- p1 + scale_color_manual(values = gen_colors, guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = gen_colors) +
  theme(axis.text.x = element_text(angle = 90))
#p1
ggsave("richness.png",p1, width = 8, height = 6)
ggsave("richness.svg",p1, width = 8, height = 6)

m1 <- aov(Richness ~ Fraction + Genotype + Experiment, data = Dat.rar$Map)
summary(m1)
m1.tuk <- TukeyHSD(m1)
m1.tuk$Genotype

