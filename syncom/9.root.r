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

library(ggplot2)
source("~/rhizogenomics/src/trunk/phosphate_code/functions.r")

Tab <- read.table("~/rhizogenomics/data/phosphate/synthetic_physiology/main.and.lateral.roots.txt",
                  sep = "\t", header = TRUE, stringsAsFactors = FALSE)
head(Tab)
Tab$Bacteria[ Tab$Bacteria == "+Bacteria" ] <- "SynCom"
Tab$Pi <- factor(Tab$Pi, levels = c("625uM","50uM"))

p1 <- ggplot(Tab, aes(x = Pi, color = Bacteria, y = Main.root)) +
  geom_boxplot(outlier.colour = NA,size = 1.5,
               position = position_dodge(width = 0.9)) +
  geom_jitter(aes(fill = Bacteria), shape = 21, color = "black",size = 2.5,
              position = position_jitterdodge(jitter.width = 0.5)) +
  scale_color_manual(values = c("#23809C","#921F3A")) +
  scale_fill_manual(values = c("#23809C","#921F3A")) +
  ylab("Main root length (cm)") +
  guides(fill = guide_legend(title = "Bacteria***"),
         color= guide_legend(title = "Bacteria***")) +
  AMOR::theme_blackbox
p1
ggsave("main.root.png",p1, width = 5, height = 4)
ggsave("main.root.svg",p1, width = 5, height = 4)

Tab$PiBac <- interaction(Tab$Pi, Tab$Bacteria)
m1 <- aov(Main.root ~ Bacteria + Rep, data = Tab)
summary(m1)
m1.tuk <- TukeyHSD(m1)
m1.tuk$PiBac
tuk <- agricolae::HSD.test(m1,trt = "Bac")
tuk

p1 <- ggplot(Tab, aes(x = Pi, color = Bacteria, y = Lateral.root)) +
  geom_boxplot(outlier.colour = NA,size = 1.5,
               position = position_dodge(width = 0.9)) +
  geom_jitter(aes(fill = Bacteria), shape = 21, color = "black",size = 2.5,
              position = position_jitterdodge(jitter.width = 0.5)) +
  scale_color_manual(values = c("#23809C","#921F3A")) +
  scale_fill_manual(values = c("#23809C","#921F3A")) +
  ylab("[# Lateral root] / [Main root (cm)]") +
  guides(fill = guide_legend(title = "Bacteria***"),
         color= guide_legend(title = "Bacteria***")) +
  
  AMOR::theme_blackbox
p1
ggsave("lateral.root.png",p1,width = 5, height = 4)
ggsave("lateral.root.svg",p1,width = 5, height = 4)

m1 <- aov(Lateral.root ~ Bacteria + Rep, data = Tab)
summary(m1)
m1.tuk <- TukeyHSD(m1)
tuk <- agricolae::HSD.test(m1,trt = "Bac")
tuk

Tab <- read.table("~/rhizogenomics/data/phosphate/synthetic_physiology/fertilization.txt",
                  sep = "\t", header = TRUE)
head(Tab)
Tab$Pi.level <- factor(Tab$Pi.level,levels = c("625uM","50uM","30uM","10uM","0uM"))
Tab$Bacteria <- factor(Tab$Bacteria, levels = c("No Bacteria", "HK10E5","HK10E6","HK10E7","SynCom"))

pal <- colorRampPalette(colors = c("#00441b","#99d8c9"))


p1 <- ggplot(subset(Tab,Figure == "A"), aes(x = Bacteria, col = Pi.level, y = Pi)) +
  geom_boxplot(outlier.colour = NA,size = 1.5,
               position = position_dodge(width = 0.9)) +
  geom_jitter(aes(fill = Pi.level), shape = 21, color = "black",size = 2.5,
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9)) +
  scale_color_manual(values = pal(5)) +
  scale_fill_manual(values = pal(5)) +
  scale_y_log10(breaks = c(1,2,5,10,12)) +
  ylab("Cellular Pi (umol/FWmg*1000)") +
  AMOR::theme_blackbox
p1
ggsave("pi.dose.syncom.png",p1,width = 5, height = 4)
ggsave("pi.dose.syncom.svg",p1,width = 5, height = 4)

Tab$Pi.level.num <- Tab$Pi.level
levels(Tab$Pi.level.num) <- c(625,50,30,10,0)
Tab$Pi.level.num <- as.numeric(as.character(Tab$Pi.level.num))
m1 <- aov(log10(Pi) ~ Pi.level.num * Bacteria, data = subset(Tab,Figure == "A"))
summary(m1)


p1 <- ggplot(subset(Tab,Figure == "B"), aes(x = Bacteria,color = Bacteria, y = Pi)) +
  geom_boxplot(outlier.colour = NA,size = 1.5,
               position = position_dodge(width = 0.9)) +
  geom_jitter(aes(fill = Bacteria), shape = 21, color = "black",size = 2.5,
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9)) +
  scale_color_manual(values = c("#23809C","#dd3497","#ae017e","#7a0177","#921F3A") ) +
  scale_fill_manual(values = c("#23809C","#dd3497","#ae017e","#7a0177","#921F3A")) +
  scale_y_log10(breaks = c(1,2,5,10,12)) +
  ylab("Cellular Pi (umol/FWmg*1000)") +
  AMOR::theme_blackbox
p1
ggsave("pi.fertilization.syncom.png",p1,width = 5, height = 4)
ggsave("pi.fertilization.syncom.svg",p1,width = 5, height = 4)
m1 <- aov(log10(Pi) ~ Bacteria, data = subset(Tab,Figure == "B"))
summary(m1)
m1.tuk <- agricolae::HSD.test(m1,trt = "Bacteria")
m1.tuk
