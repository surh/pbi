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
library(GGally)
library(AMOR)

Dat <- read.table("~/rhizogenomics/data/synthetic/wheel_phosphate/2015-07-27.master.txt",
                  sep = "\t", header = TRUE)

Dat$Condition <- NULL
Dat$Condition[ Dat$StartP == "-Pi,0.5%Suc" & Dat$EndP == "100 uM,0%Suc" ] <- "minusP100"
Dat$Condition[ Dat$StartP == "+Pi,0.5%Suc" & Dat$EndP == "100 uM,0%Suc" ] <- "plusP100"
Dat$Condition[ Dat$StartP == "-Pi,0.5%Suc" & Dat$EndP == "30 uM,0%Suc" ] <- "minusP30"
Dat$Condition[ Dat$StartP == "+Pi,0.5%Suc" & Dat$EndP == "30 uM,0%Suc" ] <- "plusP30"

dat <- data.frame(Elongation = Dat$Elongation,
                  logArea = log(Dat$shoot_area), 
                  logPi = log(Dat$Pi_content),
                  Condition = Dat$Condition)
png(filename = "correlation_phenotypes.png", width = 750, height = 750)
ggpairs(data = dat, columns = 1:3, color = "Condition", size = 3)
dev.off()

Dat.dat <- create_dataset(Tab = t(as.matrix(dat[,1:3],rownames.force = TRUE)),
                          Map = Dat)
Dat.dat <- remove_samples(Dat.dat, samples = colnames(Dat.dat$Tab)[ apply(is.na(Dat.dat$Tab),2,any) ])
Dat.pca <- PCA(Dat = Dat.dat, cor = TRUE)
p1 <- plotgg(Dat.pca, col = "Condition", point_size = 5)
ggsave("pca_3phenotypes.png",p1,width = 4,height = 4)
100*Dat.pca$sdev^2 / Dat.pca$totvar

