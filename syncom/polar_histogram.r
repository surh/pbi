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
library(reshape2)
library(phenotypicForest)
library(GGally)
#library(lme4)
library(lmerTest)



randomWord<-function(n, nLetters = 5)
  replicate(n,paste(sample(letters, nLetters, replace = TRUE),sep = '', collapse=''))

set.seed(42)
nFamily <- 20
nItemPerFamily <- sample(1:6, nFamily,replace = TRUE)
nValues <- 3

df <- data.frame(
  family = rep( randomWord(nFamily), times = nValues * nItemPerFamily),
  item   = rep( randomWord(sum(nItemPerFamily), 3), each = nValues ),
  score  = rep( paste0("V",1:nValues), times = sum(nItemPerFamily)),
  value  = round(500 * runif( sum(nItemPerFamily * nValues)),2))

# item is sample
# family is group
# score is taxa
# value is count
polarHistogram(df, circleProportion = 0.9
               familyLabels = TRUE,direction = "outwards")


polarHistogram(df, circleProportion = 0.9, normalised = FALSE,
               familyLabels = TRUE,direction = "outwards")

df.neg <- df
df.neg$value <- df.neg$value * sample(c(1,-1),size = nrow(df.neg), replace = TRUE)

polarHistogram(df.neg, circleProportion = 0.9, normalised = FALSE,
               familyLabels = TRUE,direction = "outwards")



########## TEST ELONGATION #############

Dat <- read.table("ellongation_full.txt",sep="\t", header = TRUE)
Dat$bacteria <- factor(Dat$bacteria, levels = c("none","G1G2","G2G3","G3N1","N1N2","N2N3","N3B1","B1B2","B2B3","B3G1"))

# Plot experimental reproducibility
p1 <- ggplot(Dat,aes(x = bacteria, y = Elongation, col = Experiment)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x  = element_text(size = 16, angle = 90))
p1
ggsave("elongation_experiment.png",p1)

# Plot overal distribution
p1 <- ggplot(Dat,aes(x = Elongation,fill = bacteria)) +
  geom_density(alpha = 0.3) +
  theme_classic()
p1
ggsave("elongation_distributions.png",p1)

# Fit Mixed model
m1 <- lmer(Elongation ~ start + end + bacteria +
             start*end + start*bacteria + end*bacteria +
             Plate + Experiment +
             (1|Plate) + (1|Experiment),
           data = Dat)

m2 <- lmer(Elongation ~ start + end + bacteria +
             start*end + start*bacteria + end*bacteria +
             Plate + Experiment +
             (1|Experiment),
           data = Dat)

m3 <- lmer(Elongation ~ start + end + bacteria +
             start*end + start*bacteria + end*bacteria +
             Plate + Experiment +
             (1|Plate),
           data = Dat)

anova(m1,m2)
anova(m1,m3)
AIC(m1,m2,m3)
# Take m3, Experiment does not explain as much


m4 <- lmer(Elongation ~ start + end + bacteria +
             start*end + start*bacteria + end*bacteria +
             Plate +
             (1|Plate),
           data = Dat)

m5 <- lmer(Elongation ~ start + end + bacteria +
             start*end + start*bacteria + end*bacteria +
             Experiment +
             (1|Plate),
           data = Dat)

# Keep m5
AIC(m1,m2,m3,m4,m5)
anova(m3,m5)

m6 <- lmer(Elongation ~ start + end + bacteria +
             start*end + start*bacteria + end*bacteria +
             (1|Plate),
           data = Dat)

m7 <- lmer(Elongation ~ start + end + bacteria +
             start*bacteria + end*bacteria +
             Experiment +
             (1|Plate),
           data = Dat)

m8 <- lmer(Elongation ~ start + end + bacteria +
             start*end + end*bacteria +
             Experiment +
             (1|Plate),
           data = Dat)


m9 <- lmer(Elongation ~ start + end + bacteria +
             start*end + start*bacteria +
             Experiment +
             (1|Plate),
           data = Dat)

# take m7
AIC(m1,m2,m3,m4,m5,m6,m7,m8,m9)
anova(m5,m7)

m10 <- lmer(Elongation ~ start + end + bacteria +
             start*bacteria + end*bacteria +
             (1|Plate),
           data = Dat)
m11 <- lmer(Elongation ~ start + end + bacteria +
             end*bacteria +
             Experiment +
             (1|Plate),
           data = Dat)
m12 <- lmer(Elongation ~ start + end + bacteria +
             start*bacteria +
             Experiment +
             (1|Plate),
           data = Dat)
# keep m7
AIC(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12)

m13 <- lmer(Elongation ~ start + end + bacteria +
             start*bacteria + end*bacteria +
             (1|Experiment) +
             (1|Plate),
           data = Dat)

summary(m7)
summary(m13)
AIC(m7,m13)


######################



Dat <- read.table("~/rhizogenomics/data/synthetic/wheel_phosphate/2015-07-03.master.txt",sep="\t",header = TRUE, stringsAsFactors = FALSE)

Dat$bacteria[ Dat$bacteria == "syncom1" ] <- "G1G2"
Dat$bacteria[ Dat$bacteria == "syncom2" ] <- "G2G3"
Dat$bacteria[ Dat$bacteria == "syncom3" ] <- "G3N1"
Dat$bacteria[ Dat$bacteria == "syncom4" ] <- "N1N2"
Dat$bacteria[ Dat$bacteria == "syncom5" ] <- "N2N3"
Dat$bacteria[ Dat$bacteria == "syncom6" ] <- "N3B1"
Dat$bacteria[ Dat$bacteria == "syncom7" ] <- "B1B2"
Dat$bacteria[ Dat$bacteria == "syncom8" ] <- "B2B3"
Dat$bacteria[ Dat$bacteria == "syncom9" ] <- "B3G1"
Dat$bacteria <- factor(Dat$bacteria, levels = c("none","G1G2","G2G3","G3N1","N1N2","N2N3","N3B1","B1B2","B2B3","B3G1"))

Dat$start[ Dat$start == "highP" ] <- "fullP"
Dat$condition <- paste(Dat$start,Dat$end,sep=".")

Dat$item <- 1:nrow(Dat)
Dat$family <- Dat$bacteria

# Dat$score <- "Pi"
# Dat$value <- Dat$umol.mgFW.1000
guides <- c(3, 5, 7, 9, 11, 13, 15)

Dat$score <- "size"
Dat$value <- Dat$Healthy
Dat <- Dat[ !is.na(Dat$value),]
guides <- c(0.5, 1, 1.5)

polarHistogram(df = Dat, guides = guides, direction = "outwards",spaceFamily = 3,
               circleProportion = 0.9, familyLabels = TRUE,normalised = FALSE) + ggtitle(label = "All")

polarHistogram(df = droplevels(subset(Dat, start == "lowP")), guides = guides,
               direction = "outwards",spaceFamily = 3,
               circleProportion = 0.9, familyLabels = TRUE,normalised = FALSE) + ggtitle(label = "lowP")

polarHistogram(df = droplevels(subset(Dat, start == "fullP")), guides = guides,
               direction = "outwards",spaceFamily = 3,
               circleProportion = 0.9, familyLabels = TRUE,normalised = FALSE) + ggtitle(label = "fullP")

polarHistogram(df = droplevels(subset(Dat, end == "30pi")), guides = guides,
               direction = "outwards",spaceFamily = 3,
               circleProportion = 0.9, familyLabels = TRUE,normalised = FALSE) + ggtitle(label = "30pi")

polarHistogram(df = droplevels(subset(Dat, end == "100pi")), guides = guides,
               direction = "outwards",spaceFamily = 3,
               circleProportion = 0.9, familyLabels = TRUE,normalised = FALSE) + ggtitle(label = "100pi")


boxplot(value ~ start+end+bacteria,data = Dat,las = 2,subset = start == "lowP" & end == "30pi")
boxplot(value ~ start+end+bacteria,data = Dat,las = 2,subset = start == "lowP" & end == "100pi")
boxplot(value ~ start+end+bacteria,data = Dat,las = 2,subset = start == "fullP" & end == "30pi")
boxplot(value ~ start+end+bacteria,data = Dat,las = 2,subset = start == "fullP" & end == "100pi")


ggpairs(data = data.frame(Pi = log2(Dat$umol.mgFW.1000), size = log2(Dat$Healthy), condition = Dat$condition),
        columns = 1:2, color = "condition")

