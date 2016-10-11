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

########## Process elongation ############
# Read full elongation
Tab <- read.table("~/rhizogenomics/data/synthetic/wheel_phosphate/2015-07-26.elongation.txt",
                  sep = "\t", header = TRUE, stringsAsFactors = FALSE)

names(table(Tab$Picture))[ table(Tab$Picture) == 6 ]
table(Tab$Picture)
sort(table(Tab$Picture))


# Split pictures into plates. Each picture has two plates.
# If there are two different conditions in a picture, each condition
# corresponds to a plate. If there is only one condition in a plate, split
# in half. (It will make some errors, but they should be minimal)
i <- 1
Res <- NULL
for(pic in unique(Tab$Picture)){
  #pic <- unique(Tab$Picture)[2]
  
  # Select rows from a given picture
  dat <- droplevels(subset(Tab, Picture == pic))
  exp <- unique(dat$Experiment)
  
  if(length(unique(dat$Treatment)) == 2){
    # If two conditions each gets a different plate
    dat$Plate <- NULL
    dat$Plate[ dat$Treatment == unique(dat$Treatment)[1] ] <- i
    i <- i+1
    dat$Plate[ dat$Treatment == unique(dat$Treatment)[2] ] <- i
  }else if(length(unique(dat$Treatment)) == 1){
    # If one condition and more than 10 plants, split in two, otherwise assume one plate.
    # split treatment
    if(nrow(dat) > 10){
      num1 <- round(nrow(dat) / 2)
      dat$Plate <- c(rep(i,num1),rep(i <- i+1,nrow(dat)-num1))
    }else{
      dat$Plate <- i
    }
  }else{
    # if different than one or two conditions print warning
    warning(pic)
    stop("ERROR",call. = FALSE) 
  }
  
  Res <- rbind(Res,dat)
  i <- i+1
}

#dim(aggregate(Main.root.elongation..cm. ~ Picture + Treatment,data = Res, FUN = median))

# Aggregate plates
Area <- aggregate(Elongation ~ Picture + Treatment + Experiment + Plate ,data = Res, FUN = median)

# Name conditions
Area.full <- Res
Area.full$StartP <- NA
Area.full$StartP[ grep("+P",Area.full$Treatment) ] <- "+Pi,0.5%Suc"
Area.full$StartP[ grep("-P",Area.full$Treatment) ] <- "-Pi,0.5%Suc"
Area.full$EndP <- NA
Area.full$EndP[ grep("30P",Area.full$Treatment) ] <- "30 uM,0%Suc"
Area.full$EndP[ grep("100P",Area.full$Treatment) ] <- "100 uM,0%Suc"
Area.full$Bacteria <- "none"
Area.full$Bacteria[ grep("+SC1$",Area.full$Treatment) ] <- "G1G2"
Area.full$Bacteria[ grep("+SC2$",Area.full$Treatment) ] <- "G2G3"
Area.full$Bacteria[ grep("+SC3$",Area.full$Treatment) ] <- "G3N1"
Area.full$Bacteria[ grep("+SC4$",Area.full$Treatment) ] <- "N1N2"
Area.full$Bacteria[ grep("+SC5$",Area.full$Treatment) ] <- "N2N3"
Area.full$Bacteria[ grep("+SC6$",Area.full$Treatment) ] <- "N3B1"
Area.full$Bacteria[ grep("+SC7$",Area.full$Treatment) ] <- "B1B2"
Area.full$Bacteria[ grep("+SC8$",Area.full$Treatment) ] <- "B2B3"
Area.full$Bacteria[ grep("+SC9$",Area.full$Treatment) ] <- "G1B3"
Area.full$Bacteria[ grep("+SC10$",Area.full$Treatment) ] <- "G1N1"
Area.full$Bacteria[ grep("+SC11$",Area.full$Treatment) ] <- "G2N1"
Area.full$Bacteria[ grep("+SC12$",Area.full$Treatment) ] <- "G1G3"
Area.full$Bacteria[ grep("+SC13$",Area.full$Treatment) ] <- "G2B3"
Area.full$Bacteria[ grep("+SC14$",Area.full$Treatment) ] <- "G3B3"

Area$StartP <- NA
Area$StartP[ grep("+P",Area$Treatment) ] <- "+Pi,0.5%Suc"
Area$StartP[ grep("-P",Area$Treatment) ] <- "-Pi,0.5%Suc"
Area$EndP <- NA
Area$EndP[ grep("30P",Area$Treatment) ] <- "30 uM,0%Suc"
Area$EndP[ grep("100P",Area$Treatment) ] <- "100 uM,0%Suc"
Area$Bacteria <- "none"
Area$Bacteria[ grep("+SC1$",Area$Treatment) ] <- "G1G2"
Area$Bacteria[ grep("+SC2$",Area$Treatment) ] <- "G2G3"
Area$Bacteria[ grep("+SC3$",Area$Treatment) ] <- "G3N1"
Area$Bacteria[ grep("+SC4$",Area$Treatment) ] <- "N1N2"
Area$Bacteria[ grep("+SC5$",Area$Treatment) ] <- "N2N3"
Area$Bacteria[ grep("+SC6$",Area$Treatment) ] <- "N3B1"
Area$Bacteria[ grep("+SC7$",Area$Treatment) ] <- "B1B2"
Area$Bacteria[ grep("+SC8$",Area$Treatment) ] <- "B2B3"
Area$Bacteria[ grep("+SC9$",Area$Treatment) ] <- "G1B3"
Area$Bacteria[ grep("+SC10$",Area$Treatment) ] <- "G1N1"
Area$Bacteria[ grep("+SC11$",Area$Treatment) ] <- "G2N1"
Area$Bacteria[ grep("+SC12$",Area$Treatment) ] <- "G1G3"
Area$Bacteria[ grep("+SC13$",Area$Treatment) ] <- "G2B3"
Area$Bacteria[ grep("+SC14$",Area$Treatment) ] <- "G3B3"

write.table(Area,"elongation.txt",sep = "\t", quote = FALSE, row.names = FALSE)
write.table(Area.full,"elongation_full.txt",sep = "\t", quote = FALSE, row.names = FALSE)

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

