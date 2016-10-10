library(ggplot2)
library(coefplot)

#library(reshape2)
#library(phenotypicForest)
#library(GGally)
#library(lme4)
#library(lmerTest)

############### FUNCTIONS ################
plotgg_heatmap_syncom_effects <- function(Res, cond1, cond2){
  #cond1 <- 2
  #cond2 <- 2
  
  singlecoms <- paste(rep(c("G","N","B"),each = 3),rep(1:3,times = 3),sep="")
  dat <- subset(Res, StartP == levels(Res$StartP)[cond1] & EndP == levels(Res$EndP)[cond2])
  dat$SynCom1 <- substring(dat$SynCom,1,2)
  dat$SynCom2 <- substring(dat$SynCom,3,4)
  dat$SynCom1 <- factor(dat$SynCom1 , levels = singlecoms)
  dat$SynCom2 <- factor(dat$SynCom2 , levels = rev(singlecoms))
  
  p1 <- ggplot(dat, aes(x = SynCom1, y = SynCom2)) +
    geom_tile(aes(fill = Estimate, col = p.value < 0.05), size = 3, space = 2) +
    scale_fill_gradient2(low = "#d01c8b",mid = "white", high = "#4dac26",midpoint = 0,na.value = "#404040") +  
    scale_color_manual(values = c("#8c510a","#01665e")) +
    ggtitle(paste(levels(dat$StartP)[cond1],"=>",levels(dat$EndP)[cond2])) +
    theme(panel.background = element_blank(),
          axis.title = element_text(face = "bold", size = 12),
          axis.text = element_text(color = "black", size = 10))
  
  p1
}


# Model all single communities
model_all_communities <- function(Dat,cond1,cond2){
  
  singlecoms <- paste(rep(c("G","N","B"),each = 3),rep(1:3,times = 3),sep="")
  X <- matrix(0, nrow = nrow(Dat), ncol = 10 )
  colnames(X) <- c("Inoculated", singlecoms)
  X[,"Inoculated"] <- 1*(Dat$Bacteria != "none")
  X[ grep(pattern = "G1",x = Dat$Bacteria), "G1"] <- 1
  X[ grep(pattern = "G2",x = Dat$Bacteria), "G2"] <- 1
  X[ grep(pattern = "G3",x = Dat$Bacteria), "G3"] <- 1
  X[ grep(pattern = "N1",x = Dat$Bacteria), "N1"] <- 1
  X[ grep(pattern = "N2",x = Dat$Bacteria), "N2"] <- 1
  X[ grep(pattern = "N3",x = Dat$Bacteria), "N3"] <- 1
  X[ grep(pattern = "B1",x = Dat$Bacteria), "B1"] <- 1
  X[ grep(pattern = "B2",x = Dat$Bacteria), "B2"] <- 1
  X[ grep(pattern = "B3",x = Dat$Bacteria), "B3"] <- 1
  
  Dat <- cbind(X,Dat)
  Dat$Picture <- NULL
  Dat$Treatment <- NULL
  Dat$Bacteria <- NULL
  Dat$Plate <- factor(Dat$Plate)
  
  dat <- subset(Dat, StartP == levels(Dat$StartP)[cond1] & EndP == levels(Dat$EndP)[cond2])
  dat$StartP <- NULL
  dat$EndP <- NULL
  dat <- droplevels(dat)
  dat$Inoculated <- NULL
  #dat$Experiment <- NULL
  m1 <- lm(Elongation ~ . , data = dat )
  #summary(m1)
  
  return(m1)
}

########## TEST ELONGATION #############

Dat <- read.table("elongation_full.txt",sep="\t", header = TRUE)
Dat$Bacteria <- factor(Dat$Bacteria, levels = c("none","G1G2","G2G3","G3N1","N1N2","N2N3","N3B1","B1B2","B2B3","G1B3","G1N1","G2N1","G1G3","G2B3","G3B3"))
#Dat$Bacteria

# Plot experimental reproducibility
p1 <- ggplot(Dat,aes(x = Bacteria, y = Elongation, col = Experiment)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x  = element_text(size = 16, angle = 90))
#p1
ggsave("elongation_experiment.png",p1)

# Plot overal distribution
p1 <- ggplot(Dat,aes(x = Elongation,fill = Bacteria)) +
  geom_density(alpha = 0.3) +
  theme_classic()
#p1
ggsave("elongation_distributions.png",p1)

## Fit one at a time

dir <- "elongation_images/"
communities <- levels(Dat$Bacteria)
communities <- communities[ communities != "none" ]
Res <- NULL
for(syncom in communities){
  #syncom <- communities[1]
  
  # Get data
  experiments <- unique(Dat$Experiment[ Dat$Bacteria == syncom ])
  dat <- subset(Dat, Experiment %in% experiments & Bacteria %in% c(syncom,"none"))
  dat$Plate <- factor(dat$Plate)
  
  p1 <- ggplot(dat,aes(x = Elongation,fill = Bacteria)) +
    facet_grid(StartP ~ EndP) +
    geom_density(alpha = 0.3) +
    ggtitle(syncom) +
    theme_classic()
  #p1
  filename <- paste(dir,"/",syncom,"_density.png",sep = "")
  ggsave(filename,p1, width = 4, height = 4)
  
  #m1 <- lm(Elongation ~ StartP + EndP + Bacteria + Experiment + Plate, data = dat)
  #summary(m1)
  
  m1 <- lm(Elongation ~ Bacteria + Experiment + Plate,
           data = subset(dat, StartP == "-Pi,0.5%Suc" & EndP == "100 uM,0%Suc"))
  m1.sum <- summary(m1)
  res <- data.frame(SynCom = syncom,StartP = levels(dat$StartP)[1], EndP = levels(dat$EndP)[1],
                    Estimate = m1.sum$coefficients[2,1], SE = m1.sum$coefficients[2,2],
                    t.value = m1.sum$coefficients[2,3], p.value = m1.sum$coefficients[2,4])
  Res <- rbind(Res,res)
  
  m1 <- lm(Elongation ~ Bacteria + Experiment + Plate,
           data = subset(dat, StartP == "+Pi,0.5%Suc" & EndP == "100 uM,0%Suc"))
  m1.sum <- summary(m1)
  res <- data.frame(SynCom = syncom,StartP = levels(dat$StartP)[2], EndP = levels(dat$EndP)[1],
                    Estimate = m1.sum$coefficients[2,1], SE = m1.sum$coefficients[2,2],
                    t.value = m1.sum$coefficients[2,3], p.value = m1.sum$coefficients[2,4])
  Res <- rbind(Res,res)
  m1 <- lm(Elongation ~ Bacteria + Experiment + Plate,
           data = subset(dat, StartP == "-Pi,0.5%Suc" & EndP == "30 uM,0%Suc"))
  m1.sum <- summary(m1)
  res <- data.frame(SynCom = syncom,StartP = levels(dat$StartP)[1], EndP = levels(dat$EndP)[2],
                    Estimate = m1.sum$coefficients[2,1], SE = m1.sum$coefficients[2,2],
                    t.value = m1.sum$coefficients[2,3], p.value = m1.sum$coefficients[2,4])
  Res <- rbind(Res,res)
  m1 <- lm(Elongation ~ Bacteria + Experiment + Plate,
           data = subset(dat, StartP == "+Pi,0.5%Suc" & EndP == "30 uM,0%Suc"))
  m1.sum <- summary(m1)
  res <- data.frame(SynCom = syncom,StartP = levels(dat$StartP)[2], EndP = levels(dat$EndP)[2],
                    Estimate = m1.sum$coefficients[2,1], SE = m1.sum$coefficients[2,2],
                    t.value = m1.sum$coefficients[2,3], p.value = m1.sum$coefficients[2,4])
  Res <- rbind(Res,res)
}
write.table(Res,"elongation_single_community_test.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Plot results of single community analysis
singlecoms <- paste(rep(c("G","N","B"),each = 3),rep(1:3,times = 3),sep="")
#dat <- subset(Res, StartP == levels(Res$StartP)[cond1] & EndP == levels(Res$EndP)[cond2])
dat <- Res
dat$SynCom1 <- substring(dat$SynCom,1,2)
dat$SynCom2 <- substring(dat$SynCom,3,4)
dat$SynCom1 <- factor(dat$SynCom1 , levels = singlecoms)
dat$SynCom2 <- factor(dat$SynCom2 , levels = rev(singlecoms))

p1 <- ggplot(dat, aes(x = SynCom1, y = SynCom2)) +
  facet_grid(StartP ~ EndP) +
  geom_tile(aes(fill = Estimate, col = p.value < 0.05), size = 3,
            space = 2) +
  scale_fill_gradient2(low = "#d01c8b",mid = "white",
                       high = "#4dac26",midpoint = 0,
                       na.value = "#404040",
                       guide = guide_colorbar(title = "Estimate (cm)")) +
  scale_color_manual(values = c("#8c510a","#01665e")) +
  ggtitle("Elongation") +
  theme_classic()
p1
ggsave("elongation_images/heatmap_mix_estimate.png",p1, width = 8, height = 7)
rm(singlecoms,dat)

# p1 <- plotgg_heatmap_syncom_effects(Res = Res, cond1 = 1, cond2 = 1)
# p1
# ggsave("elongation_images/heatmap_mix_estimate_minusP100.png",p1, width = 6, height = 5)
# 
# p1 <- plotgg_heatmap_syncom_effects(Res = Res, cond1 = 2, cond2 = 1)
# p1
# ggsave("elongation_images/heatmap_mix_estimate_plusP100.png",p1, width = 6, height = 5)
# 
# p1 <- plotgg_heatmap_syncom_effects(Res = Res, cond1 = 1, cond2 = 2)
# p1
# ggsave("elongation_images/heatmap_mix_estimate_minusP30.png",p1, width = 6, height = 5)
# 
# p1 <- plotgg_heatmap_syncom_effects(Res = Res, cond1 = 2, cond2 = 2)
# p1
# ggsave("elongation_images/heatmap_mix_estimate_plusP30.png",p1, width = 6, height = 5)

## Model all commubities together
singlecoms <- paste(rep(c("G","N","B"),each = 3),rep(1:3,times = 3),sep="")

Res2 <- NULL
for(i in 1:2){
  for(j in 1:2){
    #i <- 1
    #j <- 1
    m1 <- model_all_communities(Dat = Dat, cond1 = i, cond2 = j)
    m1.sum <- summary(m1)
    
    res <- data.frame(SynCom = singlecoms, StartP = levels(Dat$StartP)[i],
               EndP = levels(Dat$EndP)[j],
               Estimate = m1.sum$coefficients[ singlecoms, 1 ],
               SE = m1.sum$coefficients[ singlecoms, 2 ],
               t.value = m1.sum$coefficients[ singlecoms, 3 ],
               p.value = m1.sum$coefficients[ singlecoms, 4 ])
    Res2 <- rbind(Res2,res)
    
    
    p1 <- coefplot(m1, intercept = FALSE,innerCI = 1, outerCI = 2, coefficients = singlecoms,
                   lwdOuter = 0.5, lwdInner = 2.5, pointSize = 4,
                   color = "black", zeroColor = "red", zeroType = 1)
    p1 <- p1 + theme(axis.text.y = element_text(color = "black"),
                     panel.border = element_rect(fill = NA, color = "black", size = 2),
                     panel.background = element_rect(fill = "white")) +
      ggtitle(paste(levels(Dat$StartP)[i],"=>",levels(Dat$EndP)[j]))
    
    filename <- paste(dir,"/coefplot_",i,j,".png",sep = "")
    ggsave(filename = filename, p1 , width = 3.5, height = 5)  
  }
}



dat <- Res
dat$SynCom1 <- substring(dat$SynCom,1,2)
dat$SynCom2 <- substring(dat$SynCom,3,4)
dat$SynCom1 <- factor(dat$SynCom1 , levels = singlecoms)
dat$SynCom2 <- factor(dat$SynCom2 , levels = rev(singlecoms))


Pred <- NULL
for(i in 1:nrow(Res)){
  #i <- 1
  index1 <- Res2$StartP == dat$StartP[i] & Res2$EndP == dat$EndP[i] & Res2$SynCom == dat$SynCom1[i]
  index2 <- Res2$StartP == dat$StartP[i] & Res2$EndP == dat$EndP[i] & Res2$SynCom == dat$SynCom2[i]
  additiveguess <- Res2$Estimate[index1] + Res2$Estimate[index2]
  
  res <- data.frame(SynCom = paste(Res2$SynCom[index1],Res2$SynCom[index2],sep = ""),
             StartP = dat$StartP[i], EndP = dat$EndP[i],
             Estimate = additiveguess, SE = NA, t.value = NA,
             p.value = NA, SynCom1 = Res2$SynCom[index2],
             SynCom2 = Res2$SynCom[index1])
  Pred <- rbind(Pred,res)
}

dat$Type <- "Measured"
Pred$Type <- "Predicted"
Full <- rbind(dat,Pred)


p1 <- ggplot(Full, aes(x = SynCom1, y = SynCom2)) +
  facet_grid(StartP ~ EndP) +
  geom_tile(aes(fill = Estimate)) +
  geom_tile(aes(col = p.value < 0.05), size = 3, alpha = 0) +
  scale_fill_gradient2(low = "#d01c8b",mid = "white",
                       high = "#4dac26",midpoint = 0,
                       na.value = "#404040",
                       guide = guide_colorbar(title = "Estimate (cm)")) +
  scale_color_manual(values = c("#8c510a","#01665e")) +
  ggtitle("Elongation") +
  theme_classic()
p1
ggsave("elongation_images/heatmap_mix_full.png",p1, width = 8, height = 7)



dat$Predicted <- Pred$Estimate
p1 <- ggplot(dat,aes(x = Estimate, y = Predicted)) +
  facet_grid(StartP ~ EndP) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_classic()
p1
ggsave("elongation_images/estimate_vs_pred.png",p1, width = 7, height = 7)




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

