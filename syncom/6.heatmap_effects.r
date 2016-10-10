library(reshape2)
library(gplots)

add_condition <- function(Dat){
  
  Dat$Condition <- NULL
  Dat$Condition[ Dat$StartP == "-Pi,0.5%Suc" & Dat$EndP == "100 uM,0%Suc" ] <- "minusP_100uM"
  Dat$Condition[ Dat$StartP == "+Pi,0.5%Suc" & Dat$EndP == "100 uM,0%Suc" ] <- "plusP_100uM"
  Dat$Condition[ Dat$StartP == "-Pi,0.5%Suc" & Dat$EndP == "30 uM,0%Suc" ] <- "minusP_30uM"
  Dat$Condition[ Dat$StartP == "+Pi,0.5%Suc" & Dat$EndP == "30 uM,0%Suc" ] <- "plusP_30uM"
  
  return(Dat)
}



Elongation <- read.table("../today2/elongation_single_community_test.txt", header = TRUE, sep = "\t")
Area <- read.table("../today2/area_single_community_test.txt", header = TRUE, sep = "\t")
Pi <- read.table("../today2/pi_single_community_test.txt", header = TRUE, sep = "\t")


Elongation <- add_condition(Elongation)
Area <- add_condition(Area)
Pi <- add_condition(Pi)



Dat <- cbind(acast(data = Pi, formula = SynCom ~ Condition, value.var = "Estimate"),
             acast(data = Elongation, formula = SynCom ~ Condition, value.var = "Estimate"),
             acast(data = Area, formula = SynCom ~ Condition, value.var = "Estimate"))
colnames(Dat) <- paste( rep(c("Pi","Elongation","Area"),each = 4),colnames(Dat),sep = ".")


pal <- colorRampPalette(colors = c("#d01c8b","#ffffff", "#4dac26"))
tiff("heatmap_syncom_effects.tif",width = 2000, height = 1500, res = 250)
heatmap.2(x = scale(Dat, center = FALSE, scale = TRUE), main = "Scaled SynCom effects\non plant phenotypes",
          dendrogram = "row", trace = "none", col = pal(50),
          margins = c(12,5), Colv = FALSE)
dev.off()


