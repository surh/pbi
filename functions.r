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

gen_colors <- c("#A6CEE3",
                "#66BD63",
                "#386cb0",
                "#E5C494",
                scales::brewer_pal(type = "qual", pal = "Set3")(7)[4],
                "#DE77AE",
                scales::brewer_pal(type = "qual", pal = "Set3")(7)[6:7],
                "#000000")

#' Plot a taxon vs a variable
plotgg_taxon_vs_var <- function(Dat,taxon,var,log.transform = FALSE){
#   var <- "Pi"
#   taxon <- as.character(res$Taxon[1])
  
  dat <- Dat$Map
  dat$Taxon <- Dat$Tab[taxon,]
  dat <- dat[,c("Genotype",var,"Taxon")]
  dat <- dat[ !is.na(dat[,var]), ]
  
  if(log.transform){
    dat$Taxon <- log2(dat$Taxon + 1)
  }
  
  temp <- dat
  temp$Genotype <- "All"
  dat <- rbind(dat,temp)
  p1 <- ggplot(dat,aes_string(x = "Taxon", y = var,col = "Genotype")) +
    facet_wrap(~ Genotype,scales = "free_x") +
    geom_point(size = 4) +
    geom_smooth(method = "lm") +
    scale_color_manual(values = gen_colors) +
    scale_y_log10() +
    ggtitle(label = taxon) +
    theme_blackbox
  p1
}


#' Find Taxa that associate with a measurment
associate_measurment <- function(Dat,measurment){
  Res <- NULL
  for(otu in row.names(Dat$Tab)){
    #otu <- row.names(Dat$Tab)[1]
    
    abuns <- Dat$Tab[otu,]
    if(!all(names(abuns) == row.names(Dat$Map))){
      stop("ERROR",call. = TRUE)
    }
    
    Dat$Map$taxon <- abuns
    
    f1 <- paste("log(",measurment,") ~ Genotype + taxon + Experiment",sep = "")
    f1 <- formula(f1)
    
    m1 <- lm(f1,data = Dat$Map)
    m1.sum <- summary(m1)
    res <- m1.sum$coefficients["taxon",]
    
    res <- data.frame(Taxon = otu, Measurment = measurment, Estimate = res[1], pval = res[4])
    Dat$Map$taxon <- NULL
    
    Res <- rbind(Res,res)
  }
  return(Res)
}

#' Add linear model to a panel
panel.lm <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  #   usr <- par("usr"); on.exit(par(usr))
  #   par(usr = c(0, 1, 0, 1))
  
  # remove NAs
  mat <- cbind(x,y)
  mat <- mat[!apply(is.na(mat),1,any), ]
  x <- mat[,1]
  y <- mat[,2]
  
  m1 <- lm(y ~ x)
  points(x = x,y = y)
  abline(a = coef(m1)[1], b = coef(m1)[2],lwd=3,col = "red")
}



#' Plot correlations on the upper panels
#' 
#' With size proportional to the (absolute) correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, method = "spearman", ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  
  # remove NAs
  mat <- cbind(x,y)
  mat <- mat[!apply(is.na(mat),1,any), ]
  x <- mat[,1]
  y <- mat[,2]
  
  r <- cor(x, y,method = method)
  r.abs <- abs(r)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  #   text(0.5, 0.5, txt, cex = cex.cor * r.abs)
  text(0.5, 0.5, txt, cex = cex.cor * sqrt(r.abs))
}


#' Fit 4 models for phosphate
#' 
#' Based on the prr_propep immune function
fit_4_models <- function(f1,Dat){
  
  m1 <- matrix_glm(x = Dat, formula = f1, family = poisson(link = "log"),
                   verbose = TRUE)
  m2 <- matrix_glmNB(x = Dat, formula = f1,
                     verbose = TRUE)
  m3 <- matrix_zeroinfl(x = Dat, formula = f1, family = "poisson",
                        link = "log", verbose = TRUE)
  m4 <- matrix_zeroinfl(x = Dat, formula = f1, family = "negbin",
                        link = "log",verbose = TRUE)
  
  ###### Process models
  # Compare AIC
  aic <- data.frame(glm = m1$AIC, glm.NB = m2$AIC, zip = m3$AIC, zinb = m4$AIC)
  best.index <- apply(aic,1,which.min)
  to_remove <- sapply(best.index,length)
  to_remove <- which(to_remove  == 0)
  best.index[ to_remove ] <- 5
  best.index <- unlist(best.index)
  
  aic$best <-  colnames(aic)[ best.index ]
  table(aic$best)
  
  models <- list(m1 = m1, m2 = m2, m3 = m3, m4 = m4, aic = aic)
  
  return(models)
}

#' Take summary of all four models and combine them into one
summary_4models <- function(models){
  models.sum <- list(m1 = summary(models$m1),
                     m2 = summary(models$m2),
                     m3 = summary(models$m3),
                     m4 = summary(models$m4),
                     aic = models$aic)
  Full <- combine_models(models.sum)
  
  return(Full)
}

