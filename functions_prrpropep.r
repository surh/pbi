monte_carlo_test <- function(dat, wt = "Col-0",dist.method = "manhattan",
                             nperm = 1000){
  dat$Enrichment <- as.numeric(as.character(dat$Enrichment))
  dat <- subset(dat, Genotype != wt)
  dat$Genotype <- droplevels(dat$Genotype)
  mat <- acast(data = dat, formula = OTU ~ Genotype,
               value.var = "Enrichment")
  
  all <- mat
  #all <- all[ , -which(colnames(all) == wt) ]
  
  all.d <- dist(t(all),method=dist.method)
  all.d <- as.matrix(all.d)
  D.reps <- array(dim=c(nrow(all.d),nrow(all.d),nperm),
                  dimnames = list(colnames(all),
                                  colnames(all),
                                  1:nperm))
  for (i in 1:nperm){
    #i <- 1
    all.rep <- apply(all,2,sample)
    all.rep.d <- dist(t(all.rep),method=dist.method)
    all.rep.d <- as.matrix(all.rep.d)
    D.reps[,,i] <- all.d < all.rep.d
  }
  pval <- as.matrix(as.dist(1 - apply(D.reps,MARGIN=c(1,2),sum) / nperm))
  
}

doughnut_plot <- function(dat,top_phyla, mut_order, text_size = 10){
  count <- ftable(Enrichment~Genotype,dat)[-1,-3]
  row.names(count) <- mut_order[c(-1,-length(mut_order))]
  colnames(count) <- c("Depletion", "Enrichment")
  count <- melt(count,varnames = c("Genotype","Type"),value.name = "Total")
  count$Phylum <- NA
  count$ymin <- NA
  count$ymax <- NA

  enrich <- dat[ dat$Enrichment == 1, ]
  enrich <- calculate_enrich(enrich = enrich,top_phyla = top_phyla,
                             mut_order = mut_order)
  enrich$Type <- "Enrichment"
  deplete <- dat[ dat$Enrichment == -1, ]
  deplete <- calculate_enrich(enrich = deplete,top_phyla = top_phyla,
                              mut_order = mut_order)
  deplete$Type <- "Depletion"
  enrich <- rbind(enrich,deplete)
  enrich$Type <- factor(enrich$Type, levels = c("Enrichment","Depletion"))
  
  enrich$Total <- NA
  enrich <- rbind(enrich,count)
  
  p1 <- ggplot(enrich, aes(xmax = 4, xmin = 3)) +
    facet_grid(Type ~ Genotype) +
    geom_rect(data = subset(enrich, !is.na(ymin)), aes(ymax = ymax, ymin = ymin, fill = Phylum), color = NA) +
    xlim(c(0,4)) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = phyla_colors) +
    geom_text(data = subset(enrich, is.na(ymin)),
              aes(label = Total), x = 0, y = 0, size = text_size,
              fontface = "bold") +
    theme_blackbox +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
  #p1
  return(p1)
}



#' Plot fracgen heatmap
#' 
#' First selects the variables that match the egular expression and then
#' finds the significant ones using sig_thres. Requires columns with the
#' appropriate names in Full, meaning Taxon, Variable and q.value
#' 
#' The it processes the dataset object. Rescales the abundance to percent,
#' requires a Depth variable in the Map, and removes any taxa that are never
#' significant. Removes Soil and R fraction samples (requires a Fraction
#' variable with R and Soil values in it), and pools samples by genotype
#' (requires a Genotype variable in the Map)
#' 
#' 
#' @param Full output from combine_models
#' @param Dat Dataset object
#' @param group character vector indicating experiments to keep.
#' @param variable regular expression used to extract variable of interest
#' @param sig_thres significance threshold
#' @param pool.FUN function to use when pooling observations for the heatmap
#' @param phyl.level level at which phylum is found in Tax
#' @param mut_order
plot_fragen_heatmap <- function(Full,Dat,group, variable = "FracgenEC:",
                                sig_thres = 0.05, pool.FUN = mean, phyl.level = 3,
                                mut_order){
  #group <- c("SL15")
  #variable <- "FracgenEC:"
  #sig_thres <- 0.05
  #Dat <- Dat.fam
  #Full <- Full
  #pool.FUN <- mean
  #phyl.level <- 3
  #mut_order <- pepr_order
  
  # Get significant
  Gen <- Full[ grep(variable,Full$Variable), ]
  fracgen_effect <- as.character(unique(Gen$Taxon[ Gen$q.value < sig_thres ]))
  Gen <- Gen[ Gen$Taxon %in% fracgen_effect, ]
  Gen <- droplevels(Gen)
  Gen$Taxon <- factor(Gen$Taxon)
  Gen$significant <- Gen$q.value < sig_thres
  
  
  # Process Dat
  pool.FUN <- match.fun(pool.FUN)
  # Get samples in group and rescale
  to_remove <- row.names(Dat$Map)[ !(Dat$Map$Experiment %in% group) ]
  Dat.group <- remove_samples(Dat, samples = to_remove)
#   Dat.group <- create_dataset(Tab = t(100 * t(Dat.group$Tab) / Dat.group$Map$Depth),
#                               Map = Dat.group$Map,
#                               Tax = data.frame(ID = row.names(Dat.group$Tab), 
#                                                Taxonomy = row.names(Dat.group$Tab),
#                                                row.names = row.names(Dat.group$Tab)))
  
  Dat.group <- create_dataset(Tab = t(100 * t(Dat.group$Tab) / Dat.group$Map$Depth),
                              Map = Dat.group$Map,
                              Tax = Dat.group$Tax)
  
  
  Dat.group$Map$Genotype <- factor(Dat.group$Map$Genotype, levels = mut_order)
  # Remove non significant taxa
  Dat.group <- remove_taxons(Dat = Dat.group,
                             taxons = row.names(Dat.group$Tab)[ !(row.names(Dat.group$Tab) %in% levels(Gen$Taxon)) ])
  # Remove Soil and R samples
  Dat.group <- remove_samples(Dat = Dat.group,
                              samples = colnames(Dat.group$Tab)[ Dat.group$Map$Fraction %in% c("R","Soil") ])
  # Pool by genotype
  Dat.group <- pool_samples(Dat = Dat.group, groups = "Genotype", FUN = mean)
  
  # Create plotting object
  dat <- as.data.frame(Dat.group$Tab)
  dat$OTU <- row.names(dat)
  dat <- melt(dat,id.vars=c("OTU"),value.name="Abundance",variable.name="Genotype")
  dat$Abundance <- log2(dat$Abundance + 1)
  dat$Enrichment <- 0
  dat$Enrichment[ paste(dat$OTU,variable,dat$Genotype,sep = "") %in%
                    paste(Gen$Taxon,Gen$Variable,sep="")[ Gen$significant & Gen$Estimate > 0 ] ] <- 1
  dat$Enrichment[ paste(dat$OTU,variable,dat$Genotype,sep = "") %in%
                    paste(Gen$Taxon,Gen$Variable,sep="")[ Gen$significant & Gen$Estimate < 0 ] ] <- -1
  dat$Phylum <- get_tax_level(Dat.group$Tax[as.character(dat$OTU),],level = phyl.level)
  dat$Type <- "Mutant"
  dat$Type[ dat$Genotype == "Col-0" ] <- "WT"
  dat <- dat[order(dat$Phylum,dat$Enrichment),]
  dat$OTU <- factor(dat$OTU,levels = unique(dat$OTU))
  dat$Type <- factor(dat$Type,levels=c("WT","Mutant"))
  dat$Enrichment <- factor(dat$Enrichment,levels=c("-1","1","0"))
  dat <- dat[order(dat$Enrichment),]
  
  p1 <- ggplot(dat,aes(x=Genotype,y=OTU)) +
    #   geom_tile(aes(fill=Abundance,col=Enrichment),size=2) +
    geom_tile(aes(fill=Abundance)) +
    geom_tile(aes(col=Enrichment),size=2,alpha=1,fill = NA) +
    facet_grid(Phylum ~ Type,scales="free",space="free") +
    scale_fill_continuous(low="white",high="black") +
    scale_colour_manual(values = c("blue","red",NA)) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(face="bold",colour="black",
                                     angle = 90, size = 20),
          panel.background = element_blank(),
          strip.text.y = element_text(angle = 360))
  #p1
  
  return(p1)
}

#' Take a set of enrichments or depletions and make
#' The object for geom_rect
calculate_enrich <- function(enrich, top_phyla = top_phyla, mut_order){
  #  mut_order <- pepr_order
  
  enrich$Enrichment <- 1
  enrich$Phylum[ !(enrich$Phylum %in% top_phyla) ] <- "Low Abundance"
  enrich$Phylum <- factor(enrich$Phylum, levels = top_phyla )
  enrich$Genotype <- factor(enrich$Genotype, levels = mut_order[c(-1, -length(mut_order))])
  enrich <- acast(data = enrich, formula = Phylum ~ Genotype, fun.aggregate = sum, value.var = "Enrichment", drop = FALSE)
  
  enrich <- t(t(enrich) / colSums(enrich))
  enrich[ is.na(enrich) ] <- 0
  enrich.ymax <- apply(enrich,2,cumsum)
  enrich.ymin <- apply(enrich,2,function(vec) c(0,cumsum(vec)[-length(vec)]))
  row.names(enrich.ymin) <- row.names(enrich)
  
  enrich.ymin <- melt(enrich.ymin,varnames = c("Phylum","Genotype"),value.name = "ymin")
  enrich.ymax <- melt(enrich.ymax,varnames = c("Phylum","Genotype"),value.name = "ymax")
  enrich <- enrich.ymin
  enrich$ymax <- enrich.ymax$ymax
  
  return(enrich)
}

#' Fit four ZINB models
#' 
#' Fits the 4 zinb models, and computes the best fit
#' for each taxon based on AIC
#' 
#' First it removes any sample that is not in the specified experiments. Requires
#' a Experiment column in the Map
#' 
#' Then it rmeoves OTUs under the measurable threshold (default 25in5)
#' 
#' It reorders the Genotype variable. It also creates Fracgen and LDepth.
#' It requires Genotype, Fraction and Depth columns to be present.
#' 
#' # Fits four ZINB models and computes the best for each taxa according to AIC
#' 
#' @param Dat a Dataset object
#' @param f1 right hand side formula for the model
#' @param group character vector of experiments to include
#' @param mutant_order Order if variables in Genotype column of the Map
#' @param min_reads_otu Minimum count for an otu in a sample for that sample to
#' be included in the measurable threshold
#' @param min_samples_otu Minimum samples that pass the min_reads_otu threshold for
#' a given otu to be counted as measurable.
#' 
fit_4_models <- function(Dat, f1, group, mutant_order, min_reads_otu = 25,
                         min_samples_otu = 5){
  # Dat <- Dat.fam
  # f1 <- f1
  # group <- c("SL15","SL25","SL29")
  # mutant_order <- mutant_order
  # min_reads_otu <- 25
  # min_samples_otu <- 5
  
  # remove samples in different experiments
  # to_remove <- row.names(Dat$Map)[ !(Dat$Map$Experiment %in% group) ]
  # Dat.group <- remove_samples(Dat, samples = to_remove)
  Dat.group <- subset(Dat, Experiment %in% group, clean = TRUE, drop = TRUE)
  
  # Remove otus under the  measurable threshold
  # to_remove <- !findGoodOTUs(OTUtab = Dat.group$Tab,min_reads_otu = min_reads_otu,
  #                            min_samples_otu = min_samples_otu)
  # to_remove <- row.names(Dat.group$Tab)[ to_remove ]
  # Dat.group <- remove_taxons(Dat = Dat.group, to_remove)
  # Dat.group <- clean(Dat.group)
  Dat.group <- measurable_taxa(Dat,min_reads_otu = min_reads_otu,
                               min_samples_otu = min_samples_otu,
                               method = "absolute",table = TRUE, clean = TRUE)
  
  
  # Provess mappinfg file
  Dat.group$Map$Genotype <- factor(Dat.group$Map$Genotype, levels = mutant_order)
  Dat.group$Map$Fracgen <- factor(Dat.group$Map$Fraction, levels = c("EC","R","Soil"))
  Dat.group$Map$Fracgen <- droplevels(Dat.group$Map$Fracgen:Dat.group$Map$Genotype)
  Dat.group$Map$LDepth <- log10(Dat.group$Map$Depth)
  
  # Fit 4 models
  m1 <- matrix_glm(x = Dat.group, formula = f1, family = poisson(link = "log"),
                   verbose = TRUE)
  m2 <- matrix_glmNB(x = Dat.group, formula = f1,
                     verbose = TRUE)
  m3 <- matrix_zeroinfl(x = Dat.group, formula = f1, family = "poisson",
                        link = "log", verbose = TRUE)
  m4 <- matrix_zeroinfl(x = Dat.group, formula = f1, family = "negbin",
                        link = "log",verbose = TRUE)
  
  # Compare AIC
  #m1 <- models$m1
  #m2 <- models$m2
  #m3 <- models$m3
  #m4 <- models$m4
  #aic <- models$aic
  
  aic <- data.frame(glm = m1$AIC, glm.NB = m2$AIC, zip = m3$AIC, zinb = m4$AIC)
  best.index <- apply(aic,1,which.min)
  to_remove <- sapply(best.index,length)
  to_remove <- which(to_remove  == 0)
  best.index[ to_remove ] <- 5
  best.index <- unlist(best.index)
  
  aic$best <-  colnames(aic)[ best.index ]
  
  #models <- list(m1 = m1, m2 = m2, m3 = m3, m4 = m4, aic = aic)
  return(list(m1 = m1, m2 = m2, m3 = m3, m4 = m4, aic = aic))
}


#' Merges fur zinb models
#' 
#' Removes undefined variables and perfroms multiple testing
#' correction
#' 
#' @param models.sum list of 5 elements m1-m4 for the four zinb models
#' and aic which is a data.frame with column best that states which is
#' the best model.
combine_models <- function(models.sum){
  m1.glm <- models.sum$m1$coefficients
  m1.glmNB <- models.sum$m2$coefficients
  m1.zip <- models.sum$m3$coefficients
  m1.zinb <- models.sum$m4$coefficients
  models.sum$aic$Taxon <- row.names(models.sum$aic)
  
  Full <- NULL
  for(taxa in unique(m1.glm$Taxon)){
    #taxa <- levels(m1.glm$Taxon)[1]
    cat(taxa,"\n")
    model <- models.sum$aic$best[ models.sum$aic$Taxon == taxa ]
    if(is.na(model)) next
    if(model == "glm"){
      m1 <- m1.glm
    }else if(model == "glm.NB"){
      m1 <- m1.glmNB
    }else if(model == "zip"){
      m1 <- m1.zip
    }else if(model == "zinb"){
      m1 <- m1.zinb
    }
    m1 <- m1[m1$Taxon == taxa,]
    m1$model <- model
    Full <- rbind(Full,m1)
  }
  Full$q.value <- p.adjust(Full$p.value,method='fdr')
  
  Full <- Full[ !apply(apply(Full [ , c(3:6,8) ], 1, is.na),2,all), ]
  
  return(Full)
}

#' Fracgen Site Diversity
#' 
#' Perfoms site diversity analyis with the compare_site_diversity function
#' from the AMOR package, and then estimates a 95% confidence interval for
#' the mean site diversity, as and creates a plot
#' 
#' @param Dat A Dataset object.
#' @param nperm Number of permutations to perform fore Site Diversity
#' @param add_colors Logical. Indicates whether to use default coloring for
#' different fractions.
fracgen_sitediv <- function(Dat,nperm = 20, add_colors = TRUE){
  sitediv <- compare_site_diversity(Dat = Dat, factor = "Fracgen",
                                    divfun = total_richness, nperm = nperm)
  
  #### Plot
  # Prepare
  sitediv$Fraction <- do.call(rbind,strsplit(x = as.character(sitediv$group),
                                             split = "_"))[,1]
  sitediv$Genotype <- do.call(rbind,strsplit(x = as.character(sitediv$group),
                                             split = "_"))[,2]
  confint <- sitediv$mean + matrix(qnorm(p = c(0.025, 0.975)), 
                                   ncol = 2, nrow = nrow(sitediv), byrow = TRUE) * sitediv$sd
  colnames(confint) <- c("lower", "upper")
  sitediv <- cbind(sitediv, confint)
  
  # Plot
  p1 <- ggplot(sitediv, aes(x = nsites, y = mean,
                            col = Fraction, 
                            group = group,
                            fill = Fraction)) +
    geom_line()
  p1 <- p1 + geom_ribbon(aes(ymin = lower, ymax = upper), 
                         alpha = 0.2)
  if(add_colors){
    p1 <- p1 + scale_fill_manual(values = c("lawngreen","blue2","grey35")) +
      scale_color_manual(values = c("lawngreen","blue2","grey35"))
  }
  p1 <- p1 + ylab(label = "Richness (S)") +
    theme_blackbox
  #p1
  
  return(p1)
}

group_fracgen_sitediv <- function(Dat, group, nperm = 20, add_colors = TRUE){
  Dat.group <- remove_samples(Dat = Dat,
                              samples = row.names(Dat$Map)[ !(Dat$Map$Experiment %in% group) ])
  p1 <- fracgen_sitediv(Dat = Dat.group, nperm = nperm, add_colors = add_colors)
  return(p1)
}

plot_pco <- function(dis.pco){
  totvar <- 100 * (dis.pco$eig[ dis.pco$eig >= 0 ]) / sum(dis.pco$eig[ dis.pco$eig >= 0 ])
  totvar <- round(totvar,1)
  p1 <- plotgg(dis.pco, shape = "Fraction", col = "Fraction", point_size = 4,
               components = c("PCo1","PCo2"))
  p1 <- p1 + scale_color_manual(values = c("lawngreen","blue2","grey35")) +
    coord_fixed() +
    xlab( label = paste("PCo1 (", totvar[1],"%)",sep = "")) +
    ylab( label = paste("PCo2 (", totvar[2],"%)",sep = "")) +
    theme(axis.text = element_blank())
  #p1
  
  # All components 1,3
  p2 <- plotgg(dis.pco, shape = "Fraction", col = "Fraction", point_size = 4,
               components = c("PCo1","PCo3"))
  p2 <- p2 + scale_color_manual(values = c("lawngreen","blue2","grey35")) +
    coord_fixed() +
    xlab( label = paste("PCo1 (", totvar[1],"%)",sep = "")) +
    ylab( label = paste("PCo3 (", totvar[3],"%)",sep = "")) +
    theme(axis.text = element_blank())
  #p1
  
  return(list(p1,p2))
}

plot_group_pco <- function(Dat, dis, group){
  Dat.group <- remove_samples(Dat = Dat,
                              samples = row.names(Dat$Map)[ !(Dat$Map$Experiment %in% group) ])
  dis.group <- dis[colnames(Dat.group$Tab),colnames(Dat.group$Tab)]
  dis.pco <- PCO(Dat = Dat.group,dim = 3,distfun = function(x) as.dist(dis.group))
  res <- plot_pco(dis.pco = dis.pco)
  
  res[[1]] <- res[[1]] + ggtitle(label = paste(c("group",group), collapse = " "))
  res[[2]] <- res[[2]] + ggtitle(label = paste(c("group",group), collapse = " "))
  
  return(res)
}

process_unifrac_table <- function(file, unifrac, dir, prefix){
  load(file)
  dis <- read.table(file = unifrac)
  if( all(colnames(dis) == colnames(Dat$Tab))){
    dis.pco <- PCO(Dat = Dat,dim = 3,distfun = function(x) as.dist(dis))
    
    # All componets 1,3
    res <- plot_pco(dis.pco = dis.pco)
    filename <- paste(dir,"/",prefix,"_full_PCo1PCo2.png",sep="")
    ggsave(filename,res[[1]],width = 6,height = 6)
    filename <- paste(dir,"/",prefix,"_full_PCo1PCo3.png",sep="")
    ggsave(filename,res[[2]],width = 6,height = 6)
    
    #### For group
    group <- c("SL25","SL29")
    res <- plot_group_pco(Dat = Dat,dis = dis,group = group)
    filename <- paste(dir,"/",prefix,"_PRR_PCo1PCo2.png",sep="")
    ggsave(filename,res[[1]],width = 6,height = 6)
    filename <- paste(dir,"/",prefix,"_PRR_PCo1PCo3.png",sep="")
    ggsave(filename,res[[2]],width = 6,height = 6)
    
    group <- c("SL15","SL26")
    res <- plot_group_pco(Dat = Dat,dis = dis,group = group)
    filename <- paste(dir,"/",prefix,"_Propep_PCo1PCo2.png",sep="")
    ggsave(filename,res[[1]],width = 6,height = 6)
    filename <- paste(dir,"/",prefix,"_Propep_PCo1PCo3.png",sep="")
    ggsave(filename,res[[2]],width = 6,height = 6)
    
    ### For experiment
    group <- c("SL15")
    res <- plot_group_pco(Dat = Dat,dis = dis,group = group)
    filename <- paste(dir,"/",prefix,"_SL15_PCo1PCo2.png",sep="")
    ggsave(filename,res[[1]],width = 6,height = 6)
    filename <- paste(dir,"/",prefix,"_SL15_PCo1PCo3.png",sep="")
    ggsave(filename,res[[2]],width = 6,height = 6)
    
    group <- c("SL25")
    res <- plot_group_pco(Dat = Dat,dis = dis,group = group)
    filename <- paste(dir,"/",prefix,"_SL25_PCo1PCo2.png",sep="")
    ggsave(filename,res[[1]],width = 6,height = 6)
    filename <- paste(dir,"/",prefix,"_SL25_PCo1PCo3.png",sep="")
    ggsave(filename,res[[2]],width = 6,height = 6)
    
    group <- c("SL26")
    res <- plot_group_pco(Dat = Dat,dis = dis,group = group)
    filename <- paste(dir,"/",prefix,"_SL26_PCo1PCo2.png",sep="")
    ggsave(filename,res[[1]],width = 6,height = 6)
    filename <- paste(dir,"/",prefix,"_SL26_PCo1PCo3.png",sep="")
    ggsave(filename,res[[2]],width = 6,height = 6)
    
    group <- c("SL29")
    res <- plot_group_pco(Dat = Dat,dis = dis,group = group)
    filename <- paste(dir,"/",prefix,"_SL29_PCo1PCo2.png",sep="")
    ggsave(filename,res[[1]],width = 6,height = 6)
    filename <- paste(dir,"/",prefix,"_SL29_PCo1PCo3.png",sep="")
    ggsave(filename,res[[2]],width = 6,height = 6)
  }
}

#' Make UniFrac plot s on mutant group
#' 
#' Trying to update it so I can use it for beta diversity in general
#' 
#' @param Dat A Dataset object.
#' @param group A vector indicating which levels from the Experiment variable to include
#' @param unifrac A file path to a unifrac distance matrix produced from QIIME. If NULL,
#' it will use distfun to create a distance matrix.
#' @param prefix prefix for output filenames.
#' @param dir Output directory
#' @param mutant_colors Vector specifying the colors to use for Genotype level.
#' @param postfix postfx for output files
#' @param ordermut If specified. Levels of Genotype variable will be re-ordered
#' according to it.
#' @param method "PCA" or "PCO".
#' @param pca.cor Logical. Whether to use correlation matrix instead of covariance
#' in PCA scaling.
#' @param distfun. A function  that produces a dist object from the abundance matrix in Dat.
make_mutant_group_unifrac <- function(Dat, group,unifrac, prefix, dir,
                                      mutant_colors, postfix, ordermut = NULL,
                                      method = "PCO", pca.cor = TRUE,
                                      distfun = function(x) vegan::vegdist(t(x), method = "cao")){
#   Dat <- Dat
#   Dat <- Dat.bin
#   group <- c("SL15")
#   group <- c("SL25","SL29") # missing SL26
#   unifrac <- NULL
#   prefix <- "cao"
#   dir <- dir
#   mutant_colors <- pepr_colors
#   postfix <- "peprpropep"
#   ordermut <- NULL
#   method <- "PCA"
#   distfun <- function(x) vegan::vegdist(t(x),method = "bray")
#   distfun <- function(x) dist(t(x))
#   pca.cor <- TRUE
  
  if(is.null(unifrac)){
    distfun <- match.fun(distfun)
    dis <- distfun(Dat$Tab)
    dis <- as.matrix(dis)
  }else{
    dis <- read.table(file = unifrac)
    dis <- as.matrix(dis)
  }
  
  Dat.group <- remove_samples(Dat = Dat,
                              samples = row.names(Dat$Map)[ !(Dat$Map$Experiment %in% group) ])
  Dat.group <- clean(Dat.group)
  dis.group <- dis[colnames(Dat.group$Tab), colnames(Dat.group$Tab)]
  if(!is.null(ordermut)){
    Dat.group$Map$Genotype <- factor(Dat.group$Map$Genotype, levels = ordermut) 
  }
  
  if(method == "PCO"){
    dis.pco <- PCO(Dat = Dat.group,dim = 3,
                   distfun = function(x) as.dist(dis.group))
    p1 <- plotgg(dis.pco, shape = "Fraction", col = "Genotype",
                 point_size = 4)
    dat <- p1$data
    
    totvar <- 100 * (dis.pco$eig[ dis.pco$eig >= 0 ]) 
    totvar <- totvar / sum(dis.pco$eig[ dis.pco$eig >= 0 ])
    totvar <- round(totvar,1)
    
    xlab <- paste("PCo1 (", totvar[1],"%)",sep = "")
    ylab <- paste("PCo2 (", totvar[2],"%)",sep = "")
    
    components <- c("PCo1","PCo2")
  }else if(method == "PCA"){
    #dis.pca <- create_dataset(Tab = as.matrix(dis.group),
    #                          Map = Dat.group$Map)
    #dat.pca <- PCA(Dat = dis.pca, cor = TRUE)
    
    # For PCA need to remove all rows that are all equeal
    Dat.group <- remove_taxons(Dat = Dat.group,
                               taxons = row.names(Dat.group$Tab)[ apply(Dat.group$Tab,1,function(x) length(unique(x))) == 1 ])
    
    dat.pca <- PCA(Dat = Dat.group, cor = pca.cor)
    p1 <- plotgg(dat.pca, shape = "Fraction", col = "Genotype",
                 point_size = 4)
    #p1
    dat <- p1$data
    
    totvar <- 100*dat.pca$sdev^2 / dat.pca$totvar
    totvar <- round(totvar,1)
    
    xlab <- paste("PC1 (", totvar[1],"%)",sep = "")
    ylab <- paste("PC2 (", totvar[2],"%)",sep = "")
    
    components <- c("PC1","PC2")
    #return(p1)
  }
  
  p1 <- ggplot(dat, aes_string(x = components[1], y = components[2])) +
    geom_point(aes(shape = Fraction, col = Fraction,
                   fill = Genotype), size = 4) +
    scale_color_manual(values = rev(frac_colors)) +
    scale_fill_manual(values = mutant_colors) +
    scale_shape_manual(values = c(22,25,21)) +
    xlab( label = xlab) +
    ylab( label = ylab) +
    guides(fill = guide_legend(override.aes = list(shape = c(21)))) +
    theme_blackbox +
    theme(axis.text = element_blank())
  #p1
  
  filename <- paste(dir,"/",prefix,"_",postfix,".png",sep="")
  ggsave(filename,p1,width = 6,height = 5)
  
  return(p1)
}

############# CONSTANTS #############
pepr_colors <- c(scales::brewer_pal(type = "qual",pal = 3)(12)[1],
                 scales::brewer_pal(type = "qual",pal = 6)(9)[c(1,3,4,5,6)],
                 "#000000")
prr_colors <- c(scales::brewer_pal(type = "qual",pal = 3)(12)[1],
                scales::brewer_pal(type = "qual",pal = 2)(8)[1:7],
                "#000000")
frac_colors <- c("lawngreen","blue2","grey35")
prr_order <- c("Col-0","bak1","bak1bkk1",
  "bak1bkk1cerk1",
  "cerk1efr1fls2","eds1",
  "pad4rar1","rar1","Soil")
top_phyla <- c("Root; k__Bacteria; p__Acidobacteria",
               "Root; k__Bacteria; p__Actinobacteria",
               "Root; k__Bacteria; p__Bacteroidetes",
               "Root; k__Bacteria; p__Chloroflexi",
               "Root; k__Bacteria; p__Firmicutes",
               "Root; k__Bacteria; p__Proteobacteria",
               "Root; k__Bacteria; p__Verrucomicrobia",
               "Low Abundance")
# phyla_colors <- c("aquamarine1","blue2","yellow",
#                   "chartreuse4","darkviolet","darkorange1",
#                   "darkorchid1","black")
phyla_colors <- c("#41F0AC","#0000C0","#FFFF00","#008000",
                  "#B856D7","#FF8000","#8D8DFF","#000000")


pepr_order <- c("Col-0","pepr1pepr2","Propep1-1","Propep1-2",
                "Propep1-4","Propep1-6","Soil")

