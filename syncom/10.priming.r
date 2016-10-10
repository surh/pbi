library(ggplot2)

Tab <- read.table("~/rhizogenomics/data/phosphate/synthetic_physiology/2016-05-23.fertilization2.txt",
                  sep = "\t", header = TRUE)

agg_f1 <- Pi~Pre.Pi + Bacteria + Experiment + day
dat <- aggregate(agg_f1,data = Tab,FUN = mean)
dat$Pi.sd <- aggregate(agg_f1,data = Tab,FUN = sd)$Pi
row.names(dat) <- paste(dat$Pre.Pi,dat$Bacteria,dat$Experiment,dat$day, sep = ".")
Tab$day0 <- dat[ paste(Tab$Pre.Pi,Tab$Bacteria,Tab$Experiment,"0",sep = "."),"Pi"]
Tab$Pi.norm <- (Tab$Pi - Tab$day0) / Tab$day0
Tab$day <- factor(Tab$day)

pal <- colorRampPalette(colors = rev(c("#00441b","#99d8c9")))
p1 <- ggplot(subset(Tab, Pre.Pi %in% c("0uM","50uM","625uM")),
             aes(x = day, y = Pi.norm, color = Pre.Pi)) +
  facet_grid(~ Bacteria) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = pal(5)[c(1,4,5)]) +
  ylab("Relative Pi increase") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold",size = 12),
        strip.background = element_blank(),
        axis.line.y = element_line(color = "black",size = 2),
        axis.line.x = element_line(color = "black",size = 2))
p1
ggsave("priming_a.png",p1,width = 6, height = 4)
ggsave("priming_a.svg",p1,width = 6, height = 4)

p1 <- ggplot(subset(Tab, Pre.Pi %in% c("50uM","625uM")),
             aes(x = day, y = Pi.norm, color = Pre.Pi)) +
  facet_grid(~ Bacteria) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = pal(5)[c(4,5)]) +
  ylab("Relative Pi increase") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold",size = 12),
        strip.background = element_blank(),
        axis.line.y = element_line(color = "black",size = 2),
        axis.line.x = element_line(color = "black",size = 2))
p1
ggsave("priming_b.png",p1,width = 6, height = 4)
ggsave("priming_b.svg",p1,width = 6, height = 4)

