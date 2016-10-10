library(ggplot2)
dat <- read.table("~/rhizogenomics/data/phosphate/colonization_synthetic/Pi",
                  sep = "\t", header = TRUE)

dat$Treatment <- factor(dat$Treatment, levels = c("HighPiNone","HighPi35",
                                                  "LowPiNone","LowPi35"))
levels(dat$Treatment) <- c("625uM","625uM + SynCom","50uM","50uM + SynCom")

p1 <- ggplot(dat,aes(x = Genotype, y = Pi.content)) +
  geom_boxplot(aes(col = Treatment), outlier.colour = NA,
               position = position_dodge(width = 0.9)) +
  geom_point(aes(col = Treatment),
             position = position_jitterdodge(dodge.width = 0.9)) +
  ylab("Cellular Pi (umol/FWmg*1000)") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", size = 3, fill = NA),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(color = "black", angle = 90))
p1
#ggsave("syncomP_pi.png",p1,width = 5, height = 4)
ggsave("syncomP_pi.svg",p1,width = 7, height = 4)

dat$gentreat <- with(dat, interaction(Genotype, Treatment, sep = "."))

m1 <- aov(log(Pi.content) ~ Treatment + Genotype + gentreat, data = dat)
summary(m1)

tuk <- agricolae::HSD.test(m1,"gentreat")
tuk

