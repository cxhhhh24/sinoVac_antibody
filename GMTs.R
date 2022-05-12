library(haven)
library(stringr)
library(dplyr)
library(mgcv)
library(DescTools)
library(stringr)
library(reshape2)
library(patchwork)
library(ggplot2)


#### Young adults ####

#-- data cleaning
dat <- read.csv("cleaned_data/dat_titer.csv")
colnames(dat)
dat_new <- dat[, c("SUBJID", "titer_1", "titer_2", "titer_3", "titer_4", "titer_5", "titer_6")]
dat_new <- melt(dat_new, id = c("SUBJID"), value.name = "titer")
colnames(dat_new)[2] <- "period"
dat <- left_join(dat[, 2:5], dat_new)

dat$dose[dat$ARMCD %in% c("A-1", "A-2")] <- "3ug"
dat$dose[dat$ARMCD %in% c("B-1", "B-2")] <- "6ug"
dat$dose[dat$ARMCD %in% c("C-1", "C-2")] <- "placebo"

dat$group[substr(dat$SUBJID,0,1) == "C" & as.numeric(substr(dat$SUBJID,2,4))<=150] <- "0_14_42"
dat$group[substr(dat$SUBJID,0,1) == "C" & as.numeric(substr(dat$SUBJID,2,4))>150] <- "0_14_8m"
dat$group[substr(dat$SUBJID,0,1) == "D" & as.numeric(substr(dat$SUBJID,2,4))<=150] <- "0_28_56"
dat$group[substr(dat$SUBJID,0,1) == "D" & as.numeric(substr(dat$SUBJID,2,4))>150] <- "0_28_8m"

dat$period <- as.character(dat$period)
dat$period[dat$period == "titer_1"] <- "Baseline"
dat$period[dat$period == "titer_2" & dat$group %in%  c("0_14_42", "0_14_8m")] <- "Day 14 after dose 2"
dat$period[dat$period == "titer_2" & dat$group %in%  c("0_28_56", "0_28_8m")] <- "Day 28 after dose 2"
dat$period[dat$period == "titer_3" & dat$group %in%  c("0_14_42", "0_14_8m")] <- "Day 28 after dose 2"
dat$period[dat$period == "titer_3" & dat$group %in%  c("0_28_56")] <- "Day 28 after dose 3"
dat$period[dat$period == "titer_3" & dat$group %in%  c("0_28_8m")] <- "Day 180 after dose 2"
dat$period[dat$period == "titer_4" & dat$group %in%  c("0_14_42")] <- "Day 28 after dose 3"
dat$period[dat$period == "titer_4" & dat$group %in%  c("0_14_8m")] <- "Day 180 after dose 2"
dat$period[dat$period == "titer_4" & dat$group %in%  c("0_28_56")] <- "Day 180 after dose 3"
dat$period[dat$period == "titer_4" & dat$group %in%  c("0_28_8m")] <- "Day 28 after dose 3"
dat$period[dat$period == "titer_5" & dat$group %in%  c("0_14_42")] <- "Day 180 after dose 3"
dat$period[dat$period == "titer_5" & dat$group %in%  c("0_14_8m")] <- "Day 14 after dose 3"
dat$period[dat$period == "titer_5" & dat$group %in%  c("0_28_56")] <- "1 year after dose 3"
dat$period[dat$period == "titer_5" & dat$group %in%  c("0_28_8m")] <- "Day 180 after dose 3"
dat$period[dat$period == "titer_6" & dat$group %in%  c("0_14_42")] <- "1 year after dose 3"
dat$period[dat$period == "titer_6" & dat$group %in%  c("0_14_8m")] <- "Day 180 after dose 3"

dat$titer[dat$titer == "<4"] <- 2
dat$titer <- as.numeric(dat$titer)
dat <- dat[!is.na(dat$titer), ] #3070


#- set func to calculate GMT
dat_gmt <- dat %>% group_by(group, period, dose) %>%  summarise(gmt = exp(mean(log(titer))), n = n()) 
CI <- function(x){
  dat_gmt_ci <- c()
  for (h in unique(x$group)[1:length(unique(x$group))]) {
    tmp <- x  %>% filter(group == h)
  for (i in unique(tmp$period)[1:length(unique(tmp$period))]){
    tmp1 <- tmp %>% filter(period == i)
    for (j in unique(tmp1$dose)[1:length(unique(tmp1$dose))]){
      tmp2 <- tmp1 %>% filter(dose == j)  
      if(all(diff(tmp2$titer) == 0)){
        lci <- "-"
        uci <- "-"
      }else{
        lci <- exp(t.test(log(tmp2$titer))$conf[1])
        uci <- exp(t.test(log(tmp2$titer))$conf[2])
      }
      tmp3 <- data.frame(group = h, period = i, dose = j, lci = lci, uci = uci)
      tmp3$group <- as.character(tmp3$group)
      tmp3$period <- as.character(tmp3$period)
      tmp3$dose <- as.character(tmp3$dose)
      tmp3$lci <- as.character(tmp3$lci)
      tmp3$uci <- as.character(tmp3$uci)
      dat_gmt_ci <- rbind(dat_gmt_ci , tmp3)
     } 
    }
  }
  return(dat_gmt_ci)
}
dat_gmt_ci <- CI(dat)
dat_gmt <- left_join(dat_gmt, dat_gmt_ci)
dat_gmt$lci <- as.numeric(dat_gmt$lci)
dat_gmt$uci <- as.numeric(dat_gmt$uci)


#- calcluated gmt for younger group 
dat_gmt_2b <- dat_gmt %>% filter(group == "0_28_8m")
unique(dat_gmt_2b$period)
dat_gmt_2b$period <- factor(dat_gmt_2b$period, levels = c("Baseline", "Day 28 after dose 2","Day 180 after dose 2", 
                                                          "Day 28 after dose 3","Day 180 after dose 3"))
dat_gmt_2b$dose <- factor(dat_gmt_2b$dose, levels = c("placebo", "3ug", "6ug"))

dat_gmt_2b <- dat_gmt_2b[!is.na(dat_gmt_2b$period), ]
dat_gmt_2b$uci <- as.numeric(dat_gmt_2b$uci)
dat_gmt_2b$lci <- as.numeric(dat_gmt_2b$lci)
dat_gmt_2b$uci[is.na(dat_gmt_2b$uci)] <- dat_gmt_2b$gmt[is.na(dat_gmt_2b$uci)]
dat_gmt_2b$lci[is.na(dat_gmt_2b$lci)] <- dat_gmt_2b$gmt[is.na(dat_gmt_2b$lci)]
dat_gmt_2b <- dat_gmt_2b[order(dat_gmt_2b$period, dat_gmt_2b$dose), ]
dat_gmt_2b$x_new <- c(0.73, 1.0, 1.27, 1.73,2.0,2.27, 2.73,3.0,3.27, 3.73,4.0,4.27, 4.73,5.0,5.27)

dat_2b <- left_join(dat, dat_gmt_2b[, c("group", "period", "dose", "x_new")])
dat_2b$dose <- factor(dat_2b$dose, levels = c("placebo", "3ug", "6ug"))
dat_2b$x_new <- jitter(dat_2b$x_new, factor = 1.5)
dat_2b <- dat_2b %>% filter(group == "0_28_8m")

ggplot(data= dat_gmt_2b) +
  geom_jitter(data = dat_2b %>% filter(group == "0_28_8m" & period != "Day 14 after dose 2"),
              aes(x = x_new, y = log(titer, 2), color = dose), alpha = 0.4, show.legend = F) +

  geom_bar(aes(x = period, y = log(gmt,2), fill = dose), color = "black", stat = "identity", position = position_dodge(0.8), width = 0.7,size= 0.75,alpha = 0.5)+
  geom_errorbar(aes(x = period, ymin = log(lci,2), ymax = log(uci,2), group = dose), position = position_dodge(0.8), width = 0.2, size= 0.75) +
  scale_y_continuous(name = "Neutralization titer", expand = c(0,0), breaks = seq(1, 13, 2), labels = 2^seq(1, 13, 2), limits = c(-15, max(seq(1, 13, 2))+5))+
  scale_x_discrete(name = "", expand = c(0,0.8), labels = c("Baseline\n", "Day 28 after \ndose 2","Day 180 after \ndose 2", 
                                                            "Day 28 after \ndose 3","Day 180 after \ndose 3"))+
  scale_fill_manual(values = c("#E1ADA7","#66CC99",  "#1f78b4"), labels = c("Placebo", "3 μg", "6 μg")) +
  coord_cartesian(ylim = c(0,max(seq(1, 13, 2))), clip = "off") +
  # geom_hline(yintercept = log(8,2), linetype = "dashed",  color = "black", size = 1) + 
  geom_hline(yintercept = log(33,2), linetype = "dashed",  color = "black", size = 0.75) + 
  annotate("text", x = dat_gmt_2b$x_new, y = log(dat_gmt_2b$uci,2)+0.3, label = sprintf("%.1f",dat_gmt_2b$gmt), size=4,angle = 0,adj=0.5,
           color="black", fontface = "bold")+
  annotate("text", x = 0.1, y = -0.6, label = "No.participants", size=4,angle = 0,adj=0.5,color="black", fontface = "bold")+
  annotate("rect", xmin =  seq(0.61,4.61,1), xmax = seq(0.84,4.84,1), ymin =rep(-0.8,5), ymax= rep(-0.4,5),
           fill = c("#E1ADA7"), alpha = 0.25, color = "black", size = 0.75) +
  annotate("rect", xmin =  seq(0.88,4.91,1), xmax = seq(1.11,5.11,1), ymin =rep(-0.8,5), ymax= rep(-0.4,5),
           fill = c("#66CC99"), alpha = 0.25, color = "black", size = 0.75) +
  annotate("rect", xmin =  seq(1.15,5.15,1), xmax = seq(1.38,5.38,1), ymin =rep(-0.8,5), ymax= rep(-0.4,5),
           fill = c("#1f78b4"), alpha = 0.25, color = "black", size = 0.75) +
  annotate("text", x = dat_gmt_2b$x_new, y = -0.6, label = dat_gmt_2b$n, size=4,angle = 0,adj=0.5,color="black", fontface = "bold") +
  # annotate("text", x = 0.6, y = 13.8, label = "Cohort 2b-28d-8m", size=5,angle = 0,adj=0.5,color="black", fontface = "bold") +
  
  annotate("segment", x = 2, xend = 2.27, y = 10, yend = 10, size=0.75, color="black")+
  annotate("segment", x = 2, xend = 2, y = 9.8, yend = 10, size=0.75, color="black")+
  annotate("segment", x = 2.27, xend = 2.27, y = 9.8, yend = 10, size=0.75, color="black")+
  annotate("segment", x = 3, xend = 3.27, y = 10, yend = 10, size=0.75, color="black")+
  annotate("segment", x = 3, xend = 3, y = 9.8, yend = 10, size=0.75, color="black")+
  annotate("segment", x = 3.27, xend = 3.27, y = 9.8, yend = 10, size=0.75, color="black")+
  annotate("segment", x = 4, xend = 4.27, y = 10, yend = 10, size=0.75, color="black")+
  annotate("segment", x = 4, xend = 4, y = 9.8, yend = 10, size=0.75, color="black")+
  annotate("segment", x = 4.27, xend = 4.27, y = 9.8, yend = 10, size=0.75, color="black")+
  annotate("segment", x = 5, xend = 5.27, y = 10, yend = 10, size=0.75, color="black")+
  annotate("segment", x = 5, xend = 5, y = 9.8, yend = 10, size=0.75, color="black")+
  annotate("segment", x = 5.27, xend = 5.27, y = 9.8, yend = 10, size=0.75, color="black")+


  annotate("text", x = (2+2.27)/2, y = 10+0.3, label = "p<0.01*", size=4,angle = 0,adj=0.5,
           color="black", fontface = "bold")+
  annotate("text", x = (3+3.27)/2, y = 10+0.3, label = "p=0.80", size=4,angle = 0,adj=0.5,
           color="black", fontface = "bold")+
  annotate("text", x = (4+4.27)/2, y = 10+0.3, label = "p=0.01*", size=4,angle = 0,adj=0.5,
           color="black", fontface = "bold")+
  annotate("text", x =(5+5.27)/2, y = 10+0.3, label = "p<0.01*", size=4,angle = 0,adj=0.5,
           color="black", fontface = "bold") +
  
  theme(legend.position = c(0.02, 0.85),
        legend.justification = c(0,0),
        plot.margin = margin(1.5,0.5,1.5,0.1,unit="cm"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.ticks = element_line(colour="black",size=0.4),
        axis.line = element_line(colour="black",size=1 ),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black",size=16,margin = margin(0, 4, 0, 0),face="bold", vjust = -1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.8, "lines"), # adjust spacing between two faceted plots
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.text.x  = element_text(color="black",size=12,margin = margin(4, 0, 0, 0), vjust = -10,face="bold"),
        axis.text.y  = element_text(color="black",size=12,margin = margin(0, 4, 0, 0),face="bold"),
        panel.grid = element_blank(),
        panel.border=element_blank(),
        strip.background = element_rect(colour = "black",fill = "#6996AA"),
        strip.text =  element_text(colour = "white",size=10,face="bold")) -> p1


# pdf("fig1.pdf", height = 7, width = 11)
# 
# p1
# 
# dev.off()






#### Elder ####
#- data cleaning
dat <- read.csv("cleaned_data/dat_titer_elder.csv")
colnames(dat)
dat_new <- dat[, c("SUBJID", "titer_1", "titer_2", "titer_3", "titer_4", "titer_5", "titer_6")]
dat_new <- melt(dat_new, id = c("SUBJID"), value.name = "titer")
colnames(dat_new)[2] <- "period"
dat <- left_join(dat[, 2:5], dat_new)

dat$dose[dat$ARMCD %in% c("A")] <- "1.5ug"
dat$dose[dat$ARMCD %in% c("B")] <- "3ug"
dat$dose[dat$ARMCD %in% c("C")] <- "6ug"
dat$dose[dat$ARMCD %in% c("D")] <- "placebo"

dat$group <- "0_28_8m"

dat$period <- as.character(dat$period)
dat$period[dat$period == "titer_1"] <- "Baseline"
dat$period[dat$period == "titer_2"] <- "Day 28 after dose 2"
dat$period[dat$period == "titer_3"] <- "Day 180 after dose 2"
dat$period[dat$period == "titer_4"] <- "Before dose 3"
dat$period[dat$period == "titer_5"] <- "Day 28 after dose 3"
dat$period[dat$period == "titer_6"] <- "Day 180 after dose 3"

dat$titer[dat$titer == "<4"] <- 2
dat$titer <- as.numeric(dat$titer)

dat <- dat[!is.na(dat$titer), ]

#- set func to calculate GMT
dat_gmt <- dat %>% group_by(group, period, dose) %>%  summarise(gmt = exp(mean(log(titer))), n = n()) 
CI <- function(x){
  dat_gmt_ci <- c()
  for (h in unique(x$group)[1:length(unique(x$group))]) {
    tmp <- x  %>% filter(group == h)
    for (i in unique(tmp$period)[1:length(unique(tmp$period))]){
      tmp1 <- tmp %>% filter(period == i)
      for (j in unique(tmp1$dose)[1:length(unique(tmp1$dose))]){
        tmp2 <- tmp1 %>% filter(dose == j)  
        if(all(diff(tmp2$titer) == 0)){
          lci <- "-"
          uci <- "-"
        }else{
          lci <- exp(t.test(log(tmp2$titer))$conf[1])
          uci <- exp(t.test(log(tmp2$titer))$conf[2])
        }
        tmp3 <- data.frame(group = h, period = i, dose = j, lci = lci, uci = uci)
        dat_gmt_ci <- rbind(dat_gmt_ci , tmp3)
      } 
    }
  }
  return(dat_gmt_ci)
}
dat_gmt_ci <- CI(dat)
dat_gmt <- left_join(dat_gmt, dat_gmt_ci)


#- calcluated gmt for younger group 
dat_gmt_3 <- dat_gmt %>% filter(group == "0_28_8m")
unique(dat_gmt_3$period)
dat_gmt_3$period <- factor(dat_gmt_3$period, levels = c("Baseline", "Day 28 after dose 2","Day 180 after dose 2",
                                                          "Day 28 after dose 3","Day 180 after dose 3"))
dat_gmt_3$dose <- factor(dat_gmt_3$dose, levels = c("placebo", "1.5ug","3ug", "6ug"))

dat_gmt_3 <- dat_gmt_3[!is.na(dat_gmt_3$period), ]
dat_gmt_3$uci <- as.numeric(dat_gmt_3$uci)
dat_gmt_3$lci <- as.numeric(dat_gmt_3$lci)
dat_gmt_3$uci[is.na(dat_gmt_3$uci)] <- dat_gmt_3$gmt[is.na(dat_gmt_3$uci)]
dat_gmt_3$lci[is.na(dat_gmt_3$lci)] <- dat_gmt_3$gmt[is.na(dat_gmt_3$lci)]
dat_gmt_3 <- dat_gmt_3[order(dat_gmt_3$period, dat_gmt_3$dose), ]
dat_gmt_3$x_new <- c(0.7, 0.9, 1.1, 1.3,  1.7, 1.9, 2.1, 2.3,  2.7, 2.9, 3.1, 3.3,   3.7, 3.9, 4.1, 4.3,  4.7, 4.9, 5.1, 5.3)

dat_3 <- left_join(dat, dat_gmt_3[, c("group", "period", "dose", "x_new")])
dat_3$dose <- factor(dat_3$dose, levels = c("placebo", "1.5ug","3ug", "6ug"))
dat_3$x_new <- jitter(dat_3$x_new, factor = 1.2)
dat_3 <- dat_3[!is.na(dat_3$x_new), ]


ggplot(data= dat_gmt_3) +
  geom_jitter(data = dat_3 %>% filter(group == "0_28_8m" ),aes(x = x_new, y = log(titer, 2), color = dose), alpha = 0.4, show.legend = F) +
  geom_bar(aes(x = period, y = log(gmt,2), fill = dose), color = "black", stat = "identity", position = position_dodge(0.8),
           width = 0.7,size= 0.75,alpha = 0.5) +
  geom_hline(yintercept = log(33,2), linetype = "dashed",  color = "black", size = 0.75) + 
  geom_errorbar(aes(x = period, ymin = log(lci,2), ymax = log(uci,2), group = dose), position = position_dodge(0.8), width = 0.2, size= 0.75) +
  scale_y_continuous(name = "Neutralization titer", expand = c(0,0), breaks = seq(1, 13, 2), labels = 2^seq(1, 13, 2), limits = c(-15, max(seq(1, 13, 2))+5))+
  scale_x_discrete(name = "", expand = c(0,0.8), labels = c("Baseline\n", "Day 28 after \ndose 2","Day 180 after \ndose 2", 
                                                            "Day 28 after \ndose 3","Day 180 after \ndose 3"))+
  scale_fill_manual(values = c("#E1ADA7", "#FFD966", "#66CC99",  "#1f78b4"), labels = c("Placebo", "1.5μg","3μg", "6μg")) +
  scale_color_manual(values = c("#E1ADA7", "#FFD966", "#66CC99",  "#1f78b4"))+
  coord_cartesian(ylim = c(0,max(seq(1, 13, 2))), clip = "off") +
  # geom_hline(yintercept = log(8,2), linetype = "dashed",  color = "black", size = 1) + 
  
  annotate("text", x = dat_gmt_3$x_new, y = log(dat_gmt_3$uci,2)+0.3, label = sprintf("%.1f",dat_gmt_3$gmt), size=4,angle = 0,adj=0.5,
           color="black", fontface = "bold")+
  annotate("text", x = 0.1, y = -0.6, label = "No.participants", size=4,angle = 0,adj=0.5,color="black", fontface = "bold")+
  annotate("text", x = dat_gmt_3$x_new, y = -0.6, label = dat_gmt_3$n, size=4,angle = 0,adj=0.5,color="black", fontface = "bold") +
  # annotate("text", x = 0.6, y = 13.8, label = "Cohort 3-28d-8m", size=5,angle = 0,adj=0.5,color="black", fontface = "bold") +
  annotate("rect", xmin =  seq(0.62,4.62,1), xmax = seq(0.79,4.79,1), ymin =rep(-0.8,5), ymax= rep(-0.4,5),
           fill = c("#E1ADA7"), alpha = 0.25, color = "black", size = 0.75) +
  annotate("rect", xmin =  seq(0.82,4.82,1), xmax = seq(0.99,4.99,1), ymin =rep(-0.8,5), ymax= rep(-0.4,5),
           fill = c("#FFD966"), alpha = 0.25, color = "black", size = 0.75) +
  annotate("rect", xmin =  seq(1.02,5.02,1), xmax = seq(1.19,5.19,1), ymin =rep(-0.8,5), ymax= rep(-0.4,5),
           fill = c("#66CC99"), alpha = 0.25, color = "black", size = 0.75) +
  annotate("rect", xmin =  seq(1.22,5.22,1), xmax = seq(1.39,5.39,1), ymin =rep(-0.8,5), ymax= rep(-0.4,5),
           fill = c("#1f78b4"), alpha = 0.25, color = "black", size = 0.75) +
  
  
  annotate("segment", x = 1.85, xend = 2.08, y = 10, yend = 10, size=0.75, color="black")+
  annotate("segment", x = 1.85, xend = 1.85, y = 9.8, yend = 10, size=0.75, color="black")+
  annotate("segment", x = 2.08, xend = 2.08, y = 9.8, yend = 10, size=0.75, color="black")+

  annotate("segment", x = 1.85, xend = 2.3, y = 10.7, yend = 10.7, size=0.75, color="black")+
  annotate("segment", x = 1.85, xend = 1.85, y = 10.5, yend = 10.7, size=0.75, color="black")+
  annotate("segment", x = 2.3, xend = 2.3, y = 10.5, yend = 10.7, size=0.75, color="black")+


  annotate("segment", x = 4.85, xend = 5.08, y = 10, yend = 10, size=0.75, color="black")+
  annotate("segment", x = 4.85, xend = 4.85, y = 9.8, yend = 10, size=0.75, color="black")+
  annotate("segment", x = 5.08, xend = 5.08, y = 9.8, yend = 10, size=0.75, color="black")+

  annotate("segment", x = 5.10, xend = 5.35, y = 10, yend = 10, size=0.75, color="black")+
  annotate("segment", x = 5.35, xend = 5.35, y = 9.8, yend = 10, size=0.75, color="black")+
  annotate("segment", x = 5.10, xend = 5.10, y = 9.8, yend = 10, size=0.75, color="black")+

  annotate("segment", x = 4.85, xend = 5.35, y = 10.7, yend = 10.7, size=0.75, color="black")+
  annotate("segment", x = 4.85, xend = 4.85, y = 10.5, yend = 10.7, size=0.75, color="black")+
  annotate("segment", x = 5.35, xend = 5.35, y = 10.5, yend = 10.7, size=0.75, color="black")+


  annotate("text", x = (1.85+2.08)/2, y = 10.1+0.2, label = "p<0.01", size=4,angle = 0,adj=0.5,
           color="black", fontface = "bold")+
  annotate("text", x = (1.85+2.3)/2, y = 10.7+0.2, label = "p<0.01", size=4,angle = 0,adj=0.5,
           color="black", fontface = "bold")+

  annotate("text", x = (4.85+5.08)/2, y = 10.1+0.2, label = "p<0.01", size=4,angle = 0,adj=0.5,
           color="black", fontface = "bold")+
  annotate("text", x = (5.1+5.35)/2, y = 10.1+0.2, label = "p<0.01", size=4,angle = 0,adj=0.5,
           color="black", fontface = "bold")+
  annotate("text", x = (4.85+5.35)/2, y = 10.7+0.2, label = "p<0.01", size=4,angle = 0,adj=0.5,
           color="black", fontface = "bold")+
  
 
  theme(legend.position = "none",
        legend.justification = c(0,0),
        plot.margin = margin(1.5,0.5,1.5,0.1,unit="cm"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.ticks = element_line(colour="black",size=0.4),
        axis.line = element_line(colour="black",size=1 ),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black",size=16,margin = margin(0, 4, 0, 0),face="bold", vjust = -1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.8, "lines"), # adjust spacing between two faceted plots
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.text.x  = element_text(color="black",size=12,margin = margin(4, 0, 0, 0), vjust = -11,face="bold"),
        axis.text.y  = element_text(color="black",size=12,margin = margin(0, 4, 0, 0),face="bold"),
        panel.grid = element_blank(),
        panel.border=element_blank(),
        strip.background = element_rect(colour = "black",fill = "#6996AA"),
        strip.text =  element_text(colour = "white",size=10,face="bold")) -> p2


# pdf("fig2.pdf", height = 7, width = 11)
# p2
# dev.off()




