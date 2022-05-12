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





