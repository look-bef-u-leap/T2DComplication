# Description: This script was used to analyze the risk factors for integCompAdd
# complications in patients with t2d

library(dplyr)
library(ranger)
library(survival)
library(car)

integComp2ndEd <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\integComp2ndEd.rds" %>%
  readRDS(.)
trait2 <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\preprocess\\trait2.rds" %>%
  readRDS(.)
PIDgroup <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\grouping\\PIDgroup.rds" %>%
  readRDS(.)

integComp2ndEd <- merge(integComp2ndEd, PIDgroup[PIDgroup$group=="train_set", ], by.x = "Participant.ID", by.y = "PID")
integComp2ndEd <- integComp2ndEd[, -8]

trait3 <- trait2 %>%
  mutate(
    .,
    Ethnic.Continent = case_when(
      `Ethnic.background...Instance.0` %in% c("Irish", "British", "White",
                                              "Any other white background") ~ "Europe",
      `Ethnic.background...Instance.0` %in% c("Caribbean", "African",
                                              "Black or Black British",
                                              "Any other Black background") ~ "Africa",
      `Ethnic.background...Instance.0` %in% c("Chinese", "Pakistani",
                                              "Bangladeshi", "Indian",
                                              "Asian or Asian British",
                                              "Any other Asian background") ~ "Asia",
      `Ethnic.background...Instance.0` %in% c("White and Black African", "Mixed",
                                              "White and Asian","White and Black Caribbean",
                                              "Any other mixed background") ~ "Mixed",
      TRUE ~ "Other"
    ))
trait3$Ethnic.Continent <- as.factor(trait3$Ethnic.Continent)

trait4 <- trait3[, c(7:(ncol(trait3)-1))] %>%
  as.matrix(.) %>%
  scale(., center = TRUE, scale = TRUE)
trait5 <- cbind(trait3[, c(1, 2, ncol(trait3), 4, 5, 6)], trait4)

{ # ~ Rf and vif were used to calculate covariates important for complications ------
  comp <- integComp2ndEd
  colnames(comp)[4:5] <- c("comp", "comp.s")
  #The blood collection time is baseline, and the diagnosis time of complications
  #is after baseline. We hope to use blood components to reflect the progress of
  #complications, so the course of complications is calculated from baseline
  #comp$t2dToComp <- (comp$comp - comp$E11)
  comp$t2dToComp <- comp$comp
  comp2 <- subset(comp, E11.s == 1)
  comp3 <- comp2[, c("Participant.ID", "t2dToComp", "comp.s")]
  
  compTra <- merge(comp3, trait5, by = "Participant.ID")
  colnames(compTra)[2:3] <- c("time", "status")
  
  set.seed(123)
  rfModelComp <- ranger(Surv(time, status) ~ ., data = compTra[, c(-1)], num.trees = 50, importance = "permutation")
  imporScoresComp <- importance(rfModelComp)
  imporScoresComp <- as.data.frame(imporScoresComp)
  imporScoresComp$name <- rownames(imporScoresComp)
  
  imporScoresComp2 <- imporScoresComp[which(imporScoresComp$imporScores > 0), ]
  thrComp <- order(imporScoresComp2$imporScores, decreasing = TRUE)[1:10]
  compTra4 <- compTra[, c(1:3, thrComp+3)]
  
  compTraVif <- coxph(Surv(time, status) ~ ., compTra4[, c(-1)]) %>% vif(.)
  compTraVifDf <- as.data.frame(compTraVif)
  # Gets the column position in trait2 for variables with vif less than 4
  indComp <- which(compTraVifDf[, 1] < 4) + 3
  
  compTra5 <- compTra4[, c(1:3, indComp)]
}

{ # ~ Rf and vif were used to calculate covariates important for t2d ------
  t2d <- integComp2ndEd
  
  t2dTra <- merge(t2d[, c(1:3)], trait5, by = "Participant.ID")
  colnames(t2dTra)[2:3] <- c("time", "status")
  
  set.seed(123)
  rfModelT2d <- ranger(Surv(time, status) ~ ., data = t2dTra[, c(-1)], num.trees = 50, importance = "permutation")
  imporScoresT2d <- importance(rfModelT2d)
  imporScoresT2d <- as.data.frame(imporScoresT2d)
  imporScoresT2d$name <- rownames(imporScoresT2d)
  
  imporScoresT2d2 <- imporScoresT2d[which(imporScoresT2d$imporScores > 0), ]
  thrT2d <- order(imporScoresT2d2$imporScores, decreasing = TRUE)[1:10]
  t2dTra4 <- t2dTra[, c(1:3, thrT2d+3)]
  
  t2dTraVif <- coxph(Surv(time, status) ~ ., t2dTra4[, c(-1)]) %>% vif(.)
  t2dTraVifDf <- as.data.frame(t2dTraVif)
  # Gets the column position in trait2 for variables with vif less than 4
  indT2d <- which(t2dTraVifDf[, 1] < 4) + 3
  
  t2dTra5 <- t2dTra4[, c(1:3, indT2d)]
}

# Removetime time and status
varFiltrated <- c(colnames(t2dTra5), colnames(compTra5)) %>% unique(.)
ind2 <- match(varFiltrated[c(-2, -3)], colnames(trait5))
traRfVif2ndEd <- trait5[, ind2]

saveRDS(
  traRfVif2ndEd,
  file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\traRfVif2ndEd.rds"
)
write.csv(
  traRfVif2ndEd, 
  file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\traRfVif2ndEd.csv",
  row.names = TRUE
)

saveRDS(
  trait5,
  file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\trait5_2ndEd.rds"
)
write.csv(
  trait5, 
  file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\trait5_2ndEd.csv",
  row.names = TRUE
)

saveRDS(
  rfModelComp,
  file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\rfModelComp2ndEd.rds"
)
saveRDS(
  compTraVif,
  file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\compTraVif2ndEd.rds"
)

saveRDS(
  rfModelT2d,
  file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\rfModelT2d2ndEd.rds"
)
saveRDS(
  t2dTraVif,
  file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\t2dTraVif2ndEd.rds"
)
