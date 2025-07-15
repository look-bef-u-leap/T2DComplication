library(dplyr)
library(survival)
library(car)
library(glmnet)
library(stringr)

all_file_paths <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\LassoCox" %>%
  list.files(path = ., pattern = ".*", full.names = TRUE, recursive = TRUE)
compIndex <- which(str_detect(all_file_paths, "Comp.rds"))
comp_path <- all_file_paths[compIndex]

PIDgroup <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\grouping\\PIDgroup.rds" %>%
  readRDS(.)

trait3 <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\preprocess\\trait2_abbr_ethnic.rds" %>%
  readRDS(.)

trait3 <- merge(PIDgroup[PIDgroup$group=="train_set", ], trait3, by= "PID")
trait3 <- trait3[, -2]

for (i in 1:length(comp_path)) {
  comp <- comp_path[i] %>% readRDS(.)
  colnames(comp)[4:5] <- c("comp", "comp.s")
  # Because the covariates were collected at baseline, the course of complications
  # was calculated from baseline
  comp3 <- comp[, c("Participant.ID", "comp", "comp.s")]
  colnames(comp3) <- c("Participant.ID", "bslToComp", "comp.s")
  
  compTra <- merge(comp3, trait3, by.x = "Participant.ID", by.y = "PID")
  colnames(compTra)[2:3] <- c("time", "status")
  
  resBFpDf <- paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\univCox2ndEd\\resBFpDf_", str_sub(comp_path[i], 57, 59), "Comp.rds") %>%
    readRDS(.)
  resBFpDf2 <- resBFpDf[which(resBFpDf$BFp < 0.05), ]
  
  compTra3 <- compTra[, c("Participant.ID", "time", "status", rownames(resBFpDf2))]
  compTra3[, -c(1:3)] <- scale(compTra3[, -c(1:3)])
  
  { # lasso ------
    compTra3$time <- as.double(compTra3$time)
    compTra3$status <- as.double(compTra3$status)
    
    outcome <- data.matrix(Surv(time = compTra3$time, event = compTra3$status))
    tra <- compTra3[, -c(1:3)]
    tra2 <- model.matrix(~ .,tra)[,-1] %>%
      as.matrix(.)
    
    set.seed(1234)
    fitcv <- cv.glmnet(
      tra2,
      outcome,
      family="cox",
      alpha=1,
      nfolds=10,
      type.measure = "deviance",
      maxit= 10000000
    )
    
    coefficient <- coef(fitcv,s= fitcv$lambda.min)
    Active.Index <- which(as.numeric(coefficient)!=0)
    active.coefficients <- as.numeric(coefficient)[Active.Index]
    sig_multi_cox2 <- rownames(coefficient)[Active.Index]
  }
  
  { # cox ------
    coxTemp2 <- coxph(
      formula = Surv(time, status) ~ .,
      data = compTra3[, c("time", "status",sig_multi_cox2)]
    )
    
    cox_summ <- summary(coxTemp2)
  }
  
  saveRDS(
    fitcv,
    file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\LassoCox\\", str_sub(comp_path[i], 57, 59), "fitcv.rds")
  )
  saveRDS(
    coxTemp2,
    file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\LassoCox\\", str_sub(comp_path[i], 57, 59), "CoxTemp2.rds")
  )
  saveRDS(
    cox_summ,
    file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\LassoCox\\", str_sub(comp_path[i], 57, 59), "_cox_summ.rds")
  )
}