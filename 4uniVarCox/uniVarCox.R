library(dplyr)
library(survival)
library(stringr)

all_file_paths <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\vifLassoCox" %>%
  list.files(path = ., pattern = ".*", full.names = TRUE, recursive = TRUE)
compIndex <- which(str_detect(all_file_paths, "Comp.rds"))
comp_path <- all_file_paths[compIndex]

PIDgroup <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\grouping\\PIDgroup.rds" %>%
  readRDS(.)

trait3 <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\preprocess\\trait2_abbr_ethnic.rds" %>%
  readRDS(.)

trait3 <- merge(PIDgroup[PIDgroup$group=="train_set", ], trait3, by= "PID")
trait3 <- trait3[, -2]

adjVar <- c(
  "Sex","Ethnic.Continent","CTS","ADS","Insomnia",
  "Glc.BBC","Glc.NMR","HbA1c","Chol","LDL.Direct",
  "TG","HDL.C.BBC","HDL.C.NMR","ApoB.BBC","ApoB.NMR",
  "AAAC","BMI","BMR","BFP"
)

for (i in 1:length(comp_path)) {
  comp <- comp_path[i] %>% readRDS(.)
  colnames(comp)[4:5] <- c("comp", "comp.s")
  # Because the covariates were collected at baseline, the course of complications
  # was calculated from baseline
  comp3 <- comp[, c("Participant.ID", "comp", "comp.s")]
  colnames(comp3) <- c("Participant.ID", "bslToComp", "comp.s")
  
  compTra <- merge(comp3, trait3, by.x = "Participant.ID", by.y = "PID")
  
  colnames(compTra)[2:3] <- c("time", "status")
  
  m1 <- match(adjVar, colnames(compTra))
  covariates <- colnames(compTra)[-c(1:3, m1)]
  formHead <- paste0(adjVar, " + ") %>%
    paste(., collapse = "") %>%
    paste0("Surv(time, status) ~ ", .)
  univ_formulas <- sapply(
    covariates,
    function(x) as.formula(
      paste0(formHead, x)
      )
    )
  compTra[, -c(1:8)] <- scale(compTra[, -c(1:8)])
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = compTra)})
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           coef <- x[["coefficients"]][29, ] %>% t(.) %>% as.data.frame(.)
                           rownames(coef) <- rownames(x[["coefficients"]])[29]
                           return(coef)
                         })
  res <- dplyr::bind_rows(univ_results)
  res$BFp <- p.adjust(res$`Pr(>|z|)`, method = "bonferroni")
  
  saveRDS(
    univ_models,
    file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\univCox2ndEd\\univ_models_", str_sub(comp_path[i], 60, 66), ".rds")
  )
  saveRDS(
    res,
    file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\univCox2ndEd\\resBFpDf_", str_sub(comp_path[i], 60, 66), ".rds")
  )
  write.csv(
    res,
    file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\univCox2ndEd\\resBFpDf_", str_sub(comp_path[i], 60, 66), ".csv"),
    row.names = TRUE
  )
  
}