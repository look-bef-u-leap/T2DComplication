library(dplyr)
library(openxlsx)
library(stringr)

t2dCompMststy2ndEd <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\t2dCompMststy2ndEd.rds" %>%
  readRDS(.)
diseaseNameVsIcd10_2nd_ed <- "D:\\work2\\t2dAndComplicationsUkb\\data\\vs\\diseaseNameVsIcd10_2nd_ed.xlsx" %>% 
  read.xlsx(., sheet = 1)

t2dCompMststy <- t2dCompMststy2ndEd
colnames(t2dCompMststy)[38:43] <- c("Field42006", "Field42006.s", "Field42008", "Field42008.s", "Field42000", "Field42000.s")
diseaseNameVsIcd10 <- diseaseNameVsIcd10_2nd_ed

# Some of the study subjects had srv of 0 because there was no survival date
# data, but these subjects did not meet the reality and model requirements, so added
srv0 <- which(t2dCompMststy$srv == 0)
n00 <- which(colnames(t2dCompMststy) == "srv")
t2dCompMststy[srv0, seq(2, n00, by = 2)] <- 1

CompClass <- c(
  "microvascular",
  "macrovascular",
  "nervous system damage (neuropathy)",
  "renal system damage (nephropathy)",
  "eye damage(retinopathy)",
  "cardiovascular disease",
  "peripheral vascular disease",
  "metabolic disorder"
  )
for (i in 1:length(CompClass)) {
  t2dCompMststy$temp <- NA
  t2dCompMststy$temp.s <- NA
  
  if (i <= 2){
    SubsetTemp <- diseaseNameVsIcd10[which(diseaseNameVsIcd10$Category == CompClass[i]), ]
  }else{
    SubsetTemp <- diseaseNameVsIcd10[which(diseaseNameVsIcd10$System == CompClass[i]), ]
  }
  
  ind <- SubsetTemp[, c("ICD10", "OPSC4", "Other.Code")] %>%
    unlist(.) %>%
    na.omit(.) %>%
    as.vector(.) %>%
    match(., colnames(t2dCompMststy)) %>%
    na.omit(.) %>%
    as.vector(.)
  ind2 <- c(ind, ind + 1)
  ind3 <- ind2[order(ind2)]
  t2dCompMststy2 <- t2dCompMststy[, c(1, 2, 3, ind3, 136:139)]
  
  # Locating complication status
  m <- which(colnames(t2dCompMststy2) == "srv") - 1
  for (j in 1:nrow(t2dCompMststy2)) {
    if (sum(t2dCompMststy2[j, seq(5, m, by = 2)]) == 0) {
      t2dCompMststy2[j, "temp.s"] <- 0
      t2dCompMststy2[j, "temp"] <- t2dCompMststy2[j, "srv"]
    } else {
      t2dCompMststy2[j, "temp.s"] <- 1
      n <- which(t2dCompMststy2[j, seq(5, m, by = 2)] == 1)
      # This is the number of columns for status, with -1 indicating the location of the outcome
      n2 <- (seq(5, m, by = 2)[n] - 1)
      e <- order(as.matrix(t2dCompMststy2[j, n2]))[1]
      t2dCompMststy2[j, "temp"] <- t2dCompMststy2[j, n2][e]
    }
  }
  
  tempComp <- t2dCompMststy2[, c("Participant.ID", "E11", "E11.s", "temp", "temp.s", "srv", "srv.s")]
  colnames(tempComp)[c(4, 5)] <- c(str_sub(CompClass[i], 1, 3), paste0(str_sub(CompClass[i], 1, 3), ".s"))
  
  write.csv(
    tempComp,
    file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\lassoCox2ndEd\\", str_sub(CompClass[i], 1, 3), "Comp.csv"),
    row.names = TRUE
  )
  saveRDS(
    tempComp,
    file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\lassoCox2ndEd\\", str_sub(CompClass[i], 1, 3), "Comp.rds")
  )
}
