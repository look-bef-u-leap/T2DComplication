# Description: This script uses the mstate package to analyze the risk of 
# type 2 diabetes and overall complications

# Shared head ----
library(mstate)
library(dplyr)

integComp2ndEd <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\integComp2ndEd.rds" %>%
  readRDS(.)
traLab2ndEd <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\traLab2ndEd.rds" %>%
  readRDS(.)

PIDgroup <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\grouping\\PIDgroup.rds" %>%
  readRDS(.)
integComp2ndEd <- merge(integComp2ndEd, PIDgroup[PIDgroup$group=="train_set", ], by.x = "Participant.ID", by.y = "PID")
integComp2ndEd <- integComp2ndEd[, -8]

t2dCompAddTra <- merge(integComp2ndEd, traLab2ndEd, by = "Participant.ID")

tmat <- transMat(
  x = list(c(2, 4), c(3, 4), c(4), c()),
  names = c("baseline", "E11", "Comp", "srv")
)

msT2dComp <- msprep(
  data = t2dCompAddTra,
  trans = tmat,
  time = c(NA, "E11", "Comp", "srv"),
  status = c(NA, "E11.s", "Comp.s", "srv.s"),
  id = t2dCompAddTra$Participant.ID,
  keep = colnames(t2dCompAddTra)[8:ncol(t2dCompAddTra)]
)
n <- ncol(msT2dComp)
# Some of the study subjects were transferred from start time to end time,
# which is because there is not enough detailed survival data, but these study
# subjects do not meet the reality and model requirements, so supplement
equ <- which(msT2dComp$Tstart == msT2dComp$Tstop)
msT2dComp[equ, "Tstop"] <- (msT2dComp[equ, "Tstop"]+1)
msT2dComp[, c("Tstart", "Tstop", "time")] <- msT2dComp[, c("Tstart", "Tstop", "time")] / 365.25

events(msT2dComp)
covs <- colnames(t2dCompAddTra)[8:ncol(t2dCompAddTra)]
msT2dComp <- expand.covs(msT2dComp, covs, longnames = FALSE)

# Semi-parametric model ----
dummyVar <- colnames(msT2dComp)[(n+1):ncol(msT2dComp)]
dummyVar2 <- paste0(" ", dummyVar, " +") %>% paste(., collapse ="")
dummyVarStr <- paste0("Surv(Tstart, Tstop, status) ~ ",  dummyVar2, " strata(trans)")

options(expressions = 500000)
cfull <- coxph(
  as.formula(dummyVarStr),
  data = msT2dComp,
  method = "breslow"
)
options(expressions = 5000)

{ # ~ Back up key objects cfull ------
  saveRDS(
    cfull,
    file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\cfull2ndEd.rds"
  )
  saveRDS(
    msT2dComp,
    file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\msT2dComp2ndEd.rds"
  )
  saveRDS(
    tmat,
    file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\tmat2ndEd.rds"
  )
}
