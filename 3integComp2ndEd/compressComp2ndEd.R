# Description: This script is used to compress the time and status of all
# complications into a complication whole

library(dplyr)

t2dCompMststy2ndEd <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\t2dCompMststy2ndEd.rds" %>%
  readRDS(.)
t2dCompMststy <- t2dCompMststy2ndEd

# Some of the study subjects had srv of 0 because there was no survival date
# data, but these subjects did not meet the reality and model requirements, so added
srv0 <- which(t2dCompMststy$srv == 0)
n00 <- which(colnames(t2dCompMststy) == "srv")
t2dCompMststy[srv0, seq(2, n00, by = 2)] <- 1

t2dCompMststy$Comp <- NA
t2dCompMststy$Comp.s <- NA
for (i in 1:nrow(t2dCompMststy)) {
  if (sum(t2dCompMststy[i, seq(5, 135, by = 2)]) == 0) {
    t2dCompMststy[i, "Comp.s"] <- 0
    t2dCompMststy[i, "Comp"] <- t2dCompMststy[i, "srv"]
  } else {
    t2dCompMststy[i, "Comp.s"] <- 1
    n <- which(t2dCompMststy[i, seq(5, 135, by = 2)] == 1)
    # This is the number of columns for status, with -1 indicating the location of the outcome
    n2 <- (seq(5, 135, by = 2)[n] - 1)
    e <- order(as.matrix(t2dCompMststy[i, n2]))[1]
    t2dCompMststy[i, "Comp"] <- t2dCompMststy[i, n2][e]
  }
}

{ # Manually check the accuracy of the cycle processing pattern -----
  e11s <- c(0, 1)
  srvs <- c(0, 1)
  comps <- c(0, 1, 3)
  
  for (i in 1:length(e11s)) {
    for (j in 1:length(srvs)) {
      for (q in 1:length(comps)) {
        if (e11s[i] == 0 && comps[q] != 0) {
          r <-
            {
              (t2dCompMststy[, "E11.s"] == e11s[i]) &
                (t2dCompMststy[, "srv.s"] == srvs[j]) &
                (apply(t2dCompMststy[, seq(5, 135, by = 2)], MARGIN = 1, sum) == comps[q])
            } %>%
            which(.)
          print(paste0("nrow=", r, "---E11.s=", e11s[i], "---srv.s=", srvs[j], "---Comp.s=", comps[q]))
        } else {
          r <-
            {
              (t2dCompMststy[, "E11.s"] == e11s[i]) &
                (t2dCompMststy[, "srv.s"] == srvs[j]) &
                (apply(t2dCompMststy[, seq(5, 135, by = 2)], MARGIN = 1, sum) == comps[q])
            } %>%
            which(.) %>%
            sample(., 1)
          print(paste0("nrow=", r, "---E11.s=", e11s[i], "---srv.s=", srvs[j], "---Comp.s=", comps[q]))
          print(t2dCompMststy[r, ])
        }
      }
    }
  }
}

integComp2ndEd <- t2dCompMststy[, c("Participant.ID", "E11", "E11.s", "Comp", "Comp.s", "srv", "srv.s")]

write.csv(
  integComp2ndEd,
  file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\integComp2ndEd.csv",
  row.names = TRUE
)
saveRDS(
  integComp2ndEd,
  file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\integComp2ndEd.rds"
)
