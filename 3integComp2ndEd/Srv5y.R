# Description: Plot a prognostic curve under the dynamic progression of the
# disease in reference patients

library(dplyr)
library(mstate)
library(ggplot2)
library(patchwork)

msT2dComp2ndEd <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\msT2dComp2ndEd.rds" %>%
  readRDS(.)
cfull2ndEd <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\cfull2ndEd.rds" %>%
  readRDS(.)
tmat2ndEd <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\tmat2ndEd.rds" %>%
  readRDS(.)

msT2dComp <- msT2dComp2ndEd
cfull <- cfull2ndEd
tmat <- tmat2ndEd

# Quickly get the structure of an object
pat <- msT2dComp[1:5, 9:28]
# We are interested in the effects of different levels of metabolites and proteins,
# so we control the two variables of drinking and insomnia here.
pat[, "Alcohol.drinker.status...Instance.0"] <- msT2dComp[1, "Alcohol.drinker.status...Instance.0"]
pat[, "Sleeplessness...insomnia...Instance.0"] <-  msT2dComp[1, "Sleeplessness...insomnia...Instance.0"]

# Patient Male And High HR Q1,T2D at day 60 and complicationat day 180
patQ1 <- pat
for (i in 1:(ncol(patQ1)-2)) {
  patQ1[, i] <- sort(unique(msT2dComp[, 9]))[1]
}

patQ2 <- pat
for (i in 1:(ncol(patQ2)-2)) {
  patQ2[, i] <- sort(unique(msT2dComp[, 9]))[2]
}

patQ3 <- pat
for (i in 1:(ncol(patQ3)-2)) {
  patQ3[, i] <- sort(unique(msT2dComp[, 9]))[3]
}

patQ4 <- pat
for (i in 1:(ncol(patQ4)-2)) {
  patQ4[, i] <- sort(unique(msT2dComp[, 9]))[4]
}
  
  patList <- list(patQ1, patQ2, patQ3, patQ4)
  
  for (i in 1:length(patList)) {
    patTemp <- patList[[i]]
    patTemp$trans <- 1:5
    attr(patTemp, "trans") <- tmat
    covs <- colnames(pat)
    patTemp <- expand.covs(patTemp, covs, longnames = FALSE)
    patTemp$strata <- patTemp$trans
    msTemp <- msfit(cfull, patTemp, trans = tmat)
    
    # Predicted probability of 5-year survival
    ptTemp.5yrs <- probtrans(msTemp, predt=5, direction="fixedhorizon")
    
    pt1b <- ptTemp.5yrs[[1]]
    pt2b <- ptTemp.5yrs[[2]]
    pt3b <- ptTemp.5yrs[[3]]
    ind1 <- which(pt1b$time < 60/365.25)
    ind2 <- which(pt1b$time >= 60/365.25 & pt1b$time < 180/365.25)
    ind3 <- which(pt1b$time >= 180/365.25)
    
    tttime <- pt1b$time
    
    df <- data.frame(
      time = c(
        tttime[ind1],
        tttime[-ind1],
        tttime[ind2],
        tttime[ind3],
        tttime[ind3]
      ),
      survival_prob = c(
        1-pt1b$pstate4[ind1],
        1-pt1b$pstate4[-ind1],
        1-pt2b$pstate4[ind2],
        1-pt2b$pstate4[ind3],
        1-pt3b$pstate4[ind3]
      ),
      group = rep(
        c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5"), 
        sapply(
          list(tttime[ind1], tttime[-ind1], tttime[ind2], tttime[ind3], tttime[ind3]),
          length
        )
      )
    )
    
    plot5yTemp <- ggplot(df, aes(x = time, y = survival_prob, color = group, linetype = group)) +
      geom_line(linewidth=1.2) +
      ylim(0.6, 1) +
      scale_color_manual(
        values = c(
          "Group 1" = "green",
          "Group 2" = "green",
          "Group 3" = "#c06c84",
          "Group 4" = "#c06c84",
          "Group 5" = "#e84545"
        ),
        labels = c(
          "Group 1" = "Baseline",
          "Group 2" = "Baseline without transition",
          "Group 3" = "T2D",
          "Group 4" = "T2D without transition",
          "Group 5" = "T2D + Complication")
      ) +
      scale_linetype_manual(
        values = c(
          "Group 1" = "solid",
          "Group 2" = "dashed",
          "Group 3" = "solid",
          "Group 4" = "dashed",
          "Group 5" = "solid"
        ),labels = c(
          "Group 1" = "Baseline",
          "Group 2" = "Baseline without transition",
          "Group 3" = "T2D",
          "Group 4" = "T2D without transition",
          "Group 5" = "T2D + Complication")
      ) +
      annotate("segment", 
               x = tttime[max(ind1)], xend = tttime[max(ind1)], 
               y = 1 - pt1b$pstate4[max(ind1)], yend = 1 - pt2b$pstate4[max(ind1)], 
               color = "#c06c84", linetype = "dotted", size = 1) +
      annotate("segment", 
               x = tttime[max(ind2)], xend = tttime[max(ind2)], 
               y = 1-pt2b$pstate4[max(ind2)], yend = 1-pt3b$pstate4[max(ind2)], 
               color = "#e84545", linetype = "dotted", size = 1) +
      labs(x = "Years since baseline(years)", y = "Probability of 5-year survival", color = "", linetype = "") +
      ggtitle(paste0("Patient with Q", i, " levels")) +
      guides(
        color = guide_legend(keywidth = 2),
        linetype = guide_legend(keywidth = 2)
      ) +
      theme(
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "#dbe2ef", linewidth = 0.5),
        panel.grid.minor = element_line(color = "#dbe2ef", linewidth = 0.5),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      )
    
    print(plot5yTemp)
    assign(paste0("Srv5yQ", i), plot5yTemp)
  }

  finalProb <- (Srv5yQ1 + Srv5yQ2 + Srv5yQ3 + Srv5yQ4) +
    plot_layout(ncol = 4, byrow = TRUE, widths = c(1, 1), guides = 'collect')
  print(finalProb)
  