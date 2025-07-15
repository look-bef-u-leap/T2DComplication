# Description: This script is used to plot the incidence probability of a patient.
# In the selection of patients, we selected four levels effective value for comparison

library(dplyr)
library(mstate)
library(patchwork)
library(ggplot2)

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

for (i in 1:(ncol(pat)-2)) {
  patQ1 <- pat
  for (j in 1:(ncol(patQ1)-2)) {
    patQ1[, j] <- sort(unique(msT2dComp[, 9]))[1]
  }
  
  patQ2 <- patQ1
  patQ2[, i] <- sort(unique(msT2dComp[, 9]))[2]
  
  patQ3 <- patQ1
  patQ3[, i] <- sort(unique(msT2dComp[, 9]))[3]
  
  patQ4 <- patQ1
  patQ4[, i] <- sort(unique(msT2dComp[, 9]))[4]
  
  patList <- list(patQ1, patQ2, patQ3, patQ4)
  for (n in 1:length(patList)) {
    patTemp <- patList[[n]]
    patTemp$trans <- 1:5
    attr(patTemp, "trans") <- tmat
    covs <- colnames(pat)
    patTemp <- expand.covs(patTemp, covs, longnames = FALSE)
    patTemp$strata <- patTemp$trans
    msTemp <- msfit(cfull, patTemp, trans = tmat)
    
    { # ~ plot for cumulative hazard -------
      cumHazTemp <- plot(
        msTemp,
        xlab = "Years since baseline(years)",
        ylab = "Cumulative hazard",
        ylim = c(0, 2),
        xlim = c(0, 15),
        conf.type = "log", 
        conf.int = 0.95,
        cols = c("#FF6600", "#8e063b", "#1BC727", "#1F91DC", "#032263"),
        legend = c("Healthy to T2D", "Healthy to Death", "T2d to Complication", "T2D to Death", "Complication to Death"),
        use.ggplot = TRUE
      )+
        ggtitle(paste0("Patient with Q", n, " levels")) +
        theme(
          plot.title = element_text(size = 13, face = "bold", hjust = 0.5, color = "black"),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "#dbe2ef", linewidth = 0.5),
          panel.grid.minor = element_line(color = "#dbe2ef", linewidth = 0.5),
          legend.key = element_rect(fill = "white", color = "white")
        )
      
      print(cumHazTemp)
      assign(paste0("cumHazQ", n), cumHazTemp)
    }
    
    {#~~ Risk prediction and mapping after 0 days -----
      ptTemp <- probtrans(msTemp, predt = 0)
      
      par(mfrow = c(1, 3))
      
      dat_plot <- plot(x = ptTemp, use.ggplot = TRUE, type = "single")$data
      probTemp <- ggplot(
        data = dat_plot,
        aes(
          x = time, 
          y = prob, 
          ymin = CI_low, 
          ymax = CI_upp,
          group = state,
          col = state
        )
      ) +
        geom_ribbon(col = NA, fill = "gray", alpha = 0.5) +
        geom_line(linewidth = 1.2) +
        scale_color_manual(
          values = c(
            "srv" = "#522546",
            "Comp" = "#e84545",
            "E11" = "#c06c84",
            "baseline" = "green"
          ),
          labels = c(
            "srv" = "Death",
            "Comp" = "Complication",
            "E11" = "T2D",
            "baseline" = "Healthy"
          )
        ) +
        labs(x = "Years since baseline(years)", y = "Probability", color = "State") +
        coord_cartesian(ylim = c(0, 1), xlim = c(0, 16), expand = 0) +
        ggtitle(paste0("Patient with Q", n, " levels")) +
        theme(
          plot.title = element_text(size = 13, face = "bold", hjust = 0.5, color = "black"),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "#dbe2ef", linewidth = 0.5),
          panel.grid.minor = element_line(color = "#dbe2ef", linewidth = 0.5),
          legend.key = element_rect(fill = "white", color = "white")
        )
      
      print(probTemp)
      assign(paste0("probQ", n), probTemp)
    }
  }
  
  finalProb <- (probQ1 + probQ2 + probQ3 + probQ4) +
    plot_layout(ncol = 4, byrow = TRUE, widths = c(1, 1), guides = 'collect')
  
  pdf(
    paste0("D:\\work2\\t2dAndComplicationsUkb\\plotDemo\\integComp2ndEd\\state\\", colnames(pat)[i], "_state.pdf"),
    width = 10,
    height = 3,
    compress = TRUE,
    onefile = TRUE,
    paper = "special"
  )
  print(finalProb)
  dev.off()
  gc()
  
  finalCumHaz <- (cumHazQ1 + cumHazQ2 + cumHazQ3 + cumHazQ4) +
    plot_layout(ncol = 4, byrow = TRUE, widths = c(1, 1), guides = 'collect')
  
  pdf(
    paste0("D:\\work2\\t2dAndComplicationsUkb\\plotDemo\\integComp2ndEd\\transition\\", colnames(pat)[i], "_transition.pdf"),
    width = 10,
    height = 3,
    compress = TRUE,
    onefile = TRUE,
    paper = "special"
  )
  print(finalCumHaz)
  dev.off()
  gc()
}
