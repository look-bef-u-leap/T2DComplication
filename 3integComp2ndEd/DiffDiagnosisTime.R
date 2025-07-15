# Description: This script is used to plot the cumulative probability of Patient
# with Q4 levela confirming T2D at four different time points after baseline

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

{#~ Create Patient with Q4 levels ----
  patQ4 <- pat
  for (i in 1:(ncol(patQ4)-2)) {
    patQ4[, i] <- sort(unique(msT2dComp[, 9]))[4]
  }
  
  patQ4$trans <- 1:5
  attr(patQ4, "trans") <- tmat
  covs <- colnames(pat)
  patQ4 <- expand.covs(patQ4, covs, longnames = FALSE)
  patQ4$strata <- patQ4$trans
  msQ4 <- msfit(cfull, patQ4, trans = tmat)
}

n <- c(180, 360, 360*3, 360*5)
for (i in 1:length(n)) {
  ptQ4Temp <- probtrans(msQ4, predt = n[i]/365.25)
  dat_plot <- plot(x = ptQ4Temp, use.ggplot = TRUE, type = "single", from = 2)$data # State 2 indicates T2D
  probQ4Temp <- ggplot(
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
        "E11" = "#c06c84"
      ),
      labels = c(
        "srv" = "Death",
        "Comp" = "Complication",
        "E11" = "T2D"
      )
    ) +
    labs(x = "Years since baseline(years)", y = "Probability", color = "State") +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 16), expand = 0) +
    ggtitle(paste0("T2D diagnosed in ", n[i], " days")) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, color = "black"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "#dbe2ef", linewidth = 0.5),
      panel.grid.minor = element_line(color = "#dbe2ef", linewidth = 0.5),
      legend.key = element_rect(fill = "white", color = "white")
    )
  
  print(probQ4Temp)
  assign(paste0("probQ4d", n[i]), probQ4Temp)
}
  
finalDiffTime <- (probQ4d180 + probQ4d360 + probQ4d1080 + probQ4d1800) +
  plot_layout(ncol = 4,guides = 'collect')
print(finalDiffTime)
