library(dplyr)
library(pROC)
library(ROCR)
library(ggplot2)
library(stringr)
library(paletteer)
library(patchwork)

all_file_paths <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\lgbmNonMlr3_2ndEd" %>%
  list.files(path = ., pattern = ".*", full.names = TRUE, recursive = TRUE)
aucPIndex <- which(str_detect(all_file_paths, "_aucP.rds"))
aucP_path <- all_file_paths[aucPIndex]

PIDgroup <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\grouping\\PIDgroup.rds" %>%
  readRDS(.)
trait3 <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\preprocess\\trait2_abbr_ethnic.rds" %>%
  readRDS(.)
trait3 <- merge(PIDgroup[PIDgroup$group=="test_set", ], trait3, by= "PID")
trait3 <- trait3[, -2]

compMap <- data.frame(
  compAbbr = c("int", "ren", "per", "ner", "mic", "met", "mac", "eye", "car"),
  compLab = c("Total complications", "DKD", "DPVD", "DN", "Microvascular", "Metabolic Disorder", "Macrovascular", "DR", "DCVD")
)
timeMap <- data.frame(
  timeList = c(5, 10, 20),
  timeLab = c("5-year", "10-year", "All")
)

adjVar <- c(
  "Glc.BBC","Glc.NMR","HbA1c","Chol","LDL.Direct",
  "TG","HDL.C.BBC","HDL.C.NMR","ApoB.BBC","ApoB.NMR",
  "AAAC","BMI","BMR","BFP"
)

for (i in aucP_path) {
  abbr <- str_sub(i, 66, 68)
  time <- str_sub(i, 66, 72) %>%
    str_extract_all(., "\\d") %>%
    unlist(.) %>%
    str_c(.,  collapse = "") %>%
    as.numeric(.)
  
  comp <- paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\LassoCox\\", abbr, "Comp.rds") %>%
    readRDS(.)
  colnames(comp)[4:5] <- c("comp", "comp.s")
  # Because the covariates were collected at baseline, the course of complications
  # was calculated from baseline
  comp3 <- comp[, c("Participant.ID", "comp", "comp.s")]
  colnames(comp3) <- c("Participant.ID", "bslToComp", "comp.s")
  
  comp4 <- comp3 %>%
    mutate(
      .,
      ObservedIllness = case_when(
        bslToComp >= 365.25*time ~ 0,
        bslToComp < 365.25*time & comp.s == 1 ~ 1,
        TRUE ~ 0
      )
    )
  
  comp5 <- merge(comp4, trait3, by.x = "Participant.ID", by.y = "PID")
  comp6 <- comp5[, -c(1:3)]
  n <- ncol(comp6)
  df <- comp6[, c(2:n, 1)]
  df <- df %>%
    mutate(., ObservedIllness = as.factor(ObservedIllness))
  
  aucP <- readRDS(i) %>%
    subset(., sigVar==1)
  
  cox_summ1 <- paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\LassoCox\\", abbr, "_cox_summ.rds") %>%
    readRDS(.)
  cox_summ2 <- cox_summ1[["coefficients"]] %>%
    as.data.frame(.) %>%
    subset(., `Pr(>|z|)`<0.05)
  
  df2 <- df[, c(adjVar, intersect(aucP$abbr, rownames(cox_summ2)), "ObservedIllness")]
  
  {#~ Calculate auc  ------
    roc_HbA1c <- roc(
      df2$ObservedIllness,
      df2$HbA1c
    )
    auc_HbA1c <- auc(roc_HbA1c)
    
    aucDelong <- data.frame(var = NA, auc = NA,delongPvalueToHbA1c = NA)
    for (k in 1:(ncol(df2)-1)) {
      roc1 <- roc(df2$ObservedIllness, df2[, k])
      auc1 <- auc(roc1)
      delong1 <- roc.test(roc1, roc_HbA1c, method="delong")
      aucDelong <- rbind(aucDelong, c(colnames(df2)[k], auc1, delong1$p.value))
    }
    
    aucDelong$var <- as.character(aucDelong$var)
    aucDelong$auc <- as.numeric(aucDelong$auc)
    aucDelong$delongPvalueToHbA1c <- as.numeric(aucDelong$delongPvalueToHbA1c)
    aucDelong <- aucDelong[-1, ]
    aucDelong$delongPvalueToHbA1c <- p.adjust(aucDelong$delongPvalueToHbA1c, method = "bonferroni")
    
    aucDelong <- aucDelong %>%
      mutate(
        .,
        sigP = case_when(
          delongPvalueToHbA1c <= 0.01 ~ "**",
          delongPvalueToHbA1c > 0.01 & delongPvalueToHbA1c <= 0.05 ~ "*",
          TRUE ~ ""
        )
      )
    
    aucDelong2 <- arrange(aucDelong, desc(auc))
  }
  
  rocAll <- roc(ObservedIllness ~ ., data = df2)
  rocAll <- rocAll[aucDelong2$var]
  ggroc1 <- ggroc(rocAll, legacy.axes = TRUE) + 
    geom_abline(slope = 1, intercept = 0, color = "gray", linetype = "dashed", linewidth = 1) +
    geom_line(linewidth = 1.1) +
    theme(
      panel.background = element_rect(fill = "white", color = "white"),
      plot.background = element_rect(fill = "white", color = "white"),
      axis.text = element_text(color = "black", size = 13),
      axis.title = element_text(color = "black", size = 17),
      axis.line = element_line(color = "black", linewidth = 1),
      legend.position = "none"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "1-Specificity", y = "Sensitivity") +
    scale_colour_manual(values = paletteer_d("khroma::smoothrainbow", ncol(df2)-1))
  print(ggroc1)
  
  aucDelong3 <- aucDelong2 %>% 
    arrange(desc(auc)) %>%
    mutate(var = factor(var, levels = var))
  barPlot <- ggplot(aucDelong3, aes(var,auc))+
    coord_flip() +
    geom_col(aes(fill=var))+
    scale_fill_manual(values = paletteer_d("khroma::smoothrainbow", nrow(aucDelong))) +
    theme(legend.position = "none") +
    labs(y = "AUC") +
    theme(
      panel.background = element_rect(fill = "white", color = "white"),
      plot.background = element_rect(fill = "white", color = "white"),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y=element_blank(),
      axis.title = element_text(color = "black", size = 17),
      axis.text = element_text(color = "black", size = 13)
      ) + 
    geom_text(
      aes(label = sprintf("%.2f", auc)),
      hjust = -0.2,
      size = 3.5,
      color = "black"
    )+
    scale_y_continuous(limits = c(0, max(aucDelong3$auc)*1.3), expand = c(0, 0)) +
    annotate("rect", 
             xmin = -Inf, xmax = Inf, 
             ymin = -Inf, ymax = Inf, 
             color = "black", fill = NA, linewidth = 1.3) +
    geom_text(
      aes(y=0.01,label = var),
        size = 5,
        hjust = 0
      )
  print(barPlot)
  
  fullplot <- ggroc1 + barPlot + plot_annotation(
    title = paste0(
    timeMap[which(timeMap$timeList==time), "timeLab"],
    " incident ",
    compMap[which(compMap$compAbbr==abbr), "compLab"]
  ),
  theme = theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  )) +
  plot_layout(widths = c(7.5, 2.5))
  
  print(fullplot)
  
  ggsave(
    paste0("D:\\work2\\t2dAndComplicationsUkb\\plotDemo\\ImportantProteinRocAuc2ndEd\\", abbr, time, "Roc.pdf"),
    plot = fullplot,
    device = cairo_pdf,
    width = 10,
    height = 10
  )
  write.csv(
    aucDelong2,
    file = paste0("D:\\work2\\t2dAndComplicationsUkb\\plotDemo\\ImportantProteinRocAuc2ndEd\\", abbr, time, "AucDelong.csv"),
    row.names = FALSE
  )
  saveRDS(
    intersect(aucP$abbr, rownames(cox_summ2)),
    file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\ImpAUCROC\\", abbr, time, "_intersectVar.rds")
  )
}
