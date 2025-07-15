library(dplyr)
library(stringr)
library(pROC)
library(survminer) 
library(survival)
library(ggplot2)
library(patchwork)

all_file_paths <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\LassoCox" %>%
  list.files(path = ., pattern = ".*", full.names = TRUE, recursive = TRUE)
compIndex <- which(str_detect(all_file_paths, "Comp.rds"))
comps_paths <- all_file_paths[compIndex]

intersect_paths <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\ImpAUCROC" %>%
  list.files(path = ., pattern = ".*", full.names = TRUE, recursive = TRUE)

PIDgroup <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\grouping\\PIDgroup.rds" %>%
  readRDS(.)
trait3 <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\preprocess\\trait2_abbr_ethnic.rds" %>%
  readRDS(.)
trait3 <- merge(PIDgroup[PIDgroup$group=="test_set", ], trait3, by= "PID")
trait3 <- trait3[, -2]

adjVar1 <- c(
  "Sex","Ethnic.Continent","CTS","ADS","Insomnia",
  "Glc.BBC","Glc.NMR","HbA1c","Chol","LDL.Direct",
  "TG","HDL.C.BBC","HDL.C.NMR","ApoB.BBC","ApoB.NMR",
  "AAAC","BMI","BMR","BFP"
)
abbrMapDf <- data.frame(
  abbr=c("car", "eye", "int", "mac", "met", "mic", "ner", "per", "ren"),
  name=c("DCVD", "DR", "Total complications", "Macrovascular", "Metabolic Disorder", "Microvascular", "DN", "DPVD", "DKD")
)

cutoffDf <- data.frame(trait=NA, comp=NA, cutoff=NA, HR=NA, low=NA, high=NA, Pvalue=NA)
for (i in 1:length(comps_paths)) {
  comp <- comps_paths[i] %>% readRDS(.)
  compTra <- merge(comp, trait3, by.x = "Participant.ID", by.y = "PID")
  colnames(compTra)[c(4, 5)] <- c("comp", "comp.s")
  compTra$comp <- compTra$comp/365.25
  
  compAbbr <- str_sub(comps_paths[i], 57, 59)
  var <- intersect_paths[str_detect(intersect_paths, paste0("/", compAbbr))] %>%
    sapply(., readRDS) %>%
    unlist(.) %>%
    unique(.)
  
  compTra2 <- compTra[, c("Participant.ID", "comp", "comp.s", var, adjVar1)]
  compTra2[, -c(1:3)] <- compTra2[, -c(1:3)] %>%
    mutate(across(where(is.numeric), ~scale(., center = TRUE, scale = TRUE) %>% as.vector()))
  
  for (j in var) {
    roc_var <- roc(
      compTra2$comp.s,
      compTra2[, j]
    )
    cutoff <- roc_var$thresholds[which.max(roc_var$sensitivities+roc_var$specificities-1)]
    
    compTra2$binaryVar <- ifelse(
      compTra2[, j] <= cutoff,
      "Low",
      "High"
    ) %>%
      factor(.,levels = c("Low", "High"))
    
    compTra3 <- compTra2 %>%
      select(., !any_of(c("Participant.ID", var)))
    res.cox <- coxph(Surv(comp,comp.s)~., data = compTra3)
    coxSummary <- summary(res.cox)
    coef <- coxSummary$coefficients %>%
      as.data.frame(.)
    coef$conf.low <- (coef$coef - coef$`se(coef)`*1.96) %>%
      exp(.)
    coef$conf.high <- (coef$coef + coef$`se(coef)`*1.96) %>%
      exp(.)
    coef$`HR(95%CI)` <- paste0(round(coef$`exp(coef)`, 3), " (", round(coef$conf.low, 3), "-", round(coef$conf.high, 3), ")")
    
    fit <- survfit(Surv(comp,comp.s) ~ binaryVar,data = compTra2)
    { # ~ km plot ------
      plabel <- gsub(
        "e",
        "*10",
        format(coef["binaryVarHigh", "Pr(>|z|)"], scientific = TRUE, digits = 3)
        )
      OS <- ggsurvplot(
        fit,
        palette = "aaas",
        risk.table = TRUE,
        risk.table.col = "strata",
        conf.int = TRUE,
        legend = c(0.17, 0.5),
        legend.title = "",
        legend.labs= c("Low", "High"),
        xlab ="Follow-up time(years)",
        ylab = paste0("The probability of non-", abbrMapDf[abbrMapDf$abbr==compAbbr, "name"]),
        xlim = c(-0.7, 17),
        conf.int.alpha = c(0.2),
        censor=FALSE,
        title = paste0(
          j,
          "-",
          abbrMapDf[abbrMapDf$abbr==compAbbr, "name"]
          ),
        break.y.by = 0.01,
        ylim=c(min(fit[["lower"]]), 1.00),
        font.legend = c(12, "plain")
      )
      
      OS$plot <- OS$plot +
        theme(
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank()
        ) +
        annotate(
          "text",
          x = 0.2,
          y = (1-min(fit[["lower"]]))*0.2+min(fit[["lower"]]),
          label= paste0(
            "HR (95%CI) = ",
            coef["binaryVarHigh", "HR(95%CI)"],
            "\nP value = ",
            plabel
          ),
          size = 4.5,
          hjust = 0
        )
      
      OS$table <- OS$table +
        theme(
          plot.title = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "none"
        )
      pdf(
        paste0("D:\\work2\\t2dAndComplicationsUkb\\plotDemo\\km\\", compAbbr, "_", j, ".pdf"),
        width = 4.68,
        height = 4.68
        )
      print(OS, newpage = FALSE)
      dev.off()
      
      cutoffDf <- rbind(
        cutoffDf,
        c(
          j,
          abbrMapDf[abbrMapDf$abbr==compAbbr, "name"],
          cutoff,
          coef["binaryVarHigh", "exp(coef)"],
          coef["binaryVarHigh", "conf.low"],
          coef["binaryVarHigh", "conf.high"],
          coef["binaryVarHigh", "Pr(>|z|)"]))
    }
  }
}

write.csv(
  cutoffDf,
  file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\km\\cutoffDf.csv"
)
saveRDS(
  cutoffDf,
  file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\km\\cutoffDf.rds"
)
