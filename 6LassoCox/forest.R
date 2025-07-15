library(dplyr)
library(stringr)
library(openxlsx)
library(forestploter)
library(grid)
library(ggplot2)

labName <- c("ner", "ren", "eye", "car", "per", "met", "mic", "mac", "int")
compAbb <- c(
  "DN","DKD","DR","DCVD","DPVD",
  "Metabolic Disorder","Microvascular",
  "Macrovascular","Total complication"
)

for (i in 1:length(labName)) {
  disAbb <- labName[i]
  data <- paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\LassoCox\\", disAbb, "_cox_summ.rds") %>%
    readRDS(.)
  data <- data[["coefficients"]] %>%
    as.data.frame(.)
  
  colnames(data) <- make.names(colnames(data))
  data$VarName <- rownames(data)
  colnames(data)[c(6, 2, 5)] <- c("VarName", "estimate", "P value")
  
  data$conf.low <- (data$coef - data$se.coef.*1.96) %>%
    exp(.)
  data$conf.high <- (data$coef + data$se.coef.*1.96) %>%
    exp(.)
  
  for (j in c(1:5, 7, 8)) {
    data[, j] <- round(data[, j], 3)
  }
  
  data$`HR(95%CI)` <- paste0(data$estimate, "(", data$conf.low, "-", data$conf.high, ")")
  
  data2 <- data %>%
    mutate(
      .,
      group = case_when(
        str_detect(VarName, "Sex") ~ "Sex",
        str_detect(VarName, "Ethnic.Continent") ~ "Ethnic",
        str_detect(VarName, "Current.tobacco.smoking...Instance.0") ~ "Smoke",
        str_detect(VarName, "Alcohol.drinker.status...Instance.0") ~ "Alcohol",
        str_detect(VarName, "Sleeplessness...insomnia...Instance.0") ~ "Sleeplessness",
        TRUE ~ "Protein and Metabolite"
      )
    ) %>%
    mutate(
      .,
      Variable = case_when(
        group == "Protein and Metabolite" ~ VarName,
        group == "Sex" ~ str_sub(VarName, 4, -1),
        group == "Ethnic" ~ str_sub(VarName, 17, -1),
        group == "Smoke" ~ str_sub(VarName, 37, -1),
        group == "Alcohol" ~ str_sub(VarName, 36, -1),
        group == "Sleeplessness" ~ str_sub(VarName, 38, -1),
        TRUE ~ "temp"
      )
    )%>%
    mutate(
      .,
      `Significance symbol` = case_when(
        `P value` <= 0.01 ~ "**",
        `P value` > 0.01 & `P value` <= 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  data2 <- subset(data2, `P value`<0.05)
  
  data2$Variable <- paste0("    ", data2$Variable)
  
  ncol_df <- ncol(data2)
  na_df <- data2$group %>%
    unique(.) %>%
    length(.) %>%
    matrix(NA, nrow = ., ncol = ncol_df) %>%
    as.data.frame(.)
  colnames(na_df) <- colnames(data2)
  na_df[, c("group", "Variable")] <- unique(data2$group)
  na_df[, "VarName"] <- "AAAAA" #Create a tag for sorting
  data3 <- rbind(data2, na_df)
  data4 <- data3[order(data3$group, data3$VarName), ]
  
  data5 <- data4[, c("Variable", "estimate", "P value", "conf.low", "conf.high", "HR(95%CI)", "Significance symbol")]
  data5[is.na(data5)] <- ""
  data5$` ` <- paste(rep(" ", nrow(data5)), collapse = " ")
  p <- forest(
    data = data5[, c("Variable", "HR(95%CI)", " ", "P value", "Significance symbol")],
    lower = as.numeric(data5$conf.low),
    upper = as.numeric(data5$conf.high),
    est = as.numeric(data5$estimate),
    ci_column = 3,
    sizes = (as.numeric(data5$estimate)+0.001)*0.4, 
    ref_line = 1, 
    xlim = c(min(data$conf.low)-0.05,max(data$conf.high)+0.05)
  )
  
  ggsave(
    paste0("D:\\work2\\t2dAndComplicationsUkb\\plotDemo\\forest2ndEd\\", disAbb, ".pdf"),
    plot = p,
    device = cairo_pdf,
    limitsize = FALSE,
    width =30,
    height =55
    ) 
}
