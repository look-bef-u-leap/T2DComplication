# This script is used to draw heat maps of potential risk factors derived from multi-state models

library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(openxlsx)
library(mstate)
library(stringr)
library(tidyr)

cfull2ndEd <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\cfull2ndEd.rds" %>%
  readRDS(.)
attrTitleVsName <- "D:\\work2\\t2dAndComplicationsUkb\\data\\vs\\attrTitleVsName.xlsx" %>%
  read.xlsx(., sheet = 1)
attrTitleVsName$fieldTitle <- make.names(attrTitleVsName$fieldTitle, unique = TRUE)

coef <- summary(cfull2ndEd)[["coefficients"]] %>%
  as.data.frame(.)
coef$tran <- coef %>%
  rownames(.) %>%
  str_sub(., -1, -1)

coef$varName <- coef %>%
  rownames(.) %>%
  str_sub(., 1, -3)

coef2 <- coef %>% mutate(
  .,
  sigP = case_when(
    `Pr(>|z|)` < 0.05 & `Pr(>|z|)` >= 0.01 ~ "*",
    `Pr(>|z|)` < 0.01 ~ "**",
    TRUE ~ ""
  ),
  lastChar = substr(varName, nchar(varName), nchar(varName)),
  mappedChar = case_when(
    lastChar == "1" ~ "Q2",
    lastChar == "2" ~ "Q3",
    lastChar == "3" ~ "Q4",
    TRUE ~ lastChar
  ),
  fieldTitle = substr(varName, 1, nchar(varName) - 1)
)

ind <- match(coef2$fieldTitle, attrTitleVsName$fieldTitle)
attrTitleVsName$AbbreviationNew <- sub("_i.*", "", attrTitleVsName$Abbreviation)
coef2$Abbr <- attrTitleVsName[ind, "AbbreviationNew"]
coef2$AbbrQ <- paste0(coef2$Abbr, "_", coef2$mappedChar)

coef3 <- cbind(summary(cfull2ndEd)[["conf.int"]][, -2], coef2[, 5:13])
coef3$`exp(coef)` <- round(coef3$`exp(coef)`, 2)
coef3$`lower .95` <- round(coef3$`lower .95`, 2)
coef3$`upper .95` <- round(coef3$`upper .95`, 2)
coef3$`HR(95%ICs)` <- paste0(
  coef3$`exp(coef)`,
  "(",
  coef3$`lower .95`,
  "-",
  coef3$`upper .95`,
  ")"
  )
# For ease of presentation, assign NA manually
coef3[c(278, 280, 290), "HR(95%ICs)"] <- NA
coef3[c(278, 280, 290), "exp(coef)"] <- NA

coef3$lable <- paste(coef3$`HR(95%ICs)`, coef3$sigP, sep = "\n")

col_fun <- colorRamp2(c(0, 1, 2), c("#358faf", "white", "#fe8e63"))

Q_c <- c("Q2", "Q3", "Q4")
{# ~ Long data is converted to wide data and divided by Q ------
for (n in 1:length(Q_c)) {
  Q_width <- spread(
    coef3[which(coef3$mappedChar == Q_c[n]), c("exp(coef)", "tran", "AbbrQ", "Abbr", "mappedChar")],
    key = "tran",
    value = "exp(coef)"
  )
  rownames(Q_width) <- Q_width$Abbr %>%
    str_remove(., ".i0$")
  Q_width2 <- as.matrix(Q_width[, -c(1:3)])
  
  Q_lable <- spread(
    coef3[which(coef3$mappedChar == Q_c[n]), c("sigP", "tran", "AbbrQ", "Abbr", "mappedChar")],
    key = "tran",
    value = "sigP"
  )
  rownames(Q_lable) <- Q_lable$Abbr
  Q_lable2 <- as.matrix(Q_lable[, -c(1:3)])
  
  Q_p <- spread(
    coef3[which(coef3$mappedChar == Q_c[n]), c("Pr(>|z|)", "tran", "AbbrQ", "Abbr", "mappedChar")],
    key = "tran",
    value = "Pr(>|z|)"
  )
  rownames(Q_p) <- Q_p$Abbr
  Q_p2 <- as.matrix(Q_p[, -c(1:3)])
  
  heatmap <- Heatmap(
    Q_width2,
    name = "HR",
    col = col_fun,
    column_title = paste0(Q_c[n], " vs Q1"),
    rect_gp = gpar(col = "white", lty = 1, lwd = 2),
    column_title_gp = gpar(fontsize = 15, fontface = "bold"),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    column_names_rot = 0,
    cell_fun = function(j, i, x, y, width, height, fill) {
      if(Q_p2[i, j] < 0.05 & !is.na(Q_p2[i, j]))
        grid.text(paste0("", Q_lable2[i, j]), x, y, gp = gpar(fontsize = 24))
    },
    border = TRUE,
    border_gp = gpar(col = "#4b888f", lty = 2, lwd = 3.5)
    )
  
  assign(paste0(Q_c[n], "heatmap"), heatmap)
  }
}

Q2heatmap + Q3heatmap + Q4heatmap
