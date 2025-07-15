# This script is used to map a volcano for single-factor cox analysis

library(dplyr)
library(ggplot2)
library(stringr)
library(openxlsx)
library(patchwork)
library(ggrepel)

attrTitleVsName <- "D:\\work2\\t2dAndComplicationsUkb\\data\\vs\\attrTitleVsName.xlsx" %>%
  read.xlsx(., sheet = 1)
attrTitleVsName$fieldTitleAdj <- make.names(attrTitleVsName$fieldTitle, unique = TRUE)

labName <- c("ner", "ren", "eye", "car", "per", "met", "mic", "mac", "int")

df <- data.frame(
  coef=NA,`exp.coef.`=NA,`se.coef.`=NA,z=NA,
  `Pr...z..`=NA,`BFp`=NA,`AbbAdj`=NA,`ID`=NA,
  label=NA,group=NA
)

for (i in 1:length(labName)) {
  disAbb <- labName[i]
  data <- paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\univCox2ndEd\\resBFpDf_", disAbb, "Comp.rds") %>%
    readRDS(.)
  
  data$AbbAdj <- rownames(data)
  data$ID <- 1:nrow(data)
  data2 <- subset(data, BFp<0.05) %>%
    arrange(., desc(abs(coef)))
  n0 <- data2[1:20, "ID"]
  data$label <- ifelse(data$ID %in% n0, data$AbbAdj, "")
  data$group <- disAbb
  colnames(data) <- make.names(colnames(data))
  
  df <- rbind(df, data)
}

df <- df[-1, ]%>%
  mutate(
    .,
    group = case_when(
      group=="ner" ~ "DN",
      group=="ren" ~ "DKD",
      group=="eye" ~ "DR",
      group=="car" ~ "DCVD",
      group=="per" ~ "DPVD",
      group=="met" ~ "Metabolic",
      group=="mic" ~ "Microvascular",
      group=="mac" ~ "Macrovascular",
      TRUE ~ "Total"
    ))

df <- df[df$group != "Total", ]

mycol <- c(
  '#b12d30', '#43b5e6', 
  "#b76f9e", '#58ac41', 
  '#f1aa41',"#6cc3b9",
  "#fc3c46","#3568a3"
)

options(ggrepel.max.overlaps = Inf)
p <- ggplot() +
  geom_point(data = df, aes(x = exp.coef., y = -log10(BFp)),size = 0.8, color = 'grey') +
  facet_grid(group~.,scales = "free") +
  #coord_flip() +
  #facet_grid(. ~ group,scales = "free") +
  #One row with multiple columns cannot display all labels and there will be no prompt,
  #so it is changed to multiple rows with one column
  geom_point(data = df[df$label != "", ], aes(x = exp.coef., y = -log10(BFp),color = group)) +
  geom_vline(xintercept = 1, linewidth = 0.5, color = "grey50", lty = 'dashed')+
  geom_hline(yintercept = -log10(0.05), linewidth = 0.5, color = "grey50", lty = 'dashed')+
  scale_color_manual(values = mycol) +
  xlab(label = "Effect size(HR)") + 
  ylab(label = "-log10 (Bonferroni P)") + 
  theme_bw()+
  theme( legend.position = 'none',
         panel.grid = element_blank(),
         axis.text = element_text(size = 10),
         axis.text.x = element_text(angle = 45, vjust = 0.8),
         strip.text.x = element_text(size = 10, face = 'bold')
  ) +
  geom_text_repel(
    data = df,
    aes(x = exp.coef., y = -log10(BFp),color = group, label = label),
    fontface = 'italic',
    seed = 233,
    size = 3,
    min.segment.length = 0,
    force = 20,
    box.padding = 0.1,
    max.overlaps = Inf,
    force_pull = 2,
    segment.alpha = 0.4,
    direction = "both",
    hjust = 0.5
  )
p
options(ggrepel.max.overlaps = 10)

gc()