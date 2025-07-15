library(dplyr)
library(ggplot2)
library(shapviz)
library(ggchicklet)
library(stringr)

all_file_paths <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\lgbmNonMlr3_2ndEd" %>%
  list.files(path = ., pattern = ".*", full.names = TRUE, recursive = TRUE)
aucpIndex <- which(str_detect(all_file_paths, "_aucP.rds"))
aucp_path <- all_file_paths[aucpIndex]

abbrMapDf <- data.frame(
  abbr=c("car", "eye", "int", "mac", "met", "mic", "ner", "per", "ren"),
  name=c("DCVD", "DR", "Total complications", "Macrovascular", "Metabolic Disorder", "Microvascular", "DN", "DPVD", "DKD")
)
timeMapDf <- data.frame(
  time=c("10", "20", "5"),
  label=c("10-year incident", "All incident", "5-year incident")
)
embedXmin <- 26

for (i in 1:length(aucp_path)) {
  nameLabel <- str_sub(aucp_path[i], 66, 80) %>%
    str_remove(., "_aucP.rds")
  compAbb <- str_sub(nameLabel, 1, 3)
  timeBreak <- str_sub(nameLabel, 4, 6)
  
  lgbmRes <- paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\lgbmNonMlr3_2ndEd\\", nameLabel, "_lgbmRes.rds") %>%
    readRDS(.)
  aucP <- paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\lgbmNonMlr3_2ndEd\\", nameLabel, "_aucP.rds") %>%
    readRDS(.)
  shp <- paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\lgbmNonMlr3_2ndEd\\", nameLabel, "_shp.rds") %>%
    readRDS(.)
  
  {#~ bar plot ------
    importance_df <- merge(lgbmRes$importance_df, aucP, by.x = "feature", by.y = "abbr") %>%
      arrange(., desc(importance))
    
    n0 <- which(importance_df$sigVar==1)
    n1 <- c(1:(length(n0)+7))
    impDf2 <- importance_df[n1, ]
    
    colorPanel <- ifelse(impDf2$sigVar==1, "#FF6666", "#3d3d3d")
    fill_colors <- colorRampPalette(c("#955cb7", "#5abca3"))(nrow(impDf2)-1)
    fill_colors2 <- c(fill_colors, "white")
    
    n3 <- max(impDf2$auc)/max(impDf2$importance)
    fiplot <- ggplot() +
      geom_chicklet(
        data = impDf2,
        aes(x = reorder(feature, importance, decreasing = TRUE), y = importance), 
        alpha = 0.6,
        fill = fill_colors2,
        radius = grid::unit(2, "mm"),
        colour = NA
      ) +
      geom_chicklet(
        data = impDf2[nrow(impDf2), ],
        aes(x = reorder(feature, importance, decreasing = TRUE), y = importance), 
        fill = fill_colors[nrow(impDf2)],
        radius = grid::unit(2, "mm"),
        colour = tail(fill_colors, 1),
        linetype = "dashed",
        size = 1.1
      ) +
      geom_point(
        data = impDf2[-nrow(impDf2), ],
        aes(x = reorder(feature, importance, decreasing = TRUE), y = auc/n3), 
        color = colorPanel[-nrow(impDf2)],
        size = 2
      ) +
      geom_point(
        data = impDf2[nrow(impDf2), ],
        aes(x = reorder(feature, importance, decreasing = TRUE), y = auc/n3), 
        shape = 21,
        color = "black",
        fill = "white",
        size = 2.7
      ) +
      geom_line(
        data = impDf2, 
        aes(x = reorder(feature, importance, decreasing = TRUE), y = auc/n3, group = 1), 
        color = c(colorPanel[-1], "#3d3d3d"),
        linewidth = 0.7
      ) +
      scale_y_continuous(
        sec.axis = sec_axis(~.*n3, name = "AUC"),
        expand = c(0,0),
        limits = c(-0.001, max(impDf2$importance)*1.03)
      ) +
      labs(
        x = "", 
        y = "Importance",
        title = paste0(
          timeMapDf[timeMapDf$time==timeBreak, "label"],
          " ",
          abbrMapDf[abbrMapDf$abbr==compAbb, "name"]
        )
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = c(colorPanel[-length(colorPanel)], "gray")),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
      )
    }
  
  { # ~ shap ------
    
    #Due to the large number of scattered points in the SHAP colony map of 5-year DKD,
    #only 5,000 participants were randomly selected to draw the SHAP colony map.
    #The simplified colony map was basically consistent with the complete one
    # set.seed(123)
    # shp[sample(1:25000, size = 5000, replace = TRUE), impNameSimp]
    
    impNameSimp <- aucP[aucP$sigVar == 1, "abbr"]
    svImp2 <- sv_importance(shp[, impNameSimp], kind = "beeswarm", max_display = Inf, sort_features = FALSE) +
      theme(
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(color = "#FF6666"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "white", linetype = "solid"),
        panel.grid.minor = element_line(color = "white", linetype = "solid"),
        plot.title = element_text(hjust = 0.5)
      )
  }
  
  embedYmin <- impDf2[embedXmin, "importance"]*1.0005
  fullplot <- fiplot +
    annotation_custom(
      grob = ggplotGrob(svImp2),
      xmin = embedXmin,
      xmax = nrow(impDf2),
      ymin = embedYmin,
      ymax = min(impDf2[embedXmin:nrow(impDf2), "auc"])/n3*0.98
    )
  
  pdf(
    paste0("D:\\work2\\t2dAndComplicationsUkb\\plotDemo\\lgbmNonMlr3_2ndEd\\", compAbb, timeBreak, ".pdf"),
    width = 10.8,
    height = 8.1,
    compress = TRUE,
    onefile = TRUE,
    paper = "special"
    )
  print(fullplot)
  dev.off()
  gc()
}
