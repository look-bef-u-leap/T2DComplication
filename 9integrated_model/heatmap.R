library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(paletteer)

integrated_model_paths <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integrated_model" %>%
  list.files(path = ., pattern = ".*", full.names = TRUE, recursive = TRUE)

abbrMapDf <- data.frame(
  abbr=c("car", "eye", "int", "mac", "met", "mic", "ner", "per", "ren"),
  name=c("DCVD", "DR", "Total complications", "Macrovascular", "Metabolic Disorder", "Microvascular", "DN", "DPVD", "DKD")
)
timeMapDf <- data.frame(
  time=c("10", "20", "5"),
  label=c("10-year incident", "All incident", "5-year incident")
)

adjVar1 <- c(
  "Sex","Ethnic.Continent","CTS","ADS","Insomnia",
  "Glc.BBC","Glc.NMR","HbA1c","Chol","LDL.Direct",
  "TG","HDL.C.BBC","HDL.C.NMR","ApoB.BBC","ApoB.NMR",
  "AAAC","BMI","BMR","BFP"
)
adjVar2 <- c(
  "Glc.BBC","HbA1c", "Chol", "LDL.Direct", "TG",
  "HDL.C.BBC", "ApoB.BBC"
)

compTimeList <- sapply(1:length(integrated_model_paths), function(x){
  underscore_positions <- str_locate_all(integrated_model_paths[x], "_")[[1]]
  n0 <- underscore_positions[3, "start"]
  compTime <- str_sub(integrated_model_paths[x], 65, n0-1)
  return(compTime)
  }
  ) %>%
  unique(.)

for (i in 1:length(compTimeList)) {
  comp_time_path <- integrated_model_paths[str_detect(integrated_model_paths, compTimeList[i])]
  abbr <- str_sub(compTimeList[i], 1, 3)
  timePoint <- str_remove(compTimeList[i], paste0(abbr, "_"))
  df <- data.frame(model=c("LightGBM", "Logistic", "Random Forest", "Naive Bayes",
                           "SVM", "Single-layer neural network", "XGBoost", "CatBoost",
                           "MLP"))
  for (j in 1:length(comp_time_path)) {
    multiMOdel <- comp_time_path[j] %>%
      readRDS(.)
    auc <- sapply(
      1:length(multiMOdel),
      function(x){
        multiMOdel[[x]][[4]] %>%
          return(.)
      }
    ) %>%
      as.data.frame(.)
    auc$model <- c("LightGBM", "Logistic", "Random Forest", "Naive Bayes",
                   "SVM", "Single-layer neural network", "XGBoost", "CatBoost",
                   "MLP")
    df <- merge(df, auc, by="model")
    var <- names(multiMOdel[["bys_res"]][["bys_model"]][["isnumeric"]])
    if(all(adjVar1 %in% var)){
      var2 <- ifelse(var %in% adjVar1, "Panel1", var) %>%
        unique(.)
    }else{
      var2 <- ifelse(var %in% adjVar2, "Panel2", var) %>%
        unique(.)
    }
    varStr <- ifelse(
      length(var2)>2 & any(str_detect(var2, "Panel")),
      paste0(var2[1], " + ", "full traits"),
      paste0(var2, collapse = " + ")
    )
    varStr2 <- ifelse(
      (!str_detect(varStr, "Panel")) & str_detect(varStr, "\\+"), 
      "full traits", 
      varStr
    )
    colnames(df)[ncol(df)] <- varStr2
  }
  rownames(df) <- df$model
  df2 <- df[, -1] %>%
    t(.) %>%
    as.data.frame(.)
  
  title1 <- paste0(
    timeMapDf[timeMapDf$time==timePoint, "label"],
    " ",
    abbrMapDf[abbrMapDf$abbr==abbr, "name"]
    )
  df2$group <- title1
  df2$var <- rownames(df2)
  
  assign(paste0(abbr, timePoint), df2)
}

carModelTime <- rbind(car5,car10, car20)
eyeModelTime <- rbind(eye10, eye20)
intModelTime <- rbind(int5,int10,int20)
macModelTime <- rbind(mac5,mac10, mac20)
metModelTime <- rbind(met20)
micModelTime <- rbind(mic5,mic10, mic20)
nerModelTime <- rbind(ner10, ner20)
perModelTime <- rbind(per10, per20)
renModelTime <- rbind(ren5,ren10, ren20)

MT <- list(
  carModelTime,
  eyeModelTime,
  intModelTime,
  macModelTime,
  metModelTime,
  micModelTime,
  nerModelTime,
  perModelTime,
  renModelTime)
names(MT) <- c("car", "eye", "int", "mac", "met", "mic", "ner", "per", "ren")

for (i in 1:length(MT)) {
  df <- MT[[i]]
  dat_matrix <- df[, -c(10,11)] %>%
    as.matrix()
  
  dat_matrix_rev <- dat_matrix[, ncol(dat_matrix):1]
  
  rownames(dat_matrix) <- df$var
  
  green_pink <- colorRamp2(
    breaks = c(min(dat_matrix), mean(dat_matrix), max(dat_matrix)),
    colors = c("#4DAC26", "#F7F7F7", "#D01C8B")
  )
  
  df$group <- factor(df$group, levels = unique(df$group))
  
  pdf(
    file = paste0("D:\\work2\\t2dAndComplicationsUkb\\plotDemo\\integrated_model\\", names(MT)[i], "_circle_heatmap.pdf"),
    width = 60,
    height = 60
    )
  
  stdg <- 45
  circos.clear()
  circos.par(
    start.degree = stdg,
    canvas.xlim = c(-2, 2),
    canvas.ylim = c(-2, 2),
    gap.after = c(rep(5, length(levels(df$group)) - 1), 90-stdg),
    track.margin = c(0, 0.01),
    cell.padding = c(0, 0, 0, 0)
  )
  
  circos.heatmap(
    dat_matrix,
    split = df$group,
    cluster = FALSE,
    bg.border = "black",
    bg.lwd = 1,
    cell.border = "white",
    cell.lwd = 0.5,
    rownames.side = "outside",
    rownames.cex = 2,
    col = green_pink,
    track.height = 0.4
  )
  
  for(sector.index in get.all.sector.index()) {
    idx = which(df$group == sector.index)
    values = dat_matrix_rev[idx, , drop = FALSE]
    xlim = get.cell.meta.data("xlim", sector.index)
    ylim = get.cell.meta.data("ylim", sector.index)
    n_row = nrow(values)
    n_col = ncol(values)
    
    x_points = seq(xlim[1], xlim[2], length.out = n_row + 1)
    y_points = seq(ylim[1], ylim[2], length.out = n_col + 1)
    
    for(j1 in 1:n_row) {
      x = (x_points[j1] + x_points[j1 + 1]) / 2
      for(w in 1:n_col) {
        y = (y_points[w] + y_points[w + 1]) / 2
        circos.text(
          x, y, sprintf("%.2f", values[j1, w]), 
          sector.index = sector.index, 
          facing = "bending.inside", 
          cex = 0.8
        )
      }
    }
  }
  
  circos.track(
    track.index = get.current.track.index(),
    bg.border = NA,
    panel.fun = function(x, y){
      if(CELL_META$sector.numeric.index == length(levels(df$group))) {
        cn <- colnames(dat_matrix_rev)
        n <- length(cn)
        cell_height <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) / n
        y_coords <- seq(CELL_META$cell.ylim[1] + cell_height / 2, CELL_META$cell.ylim[2] - cell_height / 2, length.out = n)
        
        for (j2 in 1:n) {
          circos.lines(
            c(CELL_META$cell.xlim[2], CELL_META$cell.xlim[2] + convert_x(1, "mm")),
            c(y_coords[j2], y_coords[j2]),
            col = "black",
            lwd = 2
          )
        } 
        
        circos.text(
          rep(CELL_META$cell.xlim[2], n) + convert_x(1.45, "mm"),
          y_coords,
          cex = 2,
          adj = c(0, 0.5),
          cn,
          facing = "inside"
        )
      }
    }
  )
  
  circos.track(
    ylim = c(0, 1), 
    track.height = 0.05, 
    bg.col = adjustcolor(paletteer_d("ggthemes::wsj_dem_rep", length(unique(df$group))), alpha.f = 0.3),
    panel.fun = function(x, y) {
      circos.text(
        CELL_META$xcenter,
        CELL_META$ylim[2]-0.75,
        CELL_META$sector.index, 
        facing = "bending.inside", 
        cex = 2.5, 
        adj = c(0.5, 0)
      )
    }
  )
  
  heatmap_legend <- Legend(
    col_fun = green_pink, 
    direction = 'horizontal',
    at = c(0, 0.2, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
    labels = c('0.0', '0.2','0.5', '0.6', '0.7', '0.8', '0.9', '1.0'),
    break_dist = c(1, 1, 1, 2, 2, 3, 4),
    grid_height = unit(15, "mm"),
    legend_width = unit(140, "mm"),
    labels_gp = gpar(fontsize = 30)
  )
  draw(heatmap_legend, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = c("center", "center"))
  dev.off()
}
