library(lightgbm)
library(dplyr)
library(stringr)
library(ROCR)
library(pROC)
library(ggplot2)
library(treeshap)
library(shapviz)
library(tidyverse)
library(patchwork)
library(segmented)

PIDgroup <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\grouping\\PIDgroup.rds" %>%
  readRDS(.)

trait3 <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\preprocess\\trait2_abbr_ethnic.rds" %>%
  readRDS(.)

trait3 <- merge(PIDgroup, trait3, by = "PID")

# Integcomp2ded.rds was renamed integComp.rds in advance for ease of loop operation
CompAbbr <- c("int", "ren", "per", "ner", "mic", "met", "mac", "eye", "car")
timeList <- c(5, 10, 20) # 20 represents the total follow-up time

compLab <- c("Total complications", "DKD", "DPVD", "DN", "Microvascular", "Metabolic Disorder", "Macrovascular", "DR", "DCVD")
timeLab <- c("5-year", "10-year", "All")

for (i in 1:length(CompAbbr)) {
  for (j in 1:length(timeList)) {
    comp <- paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\LassoCox\\", CompAbbr[i], "Comp.rds") %>%
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
          bslToComp >= 365.25 * timeList[j] ~ 0,
          bslToComp < 365.25 * timeList[j] & comp.s == 1 ~ 1,
          TRUE ~ 0
        )
      )

    comp5 <- merge(comp4, trait3, by.x = "Participant.ID", by.y = "PID")
    comp6 <- comp5[, -c(1:3)]
    n <- ncol(comp6)
    df <- comp6[, c(2:n, 1)]
    df <- df %>%
      mutate(., ObservedIllness = as.factor(ObservedIllness))

    resBFpDf <- paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\univCox2ndEd\\resBFpDf_", CompAbbr[i], "Comp.rds") %>%
      readRDS(.)
    resBFpDf2 <- resBFpDf[which(resBFpDf$BFp < 0.05), ]
    df2 <- df[, c("group", rownames(resBFpDf2), "ObservedIllness")]

    data_train <- subset(df2, group == "train_set")[, -1]
    data_test <- subset(df2, group == "test_set")[, -1]

    data_train[, -ncol(data_train)] <- scale(data_train[, -ncol(data_train)])
    data_test[, -ncol(data_test)] <- scale(data_test[, -ncol(data_test)])

    if (table(data_train$ObservedIllness)[["1"]] <= 20) {
      next
    } else {
      params <- list(
        learning_rate = 0.01,
        num_leaves = 20,
        max_depth = 20,
        feature_fraction = 0.7,
        bagging_fraction = 0.7,
        bagging_freq = 200,
        num_iterations = 1000,
        objective = "binary"
      )
      lgbm.fit.perf <- function(dfTra, dfTes) {
        y <- dfTra$ObservedIllness %>% as.numeric(.)
        y2 <- y - 1
        x <- data.matrix(dfTra[, -ncol(dfTra)])
        colnames(x) <- colnames(dfTra)[-ncol(dfTra)]
        dtrain <- lgb.Dataset(data = x, label = y2)

        fit <- lgb.train(
          data = dtrain,
          params = params,
          nrounds = 100L,
          verbose = -1L
        )
        predTrain <- dfTra[, -ncol(dfTra)] %>%
          as.matrix(.) %>%
          fit$predict(.)
        FCPredTrain <- prediction(
          predictions = predTrain,
          labels = dfTra$ObservedIllness,
          label.ordering = c(0, 1)
        )
        aucDtrain <- performance(FCPredTrain, "auc")@y.values[[1]]

        predTest <- dfTes[, -ncol(dfTes)] %>%
          as.matrix(.) %>%
          fit$predict(.)
        FCPredTest <- prediction(
          predictions = predTest,
          labels = dfTes$ObservedIllness,
          label.ordering = c(0, 1)
        )
        aucDtest <- performance(FCPredTest, "auc")@y.values[[1]]

        importance_values <- lgb.importance(fit)
        importance_df <- data.frame(
          feature = importance_values$Feature,
          importance = importance_values$Gain
        )
        importance_df <- importance_df[order(importance_df$importance, decreasing = TRUE), ]

        res <- list(fit, predTrain, FCPredTrain, aucDtrain, predTest, FCPredTest, aucDtest, importance_df)
        names(res) <- c("fit", "predTrain", "FCPredTrain", "aucDtrain", "predTest", "FCPredTest", "aucDtest", "importance_df")
        return(res)
      }

      lgbmRes <- lgbm.fit.perf(data_train, data_test)

      stepwiseFit <- lapply(1:nrow(lgbmRes$importance_df), function(x) {
        varname <- lgbmRes$importance_df[1:x, "feature"]
        tempTra <- data_train[, c(varname, "ObservedIllness")]
        tempTes <- data_test[, c(varname, "ObservedIllness")]

        tempRes <- lgbm.fit.perf(tempTra, tempTes)
        return(tempRes)
      })

      delongList <- lapply(2:length(stepwiseFit), function(x) {
        roc1 <- roc(unlist(stepwiseFit[[x - 1]][["FCPredTest"]]@labels), unlist(stepwiseFit[[x - 1]][["FCPredTest"]]@predictions), direction = "<")
        roc2 <- roc(unlist(stepwiseFit[[x]][["FCPredTest"]]@labels), unlist(stepwiseFit[[x]][["FCPredTest"]]@predictions), direction = "<")
        delong <- roc.test(roc1, roc2, method = "delong")
        return(delong)
      })
      
      aucP <- data.frame(
        abbr = lgbmRes$importance_df$feature,
        auc = sapply(1:length(stepwiseFit), function(x){
          stepwiseFit[[x]][["aucDtest"]] %>%
            return(.)
        }
        )
      )
      aucP$ID <- 1:nrow(aucP)
      
      model.lm <- lm(auc~ID, data = aucP)
      model.lm.seg <- segmented(model.lm, seg.Z = ~ID, psi = 10)
      n2 <- confint(model.lm.seg)[, "Est."] %>%
        floor(.)
      aucP$sigVar <- rep(0, nrow(aucP))
      aucP[1:n2, "sigVar"] <- 1

      { # ~ shap ------
        unified <- unify(lgbmRes$fit, data_test[, -ncol(data_test)])
        treeshap1 <- treeshap(unified, data_test[, -ncol(data_test)])
        shp <- shapviz(treeshap1, X_pred = as.matrix(data_test[, -ncol(data_test)]))

        impNameSimp <- aucP[aucP$sigVar == 1, "abbr"]
        svImp2 <- sv_importance(shp[, impNameSimp], kind = "beeswarm", max_display = Inf, sort_features = FALSE) +
          ggtitle(paste0(timeLab[j], " incident ", compLab[i])) +
          theme(
            axis.line.x = element_line(color = "black", linewidth = 0.5),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"),
            panel.grid.major = element_line(color = "white", linetype = "solid"),
            panel.grid.minor = element_line(color = "white", linetype = "solid"),
            plot.title = element_text(hjust = 0.5)
          ) +
          labs(y = "Features")
      }

      { # ~ Feature importance ------
        importance_df <- merge(lgbmRes$importance_df, aucP, by.x = "feature", by.y = "abbr") %>%
          arrange(., desc(importance))

        colorPanel <- ifelse(importance_df$sigVar==1, "#FF6666", "#3d3d3d")
        fiplot <- ggplot() +
          geom_bar(
            data = importance_df,
            aes(x = reorder(feature, importance, decreasing = TRUE), y = importance), 
            alpha = 0.7,
            fill = "#3399cc",
            stat = "identity"
          ) +
          geom_point(
            data = importance_df,
            aes(x = reorder(feature, importance, decreasing = TRUE), y = auc/20), 
            color = colorPanel,
            size = 2
          ) +
          geom_line(
            data = importance_df, 
            aes(x = reorder(feature, importance, decreasing = TRUE), y = auc/20, group = 1), 
            color = c(colorPanel[-1], "#3d3d3d"),
            linewidth = 0.7
          ) +
          scale_y_continuous(
            sec.axis = sec_axis(~.*20, name = "AUC"),
            expand = c(0,0),
            limits = c(0, max(importance_df$importance)*1.03)
          ) +
          labs(x = "", y = "Importance") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                           color = colorPanel),
                panel.background = element_rect(fill = "white"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.ticks = element_line(color = "black")
          )
      }
      
      full_plot <- svImp2 / fiplot + plot_annotation(tag_levels = c(LETTERS[1:2]))

      ggsave(
        paste0("D:\\work2\\t2dAndComplicationsUkb\\plotDemo\\lgbmNonMlr3_2ndEd\\", CompAbbr[i],timeList[j], ".pdf"),
        plot = full_plot,
        device = cairo_pdf,
        limitsize = FALSE,
        width =20,
        height =10
      )
      saveRDS(
        lgbmRes,
        file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\lgbmNonMlr3_2ndEd\\", CompAbbr[i], timeList[j], "_lgbmRes.rds")
      )
      saveRDS(
        stepwiseFit,
        file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\lgbmNonMlr3_2ndEd\\", CompAbbr[i], timeList[j], "_stepwiseFit.rds")
      )
      saveRDS(
        delongList,
        file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\lgbmNonMlr3_2ndEd\\", CompAbbr[i], timeList[j], "_delongList.rds")
      )
      saveRDS(
        aucP,
        file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\lgbmNonMlr3_2ndEd\\", CompAbbr[i], timeList[j], "_aucP.rds")
      )
      write.csv(
        aucP,
        file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\lgbmNonMlr3_2ndEd\\", CompAbbr[i], timeList[j], "_aucP.csv"),
        row.names = TRUE
      )
      saveRDS(
        unified,
        file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\lgbmNonMlr3_2ndEd\\", CompAbbr[i], timeList[j], "_unified.rds")
      )
      saveRDS(
        treeshap1,
        file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\lgbmNonMlr3_2ndEd\\", CompAbbr[i], timeList[j], "_treeshap1.rds")
      )
      saveRDS(
        shp,
        file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\lgbmNonMlr3_2ndEd\\", CompAbbr[i], timeList[j], "_shp.rds")
      )
      saveRDS(
        svImp2,
        file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\lgbmNonMlr3_2ndEd\\", CompAbbr[i], timeList[j], "_svImp2.rds")
      )
      saveRDS(
        fiplot,
        file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\lgbmNonMlr3_2ndEd\\", CompAbbr[i], timeList[j], "_fiplot.rds")
      )
    }
  }
}
