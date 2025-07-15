library(lightgbm)
library(Hmisc)
library(caret)
library(randomForest)
library(e1071)
library(nnet)
library(xgboost)
library(catboost)
library(adabag)
library(dplyr)
library(stringr)
library(ROCR)
library(pROC)
library(neuralnet)

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

trait3 <- merge(PIDgroup, trait3, by = "PID")

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

for (i in 1:length(intersect_paths)) {
  compAbbr <- str_sub(intersect_paths[i], 58, 60)
  time <- str_sub(intersect_paths[i], 58, 65) %>%
    str_extract_all(., "\\d") %>%
    unlist(.) %>%
    str_c(.,  collapse = "") %>%
    as.numeric(.)
  
  {#~ Variable combination ------
    intersectVar <- intersect_paths[i] %>%
      readRDS(.)
    
    samPreVar1 <- expand.grid(c1 = list(adjVar1, adjVar2, NA), c2 = c(NA, as.list(intersectVar))) %>%
      split(
        .,
        seq(nrow(.))
      ) %>%
      c(
        .,
        list(intersectVar, list(adjVar1, intersectVar),
             list(adjVar2, intersectVar))
      )
    samPreVar2 <- sapply(samPreVar1, function(x) {
      fx1 <- x %>%
        unlist(.) %>%
        na.omit(.) %>%
        unique(.)
      names(fx1) <- NULL
      return(fx1)
    }) %>%
      Filter(function(x) length(x) > 0, .) %>%
      unique(.)
  }
  
  comp <- comps_paths[which(str_detect(comps_paths, compAbbr))] %>%
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
  
  data_train <- subset(df, group == "train_set")[, -1]
  data_test <- subset(df, group == "test_set")[, -1]
  
  rm(comp5)
  rm(comp6)
  rm(df)
  
  data_train[, -c(1:5, ncol(data_train))] <- scale(data_train[, -c(1:5, ncol(data_train))])
  data_test[, -c(1:5, ncol(data_test))] <- scale(data_test[, -c(1:5, ncol(data_test))])
  
  for (x0 in 1:length(samPreVar2)) {
    data_train2 <- data_train[, c(samPreVar2[[x0]], "ObservedIllness")]
    data_test2 <- data_test[, c(samPreVar2[[x0]], "ObservedIllness")]
    
    {#~ lgbm ------
      y2 <- data_train2$ObservedIllness %>% as.numeric(.) - 1
      x <- data.matrix(data_train2[,-ncol(data_train2)])
      dtrain <- lgb.Dataset(data = x,label = y2)
      
      params <- list(
        learning_rate = 0.01,
        num_leaves = 20,
        max_depth=20,
        feature_fraction = 0.7,
        bagging_fraction = 0.7,
        bagging_freq = 200,
        num_iterations = 1000,
        objective = "binary"
      )
      
      lgb_model <- lgb.train(
        data = dtrain,
        params = params,
        nrounds = 100L,
        verbose = -1L
      )
      
      lgb_test_probability <- data_test2[,-ncol(data_test2)] %>%
        data.matrix(.) %>%
        lgb_model$predict(.)
      lgb_performance <- ROCR::prediction(
        predictions = lgb_test_probability,
        labels = data_test2$ObservedIllness,
        label.ordering = c(0,1)
      )
      lgb_auc <- performance(lgb_performance, "auc")@y.values[[1]]
      lgb_res <- list(lgb_model,lgb_test_probability,lgb_performance,lgb_auc)
      names(lgb_res) <- c("lgb_model", "lgb_test_probability", "lgb_performance", "lgb_auc")
    }
    
    {#~ logistics------
      log_model <- glm(ObservedIllness~.,data = data_train2,family = "binomial")
      log_test_probability <- predict(object =log_model,newdata=data_test2,type = "response")
      log_performance <- ROCR::prediction(
        predictions = log_test_probability,
        labels = data_test2$ObservedIllness,
        label.ordering = c(0,1)
      )
      log_auc <- performance(log_performance, "auc")@y.values[[1]]
      log_res <- list(log_model,log_test_probability,log_performance,log_auc)
      names(log_res) <- c("log_model","log_test_probability","log_performance","log_auc")
    }
    
    {#~ naiveBayes------
      bys_model <- naiveBayes(ObservedIllness~. ,data = data_train2)
      bys_test_probability <- predict(bys_model,newdata = data_test2, type = "raw") %>%
        as.data.frame(.)
      bys_performance <- ROCR::prediction(
        predictions = bys_test_probability$`1`,
        labels = data_test2$ObservedIllness,
        label.ordering = c(0,1)
      )
      bys_auc <- performance(bys_performance, "auc")@y.values[[1]]
      bys_res <- list(bys_model,bys_test_probability,bys_performance,bys_auc)
      names(bys_res) <- c("bys_model","bys_test_probability","bys_performance","bys_auc")
    }
    
    {#~ SVM------
      svm_model <- svm(
        ObservedIllness~.,
        data = data_train2,
        kernel="radial",
        gamma=1/ncol(data_train2),
        probability = TRUE
      )
      svm_test_probability <- predict(svm_model,newdata = data_test2, probability = TRUE) %>%
        attr(., "probabilities") %>%
        as.data.frame(.)
      svm_performance <- ROCR::prediction(
        predictions = svm_test_probability$`1`,
        labels = data_test2$ObservedIllness,
        label.ordering = c(0,1)
      )
      svm_auc <- performance(svm_performance, "auc")@y.values[[1]]
      svm_res <- list(svm_model,svm_test_probability,svm_performance,svm_auc)
      names(svm_res) <- c("svm_model","svm_test_probability","svm_performance","svm_auc")
    }
    
    {#~ Single-layer artificial neural network------
      slann_model <- nnet(ObservedIllness~.,data = data_train2,size=10)
      slann_test_probability <- predict(slann_model,newdata = data_test2)
      slann_performance <- ROCR::prediction(
        predictions = slann_test_probability,
        labels = data_test2$ObservedIllness,
        label.ordering = c(0,1)
      )
      slann_auc <- performance(slann_performance, "auc")@y.values[[1]]
      slann_res <- list(slann_model,slann_test_probability,slann_performance,slann_auc)
      names(slann_res) <- c("slann_model","slann_test_probability","slann_performance","slann_auc")
    }
    
    {#~ xgboost------
      y2 <- data_train2$ObservedIllness %>% as.numeric(.) - 1
      x <- data.matrix(data_train2[,-ncol(data_train2)])
      dtrain <- xgb.DMatrix(data = x,label = y2)
      
      xgb_model <- xgboost(
        data = dtrain, 
        label = y2, 
        max.depth = 20, 
        eta = 0.01, 
        nrounds = 100,
        nthread = 2, 
        objective = "binary:logistic"
      )
      xgb_test_probability <- data_test2[,-ncol(data_test2)] %>%
        data.matrix(.) %>%
        predict(xgb_model, .)
      xgb_performance <- ROCR::prediction(
        predictions = xgb_test_probability,
        labels = data_test2$ObservedIllness,
        label.ordering = c(0,1)
      )
      xgb_auc <- performance(xgb_performance, "auc")@y.values[[1]]
      xgb_res <- list(xgb_model,xgb_test_probability,xgb_performance,xgb_auc)
      names(xgb_res) <- c("xgb_model","xgb_test_probability","xgb_performance","xgb_auc")
    }
    
    {#~ catboost ------
      y2 <- data_train2$ObservedIllness %>% as.numeric(.) - 1
      x <- data_train2[,-ncol(data_train2)] %>%
        as.data.frame(.)
      colnames(x) <- colnames(data_train2)[-ncol(data_train2)]
      
      train_pool <- catboost.load_pool(data = x,label = y2)
      
      catb_model <- catboost.train(
        train_pool,
        NULL,
        params = list(
          loss_function = 'Logloss',
          iterations = 200,
          metric_period=20,
          learning_rate=0.01,
          depth=10
        ) 
      )
      x_test <- data_test2[, -ncol(data_test2)] %>%
        as.data.frame(.)
      colnames(x_test) <- colnames(data_test2)[-ncol(data_test2)]
      
      catb_test_probability <- x_test %>%
        catboost.load_pool(.) %>%
        catboost.predict(catb_model, ., prediction_type = "Probability")
      catb_performance <- ROCR::prediction(
        predictions = catb_test_probability,
        labels = data_test2$ObservedIllness,
        label.ordering = c(0,1)
      )
      catb_auc <- performance(catb_performance, "auc")@y.values[[1]]
      catb_res <- list(catb_model,catb_test_probability,catb_performance,catb_auc)
      names(catb_res) <- c("catb_model","catb_test_probability","catb_performance","catb_auc")
    }
    
    {#~ RF ------
      rf_model <- randomForest(
        ObservedIllness ~ .,
        data=data_train2,
        ntree=100,
        important=TRUE,
        proximity=TRUE
      )
      rf_test_probability <- predict(rf_model,newdata=data_test2, type = "prob") %>%
        as.data.frame(.)
      rf_performance <- ROCR::prediction(
        predictions = rf_test_probability$`1`,
        labels = data_test2$ObservedIllness,
        label.ordering = c(0,1)
      )
      rf_auc <- performance(rf_performance, "auc")@y.values[[1]]
      rf_res <- list(rf_model,rf_test_probability,rf_performance,rf_auc)
      names(rf_res) <- c("rf_model","rf_test_probability","rf_performance","rf_auc")
    }
    
    {#~ MLP ------
      dtrain <- data_train2 %>%
        data.matrix(.)
      dtrain[, "ObservedIllness"] <- dtrain[, "ObservedIllness"]-1
      
      dtest <- data_test2 %>%
        data.matrix(.)
      dtest[, "ObservedIllness"] <- dtest[, "ObservedIllness"]-1
      
      mlp_model <- neuralnet(
        ObservedIllness~.,
        data=dtrain,
        linear.output = FALSE,
        threshold = 0.05,
        stepmax = 1e6,
        rep = 1,
        err.fct = 'ce',
        hidden = c(1)
      )
      mlp_test_probability <- compute(mlp_model,covariate = dtest)
      mlp_performance <- ROCR::prediction(
        predictions = mlp_test_probability$net.result,
        labels = dtest[, "ObservedIllness"],
        label.ordering = c(0,1)
      )
      mlp_auc <- performance(mlp_performance, "auc")@y.values[[1]]
      mlp_res <- list(mlp_model,mlp_test_probability,mlp_performance,mlp_auc)
      names(mlp_res) <- c("mlp_model","mlp_test_probability","mlp_performance","mlp_auc")
    }
    
    res <- list(lgb_res,log_res,rf_res,bys_res,svm_res,slann_res,xgb_res,catb_res,mlp_res)
    names(res) <- c("lgb_res","log_res","rf_res","bys_res","svm_res","slann_res","xgb_res","catb_res","mlp_res")
    
    saveRDS(
      res,
      file = paste0("D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integrated_model\\", compAbbr,"_", time, "_", x0, "MultiComb.rds")
    )
  }
}

