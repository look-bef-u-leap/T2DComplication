library(dplyr)
library(stringr)
library(tidyverse)
library(survival)
library(survminer)

attrTitleVsName <- "D:\\work2\\t2dAndComplicationsUkb\\data\\vs\\attrTitleVsName.xlsx" %>%
  openxlsx::read.xlsx(., sheet = 1)
attrTitleVsName$fieldTitleAdj <- make.names(attrTitleVsName$fieldTitle)

trait3 <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\preprocess\\trait2.rds" %>%
  readRDS(.) %>%
  mutate(
    .,
    Ethnic.Continent = case_when(
      `Ethnic.background...Instance.0` %in% c("Irish", "British", "White",
                                              "Any other white background") ~ "Europe",
      `Ethnic.background...Instance.0` %in% c("Caribbean", "African",
                                              "Black or Black British",
                                              "Any other Black background") ~ "Africa",
      `Ethnic.background...Instance.0` %in% c("Chinese", "Pakistani",
                                              "Bangladeshi", "Indian",
                                              "Asian or Asian British",
                                              "Any other Asian background") ~ "Asia",
      `Ethnic.background...Instance.0` %in% c("White and Black African", "Mixed",
                                              "White and Asian","White and Black Caribbean",
                                              "Any other mixed background") ~ "Mixed",
      TRUE ~ "Other"
    ))
trait3$Ethnic.Continent <- as.factor(trait3$Ethnic.Continent)
trait3 <- trait3[, c(1, 2, 1782, 4:1781)]

n0 <- match(colnames(trait3), attrTitleVsName$fieldTitleAdj)
colnames(trait3) <- attrTitleVsName[n0, "Abbreviation"] %>%
  str_remove(., "_i0$")
colnames(trait3)[3] <- "Ethnic.Continent"

all_file_paths <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\vifLassoCox" %>%
  list.files(path = ., pattern = ".*", full.names = TRUE, recursive = TRUE)
all_file_paths
compIndex <- which(str_detect(all_file_paths, "Comp.rds"))
comp_paths <- all_file_paths[compIndex]
for (i in 1:length(comp_paths)) {
  tempComp <- comp_paths[i] %>%
    readRDS(.)
  
  tempComp2 <- tempComp %>%
    mutate(
      .,
      T2D5 = case_when(
        tempComp[, 2] >= 365.25*5 ~ as.factor(0),
        tempComp[, 2] < 365.25*5 & tempComp[, 3] == 1 ~ as.factor(1),
        TRUE ~ as.factor(0)
      ),
      T2D10 = case_when(
        tempComp[, 2] >= 365.25*10 ~ as.factor(0),
        tempComp[, 2] < 365.25*10 & tempComp[, 3] == 1 ~ as.factor(1),
        TRUE ~ as.factor(0)
      ),
      T2D20 = tempComp[, 3] %>%
        as.factor(.),
      comp5 = case_when(
        tempComp[, 4] >= 365.25*5 ~ as.factor(0),
        tempComp[, 4] < 365.25*5 & tempComp[, 5] == 1 ~ as.factor(1),
        TRUE ~ as.factor(0)
      ),
      comp10 = case_when(
        tempComp[, 4] >= 365.25*10 ~ as.factor(0),
        tempComp[, 4] < 365.25*10 & tempComp[, 5] == 1 ~ as.factor(1),
        TRUE ~ as.factor(0)
      ),
      comp20 = tempComp[, 5] %>%
        as.factor(.),
      srv5 = case_when(
        tempComp[, 6] >= 365.25*5 ~ as.factor(0),
        tempComp[, 6] < 365.25*5 & tempComp[, 7] == 1 ~ as.factor(1),
        TRUE ~ as.factor(0)
      ),
      srv10 = case_when(
        tempComp[, 6] >= 365.25*10 ~ as.factor(0),
        tempComp[, 6] < 365.25*10 & tempComp[, 7] == 1 ~ as.factor(1),
        TRUE ~ as.factor(0)
      ),
      srv20 = tempComp[, 7] %>%
        as.factor(.),
    )
  colnames(tempComp2)[11:13] <- paste0(colnames(tempComp2)[4], c("5", "10", "20"))
  trait3 <- merge(trait3, tempComp2[, c(1, 4, 5, 11:13)], by.x = "PID", by="Participant.ID")
}

traComp <- merge(trait3, tempComp2[, -c(4, 5, 11:13)], by.x = "PID", by="Participant.ID")
traComp$CTS <- fct_relevel(traComp$CTS,"Prefer not to answer", "No","Only occasionally","Yes, on most or all days")
traComp$ADS <- fct_relevel(traComp$ADS,"Prefer not to answer", "Never","Previous","Current")
traComp$Insomnia <- fct_relevel(traComp$Insomnia,"Prefer not to answer", "Never/rarely","Sometimes","Usually")

  {# group ------
    traComp$group <- NA
    set.seed(123)
    n <- nrow(traComp)
    indices <- sample(1:n, size = n / 2)
    traComp[indices, "group"] <- "train_set"
    traComp[-indices, "group"] <- "test_set"
    traComp$group <- as.factor(traComp$group)
  }
  
  dfCla <- sapply(traComp, class) %>%
    as.data.frame(.)
  colnames(dfCla) <- "class"
  dfCla <- arrange(dfCla, class)
  
  Chisquared <- rownames(dfCla)[c(1, 2, 6:17, 19:38)]
  fisher <- c("met5")
  Mann_Whitney_U <- rownames(dfCla)[c(3, 4, 5)]
  t_test <- colnames(traComp)[7:1781]
  
  n1 <- match(c(Chisquared, Mann_Whitney_U, t_test), rownames(dfCla))
  log_rank <- rownames(dfCla)[-n1]
  log_rank <- log_rank[-c(1:3)] %>%
    sort(.)
  
  group_diff <- data.frame(traitName=NA, p=NA,method=NA)
  {# Chi_squared_test ------
    for (i in 1:length(Chisquared)) {
      data <- as.matrix(
        rbind(
          table(traComp[indices, Chisquared[i]]),
          table(traComp[-indices, Chisquared[i]])
        )
      )
      rownames(data) <- c("train_set", "test_set")
      result <- chisq.test(data)
      group_diff <- rbind(group_diff, c(Chisquared[i], result$p.value, "Chi_squared"))
    }
  }
  
  {# fisher_test ------
    data <- as.matrix(
      rbind(
        table(traComp[indices, fisher]),
        table(traComp[-indices, fisher])
      )
    )
    rownames(data) <- c("train_set", "test_set")
    result <- fisher.test(data)
    group_diff <- rbind(group_diff, c(fisher, result$p.value, "fisher"))
  }
  
  {# Mann_Whitney_U & Wilcoxon rank sum test ------
    for (i in 1:length(Mann_Whitney_U)) {
      result <- wilcox.test(
        as.numeric(traComp[indices, Mann_Whitney_U[i]]),
        as.numeric(traComp[-indices, Mann_Whitney_U[i]]),
        paired = FALSE,
        exact = FALSE, 
        correct = FALSE
      )
      group_diff <- rbind(group_diff, c(Mann_Whitney_U[i], result$p.value, "Mann_Whitney_U"))
    }
  }
  
  {# t_test ------
    for (i in 1:length(t_test)) {
      res1 <- var.test(get(t_test[i]) ~ group, data = traComp)
      if (res1$p.value>0.05){
        result <- t.test(get(t_test[i]) ~ group, data = traComp, var.equal = T)
        group_diff <- rbind(group_diff, c(t_test[i], result$p.value, "t_test"))
      }else{
        result <- wilcox.test(get(t_test[i])~group, data=traComp, var.equal=FALSE)
        group_diff <- rbind(group_diff, c(t_test[i], result$p.value, "Mann_Whitney_U"))
      }
    }
  }
  
  {# log rank ------
    for (i in seq(1, 21, by=2)) {
      result <- survdiff(
        Surv(get(log_rank[i]),get(log_rank[i+1])) ~ group,
        traComp
      )
      group_diff <- rbind(group_diff, c(log_rank[i], result$pvalue, "log_rank"))
    }
  }
  
  group_diff <- group_diff[-1, ]
  group_diff$p <- as.numeric(group_diff$p)
  for (i in p.adjust.methods) {
    group_diff[, i] <- p.adjust(group_diff$p, method = i)
  }

  saveRDS(
    group_diff,
    file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\grouping\\group_diff.rds"
  )
  write.csv(
    group_diff,
    file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\grouping\\group_diff.csv",
    row.names = FALSE
  )
  saveRDS(
    traComp[, c("PID", "group")],
    file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\grouping\\PIDgroup.rds"
  )
  write.csv(
    traComp[, c("PID", "group")],
    file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\grouping\\PIDgroup.csv",
    row.names = FALSE
  )
  