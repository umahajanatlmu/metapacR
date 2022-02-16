require(tidyverse)
##----------------------------------------------------------------
##                          base model                          --
##----------------------------------------------------------------
runBaseModel <- function(train,
                         test,
                         features,
                         response,
                         time,
                         cutoff_metric,
                         exclude_algos,
                         save.location) {
  automl <- data.frame()
  factor.models <- c("XGBoost", "DeepLearning", "GLM")
  ## model
  model.fs <- h2o.automl(
    x = features,
    y = response,
    training_frame = train,
    leaderboard_frame = test,
    keep_cross_validation_models = TRUE,
    keep_cross_validation_predictions = TRUE,
    keep_cross_validation_fold_assignment = TRUE,
    balance_classes = TRUE,
    max_runtime_secs = time,
    exclude_algos = exclude_algos,
    seed = 1234,
    sort_metric = "logloss",
    nfolds = 10,
    verbosity = "info"
  )
  ## get model
  fs.Model <- h2o.getModel(model.fs@leader@model_id)
  if (grepl("StackedEnsemble", model.fs@leader@model_id) == TRUE) {
    perf <- h2o.performance(fs.Model, newdata = test)
    ensemble_auc_test <- h2o.auc(perf)
    print(sprintf("Ensemble Test AUC:  %s", ensemble_auc_test))
    ##---------------------------------------------------------------
    ##                         metalearner                         --
    ##---------------------------------------------------------------
    metalearner <-
      h2o.getModel(model.fs@leader@model$metalearner$name)
    modelImp <- h2o.varimp(metalearner)
    highestImpName <- modelImp[1, 1]
    ## get model
    model  <-
      h2o.getModel(model.fs@leader@model$metalearner_model@model$names[highestImpName %in%
                                                                         model.fs@leader@model$metalearner_model@model$names])
    ## variable importance
    fs.Table <- as.data.frame(h2o.varimp(model))
    ## clear feature names
    if (grepl(paste(factor.models, collapse = "|"), model@model_id) == TRUE) {
      fs.Table$variable <- gsub("^(.*)[.].*", "\\1", fs.Table$variable)
    }
    ## filter variable table by threshould
    filtered.fs.Table <- fs.Table %>%
      distinct(variable, .keep_all = TRUE) %>%
      arrange(desc(scaled_importance))
    ## if not stacked ensemble but classifier models
  } else if (grepl(paste(factor.models, collapse = "|"), model.fs@leader@model_id) == TRUE) {
    ## variable importance
    fs.Table <- as.data.frame(h2o.varimp(fs.Model))
    ## clear feature names
    fs.Table$variable <- gsub("^(.*)[.].*", "\\1", fs.Table$variable)
    ## clear feature names
    filtered.fs.Table <- fs.Table %>%
      distinct(variable, .keep_all = TRUE) %>%
      arrange(desc(scaled_importance))
    ## if not stacked ensemble and classifier models
  } else {
    fs.Table <- as.data.frame(h2o.varimp(fs.Model))
    filtered.fs.Table <-
      fs.Table %>%
      distinct(variable, .keep_all = TRUE) %>%
      arrange(desc(scaled_importance))
  }
  ##----------------------------------------------------------------
  ##                         model metric                         --
  ##----------------------------------------------------------------
  logloss  <- h2o.logloss(h2o.performance(fs.Model, test))
  cutoff <-
    h2o.find_threshold_by_max_metric(h2o.performance(fs.Model, test),
                                     cutoff_metric)
  accuracy <-
    h2o.accuracy(h2o.performance(fs.Model, newdata = test), cutoff)
  auc <- h2o.auc(h2o.performance(fs.Model, newdata = test))
  ##----------------------------------------------------------------
  ##                      model metric table                      --
  ##----------------------------------------------------------------
  automl[1, "model.name"] <-  model.fs@leader@model_id
  automl[1, "logloss"] <-  logloss
  automl[1, "number.of.features"] <-  length(features)
  automl[1, "auc"] <-  h2o.auc(h2o.performance(fs.Model, newdata = test))
  automl[1, "accuracy"] <-
    h2o.accuracy(h2o.performance(fs.Model, newdata = test), cutoff)
  automl[1, "specificity"] <-
    h2o.specificity(h2o.performance(fs.Model, newdata = test), cutoff)
  automl[1, "sensitivity"] <-
    h2o.sensitivity(h2o.performance(fs.Model, newdata = test), cutoff)
  ##----------------------------------------------------------------
  ##                          leaderbord                          --
  ##----------------------------------------------------------------
  lb <- h2o.get_leaderboard(object = model.fs, extra_columns = 'ALL')
  leaderboard <- as.data.frame(lb)
  ## save leaderboard
  write.csv(leaderboard, paste0(save.location, "leaderboard.csv"))
  ##----------------------------------------------------------------
  ##                          save model                          --
  ##----------------------------------------------------------------
  h2o.saveModel(fs.Model, path = save.location, force = TRUE)
  ##---------------------------------------------------------------
  ##                   variable importance plot                   --
  ##---------------------------------------------------------------
  filtered.fs.Table <- filtered.fs.Table %>%
    mutate(nrow = ifelse(n() < 20, "keep","discard")) %>%
    filter(nrow == "keep")
  plot <- ggplot(filtered.fs.Table,
                 aes(x=reorder(variable, scaled_importance),
                     y=scaled_importance)) +
    geom_bar(stat = "identity", color="black", alpha=0.5,
             fill = ifelse(filtered.fs.Table$scaled_importance < 0.05, "#E41A1C", "#377EB8")) +
    theme_bw() +
    theme(
      axis.line = element_line(size = 0.75),
      axis.text = element_text(
        size = 11,
        face = "bold",
        colour = "black"
      ),
      axis.title = element_text(size = 12, face = "bold")
    ) +
    labs(x="", y="scaled importance",
         title= paste("variable importance plot: top", length(filtered.fs.Table$nrow), "variable")) +
    coord_flip()

  ##----------------------------------------------------------------
  ##                  Enlist variable importance                  --
  ##----------------------------------------------------------------
  ## save table
  write.csv(filtered.fs.Table, paste0(save.location, "filtered.fs.Table.csv"))
  ##---------------------------------------------------------------
  ##                       save all models                       --
  ##---------------------------------------------------------------
  ## create directory
  all.model <- paste(save.location,"all_models", sep = "")
  lapply(all.model, function(x) if (!dir.exists(x)) dir.create(x))
  ## save models
  for(n in 1:nrow(lb)) {
    modelAutoML <- h2o.getModel(model.fs@leaderboard[n, "model_id"])
    h2o.saveModel(object = modelAutoML, path = all.model)
  }
  ## performance
  performance <- h2o.performance(fs.Model, newdata = test)
  cm <- as.data.frame(h2o.confusionMatrix(performance, cutoff))
  ### performance of model after 10-fold cross-validation
  ## empty dataset
  cv.results <- data.frame()
  ## define folds
  folds <- cut(seq(1, nrow(test)), breaks = 10, labels = FALSE)
  ##----------------------------------------------------------------
  ##            erformance on 10 fold cross validation            --
  ##----------------------------------------------------------------
  for (i in 1:10) {
    ## segement data by fold using the which() function
    Indexes <- which(folds == i, arr.ind = TRUE)
    test.cv <- test[-Indexes,]
    perf.cv <- h2o.performance(fs.Model, newdata = test.cv)
    ## enlist parametes
    cv.results[i, "fold"] <- i
    cv.results[i, "MSE"] <- perf.cv@metrics$MSE
    cv.results[i, "RMSE"] <- perf.cv@metrics$RMSE
    cv.results[i, "R2"] <- perf.cv@metrics$r2
    cv.results[i, "logloss"] <- perf.cv@metrics$logloss
    cv.results[i, "AUC"] <- perf.cv@metrics$AUC
    cv.results[i, "PRAUC"] <- perf.cv@metrics$pr_auc
    cv.results[i, "Gini"] <- perf.cv@metrics$Gini
    cv.results[i, "Mean_per_class_error"] <- perf.cv@metrics$mean_per_class_error
  }
  ## summarize datasets
  cv.results.summary <-
    cv.results[, !colnames(cv.results) %in% "fold"] %>%
    gather(factor_key = TRUE) %>%
    group_by(key) %>%
    summarise(
      mean = mean(value),
      sd = sd(value),
      max = max(value),
      min = min(value)
    )
  ## save CV data
  write.csv(cv.results.summary, paste0(save.location,
                                       "cvResultsSummary.csv"))
  return(list(
    result=automl,
    leaderboard =leaderboard,
    varImpTable = filtered.fs.Table,
    plot=plot,
    model=fs.Model,
    performance=performance,
    confusionMatrix=cm,
    cv.results=cv.results.summary,
    cutoff = cutoff
  ))
}


##----------------------------------------------------------------
##                       iterative model                       --
##----------------------------------------------------------------

runIterativeModel <- function(model, train, test, features, response, time, save.location,cutoff_metric, featureSelectionThreshold = 0.05) {
  automl <- data.frame()
  loop <- 1
  algoInclude <- gsub("(_).*","", model@model_id)
  factor.models <- c("XGBoost", "DeepLearning", "GLM")
  ## loop
  repeat {
    model.fs <- h2o.automl(
      x = features,
      y = response,
      training_frame = train,
      leaderboard_frame = test,
      keep_cross_validation_models = TRUE,
      keep_cross_validation_predictions = TRUE,
      keep_cross_validation_fold_assignment = TRUE,
      balance_classes = TRUE,
      max_runtime_secs = time,
      include_algos = algoInclude,
      seed = 1234,
      sort_metric = "logloss",
      nfolds = 10
    )
    ## model selection
    fs.Model <- h2o.getModel(model.fs@leader@model_id)
    ## if stackensemble
    if (grepl("StackedEnsemble", model.fs@leader@model_id) == TRUE) {
      perf <- h2o.performance(fs.Model, newdata = test)
      ensemble_auc_test <- h2o.auc(perf)
      print(sprintf("Ensemble Test AUC:  %s", ensemble_auc_test))
      ##---------------------------------------------------------------
      ##                         metalearner                         --
      ##---------------------------------------------------------------
      metalearner <-
        h2o.getModel(model.fs@leader@model$metalearner$name)
      modelImp <- h2o.varimp(metalearner)
      highestImpName <- modelImp[1, 1]
      ## model
      model  <-
        h2o.getModel(model.fs@leader@model$metalearner_model@model$names[highestImpName])
      ## vip table
      fs.Table <- as.data.frame(h2o.varimp(model))
      ## if not stackedensemble but classifier model
      if (grepl(paste(factor.models, collapse = "|"), model@model_id) == TRUE) {
        fs.Table$variable <- gsub("^(.*)[.].*", "\\1", fs.Table$variable)
      }
      filtered.fs.Table <-
        fs.Table[fs.Table$scaled_importance > featureSelectionThreshold, ]
      filtered.fs.Table <- filtered.fs.Table %>%
        distinct(variable, .keep_all = TRUE) %>%
        arrange(desc(scaled_importance))
      ## if not stackedensemble and classifier model
    } else  if (grepl(paste(factor.models, collapse = "|"),
                      model.fs@leader@model_id) == TRUE) {
      fs.Table <- as.data.frame(h2o.varimp(fs.Model))
      fs.Table$variable <-
        gsub("^(.*)[.].*", "\\1", fs.Table$variable)
      filtered.fs.Table <-
        fs.Table[fs.Table$scaled_importance > featureSelectionThreshold, ]
      filtered.fs.Table <- filtered.fs.Table %>%
        distinct(variable, .keep_all = TRUE) %>%
        arrange(desc(scaled_importance))
    } else {
      fs.Table <- as.data.frame(h2o.varimp(fs.Model))
      filtered.fs.Table <-
        fs.Table[fs.Table$scaled_importance > featureSelectionThreshold, ]
      filtered.fs.Table <- filtered.fs.Table %>%
        distinct(variable, .keep_all = TRUE) %>%
        arrange(desc(scaled_importance))
    }
    ## define features
    features <-
      features[features %in% filtered.fs.Table$variable]
    ## model metric
    logloss  <- h2o.logloss(h2o.performance(fs.Model, newdata = test))
    cutoff <-
      h2o.find_threshold_by_max_metric(h2o.performance(fs.Model, test), cutoff_metric)
    accuracy <-
      h2o.accuracy(h2o.performance(fs.Model, newdata = test), cutoff)
    auc <- h2o.auc(h2o.performance(fs.Model, newdata = test))
    ## model metric table
    automl[loop, "model"] <-  model.fs@leader@model_id
    automl[loop, "iteration.number"] <-  loop
    automl[loop, 'logloss'] <-  logloss
    automl[loop, "number.of.features"] <-  length(features)
    automl[loop, "auc"] <-
      h2o.auc(h2o.performance(fs.Model, newdata = test))
    automl[loop, "accuracy"] <-
      h2o.accuracy(h2o.performance(fs.Model, newdata = test), cutoff)
    automl[loop, "specificty"] <-
      h2o.specificity(h2o.performance(fs.Model, newdata = test), cutoff)
    automl[loop, "sensitivity"] <-
      h2o.sensitivity(h2o.performance(fs.Model, newdata = test), cutoff)
    automl[loop, "cutoff"] <-
      cutoff
    ## save model
    h2o.saveModel(fs.Model, path = save.location, force = TRUE)
    ## break loop
    if ((automl[loop, "auc"] > 0.98) ||
        loop > 2 && (length(features) ==
                     ifelse((automl[loop-1, "number.of.features"] == automl[loop-2, "number.of.features"])==TRUE,
                            automl[loop-1, "number.of.features"], automl[loop, 4]) ) == TRUE )
    {
      break
    }
    loop <- loop + 1

  }
  ##---------------------------------------------------------------
  ##                   variable importance plot                   --
  ##---------------------------------------------------------------
  filtered.fs.Table <- filtered.fs.Table %>%
    mutate(nrow = ifelse(n() < 20, "keep","discard")) %>%
    filter(nrow == "keep")
  plot <- ggplot(filtered.fs.Table,
                 aes(x=reorder(variable, scaled_importance),
                     y=scaled_importance)) +
    geom_bar(stat = "identity", color="black", alpha=0.5,
             fill = ifelse(filtered.fs.Table$scaled_importance < 0.05, "#E41A1C", "#377EB8")) +
    theme_bw() +
    theme(
      axis.line = element_line(size = 0.75),
      axis.text = element_text(
        size = 11,
        face = "bold",
        colour = "black"
      ),
      axis.title = element_text(size = 12, face = "bold")
    ) +
    labs(x="", y="scaled importance",
         title= paste("variable importance plot: top", length(filtered.fs.Table$nrow), "variable")) +
    coord_flip()

  ##----------------------------------------------------------------
  ##                  Enlist variable importance                  --
  ##----------------------------------------------------------------
  ## save table
  write.csv(filtered.fs.Table, paste0(save.location, "filtered.fs.Table.csv"))
  ## performance
  performance <- h2o.performance(fs.Model, newdata = test)
  cm <- as.data.frame(h2o.confusionMatrix(performance))
  ### performance of model after 10-fold cross-validation
  ## empty dataset
  cv.results <- data.frame()
  ## define folds
  folds <- cut(seq(1, nrow(test)), breaks = 10, labels = FALSE)
  ##----------------------------------------------------------------
  ##            performance on 10 fold cross validation            --
  ##----------------------------------------------------------------
  for (i in 1:10) {
    ## segement data by fold using the which() function
    Indexes <- which(folds == i, arr.ind = TRUE)
    test.cv <- test[-Indexes,]
    perf.cv <- h2o.performance(fs.Model, newdata = test.cv)
    ## enlist parametes
    cv.results[i, "fold"] <- i
    cv.results[i, "MSE"] <- perf.cv@metrics$MSE
    cv.results[i, "RMSE"] <- perf.cv@metrics$RMSE
    cv.results[i, "R2"] <- perf.cv@metrics$r2
    cv.results[i, "logloss"] <- perf.cv@metrics$logloss
    cv.results[i, "AUC"] <- perf.cv@metrics$AUC
    cv.results[i, "PRAUC"] <- perf.cv@metrics$pr_auc
    cv.results[i, "Gini"] <- perf.cv@metrics$Gini
    cv.results[i, "Mean_per_class_error"] <- perf.cv@metrics$mean_per_class_error
  }
  ## summarize datasets
  cv.results.summary <-
    cv.results[, !colnames(cv.results) %in% "fold"] %>%
    gather(factor_key = TRUE) %>%
    group_by(key) %>%
    summarise(
      mean = mean(value),
      sd = sd(value),
      max = max(value),
      min = min(value)
    )
  ## save CV data
  write.csv(cv.results.summary, paste0(save.location,
                                       "cvResultsSummary.csv"))
  return(list(
    result = automl,
    varImpTable = filtered.fs.Table,
    plot=plot,
    model=fs.Model,
    performance=performance,
    confusionMatrix=cm,
    cv.results=cv.results.summary,
    cutoff =cutoff
  ))
}

##----------------------------------------------------------------
##                       feature reduction                      --
##----------------------------------------------------------------

runFeatureReduction <- function(model, train, test, response, time, save.location, cutoff_metric, featureSelectionThreshold = 0.05) {
  automl <- data.frame()
  loop <- 1
  algoInclude <- gsub("(_).*","", model@model_id)
  featuresToUse <- as.data.frame(model@model$variable_importances)
  featuresToUse <- featuresToUse %>%
    filter(scaled_importance > featureSelectionThreshold)
  featuresToUse <- rev(featuresToUse$variable)


  for (f in featuresToUse) {
    if (length(featuresToUse) != 1) {
      featuresToUse <- featuresToUse[!featuresToUse %in% f]
    }
    featuresToUse <- featuresToUse
    model.fs <- h2o.automl(
      x = featuresToUse,
      y = response,
      training_frame = train,
      leaderboard_frame = test,
      keep_cross_validation_models = TRUE,
      keep_cross_validation_predictions = TRUE,
      keep_cross_validation_fold_assignment = TRUE,
      balance_classes = TRUE,
      max_runtime_secs = time,
      include_algos = algoInclude,
      seed = 1234,
      sort_metric = "logloss",
      nfolds = 10
    )
    ## model selection
    fs.Model <- h2o.getModel(model.fs@leader@model_id)
      ## model metric
      logloss  <- h2o.logloss(h2o.performance(fs.Model, newdata = test))
      cutoff <-
        h2o.find_threshold_by_max_metric(h2o.performance(fs.Model, test), cutoff_metric)
      accuracy <-
        h2o.accuracy(h2o.performance(fs.Model, newdata = test), cutoff)
      auc <- h2o.auc(h2o.performance(fs.Model, newdata = test))
      ## model metric table
      automl[loop, "model"] <-  model.fs@leader@model_id
      automl[loop, "metabolite"] <-  f
      automl[loop, 'logloss'] <-  logloss
      automl[loop, "number.of.features"] <-  length(featuresToUse[!featuresToUse %in% f])
      automl[loop, "auc"] <-
        h2o.auc(h2o.performance(fs.Model, newdata = test))
      automl[loop, "accuracy"] <-
        h2o.accuracy(h2o.performance(fs.Model, newdata = test), cutoff)
      automl[loop, "specificty"] <-
        h2o.specificity(h2o.performance(fs.Model, newdata = test), cutoff)
      automl[loop, "sensitivity"] <-
        h2o.sensitivity(h2o.performance(fs.Model, newdata = test), cutoff)
      ## save model
      h2o.saveModel(fs.Model, path = save.location, force = TRUE)
      loop <- loop + 1
      featuresToUse <- featuresToUse

  }
  ## plot
  dfdata <- automl %>%
    select(metabolite, auc, number.of.features) %>%
    gather("variable", "value", -metabolite, -number.of.features)

  plot <- ggplot(dfdata, aes(x=as.factor(number.of.features),
               y=value,
               group=variable,
               fill=variable),
           color="black") +
    geom_point(shape = 21,
               color = "black",
               size = 4,
               show.legend = FALSE) +
    geom_smooth(method = "lm", formula = y ~ x + I(x^2),
                aes(color = variable),
                se = TRUE,
                show.legend = FALSE,
                size=2,
                alpha=0.2) +
    theme_bw() +
    theme(
      axis.line = element_line(size = 0.75),
      axis.text = element_text(
        size = 11,
        face = "bold",
        colour = "black"
      ),
      axis.title = element_text(size = 12, face = "bold")
    ) +
    labs(x="feature reduction",
         y="Area under curve (AUC)",
         title= "accuracy following feature reduction") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_color_manual(values = RColorBrewer::brewer.pal(2, "Set1")) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(2, "Set1")) +
    ggrepel::geom_label_repel(aes(label = metabolite),
                              colour = "black",
                              fontface = "bold",
                              show.legend = FALSE,
                              min.segment.length = 0,
                              seed = 42,
                              box.padding = 0.5,
                              alpha=0.5)


  return(list(plot=plot,
              data=automl,
              cutoff=cutoff))

}

##----------------------------------------------------------------
##                          leave one out                       --
##----------------------------------------------------------------

runLeaveOneOut <- function(model, train, test, response, time, save.location, cutoff_metric, featureSelectionThreshold = 0.05) {
  automl <- data.frame()
  loop <- 1
  algoInclude <- gsub("(_).*","", model@model_id)
  featuresToUse <- as.data.frame(model@model$variable_importances)
  featuresToUse <- featuresToUse %>%
    filter(scaled_importance > featureSelectionThreshold)
  featuresToUse <- rev(featuresToUse$variable)
  for (f in featuresToUse) {
    model.fs <- h2o.automl(
      x = featuresToUse[!featuresToUse %in% f],
      y = response,
      training_frame = train,
      leaderboard_frame = test,
      keep_cross_validation_models = TRUE,
      keep_cross_validation_predictions = TRUE,
      keep_cross_validation_fold_assignment = TRUE,
      balance_classes = TRUE,
      max_runtime_secs = time,
      include_algos = algoInclude,
      seed = 1234,
      sort_metric = "logloss",
      nfolds = 10
    )
    ## model selection
    fs.Model <- h2o.getModel(model.fs@leader@model_id)
    ## vip table
    fs.Table <- as.data.frame(h2o.varimp(model))
    ## model metric
    logloss  <- h2o.logloss(h2o.performance(fs.Model, newdata = test))
    cutoff <-
      h2o.find_threshold_by_max_metric(h2o.performance(fs.Model, test), cutoff_metric)
    accuracy <-
      h2o.accuracy(h2o.performance(fs.Model, newdata = test), cutoff)
    auc <- h2o.auc(h2o.performance(fs.Model, newdata = test))
    ## model metric table
    automl[loop, "model"] <-  model.fs@leader@model_id
    automl[loop, "metabolite"] <-  featuresToUse[featuresToUse %in% f]
    automl[loop, 'logloss'] <-  logloss
    automl[loop, "number.of.features"] <-  length(featuresToUse[!featuresToUse %in% f])
    automl[loop, "auc"] <-
      h2o.auc(h2o.performance(fs.Model, newdata = test))
    automl[loop, "accuracy"] <-
      h2o.accuracy(h2o.performance(fs.Model, newdata = test), cutoff)
    automl[loop, "specificty"] <-
      h2o.specificity(h2o.performance(fs.Model, newdata = test), cutoff)
    automl[loop, "sensitivity"] <-
      h2o.sensitivity(h2o.performance(fs.Model, newdata = test), cutoff)
    ## save model
    h2o.saveModel(fs.Model, path = save.location, force = TRUE)
    loop <- loop + 1
  }
  ## plot
  plot <- automl %>%
    select(metabolite, auc, accuracy) %>%
    gather("variable", "value", -metabolite) %>%
    ggplot(aes(x=reorder(metabolite, -value),
               y=value,
               group=variable,
               fill=variable)) +
    geom_point(shape = 21,
               color = "black",
               size = 4) +
    geom_smooth(method = "lm", formula = y ~ x + I(x^2),
                aes(color = variable),
                se = FALSE,
                show.legend = FALSE,
                size=2) +
    theme_bw() +
    theme(
      axis.line = element_line(size = 0.75),
      axis.text = element_text(
        size = 11,
        face = "bold",
        colour = "black"
      ),
      axis.title = element_text(size = 12, face = "bold")
    ) +
    labs(x="feature left out",
         y="",
         title= "accuracy following dropout of the feature") +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_color_manual(values = RColorBrewer::brewer.pal(2, "Set1")) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(2, "Set1"))

  return(list(plot=plot,
              data=automl,
              cutoff=cutoff))

}

modelOdds <- function(model, test, levels) {

  model.x <- model@parameters$x
  ## training set AUC +/- sd
  auc <- model@model$training_metrics@metrics$AUC
  n_p <- model@model$training_metrics@metrics$cm$table[3,2]
  n_n <- model@model$training_metrics@metrics$cm$table[3,1]
  print(auc)
  print(se_auc(auc, n_p, n_n))


  cm <- as.data.frame(h2o.confusionMatrix(h2o.performance(model, test)))

  lvs <- sort(levels)
  truth <-
    factor(rep(lvs, times = c(cm[[lvs[1]]][1] + cm[[lvs[2]]][1], cm[[lvs[1]]][2] + cm[[lvs[2]]][2])),
           levels = rev(lvs))
  pred <- factor(c(rep(lvs, times = c(cm[[lvs[1]]][1], cm[[lvs[2]]][1])),
                   rep(lvs, times = c(cm[[lvs[1]]][2], cm[[lvs[2]]][2]))),
                 levels = rev(lvs))

  xtab <- table(pred, truth)
  or <- oddsratio(xtab)
  return(or)
}


plotVarImp <- function (model) {
  ##---------------------------------------------------------------
  ##                   variable importance plot                   --
  ##---------------------------------------------------------------
  df <- as.data.frame(h2o.varimp(model))

  plot <- df %>%
    mutate(nrow = ifelse(n() < 20, "keep","discard")) %>%
    filter(nrow == "keep")
  plot <- ggplot(plot,
                 aes(x=reorder(variable, scaled_importance),
                     y=scaled_importance)) +
    geom_bar(stat = "identity", color="black", alpha=0.5,
             fill = ifelse(df$scaled_importance < 0.05, "#E41A1C", "#377EB8")) +
    theme_bw() +
    theme(
      axis.line = element_line(size = 0.75),
      axis.text = element_text(
        size = 11,
        face = "bold",
        colour = "black"
      ),
      axis.title = element_text(size = 12, face = "bold")
    ) +
    labs(x="", y="scaled importance",
         title= paste("variable importance plot: top", length(df$nrow), "variable")) +
    coord_flip()
}
