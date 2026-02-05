# function to train random forest

train_rf<-function(Xtrain,Ytrain,n_tree=10,k_cv=5,method_control="cv",n_cores=1){
  repGrid <- expand.grid(.mtry = c(10,20,30,40,50,60,70))
  set.seed(123)
  # create the cluster for caret to use
  cl <- makePSOCKcluster(n_cores)
  registerDoParallel(cl)

  if(method_control=="cv"){
    model_out <- train(
      x=Xtrain, y=Ytrain, method = "rf",metric = "Accuracy",ntree=n_tree,preProcess = c("center", "scale"),
      trControl = trainControl(method = "cv",number = k_cv,trim = T,returnData = F),
      tuneGrid = repGrid
    )
  }else if(method_control=="oob"){
    model_out <- train(
      x=Xtrain, y=Ytrain, method = "rf",metric = "Accuracy",ntree=n_tree,preProcess = c("center", "scale"),
      trControl = trainControl(method = "oob", trim = T,returnData = F,classProbs = T),
      tuneGrid = repGrid
    )
  }
  stopCluster(cl)
  registerDoSEQ()
  return(model_out)

}

# function to train logistic regression

train_log_reg<-function(Xtrain,Ytrain,k_cv=5,n_cores=1){
  repGrid <- expand.grid(decay  = c(0.1,0.01,0.001))
  set.seed(123)
  # create the cluster for caret to use
  cl <- makePSOCKcluster(n_cores)
  registerDoParallel(cl)
  model_out <- train(
    x=Xtrain, y=Ytrain, method = "multinom",metric="Accuracy", preProcess = c("center", "scale"),
    trControl = trainControl(method="cv", number=k_cv,
                             trim = T,returnData = F,classProbs = T),
    tuneGrid = repGrid,
    MaxNWts = 10000
    )
  
  stopCluster(cl)
  registerDoSEQ()
  return(model_out)
  
}

# function to train knn

train_knn<-function(Xtrain,Ytrain,k_cv=5,n_cores=1){
  repGrid <- expand.grid(k  = c(3, 5, 10, 15, 20))
  set.seed(123)
  # create the cluster for caret to use
  cl <- makePSOCKcluster(n_cores)
  registerDoParallel(cl)
  
  model_out <- train(
    x=Xtrain, y=Ytrain, method = "knn",metric="Accuracy", preProcess = c("center", "scale"),
    trControl = trainControl(method="cv", number=k_cv,
                             trim = T,returnData = F,classProbs = T),
    tuneGrid = repGrid
  )
  
  stopCluster(cl)
  registerDoSEQ()
  return(model_out)
  
}

# function to train fnn

train_fnn<-function(Xtrain,Ytrain,k_cv=5,n_cores=1){
  repGrid <- expand.grid(size = c(1, 3, 5),  # Number of hidden units (1, 3, 5)
                           decay = c(0, 0.1, 0.5))  # Regularization parameter
  set.seed(123)
  # create the cluster for caret to use
  cl <- makePSOCKcluster(n_cores)
  registerDoParallel(cl)
  
  model_out <- train(
    x=Xtrain, y=Ytrain, method = "nnet",metric="Accuracy", preProcess = c("center", "scale"),
    trControl = trainControl(method="cv", number=k_cv,
                             trim = T,returnData = F,classProbs = T),
    tuneGrid = repGrid
  )
  
  stopCluster(cl)
  registerDoSEQ()
  return(model_out)
  
}

# function to train svm

train_svm<-function(Xtrain,Ytrain,k_cv=5,n_cores=1){
  repGrid <- expand.grid(C = c(0.1, 1, 10))
  set.seed(123)
  # create the cluster for caret to use
  cl <- makePSOCKcluster(n_cores)
  registerDoParallel(cl)
  
  model_out <- train(
    x=Xtrain, y=Ytrain, method = "svmLinear",metric="Accuracy", preProcess = c("center", "scale"),
    trControl = trainControl(method="cv", number=k_cv,
                             trim = T,returnData = F,classProbs = T),
    tuneGrid = repGrid
  )
  
  stopCluster(cl)
  registerDoSEQ()
  return(model_out)
  
}

# function to train nb

train_nb<-function(Xtrain,Ytrain,k_cv=5,n_cores=1){
  repGrid <- expand.grid(laplace = c(0, 1),
                         usekernel = c(TRUE, FALSE),adjust=c(1,2))  
  set.seed(123)
  # create the cluster for caret to use
  cl <- makePSOCKcluster(n_cores)
  registerDoParallel(cl)
  
  model_out <- train(
    x=Xtrain, y=Ytrain, method = "naive_bayes",metric="Accuracy", preProcess = c("center", "scale"),
    trControl = trainControl(method="cv", number=k_cv,
                             trim = T,returnData = F,classProbs = T),
    tuneGrid = repGrid
  )
  
  stopCluster(cl)
  registerDoSEQ()
  return(model_out)
  
}

# function to train one class svm

train_oneclass_svm<-function(Xtrain,Ytrain,k_cv=5,n_cores=1,ref_label="algae"){
  
  pre_proc <- preProcess(Xtrain, method = c("center", "scale"))
  Xtrain <- predict(pre_proc, Xtrain)
  


  # cross-validation
  folds <- createFolds(Ytrain, k = k_cv, list = TRUE, returnTrain = TRUE)
  vec_scores<-c()
  list_cross_val_model<-list()
  nu_val<-0
  for(i in 1:k_cv){
    nu_val<-nu_val+0.1
    # Create training and validation data for this fold
    train_fold <- Xtrain[folds[[i]], ]
    Y_truth_train<- Ytrain[folds[[i]]]
    
    valid_fold <- Xtrain[-folds[[i]], ]
    Y_truth_val<- Ytrain[-folds[[i]]]
    # Train a model
    inds<-which(Y_truth_train==ref_label) # select only normal data (indicated by ref label)
    train_fold_normal<-train_fold[inds,]
    
    set.seed(123)
    cl <- makePSOCKcluster(n_cores)
    registerDoParallel(cl)
    
    model_out_i <- svm(
      train_fold_normal,
      type = "one-classification",
      kernel = "radial",  # Radial basis function kernel
      nu = nu_val         # Fraction of outliers allowed
    )
    
    stopCluster(cl)
    registerDoSEQ()
    
    # Make predictions on the validation set
    y_val <- predict(model_out_i, newdata = valid_fold)
    inds<-which(y_val==T)
    y_val[inds]<-ref_label # normal
    all_lab<-unique(Ytrain)
    ind_ref<-which(all_lab %in% ref_label)
    y_val[-inds]<-all_lab[-ind_ref] # outlier
    y_val<-as.vector(y_val)
    # get accuracy on validation set
    score<-get_accuracy(y_truth = Y_truth_val,y_pred = y_val,pos_label = ref_label)
    list_cross_val_model[[as.character(i)]]<-model_out_i
    vec_scores<-c(vec_scores,score)
  }
  model_out<-list_cross_val_model
  return(list(model_out=model_out,vec_scores=vec_scores))

  
}


# function to get gmm model

train_gmm<-function(Xtrain,Ytrain,k_cv=5,n_cores=1){
  
  pre_proc <- preProcess(Xtrain, method = c("center", "scale"))
  test_df_non_char <- predict(pre_proc, Xtrain)
  
  # by default number of components is estimated through BIC
  # cross-validation
  folds <- createFolds(Ytrain, k = k_cv, list = TRUE, returnTrain = TRUE)
  vec_scores<-c()
  list_cross_val_model<-list()
  
  for(i in 1:k_cv){
    # Create training and validation data for this fold
    train_fold <- Xtrain[folds[[i]], ]
    valid_fold <- Xtrain[-folds[[i]], ]
    Y_truth_val<- Ytrain[-folds[[i]]]
    # Train a model
    set.seed(123)
    cl <- makePSOCKcluster(n_cores)
    registerDoParallel(cl)
    
    model_out_i <- Mclust(train_fold)
    
    stopCluster(cl)
    registerDoSEQ()
    # Make predictions on the validation set
    y_val <- predict(model_out_i, newdata = valid_fold)
    y_val<-y_val$classification
    score<-adjustedRandIndex(x = y_val,y = Y_truth_val)
    list_cross_val_model[[as.character(i)]]<-model_out_i
    vec_scores<-c(vec_scores,score)
  }
  model_out<-list_cross_val_model
  return(list(model_out=model_out,vec_scores=vec_scores))
}


# function to perform isolation forest
train_isoforest<-function(Xtrain,Ytrain,k_cv=5,n_cores=1,ref_label="algae"){
  
  pre_proc <- preProcess(Xtrain, method = c("center", "scale"))
  Xtrain <- predict(pre_proc, Xtrain)
  
  
  
  # cross-validation
  folds <- createFolds(Ytrain, k = k_cv, list = TRUE, returnTrain = TRUE)
  vec_scores<-c()
  list_cross_val_model<-list()
  cont_val<-0
  for(i in 1:k_cv){
    cont_val<-cont_val+0.1
    # Create training and validation data for this fold
    train_fold <- Xtrain[folds[[i]], ]
    Y_truth_train<- Ytrain[folds[[i]]]
    
    valid_fold <- Xtrain[-folds[[i]], ]
    Y_truth_val<- Ytrain[-folds[[i]]]
    # Train a model
    inds<-which(Y_truth_train==ref_label) # select only normal data (indicated by ref label)
    train_fold_normal<-train_fold[inds,]
    
    set.seed(123)
    cl <- makePSOCKcluster(n_cores)
    registerDoParallel(cl)
    model_out_i <- iso_forest_train(Xtrain = train_fold_normal,cont_val = cont_val,random_state = as.integer(42))
    stopCluster(cl)
    registerDoSEQ()

    # Make predictions on the validation set
    y_val <- iso_forest_predict(isof_model=model_out_i, test_df = valid_fold)
    inds<-which(y_val==T)
    y_val[inds]<-ref_label # normal
    all_lab<-unique(Ytrain)
    ind_ref<-which(all_lab %in% ref_label)
    y_val[-inds]<-all_lab[-ind_ref] # outlier
    y_val<-as.vector(y_val)
    # get accuracy on validation set
    score<-get_accuracy(y_truth = Y_truth_val,y_pred = y_val,pos_label = ref_label)
    list_cross_val_model[[as.character(i)]]<-model_out_i
    vec_scores<-c(vec_scores,score)
  }
  model_out<-list_cross_val_model
  return(list(model_out=model_out,vec_scores=vec_scores))
}





# function to train data using different models

data_training<-function(data, method="rf",n_cores=1,n_samples_down=2000,method_ctrl="cv",n_tree=10, n_seed=123,k_cv=5,return_data=F,
                        var_down="Species",ref_var="type",ref_label="algae"){
  start<-Sys.time()
  # ---- downsampling
  print("downsampling.....")
  set.seed(n_seed)
  data_algae<-data[data$type=="algae",]
  data_non_algae<-data[data$type=="non_algae",]
  if(n_samples_down!="none"){
    data_final_algae<-get_downsample(df = data_algae,var = var_down,n_samples = n_samples_down)
    data_final_non_algae<-get_downsample(df = data_non_algae,var = var_down,n_samples = n_samples_down)
  }else{
    data_final_algae<-data_algae
    data_final_non_algae<-data_non_algae
  }
  data_final<-rbind(data_final_algae,data_final_non_algae)
  print("Total events:")
  print(nrow(data_final))
  print("Number of events per Species:")
  print(table(data_final$Species))
  print("Number of events per gate:")
  print(table(data_final$type))
  if(return_data==T){
    return(data_final)
  }
  #----- prepare data: get X_train, Y_train
  Y_train<-data_final[[ref_var]]
  
  
  # remove characters columns if needed
  print("check characters columns...")
  vec_type_col<-sapply(data_final,function(x){
    return(class(x))
  })
  inds_rm<-grep("character",vec_type_col)
  if(length(inds_rm)!=0){
    X_train<-data_final[,-inds_rm, with=F]
  }else{
    X_train<-data_final # test data already without characters columns
  }
  
  #------ perform the training
  set.seed(123)
  if(method=="rf"){
    print("training using rf model...")
    if(length(n_tree)==1){
      model_trained<-train_rf(Xtrain = X_train,Ytrain = Y_train, n_tree = n_tree,method_control = method_ctrl,n_cores=n_cores,k_cv=k_cv)
    }else{
      list_all_models<-list()
      for(n_tree_val in n_tree){
        print(sprintf("n_tree:%s",n_tree_val))
        set.seed(123)
        model_trained_tree_val<-train_rf(Xtrain = X_train,Ytrain = Y_train, n_tree = n_tree_val,method_control = method_ctrl,n_cores=n_cores,k_cv=k_cv)
        list_all_models[[as.character(n_tree_val)]]<-model_trained_tree_val
      }
      end<-Sys.time()
      time_taken<-end-start
      print("Execution time:")
      print(time_taken)
      print("Done")
      return(list_all_models)
    }
  }else if(method=="log_reg"){

    print("training using log reg model...")
    model_trained<-train_log_reg(Xtrain = X_train,Ytrain = Y_train,n_cores=n_cores,k_cv=k_cv)
    
  }else if(method=="knn"){
    print("training using knn model...")
    model_trained<-train_knn(Xtrain = X_train,Ytrain = Y_train,n_cores=n_cores,k_cv=k_cv)
    
  }else if(method=="fnn"){
    print("training using fnn model...")
    model_trained<-train_fnn(Xtrain = X_train,Ytrain = Y_train,n_cores=n_cores,k_cv=k_cv)
    
  }else if(method=="svm"){
    print("training using svm model...")
    model_trained<-train_svm(Xtrain = X_train,Ytrain = Y_train,n_cores=n_cores,k_cv=k_cv)
    
  }else if(method=="nb"){
    print("training using nb model...")
    model_trained<-train_nb(Xtrain = X_train,Ytrain = Y_train,n_cores=n_cores,k_cv=k_cv)
    
  }else if(method=="one_class_svm"){

    print("training using one class svm model...")
    model_trained<-train_oneclass_svm(Xtrain = X_train,Ytrain = Y_train,n_cores=n_cores,k_cv=k_cv)
    
  }else if(method=="gmm"){
    print("training using gmm model...")
    model_trained<-train_gmm(Xtrain = X_train,Ytrain = Y_train,k_cv=k_cv,n_cores=n_cores) 
    
  }else if(method=="isoforest"){
    print("training using isoforest model...")
    model_trained<-train_isoforest(Xtrain = X_train,Ytrain = Y_train,k_cv=k_cv,n_cores=n_cores) 
    
  }
  end<-Sys.time()
  time_taken<-end-start
  print("Execution time:")
  print(time_taken)
  print("Done")
  return(model_trained)
}


# function to train data in multiple simulations

data_train<-function(data, method="rf",n_cores=1,n_samples_down=2000,method_ctrl="cv",n_tree=10, rep_seed=123,k_cv=5,n_rep=6,
                     var_down="Species",ref_label="algae",ref_var="type"){
  list_all_models<-list()
  for(i in 1:n_rep){
    print(sprintf("####### Cross-validation repetition %d #########",i))
    out_model_rep_i<-data_training(data = data,method=method,n_cores = n_cores,n_samples_down = n_samples_down,method_ctrl = method_ctrl,n_tree=n_tree,
                                   n_seed = rep_seed+i,k_cv=k_cv,return_data=F,var_down=var_down,ref_label=ref_label,ref_var=ref_var)
    list_all_models[[as.character(i)]]<-out_model_rep_i
  }
  print("Get best model...")
  if(method=="gmm" || method=="one_class_svm" || method=="isoforest"){
    out_best_model<-get_best_model_v2(list_models = list_all_models)
  }else{
    out_best_model<-get_best_model(list_models = list_all_models)
  }

  vec_max_values<-round(out_best_model$vec_max_values,3)
  
  max_value<-max(vec_max_values)
  
  best_model<-out_best_model$best_model
  
  return(list(vec_max_values=vec_max_values,max_value=max_value,best_model=best_model))
}


