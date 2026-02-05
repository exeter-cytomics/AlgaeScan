# if positive=NULL, the F1_score function will use the first factor level (as sorted alphabetically).
# function to perform prediction based on model and get F1 score
data_prediction<-function(input_model,test_df,normalize=F,prob=F){
  start<-Sys.time()
  #----- remove characters columns if needed
  print("check characters columns...")
  vec_type_col<-sapply(test_df,function(x){
    return(class(x))
  })
  inds_rm<-grep("character",vec_type_col)
  if(length(inds_rm)!=0){
    test_df_non_char<-test_df[,-inds_rm, with=F]
  }else{
    test_df_non_char<-test_df # test data already without characters columns
  }
  if(normalize==T){
    pre_proc <- preProcess(test_df_non_char, method = c("center", "scale"))
    test_df_non_char <- predict(pre_proc, test_df_non_char)
    
  }
  #----- prediction step
  print("prediction using provided trained model")
  if(prob==F){
    if("train" %in% class(input_model)){
      y_pred<-predict(input_model,test_df_non_char)
    }else if("Mclust" %in% class(input_model)){
      y_pred<-predict(input_model,test_df_non_char)
      y_pred<-y_pred$classification
    }else if("svm" %in% class(input_model)){
      y_pred<-predict(input_model,test_df_non_char)
      inds<-which(y_pred==T)
      y_pred[inds]<-1 # normal
      y_pred[-inds]<-0 # outlier
      y_pred<-as.vector(y_pred)
    }else if("sklearn.ensemble._iforest.IsolationForest" %in% class(input_model)){
      y_pred <- iso_forest_predict(isof_model=input_model, test_df = test_df_non_char)
      inds<-which(y_pred==1)
      y_pred[inds]<-1 # normal
      y_pred[-inds]<-0 # outlier
      y_pred<-as.vector(y_pred)
    }
    y_pred<-as.character(y_pred)
  }else{
    if("sklearn.ensemble._iforest.IsolationForest" %in% class(input_model)){
      y_pred<-iso_forest_predict_prob(isof_model=input_model, test_df = test_df_non_char)
    }else{
      y_pred<-predict(input_model,test_df_non_char,type="prob")
      # y_pred is a dataframe of probabilities
      y_pred<-as.data.table(y_pred)
      colnames(y_pred)<-paste(colnames(y_pred),"prob",sep="_")
    }
 
  }

  

  end<-Sys.time()
  time_taken<-end-start
  print("Execution time:")
  print(time_taken)
  print("Done")
  return(y_pred)
}
