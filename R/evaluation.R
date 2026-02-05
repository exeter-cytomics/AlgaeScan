# function to get F1 score
get_F1_score<-function(y_truth,y_pred,pos_label="algae", multiclass=F){
  print("get F1 score...")
  y_truth<-as.character(y_truth)
  y_pred<-as.character(y_pred)
  if(multiclass==F){
    score<-F1_Score(y_true = y_truth,y_pred = y_pred,positive = pos_label)
    score<-round(score,3)
    return(score)
  }else{
    print("perform multiclass evaluation...")
    y_truth<-as.factor(y_truth)
    y_pred<-as.factor(y_pred)
    
    print(levels(y_truth))
    print(levels(y_pred))
    check_levels<-levels(y_truth) %in% levels(y_pred)
    if(any(!check_levels)){ # check if there is any false (the ! negates check_levels)
      # Align levels
      print("aligning levels....")
      all_levels <- union(levels(y_truth), levels(y_pred))
      y_truth <- factor(y_truth, levels = all_levels)
      y_pred <- factor(y_pred, levels = all_levels)
    }

    # Create a confusion matrix
    conf_matrix <- confusionMatrix(y_pred, y_truth)
    
    # Extract precision, recall, and calculate F1 scores for each class
    precision <- conf_matrix$byClass[, "Precision"]
    recall <- conf_matrix$byClass[, "Recall"]
    f1_scores <- 2 * (precision * recall) / (precision + recall)
    macro_f1 <- mean(f1_scores, na.rm = TRUE)
    return(list(macro_score=macro_f1,all_scores=f1_scores,conf_matrix=conf_matrix))
  }
}


# function to get ARI score
get_ari_score<-function(y_truth,y_pred){
  print("get ARI score...")

  score<-adjustedRandIndex(x = y_pred,y = y_truth)
  score<-round(score,3)
  return(score)
}

# function to get accuracy, which is the proportion of correct predictions
get_accuracy<-function(y_truth,y_pred,pos_label="algae"){
  y_truth<-as.character(y_truth)
  y_pred<-as.character(y_pred)
  
  inds<-which(y_truth==pos_label)
  if(length(inds)!=0){
    y_truth[inds]<-1
    y_truth[-inds]<-0
  }else{
    stop("reference label not found in y_truth")
  }
  
  inds<-which(y_pred==pos_label)
  if(length(inds)!=0){
    y_pred[inds]<-1
    y_pred[-inds]<-0
  }else{
    stop("reference label not found in y_pred")
  }
  
  score <- sum(y_truth == y_pred) / length(y_truth)
  
  score<-round(score,3)
  return(score)
}
