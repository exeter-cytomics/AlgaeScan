# function to perform pca analysis
pca_analysis<-function(df,select_measure="all",n_pc=10){
  if(select_measure=="all"){
    # for analysis of all measures at once
    df_final<-copy(df) # to avoid inplace modifications from data.table package
  }else{
    # for analysis of single measure
    df_final<-df[df$measure==select_measure,] # select single measure
  }
  
  # remove characters colums for pca analysis
  inds_cols_to_rm<-grep("Species|measure|type|Class",colnames(df_final))
  colums_names_to_rm<-(colnames(df_final)[inds_cols_to_rm])
  
  df_final_copy<-copy(df_final) # to keep info about characters column
  
  df_final[,(colums_names_to_rm) := NULL] # inplace modifications, removal of char columns
  # perform pca 
  print("pca computation....")
  set.seed(123)
  pc <- prcomp(df_final,center = TRUE,scale. = TRUE)
  
  # extract results
  print("extracting results....")
  
  summary_pc_results<-summary(pc) # PC1 and PC2 most important components
  
  df_main_pc<-pc$rotation[,c(1,2)]
  
  df_main_pc<-data.table(df_main_pc,keep.rownames = TRUE)

  df_main_pc_sorted_pc1<-df_main_pc[order(df_main_pc$PC1,decreasing = T),]
  
  df_2d_expr<-data.table(pc$x)
  df_2d_expr<-cbind(df_2d_expr,df_final_copy$Species)
  colnames(df_2d_expr)[ncol(df_2d_expr)]<-"Species"
  df_2d_expr<-cbind(df_2d_expr,df_final_copy$measure)
  colnames(df_2d_expr)[ncol(df_2d_expr)]<-"measure"
  
  
  if("Class" %in% colnames(df_final_copy)){
    df_2d_expr<-cbind(df_2d_expr,df_final_copy$Class)
    colnames(df_2d_expr)[ncol(df_2d_expr)]<-"class"
  }
  # first 10 most important features for PC1
  first_n_feat_pc1<-df_main_pc_sorted_pc1$rn[1:n_pc]
  df_first_n_feat_pc1<-df[,first_n_feat_pc1,
                           with=F]

  # select first 10 PCA components
  vec_type_col<-sapply(df_2d_expr,function(x){
    return(class(x))
  })
  inds<-which(vec_type_col=="character")
  df_2d_expr_select_feat<-df_2d_expr[,c(1:n_pc,inds),with=F]
  
  return(list(df_expr_pc=df_2d_expr,df_expr_first_n_pc=df_2d_expr_select_feat,summary_pc_results=summary_pc_results,df_first_n_feat_pc1=df_first_n_feat_pc1))
}

# function to perform umap analysis
umap_analysis<-function(df,downs=0.05,n_neigh=5,min_dist=0.05,n_epch=500,spread_val=1,metric_select="euclidean"){
  # downsampling keeping same distribution
  print("downsampling...")
  set.seed(123)
  ind_res<-createDataPartition(y=factor(df$Species),times = 1,p = downs) 
  ind_res<-ind_res$Resample1
  df<-df[ind_res,]
  
  #----- remove characters columns if needed
  
  print("check characters columns...")
  vec_type_col<-sapply(df,function(x){
    return(class(x))
  })
  inds_rm<-grep("character",vec_type_col)
  if(length(inds_rm)!=0){
    df_non_char<-df[,-inds_rm, with=F]
  }else{
    df_non_char<-df # test data already without characters columns
  }
  # ------ remove prob columns if needed
  inds_check_prob<-grep("_prob$",colnames(df_non_char))
  if(length(inds_check_prob)!=0){
    df_non_char<-df_non_char[,-inds_check_prob,with=F]
    
  }
  #-------- umap exection
  print("umap computation..." )
  
  set.seed(123)
  umap_result <- umap(df_non_char, 
                      n_neighbors = n_neigh, 
                      min_dist = min_dist, 
                      metric = metric_select, 
                      init = "spectral", 
                      n_epochs = n_epch, 
                      spread = spread_val,
                      random_state = 123,)
  
  print("extracting results..." )
  
  umap_df <- data.table(umap_result$layout)
  colnames(umap_df)<-c("UMAP1","UMAP2")
  umap_df<-cbind(umap_df,df$Species)
  colnames(umap_df)[ncol(umap_df)]<-"Species"
  umap_df<-cbind(umap_df,df$measure)
  colnames(umap_df)[ncol(umap_df)]<-"measure"
  umap_df<-cbind(umap_df,df$type)
  colnames(umap_df)[ncol(umap_df)]<-"type"
  if("type_predicted" %in% colnames(df)){
    umap_df<-cbind(umap_df,df$type_predicted)
    colnames(umap_df)[ncol(umap_df)]<-"type_predicted"
  }
  if("species_predicted" %in% colnames(df)){
    umap_df<-cbind(umap_df,df$species_predicted)
    colnames(umap_df)[ncol(umap_df)]<-"species_predicted"
  }
  if("Species_with_noise" %in% colnames(df)){
    umap_df<-cbind(umap_df,df$Species_with_noise)
    colnames(umap_df)[ncol(umap_df)]<-"Species_with_noise"
  }
  if("species_predicted_with_noise" %in% colnames(df)){
    umap_df<-cbind(umap_df,df$species_predicted_with_noise)
    colnames(umap_df)[ncol(umap_df)]<-"species_predicted_with_noise"
  }
  ind<-grep("class|Class",colnames(df))
  col_name<-colnames(df)[ind]
  if(length(ind)!=0){
    umap_df<-cbind(umap_df,df[,col_name,with=F])
    colnames(umap_df)[ncol(umap_df)]<-"class"
  }
  inds_check_prob<-grep("_prob$",colnames(df))
  if(length(inds_check_prob)!=0){
    umap_df<-cbind(umap_df,df[,inds_check_prob,with=F])
    
  }
  
  return(umap_df)
}

  
