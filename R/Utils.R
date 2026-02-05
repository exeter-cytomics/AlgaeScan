# function to combine all original csv data in one file
get_all_data_csv<-function(paths_data){
  list_all_dfs<-list()
  for(p in paths_data){
    print("reading...")
    print(p)
    data_x<-fread(file = p)
    splitted_seq<-strsplit(p,"/")[[1]]
    name_data<-tail(splitted_seq,1)
    data_x$measure<-rep(name_data,nrow(data_x))
    if(grepl("na|non_algae", name_data)){
      data_x$type<-rep("non-algae",nrow(data_x))
    }else{
      data_x$type<-rep("algae",nrow(data_x))
    }
    list_all_dfs[[name_data]]<-data_x
  }
  df_final<-do.call(rbind,list_all_dfs)
  return(df_final)
}


# function to add classes info
add_class_info<-function(df,df_info){
  all_species<-unique(df$Species)
  df$Class<-rep("NA",nrow(df))
  for(s in all_species){
    print(s)
    class_current_s<-unique(df_info$Class[df_info$Species==s])
    print(class_current_s)
    if(length(class_current_s)>1){
      stop("more than one class for species, probable error")
    }else if(length(class_current_s)==0){
      class_current_s<-"NA"
    }
    inds<-which(df$Species==s)
    df$Class[inds]<-class_current_s
  }
  return(df)
}

# function to get train data and test for all species

partition_data<-function(df_tot){
  all_species<-unique(df_tot$Species)
  list_train_dfs<-list()
  list_test_dfs<-list()
  for(s in all_species){
    print(s)
    df_tot_species_s<-df_tot[df_tot$Species==s]
    set.seed(123)
    index_train<-createDataPartition(df_tot_species_s$Species, p = 0.8, list = T)
    index_train<-index_train$Resample1
    df_tot_species_s_train<-df_tot_species_s[index_train,]
    df_tot_species_s_test<-df_tot_species_s[-index_train,]
    list_train_dfs[[s]]<-df_tot_species_s_train
    list_test_dfs[[s]]<-df_tot_species_s_test
  }
  print("combining dfs...")
  df_test<-do.call(rbind,list_test_dfs)
  df_train<-do.call(rbind,list_train_dfs)
  
  return(list(df_test=df_test,df_train=df_train))
  
}

# function to perform downsample based on variable
get_downsample<-function(df,var="Species",n_samples=1000){
  all_unique_val<-unique(unlist(df[,..var]))
  list_all_dfs<-list()
  for(val in all_unique_val){
    print(val)
    inds_val<-which(df[,..var]==val)
    df_val<-df[inds_val,]
    if(n_samples>nrow(df_val)){
      n_sample<-nrow(df_val)
    }
    if(nrow(df_val)<n_samples){
      n_samples<-nrow(df_val)
    }
    inds_val_down<-sample.int(n=nrow(df_val),size = n_samples)
    df_val_down<-df_val[inds_val_down,]
    list_all_dfs[[val]]<-df_val_down
  }
  df_final<-do.call(rbind,list_all_dfs)
  return(df_final)
}

# function to get best models across all simulations
get_best_model<-function(list_models){
  # get max acc values across models (including nested models) 
  vec_max_values<-c()
  vec_max_values_name<-c()
  for(i in 1:length(list_models)){
    print(i)
    m<-list_models[[i]]
    if(class(m)=="list"){
      vec_max_i<-sapply(m,function(m_i){
        max_acc_m_i<-max(m_i$results$Accuracy)
        return(max_acc_m_i)
      })
      max_across_ma_i<-max(vec_max_i)
      vec_max_values<-c(vec_max_values,max_across_ma_i)
    }else if(class(m)=="train"){
      max_acc_m<-max(m$results$Accuracy)
      vec_max_values<-c(vec_max_values,max_acc_m)
    }
    vec_max_values_name<-c(vec_max_values_name,sprintf("cross_val_rep_%d",i))
    
  }
  names(vec_max_values)<-vec_max_values_name
  # get best model based on max accuracy values
  best_model<-list_models[[which.max(vec_max_values)]]
  if(class(best_model)=="list"){
    vec_max_i<-sapply(best_model,function(m_i){
      max_acc_m_i<-max(m_i$results$Accuracy)
      return(max_acc_m_i)
    })
    best_model<-best_model[[which.max(vec_max_i)]]
    
  }
  return(list(best_model=best_model,vec_max_values=vec_max_values))
}

# function to get best models across all simulations for gmm or other unsupervised approaches
get_best_model_v2<-function(list_models){
  # get max acc values across models (including nested models) 
  vec_max_values<-c()
  vec_max_values_name<-c()
  for(i in 1:length(list_models)){
    print(i)
    list_i<-list_models[[i]]
    max_score<-max(list_i$vec_scores)
    vec_max_values<-c(vec_max_values,max_score)
    vec_max_values_name<-c(vec_max_values_name,sprintf("cross_val_rep_%d",i))
  }
  ind<-which.max(vec_max_values)
  list_model_best<-list_models[[ind]]
  ind<-which.max(list_model_best$vec_scores)
  best_model<-list_model_best$model_out[[ind]]
  names(vec_max_values)<-vec_max_values_name
  return(list(best_model=best_model,vec_max_values=vec_max_values))
}


# function to get time in secs from Rprof
interpret_rprof <- function(profile_path, sampling_interval = 0.02) {
  # Summarize profile data
  profile_summary <- summaryRprof(profile_path)
  
  # Add calculated time columns
  profile_summary$by.self$time_sec <- profile_summary$by.self$self.time * sampling_interval
  profile_summary$by.self$time_min <- profile_summary$by.self$time_sec / 60
  
  profile_summary$by.total$time_sec <- profile_summary$by.total$total.time * sampling_interval
  profile_summary$by.total$time_min <- profile_summary$by.total$time_sec / 60
  
  return(profile_summary)
}

