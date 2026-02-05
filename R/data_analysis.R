
# import libraries
library(ggplot2) # for graphics
library(caret) # for ML data split and other ML algorithms
library(e1071) # for svm
library(parallel) # for parallel processing
library(pryr) # for memory usage tracking
library(data.table) # load large cs file (faster than base R)
library(stringr) # string manipulation
library(umap)  # to perform DR using umap
library(doParallel) # parallelization with caret 
library(randomForest) # to perform random forest with caret
library(nnet) # to perform multinomial logistic regression.
library(MLmetrics) # to get F1 score
library(kernlab) # to perform svm
library(naivebayes) # to perform naive bayes
library(mclust) # to perform GMM and calculate ARI
library(e1071) # to perform one class svm
library(reticulate) # to integrate Python in R
# import functions

source("/home/rstudio/marine_data/Analysis_code/Utils.R")
source("/home/rstudio/marine_data/Analysis_code/preliminary_analysis_code.R")
source("/home/rstudio/marine_data/Analysis_code/visualization.R")
source("/home/rstudio/marine_data/Analysis_code/training_code.R")
source("/home/rstudio/marine_data/Analysis_code/prediction_code.R")
source("/home/rstudio/marine_data/Analysis_code/evaluation.R")
source_python("/home/rstudio/marine_data/Analysis_code/isolationforest_code.py")

# install packages (in case of additional packages not included in the docker image)
install.packages("MLmetrics")
install.packages("kernlab")
install.packages("naivebayes")
install.packages("mclust")
######################################################### execute code ####################################

# check total physical memory system (linux system or container)
mem_info <- readLines("/proc/meminfo")
total_mem_kb <- as.numeric(gsub("\\D", "", grep("MemTotal", mem_info, value = TRUE)))
total_mem_gb <- total_mem_kb / 1024^2
cat("Total Physical Memory:", total_mem_gb, "GB")

mem_used()


################### import data and experiment information first batch ###################

# Get csv data 
all_paths_data_original<-list.files("/home/rstudio/marine_data/Aggregated data for paper/Pipeline/A02_001_dataset_02_2024/",full.names = T)
paths_all_data<-all_paths_data_original[grep("m1\\.csv|m2\\.csv|m3\\.csv|na_m1\\.csv|na_m2\\.csv|na_m3\\.csv",all_paths_data_original)]

df_total_data<-get_all_data_csv(paths_data = paths_all_data)
object_size(df_total_data) # check memory occupied by object

unique(df_total_data$type)


######################################################################################################
#####################################  Quality checking  #######################################
######################################################################################################

df_counts_species<-df_total_data[, .N, by = .(measure, Species)]
df_counts_species$measure<-str_remove(df_counts_species$measure,"A02_001_")
df_counts_species$measure<-str_remove(df_counts_species$measure,".csv")

inds<-grep("na",df_counts_species$measure)
df_counts_species_algae_only<-df_counts_species[-inds,] # algae only
df_counts_species_na_only<-df_counts_species[inds,] # non algae only

# checking number of events across measures

table(df_total_data$Species)


get_line_plot(df = df_counts_species_algae_only, x_var = "measure",y_var = "N",color_var = "Species",
              legend_title = "Species",show_legend = T,x_lab = "Measures",y_lab = "Number of events",size_title_x = 30,size_title_y = 30,size_axis_text = 25)

get_line_plot(df = df_counts_species_algae_only, x_var = "measure",y_var = "N",color_var = "Species",
              legend_title = "Species",show_legend = F,x_lab = "Measures",y_lab = "Number of events",
              size_title_x = 30,size_title_y = 30,size_axis_text = 25,col_single_species = "RCC1502")

# checking number of events across Species based on measure

df_counts_species_algae_only_m1<-df_counts_species_algae_only[df_counts_species_algae_only$measure=="m1",]

df_counts_species_algae_only_m1<-df_counts_species_algae_only_m1[order(df_counts_species_algae_only_m1$N,decreasing = T),]

df_counts_species_algae_only_m1$Species<-factor(df_counts_species_algae_only_m1$Species, level=df_counts_species_algae_only_m1$Species)

get_bar_plot(df = df_counts_species_algae_only_m1,x_var = "Species",y_var = "N",size_axis_text=18,
             size_title_x=20,size_title_y=20,x_lab = "Species",y_lab = "Cells count")

# clean dataset

df_total_data<-df_total_data[df_total_data$Species!="RCC1502",] # remove bad quality samples

df_total_data$measure<- str_remove(df_total_data$measure,"A02_001_")
df_total_data$measure<- str_remove(df_total_data$measure,".csv")

# add class for species

df_species_info<-read.csv(file = "~/marine_data/Aggregated data for paper/Pipeline/species list/species list-Roscoff_first_sheet.csv")

mapping <- setNames(df_species_info$Class,df_species_info$Roscoff.Culture.Collection.Identifier)

df_total_data$Class <- mapping[df_total_data$Species]

unique(df_total_data$Class)

# df_total_data contains algae + non algae (i.e., noise) for all measures (day1, day 2, day 7) considering all events for each sample
fwrite(x = df_total_data,file = "/home/rstudio/marine_data/results/df_total_data.csv")



######################################################################################################
#####################################  dimensional reduction analysis (first batch)  #######################################
######################################################################################################

df_total_data<-fread(file = "/home/rstudio/marine_data/results/df_total_data.csv",check.names = F)

#-------------------------- algae ------------------------------
df_total_data$Time<-NULL # remove time channel

df_total_data_alg<-df_total_data[df_total_data$type=="algae"]
df_total_data_alg_copy<-copy(df_total_data_alg) # to keep info about metadata


#--------- PCA

# single measure
out_pca<-pca_analysis(df = df_total_data_alg,select_measure = "m1") #either m1 or m2 or m3,  by default n_pc=10

df_expr_pc<-out_pca$df_expr_pc

df_expr_first_10_pc<-out_pca$df_expr_first_n_pc # first 10 components for UMAP

out_pca$summary_pc_results

out_pca<-pca_analysis(df = df_total_data_alg,select_measure = "m1",n_pc = 5) #either m1 or m2 or m3

df_expr_first_5_pc<-out_pca$df_expr_first_n_pc # first 5 components for UMAP

out_pca<-pca_analysis(df = df_total_data_alg,select_measure = "m1",n_pc = 20) #either m1 or m2 or m3

df_expr_first_20_pc<-out_pca$df_expr_first_n_pc # first 5 components for UMAP

# all measure 

out_pca<-pca_analysis(df = df_total_data_alg,select_measure = "all") # all measures

df_expr_pc<-out_pca$df_expr_pc

df_expr_first_10_pc<-out_pca$df_expr_first_10_pc # first 10 components for UMAP

out_pca$summary_pc_results

# downsample for easy visualization
set.seed(123)
ind_res<-createDataPartition(y=factor(df_expr_pc$Species),times = 1,p = 0.05) # downsampling by keeping same probability distribution
ind_res<-ind_res$Resample1
df_expr_pc_downsampled<-df_expr_pc[ind_res,]


# species layer
get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (13%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25) # m1


get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (14%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25) # m2


get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (12%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25) # m3

# color single species
get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (13%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25,
                 col_single_species = "RCC6336") 

get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (13%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25,
                 col_single_species = "RCC950") 

get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (13%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25,
                 col_single_species = "RCC1507") 

get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (13%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25,
                 col_single_species = "RCC1511") 

get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (13%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25,
                 col_single_species = "RCC539") 


# measure layers
get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "measure",shape_var = "measure",size_p = 1.5,
                 legend_title = "measure",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (13%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25)

#----------- UMAP

# testing several parameters and number of components
# umap_df<-umap_analysis(df_expr_first_10_pc,downs = 0.03,n_neigh = 10,min_dist = 0.05,n_epch = 500,metric_select = "euclidean") # 3% downsampling for single measures, 1% for all measures
# 
# umap_df<-umap_analysis(df_expr_first_10_pc,downs = 0.01,n_neigh = 10,min_dist = 0.05,n_epch = 500,metric_select = "cosine") # 3% downsampling for single measures, 1% for all measures
# 
# umap_df<-umap_analysis(df_expr_first_5_pc,downs = 0.01,n_neigh = 10,min_dist = 0.05,n_epch = 500,metric_select = "cosine") # 3% downsampling for single measures, 1% for all measures
# 
# umap_df<-umap_analysis(df_expr_first_20_pc,downs = 0.01,n_neigh = 10,min_dist = 0.05,n_epch = 500,metric_select = "cosine") # 3% downsampling for single measures, 1% for all measures

df_total_data_alg_selected_m<-df_total_data_alg[df_total_data_alg$measure=="m3",]

umap_df<-umap_analysis(df_total_data_alg,downs = 0.01,n_neigh = 15,min_dist = 0.3,n_epch = 200,metric_select = "cosine") # 3% downsampling for single measures, 1% for all measures

umap_df<-umap_analysis(df_total_data_alg_selected_m,downs = 0.03,n_neigh = 15,min_dist = 0.3,n_epch = 200,metric_select = "cosine") # 3% downsampling for single measures, 1% for all measures


#best one m1 : downs=0.03, n_neight=15,dist=0.3,epocs=200, metric=cosine
#best one all m: downs=0.01, n_neight=15,dist=0.3,epocs=200, metric=cosine

# species layer
get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = T,x_lab = "UMAP1",y_lab = "UMAP2",size_title_x = 30,size_title_y = 30,size_axis_text = 25)

# class layer
get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "class",size_p = 0.05,
                 legend_title = "Class",show_legend = T,x_lab = "UMAP1",y_lab = "UMAP2",size_title_x = 30,size_title_y = 30,size_axis_text = 25)

# color single species
get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",size_title_x = 30,size_title_y = 30,size_axis_text = 25,
                 col_single_species = "RCC1435")

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",size_title_x = 30,size_title_y = 30,size_axis_text = 25,
                 col_single_species = "RCC950")

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",size_title_x = 30,size_title_y = 30,size_axis_text = 25,
                 col_single_species = "RCC69")

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",size_title_x = 30,size_title_y = 30,size_axis_text = 25,
                 col_single_species = "RCC539")

# measure layer
get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "measure",shape_var = "measure",size_p = 1.5,
                 legend_title = "Measure",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",size_title_x = 30,size_title_y = 30,size_axis_text = 25)

# Story: We started with PCA to analyze general patterns. PCA did not identify well separated structure when considering the species layer, suggesting that data may show more complex relationships. We retained the 10 first PCA components and we applied a UMAP algorithm

# The umap algorithm showed better well separated structures  with the species layer.

# both PCA and Umap show overlapping structures with the measurement layer, suggesting that the measurement of the data does not represent important information in the data.


#--------- non algae ---------
df_total_data_na<-df_total_data[df_total_data$type=="non-algae"]
df_total_data_na_copy<-copy(df_total_data_na)

#--------- PCA

# single measure
out_pca<-pca_analysis(df = df_total_data_na,select_measure = "m1") #either m1 or m2 or m3

df_expr_pc<-out_pca$df_expr_pc

df_expr_first_10_pc<-out_pca$df_expr_first_10_pc # first 10 components for UMAP

out_pca$summary_pc_results


# all measure 

out_pca<-pca_analysis(df = df_total_data_na,select_measure = "all") # all measures

df_expr_pc<-out_pca$df_expr_pc

df_expr_first_10_pc<-out_pca$df_expr_first_10_pc # first 10 components for UMAP

out_pca$summary_pc_results


# downsample for easy visualization
ind_res<-createDataPartition(y=factor(df_expr_pc$Species),times = 1,p = 0.01) # downsampling by keeping same probability distribution
ind_res<-ind_res$Resample1
df_expr_pc_downsampled<-df_expr_pc[ind_res,]


# species layer
get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (13%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25) # m1


get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (14%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25) # m2


get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (12%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25) # m3

# measure layers
get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "measure",shape_var = "measure",size_p = 1.5,
                 legend_title = "measure",show_legend = F,x_lab = "PC1 (71%)",y_lab = "PC2 (11%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25)

#----------- UMAP

umap_df<-umap_analysis(df_total_data_na,downs = 0.002,n_neigh = 15,min_dist = 0.3,n_epch = 200,
                       metric_select = "cosine",spread_val = 1) # 0.2% for all measures of non-algae


# species layer
get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",size_title_x = 30,size_title_y = 30,size_axis_text = 25)

# measure layer
get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "measure",shape_var = "measure",size_p = 1.5,
                 legend_title = "Measure",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",size_title_x = 30,size_title_y = 30,size_axis_text = 25)


######################################################################################################
#####################################  Preprocessing for training and testing  #######################################
######################################################################################################

# separating train vs test (80% train vs 20% test) for each species and algae vs non-algae (i.e., 80% of species 1,2,3 etc....)
df_total_data<-fread(file = "/home/rstudio/marine_data/results/df_total_data.csv",check.names = F)
df_total_data$Time<-NULL # remove time channel

table(df_total_data$Species)

# some models in in caret do not like symbol "-" in class name, it converts "-" to "." 
# thus changing the class name,which becomes different than test class name. Caret shows an error, so we temporary convert to "_"

df_total_data$type<-str_replace(df_total_data$type,"non-algae","non_algae") 

unique(df_total_data$type)

# considering all measurements together (algae + non algae)
list_dfs<-partition_data(df_tot = df_total_data)

test_data<-list_dfs$df_test # 20% for each species

train_data<-list_dfs$df_train # 80% for each species

fwrite(train_data,file = "/home/rstudio/marine_data/results/train_data_all_m.csv",row.names = F)
fwrite(test_data,file = "/home/rstudio/marine_data/results/test_data_all_m.csv",row.names = F)

# separating measurements (algae + non algae)
# Note: we do not need to separate measurements, because during DR we saw that measurements data overlap.
# df_total_data_m1<-df_total_data[grep("m1|na_m2|na_m3",df_total_data$measure),] # m1 algae + non algae (m1,m2,m3)
# 
# list_dfs<-partition_data(df_tot = df_total_data_m1)
# 
# test_data<-list_dfs$df_test # 20% for each species
# 
# train_data<-list_dfs$df_train # 80% for each species
# 
# train_data[,"Species"]=="RCC10364"


print("Number of events per Species:")
print(table(test_data$Species))
print("Number of events per gate:")
print(table(test_data$type))


# The analysis proceeds as follows:
# 1) We need to generate first the model to separate algae from non-algae (two classes model, the ML class variable is the Class column)
# 2) We generate the second model that can identify the algae species (the ML class variable is Species column)

# First however, we need to test different models on subset  of the training to understand which model works best.
# So in next sections (algae vs noise analysis), we aim to generate the first model. But first we perform simulations to understand which models to use.
######################################################################################################
#####################################  Algae vs noise analysis (first model)  #######################################
######################################################################################################

###############  training and testing steps #######################################
train_data<-fread(file = "/home/rstudio/marine_data/results/train_data_all_m.csv")
test_data<-fread(file = "/home/rstudio/marine_data/results/test_data_all_m.csv")
########### supervised algorithms #################
# note: sample down for each species for type of category (algae or non-algae), example: 2000 events per type for x specie = 4000 events in total per x species

# simulations (multiple repetitions) using a small amounts of events to understand which model works best.
# A model that is stable across multiple repetitions (similar results across multiple repetitions) it is a good thing because it means using different training datasets 
# (i.e., a condition simulated using different parts of the same training set through different seed for downsampling)
# is not dramatically altering the model accuracy. So it means that even if we do not train on everything, which would consume time and memory, 
# we would still get similar results.
#---- random forest-----

# training 

out_train<-data_train(data = train_data,method = "rf",n_cores = 16,n_samples_down = 300, method_ctrl = "cv",n_tree = c(5,10,15),rep_seed = 123,k_cv = 5, n_rep = 10)

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model

class(best_model)
# prediction

y_pred<-data_prediction(input_model = best_model,test_df = test_data)

test_score<-get_F1_score(y_truth = test_data$type,y_pred = y_pred,pos_label = "algae")


# get confusion matrix
y_truth<-ifelse(test_data$type=="non_algae","unknown","algae") # Better to use the term unknown instead of non-algae (since they could be algae, just we do not know)
y_pred<-ifelse(y_pred=="non_algae","unknown","algae")
unique(y_truth)
unique(y_pred)
conf_plot<-get_confusion_matrix(y_pred = y_pred,y_true = y_truth, show_legend = T)

vec_results_rf<-c("RF",vec_max_values,max_value,test_score)


#---- multinomial logistic regression -----

# training 

out_train<-data_train(data = train_data,method = "log_reg",n_cores = 16,n_samples_down = 300, method_ctrl = "cv",rep_seed = 123,k_cv = 5, n_rep = 10)

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model

# prediction

y_pred<-data_prediction(input_model = best_model,test_df = test_data)

y_truth<-test_data$type

test_score<-get_F1_score(y_truth = y_truth,y_pred = y_pred,pos_label = "algae")

# get confusion matrix
conf_plot<-get_confusion_matrix(y_pred = y_pred,y_true = y_truth, show_legend = T)

vec_results_lr<-c("MLR",vec_max_values,max_value,test_score)

#----- KNN -------

# training 

out_train<-data_train(data = train_data,method = "knn",n_cores = 16,n_samples_down = 300, method_ctrl = "cv",rep_seed = 123,k_cv = 5, n_rep = 10)

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model


# prediction

y_pred<-data_prediction(input_model = best_model,test_df = test_data)


test_score<-get_F1_score(y_truth = test_data$type,y_pred = y_pred,pos_label = "algae")


# get confusion matrix
conf_plot<-get_confusion_matrix(y_pred = y_pred,y_true = test_data$type, show_legend = T)

vec_results_knn<-c("KNN",vec_max_values,max_value,test_score)


#---- neural network -------

# training 

out_train<-data_train(data = train_data,method = "fnn",n_cores = 16,n_samples_down = 300, method_ctrl = "cv",rep_seed = 123,k_cv = 5, n_rep = 10)

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model

# prediction

y_pred<-data_prediction(input_model = best_model,test_df = test_data)


test_score<-get_F1_score(y_truth = test_data$type,y_pred = y_pred,pos_label = "algae")

# get confusion matrix
conf_plot<-get_confusion_matrix(y_pred = y_pred,y_true = test_data$type, show_legend = T)

vec_results_fnn<-c("FNN",vec_max_values,max_value,test_score)

#---- support vector machine -------

# training 

out_train<-data_train(data = train_data,method = "svm",n_cores = 16,n_samples_down = 300, method_ctrl = "cv",rep_seed = 123,k_cv = 5, n_rep = 10)

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model

# prediction

y_pred<-data_prediction(input_model = best_model,test_df = test_data)


test_score<-get_F1_score(y_truth = test_data$type,y_pred = y_pred,pos_label = "algae")

# get confusion matrix
conf_plot<-get_confusion_matrix(y_pred = y_pred,y_true = test_data$type, show_legend = T)

vec_results_svm<-c("SVM",vec_max_values,max_value,test_score)

#---- naive bayes -------

# training 

out_train<-data_train(data = train_data,method = "nb",n_cores = 16,n_samples_down = 300, method_ctrl = "cv",rep_seed = 123,k_cv = 5, n_rep = 10)

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model

# prediction

y_pred<-data_prediction(input_model = best_model,test_df = test_data)

test_score<-get_F1_score(y_truth = test_data$type,y_pred = y_pred,pos_label = "algae")

# get confusion matrix
conf_plot<-get_confusion_matrix(y_pred = y_pred,y_true = test_data$type, show_legend = T)

vec_results_nb<-c("NB",vec_max_values,max_value,test_score)

# ------- combine all algorithms results in one table
df_all_results<-as.data.frame(rbind(vec_results_rf,vec_results_lr,vec_results_knn,vec_results_fnn,vec_results_svm,vec_results_nb))
colnames(df_all_results)<-c("Name_model",names(vec_max_values),"Best_rep","Pred_F1score")
row.names(df_all_results)<-NULL

#----- export results
write.csv(df_all_results,file = "/home/rstudio/marine_data/results/df_results_supervised_models_algae_vs_non_algae.csv",row.names = F)


########### unsupervised algorithms  #################



#---- GMM -----

# training 

out_train<-data_train(data = train_data,method = "gmm",n_cores = 16,n_samples_down = 300, method_ctrl = "cv",rep_seed = 123,k_cv = 5, n_rep = 10)

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model

# prediction

y_pred<-data_prediction(input_model = best_model,test_df = test_data)

test_score<-get_ari_score(y_truth = test_data$type,y_pred = y_pred)

# get confusion matrix
#conf_plot<-get_confusion_matrix(y_pred = out_pred$y_pred,y_true = out_pred$y_truth, show_legend = T)

vec_results_gmm<-c("GMM",vec_max_values,max_value,test_score)

#---- one class svm -----

# training 

out_train<-data_train(data = train_data,method = "one_class_svm",n_cores = 1,n_samples_down = 300, method_ctrl = "cv",rep_seed = 123,k_cv = 5, n_rep = 10,
                      ref_label = "algae",ref_var = "type")

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model

best_model$finalModel
class(best_model)
# prediction

y_pred<-data_prediction(input_model = best_model,test_df = test_data,normalize = T) # svm was trained on normalized data, so also test data needs to be normalized

y_truth<-test_data$type

y_truth<-ifelse(test_data$type=="algae","1","0")

test_score<-get_F1_score(y_truth = y_truth,y_pred = y_pred,pos_label = "1")

# get confusion matrix
conf_plot<-get_confusion_matrix(y_pred = y_pred,y_true = y_truth, show_legend = T)

vec_results_one_class_svm<-c("One_class_svm",vec_max_values,max_value,test_score)


#---- isolation forest -----

# training 

out_train<-data_train(data = train_data,method = "isoforest",n_cores = 16,n_samples_down = 300, method_ctrl = "cv",rep_seed = 123,k_cv = 5, n_rep = 10)

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model

class(best_model)
# prediction

y_pred<-data_prediction(input_model = best_model,test_df = test_data,normalize = T) # isolation forest was trained on normalized data, so also test data needs to be normalized

y_truth<-test_data$type

y_truth<-ifelse(test_data$type=="algae","1","0")

test_score<-get_F1_score(y_truth = y_truth,y_pred = y_pred,pos_label = "1")

# get confusion matrix
conf_plot<-get_confusion_matrix(y_pred = y_pred,y_true = y_truth, show_legend = T)

vec_results_isoforest<-c("isoforest",vec_max_values,max_value,test_score)

# ------- combine all algorithms results in one table
df_all_results<-as.data.frame(rbind(vec_results_gmm,vec_results_one_class_svm,vec_results_isoforest))
colnames(df_all_results)<-c("Name_model",names(vec_max_values),"Best_rep","Pred_F1score")
row.names(df_all_results)<-NULL

#----- export results
write.csv(df_all_results,file = "/home/rstudio/marine_data/results/df_results_unsupervised_models_algae_vs_non_algae.csv",row.names = F)

#----------------------- memory and speed profiling -------------------

# change the function based on the algorithm you want to test the speed and memory (1 rep)
Rprof("profile.out",memory.profiling = TRUE)
out_train<-data_train(data = train_data,method = "gmm",n_cores = 1,n_samples_down = 300, method_ctrl = "cv",rep_seed = 123,k_cv = 5, n_rep = 1)
Rprof(NULL)

df_mem<-summaryRprof("profile.out",memory = "stats")
df_mem<-as.data.frame(t(as.data.frame(df_mem$`"data_train":"data_training"`)))
tot_mem_bytes<-df_mem$max.vsize.small+df_mem$max.vsize.large
tot_mem_mb<- tot_mem_bytes/(1024^2)

Rprof("profile.out",memory.profiling = TRUE)
best_model<-out_train$best_model
y_pred<-data_prediction(input_model = best_model,test_df = test_data,normalize = F)
Rprof(NULL)

df_mem<-summaryRprof("profile.out",memory = "stats")
df_mem<-as.data.frame(t(as.data.frame(df_mem$`"data_prediction":"["`)))
tot_mem_bytes<-df_mem$max.vsize.small+df_mem$max.vsize.large
tot_mem_mb<- tot_mem_bytes/(1024^2)


interpret_rprof("profile.out")

# ---------- plots results cross validation -------------

df_all_results_supervised<-read.csv("/home/rstudio/marine_data/results/df_results_supervised_models_algae_vs_non_algae.csv")
df_all_results_unsupervised<-read.csv("/home/rstudio/marine_data/results/df_results_unsupervised_models_algae_vs_non_algae.csv")

# supervised
cross_val_cols<-colnames(df_all_results_supervised)[grep("cross",colnames(df_all_results_supervised))]
df_all_results_reshaped<-reshape2::melt(data = df_all_results_supervised,
                                        id.vars = c("Name_model","Pred_F1score","Best_rep"),
                                        variable.name = "Cross_validation_rep",value.name = "Value",measure.vars=cross_val_cols)


df_all_results_reshaped$Value<-as.numeric(df_all_results_reshaped$Value)

df_all_results_reshaped$Cross_validation_rep<-str_remove(df_all_results_reshaped$Cross_validation_rep,"cross_val_rep_")

df_all_results_reshaped$Cross_validation_rep<-as.integer(df_all_results_reshaped$Cross_validation_rep)

sorted_levels<-unique(sort(as.integer(df_all_results_reshaped$Cross_validation_rep)))

df_all_results_reshaped$Cross_validation_rep<-factor(df_all_results_reshaped$Cross_validation_rep,levels=as.character(sorted_levels))

get_line_plot_v2(df = df_all_results_reshaped,x_var = "Cross_validation_rep",y_var = "Value",color_var = "Name_model",size_axis_text = 18,
                 size_title_x = 20,size_title_y = 20,legend_title = "Model",show_legend = F,x_lab = "Repetition",y_lab = "Cross-validation accuracy")

# unsupervised
cross_val_cols<-colnames(df_all_results_unsupervised)[grep("cross",colnames(df_all_results_unsupervised))]
df_all_results_reshaped<-reshape2::melt(data = df_all_results_unsupervised,
                                        id.vars = c("Name_model","Pred_F1score","Best_rep"),
                                        variable.name = "Cross_validation_rep",value.name = "Value",measure.vars=cross_val_cols)


df_all_results_reshaped$Value<-as.numeric(df_all_results_reshaped$Value)

df_all_results_reshaped$Cross_validation_rep<-str_remove(df_all_results_reshaped$Cross_validation_rep,"cross_val_rep_")

df_all_results_reshaped$Cross_validation_rep<-as.integer(df_all_results_reshaped$Cross_validation_rep)

sorted_levels<-unique(sort(as.integer(df_all_results_reshaped$Cross_validation_rep)))

df_all_results_reshaped$Cross_validation_rep<-factor(df_all_results_reshaped$Cross_validation_rep,levels=as.character(sorted_levels))

get_line_plot_v2(df = df_all_results_reshaped,x_var = "Cross_validation_rep",y_var = "Value",color_var = "Name_model",size_axis_text = 18,
                 size_title_x = 20,size_title_y = 20,legend_title = "Model",show_legend = F,x_lab = "Repetition",y_lab = "Cross-validation accuracy")


############################## training final best model on large events first batch ###############################
# Once we understood which is model is best, we train on a very large section of the training set (or all the training set)
train_data<-fread(file = "/home/rstudio/marine_data/results/train_data_all_m.csv")
test_data<-fread(file = "/home/rstudio/marine_data/results/test_data_all_m.csv")
#---- random forest-----

# training 

# no repetitions, we are not simulating now. We still test different number of trees.  n_samples_down = 5000 or "none"
out_train<-data_train(data = train_data,method = "rf",n_cores = 16,n_samples_down = "none", method_ctrl = "cv",n_tree = c(5,10,15),rep_seed = 123,k_cv = 5, n_rep = 1)

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model

class(best_model)

saveRDS(object = best_model,file =  "/home/rstudio/marine_data/results/final_rf_model_alg_vs_noise_all_train_data.rds")


# prediction

y_pred<-data_prediction(input_model = best_model,test_df = test_data)

test_score<-get_F1_score(y_truth = test_data$type,y_pred = y_pred,pos_label = "algae")

# get confusion matrix

y_truth<-ifelse(test_data$type=="non_algae","unknown","algae") # Better to use the term unknown instead of non-algae (since they could be algae, just we do not know)
y_pred<-ifelse(y_pred=="non_algae","unknown","algae")
unique(y_truth)
unique(y_pred)
conf_plot<-get_confusion_matrix(y_pred = y_pred,y_true = y_truth, show_legend = T)

# add prediction to test data

test_data$type_predicted<-y_pred

# export results
fwrite(test_data, file = "/home/rstudio/marine_data/results/pred_algae_rf.csv",row.names = F)

#---- isolation forest -----

# training 

out_train<-data_train(data = train_data,method = "isoforest",n_cores = 16,n_samples_down = "none", method_ctrl = "cv",rep_seed = 123,k_cv = 5, n_rep = 1)

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model

class(best_model)
# make sure reticulate uses the Python environment because best_model is a scikit learn object
use_python("/usr/bin/python3", required = TRUE)  # adjust path if needed

# import joblib
joblib <- import("joblib")

# dump your Python model object to a file
joblib$dump(best_model, "/home/rstudio/marine_data/results/final_if_model_alg_vs_noise_all_train_data.pkl")


# prediction

y_pred<-data_prediction(input_model = best_model,test_df = test_data,normalize = T) # isolation forest was trained on normalized data, so also test data needs to be normalized

y_truth<-test_data$type

y_truth<-ifelse(y_truth=="algae","1","0")

test_score<-get_F1_score(y_truth = y_truth,y_pred = y_pred,pos_label = 1)

# get confusion matrix
y_pred<-ifelse(y_pred=="1","algae","unknown")
y_truth<-ifelse(y_truth=="1","algae","unknown")
unique(y_truth)
unique(y_pred)
conf_plot<-get_confusion_matrix(y_pred = y_pred,y_true = y_truth, show_legend = T)

# add prediction to test data

test_data$type_predicted<-y_pred

# export results
fwrite(test_data, file = "/home/rstudio/marine_data/results/pred_algae_if.csv",row.names = F)

#---- one class svm -----

# training 

out_train<-data_train(data = train_data,method = "one_class_svm",n_cores = 16,n_samples_down = 5000, method_ctrl = "cv",rep_seed = 123,k_cv = 5, n_rep = 1)

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model

class(best_model)
# prediction

y_pred<-data_prediction(input_model = best_model,test_df = test_data,normalize = T) # svm was trained on normalized data, so also test data needs to be normalized

y_truth<-test_data$type

y_truth<-ifelse(y_truth=="algae",1,0)

test_score<-get_F1_score(y_truth = y_truth,y_pred = y_pred,pos_label = 1)
# get confusion matrix

y_pred_conf<-ifelse(y_pred==1,"algae","noise")
y_truth<-ifelse(test_data$type=="non-algae","noise","algae")
unique(y_truth)
unique(y_pred)
conf_plot<-get_confusion_matrix(y_pred = y_pred_conf,y_true = y_truth, show_legend = T)

# add prediction to test data

test_data$type_predicted<-y_pred_conf

# export results
write.csv(test_data, file = "/home/rstudio/marine_data/results/pred_algae_ocsvm.csv",row.names = F)

####################  visualize results  #######################################

# ------- show algae vs non algae/unknown


# UMAP for visualize algae using RF

df_pred_algae_rf<-fread(file = "/home/rstudio/marine_data/results/pred_algae_rf.csv",check.names = F)

#umap_df<-umap_analysis(df_pred_algae_rf,downs = 0.01,n_neigh = 10,min_dist = 0.05,n_epch = 500) # old version
#umap_df<-umap_analysis(df_pred_algae_rf,downs = 0.01,n_neigh = 15,min_dist = 0.3,n_epch = 500,metric_select = "cosine") 
#umap_df<-umap_analysis(df_pred_algae_rf,downs = 0.005,n_neigh = 15,min_dist = 0.3,n_epch = 200,spread_val = 5,metric_select = "cosine") 

umap_df<-umap_analysis(df_pred_algae_rf,downs = 0.005,n_neigh = 15,min_dist = 0.1,n_epch = 200,spread_val = 5,metric_select = "cosine") 

umap_df$type<-ifelse(umap_df$type=="non_algae","unknown","algae")
umap_df$type_predicted<-ifelse(umap_df$type_predicted=="non_algae","unknown","algae")

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "type",size_p = 0.5,
                 legend_title = "Gate",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25) # reference plot

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "type_predicted",size_p = 0.5,
                 legend_title = "Gate",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25)

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "Species",size_p = 0.5,
                 legend_title = "Gate",show_legend = T,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25)

# UMAP for visualize algae using IF
df_pred_algae_if<-fread(file = "/home/rstudio/marine_data/results/pred_algae_if.csv",check.names = F)

umap_df<-umap_analysis(df_pred_algae_if,downs = 0.005,n_neigh = 15,min_dist = 0.1,n_epch = 200,spread_val = 5,metric_select = "cosine") 

umap_df$type<-ifelse(umap_df$type=="non_algae","unknown","algae")

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "type_predicted",size_p = 0.5,
                 legend_title = "Gate",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25)


######################################################################################################
#####################################  Algae species analysis (second model)  #######################################
######################################################################################################

# train_data_all_m_algae.csv and test_data_all_m_algae.csv generated by subsetting train_data and test data to select only algae

# in this section we need to do the same thing we did in the previous section, but in this case we are developing the second model to identify algae species (the class variable is the Species column).

# We first need to perform simulations to understand which model works best. Then we develop the final model using a large chunk of training set using the model selected.
###############  training and testing steps #######################################
train_data<-fread(file = "/home/rstudio/marine_data/results/train_data_all_m_algae.csv")
test_data<-fread(file = "/home/rstudio/marine_data/results/test_data_all_m_algae.csv")

########### supervised algorithms #################
# note: sample down for each species for type of category (algae or non-algae), example: 2000 events per type for x specie = 4000 events in total per x species

#---- random forest-----

# training 

out_train<-data_train(data = train_data,method = "rf",n_cores = 16,n_samples_down = 300, method_ctrl = "cv",n_tree = c(5,10,15),rep_seed = 123,k_cv = 5, n_rep = 10,ref_var = "Species")

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model

class(best_model)

# prediction

y_pred<-data_prediction(input_model = best_model,test_df = test_data)

test_score<-get_F1_score(y_truth = test_data$Species,y_pred = y_pred,multiclass = T)

test_score<-test_score$macro_score

vec_results_rf<-c("RF",vec_max_values,max_value,test_score)

#---- multinomial logistic regression -----

# training 

out_train<-data_train(data = train_data,method = "log_reg",n_cores = 1,n_samples_down = 300, method_ctrl = "cv",rep_seed = 123,k_cv = 5, n_rep = 10,ref_var = "Species")

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model

# prediction

y_pred<-data_prediction(input_model = best_model,test_df = test_data)

test_score<-get_F1_score(y_truth = test_data$Species,y_pred = y_pred,multiclass = T)

test_score<-test_score$macro_score


vec_results_lr<-c("MLR",vec_max_values,max_value,test_score)

#----- KNN -------

# training 

out_train<-data_train(data = train_data,method = "knn",n_cores = 16,n_samples_down = 300, method_ctrl = "cv",rep_seed = 123,k_cv = 5, n_rep = 10,ref_var = "Species")

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model


# prediction

y_pred<-data_prediction(input_model = best_model,test_df = test_data)


test_score<-get_F1_score(y_truth = test_data$Species,y_pred = y_pred,multiclass = T)

test_score<-test_score$macro_score


vec_results_knn<-c("KNN",vec_max_values,max_value,test_score)


#---- neural network -------

# training 

out_train<-data_train(data = train_data,method = "fnn",n_cores = 16,n_samples_down = 300, method_ctrl = "cv",rep_seed = 123,k_cv = 5, n_rep = 10,ref_var = "Species")

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model

# prediction

y_pred<-data_prediction(input_model = best_model,test_df = test_data)

test_score<-get_F1_score(y_truth = test_data$Species,y_pred = y_pred,multiclass = T)

test_score<-test_score$macro_score



vec_results_fnn<-c("FNN",vec_max_values,max_value,test_score)

#---- support vector machine -------

# training 

out_train<-data_train(data = train_data,method = "svm",n_cores = 16,n_samples_down = 300, method_ctrl = "cv",rep_seed = 123,k_cv = 5, n_rep = 10,ref_var = "Species")

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model

# prediction

y_pred<-data_prediction(input_model = best_model,test_df = test_data)

test_score<-get_F1_score(y_truth = test_data$Species,y_pred = y_pred,multiclass = T)

test_score<-test_score$macro_score

vec_results_svm<-c("SVM",vec_max_values,max_value,test_score)

#---- naive bayes -------

# training 

out_train<-data_train(data = train_data,method = "nb",n_cores = 16,n_samples_down = 300, method_ctrl = "cv",rep_seed = 123,k_cv = 5, n_rep = 10,ref_var = "Species")

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model

# prediction

y_pred<-data_prediction(input_model = best_model,test_df = test_data)

test_score<-get_F1_score(y_truth = test_data$Species,y_pred = y_pred,multiclass = T)

test_score<-test_score$macro_score

vec_results_nb<-c("NB",vec_max_values,max_value,test_score)

# ------- combine all algorithms results in one table -------
df_all_results<-as.data.frame(rbind(vec_results_rf,vec_results_lr,vec_results_knn,vec_results_fnn,vec_results_svm,vec_results_nb))
colnames(df_all_results)<-c("Name_model",names(vec_max_values),"Best_rep","Pred_F1score")
row.names(df_all_results)<-NULL

#----- export results --------
write.csv(df_all_results,file = "/home/rstudio/marine_data/results/df_results_supervised_models_algae_species.csv",row.names = F)


# ---------- plots results cross validation -------------

df_all_results_supervised<-read.csv("/home/rstudio/marine_data/results/df_results_supervised_models_algae_species.csv")

# supervised
cross_val_cols<-colnames(df_all_results_supervised)[grep("cross",colnames(df_all_results_supervised))]
df_all_results_reshaped<-reshape2::melt(data = df_all_results_supervised,
                                        id.vars = c("Name_model","Pred_F1score","Best_rep"),
                                        variable.name = "Cross_validation_rep",value.name = "Value",measure.vars=cross_val_cols)


df_all_results_reshaped$Value<-as.numeric(df_all_results_reshaped$Value)

df_all_results_reshaped$Cross_validation_rep<-str_remove(df_all_results_reshaped$Cross_validation_rep,"cross_val_rep_")

df_all_results_reshaped$Cross_validation_rep<-as.integer(df_all_results_reshaped$Cross_validation_rep)

sorted_levels<-unique(sort(as.integer(df_all_results_reshaped$Cross_validation_rep)))

df_all_results_reshaped$Cross_validation_rep<-factor(df_all_results_reshaped$Cross_validation_rep,levels=as.character(sorted_levels))

get_line_plot_v2(df = df_all_results_reshaped,x_var = "Cross_validation_rep",y_var = "Value",color_var = "Name_model",size_axis_text = 18,
                 size_title_x = 20,size_title_y = 20,legend_title = "Model",show_legend = T,x_lab = "Repetition",y_lab = "Cross-validation accuracy")



############################### training final best model on large events first batch ################################
train_data<-fread(file = "/home/rstudio/marine_data/results/train_data_all_m_algae.csv")
test_data<-fread(file = "/home/rstudio/marine_data/results/test_data_all_m_algae.csv")

# ------ training using MLR -----

out_train<-data_train(data = train_data,method = "log_reg",n_cores = 1,n_samples_down = 5000, method_ctrl = "cv",
                      rep_seed = 123,k_cv = 5, n_rep = 1,ref_var = "Species")

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model

saveRDS(object = best_model,file =  "/home/rstudio/marine_data/results/final_mlr_model_alg_species.rds")

# ------ training using RF -----

out_train<-data_train(data = train_data,method = "rf",n_cores = 8,n_samples_down = 5000, method_ctrl = "cv",n_tree = c(5,10,15),
                      rep_seed = 123,k_cv = 5, n_rep = 1,ref_var = "Species")

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model

best_model$finalModel
saveRDS(object = best_model,file =  "/home/rstudio/marine_data/results/final_rf_model_alg_species.rds")

#-------- prediction -------------

model_b_rf<-readRDS("/home/rstudio/marine_data/results/final_rf_model_alg_species.rds")

model_b_mlr<-readRDS("/home/rstudio/marine_data/results/final_mlr_model_alg_species.rds")


df_pred_algae_rf<-fread(file = "/home/rstudio/marine_data/results/pred_algae_rf.csv",check.names = F)

df_only_algae_rf_pred<-df_pred_algae_rf[df_pred_algae_rf$type_predicted=="algae",] # predicting on predicted algae events

y_pred_from_pred<-data_prediction(input_model = model_b_mlr,test_df = df_only_algae_rf_pred) # change variable based on chosen model

test_score_out<-get_F1_score(y_truth = df_only_algae_rf_pred$Species,y_pred = y_pred_from_pred,multiclass = T)

test_score<-test_score_out$macro_score

all_scores<-test_score_out$all_scores


# df_only_algae_ref<-df_pred_algae_rf[df_pred_algae_rf$type=="algae",] # predicting on reference algae events
# 
# y_pred_from_ref<-data_prediction(input_model = best_model,test_df = df_only_algae_ref)
# 
# test_score_out<-get_F1_score(y_truth = df_only_algae_ref$Species,y_pred = y_pred_from_ref,multiclass = T)
# 
# test_score<-test_score_out$macro_score
# 
# 
# y_pred_all<-data_prediction(input_model = best_model,test_df = df_pred_algae_rf) # predicting on all reference events (algae + noise)
# 
# test_score_out<-get_F1_score(y_truth = df_pred_algae_rf$Species,y_pred = y_pred_all,multiclass = T)
# 
# test_score<-test_score_out$macro_score

# add prediction to test data

df_only_algae_rf_pred$species_predicted<-y_pred_from_pred


# export results
fwrite(df_only_algae_ref, file = "/home/rstudio/marine_data/results/pred_algae_species_mlr_from_ref_algae.csv",row.names = F) # species predicted (using mlr model) with mlr on reference algae events

fwrite(df_only_algae_rf_pred, file = "/home/rstudio/marine_data/results/pred_algae_species_mlr_from_pred_rf_algae.csv",row.names = F) # species predicted (using mlr model) on algae events predicted with rf

fwrite(df_pred_algae_rf, file = "/home/rstudio/marine_data/results/pred_algae_species_mlr_all_events.csv",row.names = F) # species predicted (using mlr model) on all events


fwrite(df_only_algae_ref, file = "/home/rstudio/marine_data/results/pred_algae_species_rf_from_ref_algae.csv",row.names = F) # species predicted (using rf model) with rf on reference algae events

fwrite(df_only_algae_rf_pred, file = "/home/rstudio/marine_data/results/pred_algae_species_rf_from_pred_rf_algae.csv",row.names = F) # species predicted (using rf model) on algae events predicted with rf

fwrite(df_pred_algae_rf, file = "/home/rstudio/marine_data/results/pred_algae_species_rf_all_events.csv",row.names = F) # species predicted (using rf model) on all events

####################  visualize results  #######################################

# ------- show algae species

df_only_algae_rf_pred<-fread(file = "/home/rstudio/marine_data/results/pred_algae_species_mlr_from_pred_rf_algae.csv",check.names = F)

df_only_algae_rf_pred<-fread(file = "/home/rstudio/marine_data/results/pred_algae_species_rf_from_pred_rf_algae.csv",check.names = F)


# UMAP to visualize algae

umap_df<-umap_analysis(df_only_algae_rf_pred,downs = 0.009,n_neigh = 15,min_dist = 0.5,n_epch = 200,
                       spread_val = 6,metric_select = "cosine") # best tuning parameters

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "Species",size_p = 0.5,
                 legend_title = "Species",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25)


get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "species_predicted",size_p = 0.5,
                 legend_title = "Species",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25)


get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "class",size_p = 0.5,
                 legend_title = "Classes",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25)


#------- bar plots F1 scores
# mlr
species_names<-str_remove(names(all_scores),"Class: ")
f1_scores<-as.vector(all_scores)
df_scores<-as.data.frame(cbind(species_names,f1_scores),stringsAsFactors = F)
df_scores<-df_scores[order(df_scores$f1_scores,decreasing = T),]
df_scores$f1_scores<-as.numeric(df_scores$f1_scores)
df_scores$species_names<-factor(df_scores$species_names,levels = df_scores$species_names)
species_names_order_mlr<-df_scores$species_names

# rf
species_names<-str_remove(names(all_scores),"Class: ")
f1_scores<-as.vector(all_scores)
df_scores<-as.data.frame(cbind(species_names,f1_scores),stringsAsFactors = F)
df_scores<-df_scores[order(df_scores$f1_scores,decreasing = T),]
df_scores$f1_scores<-as.numeric(df_scores$f1_scores)
df_scores$species_names<-factor(df_scores$species_names,levels = species_names_order_mlr)
df_scores$species_names
get_bar_plot(df = df_scores,x_var = "species_names",y_var = "f1_scores",size_axis_text=18,size_title_x=20,size_title_y=20,
             x_lab = "Species",y_lab = "F1 score")



####################  Identification of unknown species first batch (removed species from training)  #######################################

# in this section we to develop another type of first model. The model stills aims to separate between algae and non algae (class variable is the Class column)
# but without some species in the training data.
train_data<-fread(file = "/home/rstudio/marine_data/results/train_data_all_m.csv")
test_data<-fread(file = "/home/rstudio/marine_data/results/test_data_all_m.csv")

# filter out some species
species_to_rm<-c("RCC1511","RCC1507","RCC539")
inds_to_rm<-which(train_data$Species %in% species_to_rm)
train_data_filtered<-train_data[-inds_to_rm,]

#---- random forest-----

# training 

out_train<-data_train(data = train_data_filtered,method = "rf",n_cores = 16,n_samples_down = 5000, method_ctrl = "cv",
                      n_tree = c(5,10,15),rep_seed = 123,k_cv = 5, n_rep = 1,ref_var = "type")

vec_max_values<-out_train$vec_max_values

max_value<-out_train$max_value

best_model<-out_train$best_model

class(best_model)
saveRDS(object = best_model,file =  "/home/rstudio/marine_data/results/final_rf_model_alg_vs_noise_filtered_species.rds")

# prediction

y_pred<-data_prediction(input_model = best_model,test_df = test_data)

test_score<-get_F1_score(y_truth = test_data$type,y_pred = y_pred,pos_label = "algae")

y_prob<-data_prediction(input_model = best_model,test_df = test_data,prob = T)

# get confusion matrix
y_truth<-ifelse(test_data$type=="non_algae","unknown","algae") # Better to use the term unknown instead of non-algae (since they could be algae, just we do not know)
y_pred<-ifelse(y_pred=="non_algae","unknown","algae")
unique(y_truth)
unique(y_pred)
conf_plot<-get_confusion_matrix(y_pred = y_pred,y_true = y_truth, show_legend = T)


# add prediction to test data

test_data$type_predicted<-y_pred

test_data<-cbind(test_data,y_prob)

test_data$type<-ifelse(test_data$type=="non_algae","unknown","algae")

# export results
fwrite(test_data, file = "/home/rstudio/marine_data/results/pred_algae_rf_filtered_train.csv",row.names = F)

#------- visualize results ------------
df_pred_algae_rf_filtered_train<-fread(file = "/home/rstudio/marine_data/results/pred_algae_rf_filtered_train.csv",check.names = F)

# UMAP for visualize algae

umap_df<-umap_analysis(df_pred_algae_rf_filtered_train,downs = 0.005,n_neigh = 15,min_dist = 0.1,n_epch = 200,spread_val = 5,metric_select = "cosine") 

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "type",size_p = 0.5,
                 legend_title = "Gate",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25)

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "Species",size_p = 0.5,
                 legend_title = "Gate",show_legend = T,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25,col_single_algae_species = c("RCC1511","RCC1507","RCC539"))

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "Species",size_p = 0.5,
                 legend_title = "Gate",show_legend = T,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25)

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "type_predicted",size_p = 0.5,
                 legend_title = "Gate",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25)


get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "algae_prob",size_p = 0.5,
                 legend_title = "Probability for Algae",show_legend = T,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25,color_as_prob = T)

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "non-algae_prob",size_p = 0.5,
                 legend_title = "Probability for noise",show_legend = T,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25,color_as_prob = T)


values_prob<-unlist(umap_df[umap_df$Species=="RCC539","algae_prob",with=F])
mean(values_prob)
hist(values_prob)

values_prob<-unlist(umap_df[umap_df$Species=="RCC1511","algae_prob",with=F])
mean(values_prob)
hist(values_prob)

values_prob<-unlist(umap_df[umap_df$Species=="RCC1507","algae_prob",with=F])
mean(values_prob)
hist(values_prob)

values_prob<-unlist(umap_df[umap_df$Species=="RCC539","non-algae_prob",with=F])
mean(values_prob)

######################################################################################################
################### import data and experiment information second batch #######################################
######################################################################################################

# Get csv data 
all_paths_data<-list.files("~/marine_data/Aggregated data for paper/Pipeline/second_batch_2025/data_csv_all",full.names = T,pattern = ".csv")
#paths_all_data<-all_paths_data_original[grep("m1\\.csv|m2\\.csv|m3\\.csv|na_m1\\.csv|na_m2\\.csv|na_m3\\.csv",all_paths_data_original)]

df_total_data<-get_all_data_csv(paths_data = all_paths_data)
object_size(df_total_data) # check memory occupied by object

unique(df_total_data$type)


######################################################################################################
#####################################  Quality checking second batch #######################################
######################################################################################################


df_counts_species<-df_total_data[, .N, by = .(measure, Species)]
df_counts_species$measure<-str_remove(df_counts_species$measure,"second_batch_")
df_counts_species$measure<-str_remove(df_counts_species$measure,".csv")

inds<-grep("non_algae",df_counts_species$measure)
df_counts_species_algae_only<-df_counts_species[-inds,] # algae only
df_counts_species_na_only<-df_counts_species[inds,] # non algae only

# checking number of events across measures

table(df_total_data$Species)


get_line_plot(df = df_counts_species_algae_only, x_var = "measure",y_var = "N",color_var = "Species",
              legend_title = "Species",show_legend = T,x_lab = "Measures",y_lab = "Number of events",size_title_x = 30,size_title_y = 30,size_axis_text = 25)

get_line_plot(df = df_counts_species_algae_only, x_var = "measure",y_var = "N",color_var = "Species",
              legend_title = "Species",show_legend = F,x_lab = "Measures",y_lab = "Number of events",
              size_title_x = 30,size_title_y = 30,size_axis_text = 25,col_single_species = "RCC1502")

# checking number of events across Species based on measure

df_counts_species_algae_only_m1<-df_counts_species_algae_only[df_counts_species_algae_only$measure=="algae_m1",]

df_counts_species_algae_only_m1<-df_counts_species_algae_only_m1[order(df_counts_species_algae_only_m1$N,decreasing = T),]

df_counts_species_algae_only_m1$Species<-factor(df_counts_species_algae_only_m1$Species, level=df_counts_species_algae_only_m1$Species)

get_bar_plot(df = df_counts_species_algae_only_m1,x_var = "Species",y_var = "N",size_axis_text=18,
             size_title_x=20,size_title_y=20,x_lab = "Species",y_lab = "Cells count")

# clean dataset

df_total_data$measure<- str_remove(df_total_data$measure,"second_batch_")
df_total_data$measure<- str_remove(df_total_data$measure,".csv")

# add class for species

df_species_info<-read.csv(file = "~/marine_data/Aggregated data for paper/Pipeline/species list/species list-Roscoff_first_sheet.csv")

mapping <- setNames(df_species_info$Class,df_species_info$Roscoff.Culture.Collection.Identifier)

df_total_data$Class <- mapping[df_total_data$Species]

unique(df_total_data$Class)

# df_total_data contains algae + non algae (i.e., noise) for all measures (day1, day 2, day 7) considering all events for each sample
fwrite(x = df_total_data,file = "/home/rstudio/marine_data/results/df_total_data_second_batch.csv")

######################################################################################################
#####################################  dimensional reduction analysis second batch  #######################################
######################################################################################################

df_total_data<-fread(file = "/home/rstudio/marine_data/results/df_total_data.csv",check.names = F)

df_total_data_2nd<-fread(file = "/home/rstudio/marine_data/results/df_total_data_second_batch.csv",check.names = F)

inds_check<-which(df_total_data_2nd$Species %in% df_total_data$Species)

species_not_in_first_batch<-unique(df_total_data_2nd$Species[-inds_check]) 

#-------------------------- algae ------------------------------
df_total_data_2nd$Time<-NULL # remove time channel

df_total_data_2nd_alg<-df_total_data_2nd[df_total_data_2nd$type=="algae"]
df_total_data_2nd_alg_copy<-copy(df_total_data_2nd_alg) # to keep info about metadata

#--------- PCA

# single measure
out_pca<-pca_analysis(df = df_total_data_2nd_alg,select_measure = "algae_m1") #either m1 or m2 or m3,  by default n_pc=10

df_expr_pc<-out_pca$df_expr_pc

df_expr_first_10_pc<-out_pca$df_expr_first_n_pc # first 10 components for UMAP

out_pca$summary_pc_results

# all measure 

out_pca<-pca_analysis(df = df_total_data_2nd_alg,select_measure = "all") # all measures

df_expr_pc<-out_pca$df_expr_pc

df_expr_first_10_pc<-out_pca$df_expr_first_10_pc # first 10 components for UMAP

out_pca$summary_pc_results

# downsample for easy visualization
set.seed(123)
ind_res<-createDataPartition(y=factor(df_expr_pc$Species),times = 1,p = 0.05) # downsampling by keeping same probability distribution
ind_res<-ind_res$Resample1
df_expr_pc_downsampled<-df_expr_pc[ind_res,]


# species layer
get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (13%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25) # m1


get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (14%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25) # m2


get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (12%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25) # m3

# color single species
get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (13%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25,
                 col_single_species = "RCC6336") 

get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (13%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25,
                 col_single_species = "RCC950") 

get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (13%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25,
                 col_single_species = "RCC1507") 

get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (13%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25,
                 col_single_species = "RCC1511") 

get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (13%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25,
                 col_single_species = "RCC539") 


# measure layers
get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "measure",shape_var = "measure",size_p = 1.5,
                 legend_title = "measure",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (13%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25)

#----------- UMAP


df_total_data_2nd_alg_m1<-df_total_data_2nd_alg[df_total_data_2nd_alg$measure=="algae_m1",]

umap_df<-umap_analysis(df_total_data_2nd_alg,downs = 0.01,n_neigh = 15,min_dist = 0.3,n_epch = 200,metric_select = "cosine") # 3% downsampling for single measures, 1% for all measures

umap_df<-umap_analysis(df_total_data_2nd_alg_m1,downs = 0.03,n_neigh = 15,min_dist = 0.3,n_epch = 200,metric_select = "cosine") # 3% downsampling for single measures, 1% for all measures


#best one m1 : downs=0.03, n_neight=15,dist=0.3,epocs=200, metric=cosine
#best one all m: downs=0.01, n_neight=15,dist=0.3,epocs=200, metric=cosine

# species layer
get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",size_title_x = 30,size_title_y = 30,size_axis_text = 25)

# class layer
get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "class",size_p = 0.05,
                 legend_title = "Class",show_legend = T,x_lab = "UMAP1",y_lab = "UMAP2",size_title_x = 30,size_title_y = 30,size_axis_text = 25)

# color single species
get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",size_title_x = 30,size_title_y = 30,size_axis_text = 25,
                 col_single_species = "RCC1435")

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",size_title_x = 30,size_title_y = 30,size_axis_text = 25,
                 col_single_species = "RCC950")

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",size_title_x = 30,size_title_y = 30,size_axis_text = 25,
                 col_single_species = "RCC69")

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",size_title_x = 30,size_title_y = 30,size_axis_text = 25,
                 col_single_species = "RCC539")

# measure layer
umap_df$measure<- str_remove(umap_df$measure,"algae_")
get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "measure",shape_var = "measure",size_p = 1.5,
                 legend_title = "Measure",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",size_title_x = 30,size_title_y = 30,size_axis_text = 25)

# Story: We started with PCA to analyze general patterns. PCA did not identify well separated structure when considering the species layer, suggesting that data may show more complex relationships. We retained the 10 first PCA components and we applied a UMAP algorithm

# The umap algorithm showed better well separated structures  with the species layer.

# both PCA and Umap show overlapping structures with the measurement layer, suggesting that the measurement of the data does not represent important information in the data.


#--------- non algae ---------
df_total_data_2nd_na<-df_total_data_2nd[df_total_data_2nd$type=="non-algae"]
df_total_data_2nd_na_copy<-copy(df_total_data_2nd_na)

#--------- PCA

# single measure
out_pca<-pca_analysis(df = df_total_data_2nd_na,select_measure = "m1") #either m1 or m2 or m3

df_expr_pc<-out_pca$df_expr_pc

df_expr_first_10_pc<-out_pca$df_expr_first_10_pc # first 10 components for UMAP

out_pca$summary_pc_results


# all measure 

out_pca<-pca_analysis(df = df_total_data_2nd_na,select_measure = "all") # all measures

df_expr_pc<-out_pca$df_expr_pc

df_expr_first_10_pc<-out_pca$df_expr_first_10_pc # first 10 components for UMAP

out_pca$summary_pc_results


# downsample for easy visualization
ind_res<-createDataPartition(y=factor(df_expr_pc$Species),times = 1,p = 0.01) # downsampling by keeping same probability distribution
ind_res<-ind_res$Resample1
df_expr_pc_downsampled<-df_expr_pc[ind_res,]


# species layer
get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (13%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25) # m1


get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (14%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25) # m2


get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "PC1 (70%)",y_lab = "PC2 (12%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25) # m3

# measure layers
get_scatter_plot(df = df_expr_pc_downsampled,x_var = "PC1",y_var = "PC2",color_var = "measure",shape_var = "measure",size_p = 1.5,
                 legend_title = "measure",show_legend = F,x_lab = "PC1 (71%)",y_lab = "PC2 (11%)",size_title_x = 30,size_title_y = 30,size_axis_text = 25)

#----------- UMAP

umap_df<-umap_analysis(df_expr_first_10_pc,downs = 0.01,n_neigh = 10,min_dist = 0.05,n_epch = 500) # 3% downsampling for single measures, 1% for all measures

# species layer
get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "Species",size_p = 0.05,
                 legend_title = "Species",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",size_title_x = 30,size_title_y = 30,size_axis_text = 25)

# measure layer
get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "measure",shape_var = "measure",size_p = 1.5,
                 legend_title = "Measure",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",size_title_x = 30,size_title_y = 30,size_axis_text = 25)


######################################################################################################
#####################################  Algae vs noise analysis (first model) second batch  ##########################
######################################################################################################
df_total_data<-fread(file = "/home/rstudio/marine_data/results/df_total_data.csv",check.names = F)

df_total_data_2nd<-fread(file = "/home/rstudio/marine_data/results/df_total_data_second_batch.csv",check.names = F)
test_data<-df_total_data_2nd
test_data$Time<-NULL # remove time channel
test_data$type<-str_replace(test_data$type,"non-algae","non_algae") 

#---- random forest -----
model_a<-readRDS("/home/rstudio/marine_data/results/final_rf_model_alg_vs_noise_all_train_data.rds") # random forest model

# prediction

y_pred<-data_prediction(input_model = model_a,test_df = test_data) 

y_truth<-test_data$type

test_score<-get_F1_score(y_truth = y_truth,y_pred = y_pred,pos_label = "algae")

y_prob<-data_prediction(input_model = model_a,test_df = test_data,prob = T) 

# get confusion matrix
y_truth<-ifelse(y_truth=="non_algae","unknown","algae") # Better to use the term unknown instead of non-algae (since they could be algae, just we do not know)
y_pred<-ifelse(y_pred=="non_algae","unknown","algae")
unique(y_truth)
unique(y_pred)
conf_plot<-get_confusion_matrix(y_pred = y_pred,y_true = y_truth, show_legend = T)

# add prediction to test data

test_data$type_predicted<-y_pred

test_data<-cbind(test_data,y_prob)

# export results
fwrite(test_data, file = "/home/rstudio/marine_data/results/pred_algae_rf_second_batch.csv",row.names = F)

#---- isolation forest -----

# make sure reticulate uses the Python environment to load the scikit learn model
use_python("/usr/bin/python3", required = TRUE)  # adjust path if needed

# import joblib
joblib <- import("joblib")

# dump your Python model object to a file
model_a<-joblib$load("/home/rstudio/marine_data/results/final_if_model_alg_vs_noise_all_train_data.pkl")

# prediction

y_pred<-data_prediction(input_model = model_a,test_df = test_data,normalize = T) # isolation forest was trained on normalized data, so also test data needs to be normalized

y_truth<-test_data$type

y_truth<-ifelse(y_truth=="algae","1","0")

test_score<-get_F1_score(y_truth = y_truth,y_pred = y_pred,pos_label = "1")

y_prob<-data_prediction(input_model = model_a,test_df = test_data,prob = T) 
# they are actually normalized anomaly scores:
# higher values ⇒ more normal (inlier)
# lower values ⇒ more abnormal (outlier)
# in our case, the normal ones are the algae, we trained the isolation forest on algae only data (one class,like svm).

# get confusion matrix
y_pred<-ifelse(y_pred=="1","algae","unknown")
y_truth<-ifelse(y_truth=="1","algae","unknown")
unique(y_truth)
unique(y_pred)
conf_plot<-get_confusion_matrix(y_pred = y_pred,y_true = y_truth, show_legend = T)

# add prediction to test data

test_data$type_predicted<-y_pred

test_data<-cbind(test_data,y_prob)

# export results
fwrite(test_data, file = "/home/rstudio/marine_data/results/pred_algae_if_second_batch.csv",row.names = F)

####################  visualize results  #######################################


# ------- show algae vs non algae

# UMAP for visualize algae using RF

df_pred_algae_rf_2nd<-fread(file = "/home/rstudio/marine_data/results/pred_algae_rf_second_batch.csv",check.names = F)

# #umap_df<-umap_analysis(df_pred_algae_if,downs = 0.001,n_neigh = 10,min_dist = 0.05,n_epch = 500) # 1% downsampling for single measures, 1% for all measures
# umap_df<-umap_analysis(df_pred_algae_rf,downs = 0.001,n_neigh = 15,min_dist = 0.3,n_epch = 500,metric_select = "cosine") 
# umap_df<-umap_analysis(df_pred_algae_rf,downs = 0.001,n_neigh = 15,min_dist = 0.5,n_epch = 500,metric_select = "cosine") 
# umap_df<-umap_analysis(df_pred_algae_rf,downs = 0.001,n_neigh = 15,min_dist = 0.5,n_epch = 500,spread_val = 3,metric_select = "cosine") 

# best tuning 
umap_df<-umap_analysis(df_pred_algae_rf_2nd,downs = 0.001,n_neigh = 15,min_dist = 0.3,n_epch = 200,spread_val = 5,metric_select = "cosine") 

umap_df$type<-ifelse(umap_df$type=="non_algae","unknown","algae")
umap_df$type_predicted<-ifelse(umap_df$type_predicted=="non_algae","unknown","algae")

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "type",size_p = 0.5,
                 legend_title = "Gate",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25) # reference plot

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "type_predicted",size_p = 0.5,
                 legend_title = "Gate",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25)

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "Species",size_p = 0.5,
                 legend_title = "Gate",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25,col_single_species = species_not_in_first_batch) # reference plot

# UMAP for visualize algae using IF
df_pred_algae_if_2nd<-fread(file = "/home/rstudio/marine_data/results/pred_algae_if_second_batch.csv",check.names = F)

umap_df<-umap_analysis(df_pred_algae_if_2nd,downs = 0.001,n_neigh = 15,min_dist = 0.3,n_epch = 200,spread_val = 5,metric_select = "cosine") 

umap_df$type<-ifelse(umap_df$type=="non_algae","unknown","algae")


get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "type_predicted",size_p = 0.5,
                 legend_title = "Gate",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25)


####################  identification of unknown species (probabilities) second batch  #######################################

# we check the name of the species absent in the first batch
df_total_data<-fread(file = "/home/rstudio/marine_data/results/df_total_data.csv",check.names = F) # first batch data

df_total_data_2nd<-fread(file = "/home/rstudio/marine_data/results/df_total_data_second_batch.csv",check.names = F) # second batch data

inds_check<-which(df_total_data_2nd$Species %in% df_total_data$Species)

species_not_in_first_batch<-unique(df_total_data_2nd$Species[-inds_check]) 


# check prob predictions on unknown species
df_pred_algae_rf_2nd<-fread(file = "/home/rstudio/marine_data/results/pred_algae_rf_second_batch.csv",check.names = F)

umap_df<-umap_analysis(df_pred_algae_rf_2nd,downs = 0.001,n_neigh = 15,min_dist = 0.3,n_epch = 200,spread_val = 5,metric_select = "cosine") 

umap_df$type<-ifelse(umap_df$type=="non_algae","unknown","algae")
umap_df$type_predicted<-ifelse(umap_df$type_predicted=="non_algae","unknown","algae")

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "type",size_p = 0.5,
                 legend_title = "Gate",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25) # reference plot

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "type_predicted",size_p = 0.5,
                 legend_title = "Gate",show_legend = T,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25)


get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "Species",size_p = 0.5,
                 legend_title = "Species",show_legend = T,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25,col_single_algae_species = species_not_in_first_batch)


get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "algae_prob",size_p = 0.5,
                 legend_title = "Probability for Algae",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25,
                 color_as_prob = T)

######################################################################################################
#####################################  Algae species analysis (second model) second batch  #######################################
######################################################################################################

model_b_rf<-readRDS("/home/rstudio/marine_data/results/final_rf_model_alg_species.rds")

model_b_mlr<-readRDS("/home/rstudio/marine_data/results/final_mlr_model_alg_species.rds")


df_pred_algae_rf<-fread(file = "/home/rstudio/marine_data/results/pred_algae_rf_second_batch.csv",check.names = F)

df_only_algae_rf_pred<-df_pred_algae_rf[df_pred_algae_rf$type_predicted=="algae",] # predicting on predicted algae events

ind_species_to_remove<-which(df_only_algae_rf_pred$Species %in% species_not_in_first_batch)

df_only_algae_rf_pred_filtered<-df_only_algae_rf_pred[-ind_species_to_remove,]

df_only_algae_rf_pred_removed<-df_only_algae_rf_pred[ind_species_to_remove,]

y_pred_from_pred<-data_prediction(input_model = model_b_rf,test_df = df_only_algae_rf_pred_filtered) # change variable based on chosen model

test_score_out<-get_F1_score(y_truth = df_only_algae_rf_pred_filtered$Species,y_pred = y_pred_from_pred,multiclass = T)

test_score<-test_score_out$macro_score

all_scores<-test_score_out$all_scores

# add prediction to test data

df_only_algae_rf_pred_filtered$species_predicted<-y_pred_from_pred
df_only_algae_rf_pred_removed$species_predicted<-"Species_not_trained"
df_only_algae_rf_pred_final<-rbind(df_only_algae_rf_pred_filtered,df_only_algae_rf_pred_removed)

# export results

fwrite(df_only_algae_rf_pred_final, file = "/home/rstudio/marine_data/results/pred_algae_species_mlr_from_pred_rf_algae_2nd_batch.csv",row.names = F) # species predicted (using mlr model) on algae events predicted with rf

fwrite(df_only_algae_rf_pred_final, file = "/home/rstudio/marine_data/results/pred_algae_species_rf_from_pred_rf_algae_2nd_batch.csv",row.names = F) # species predicted (using mlr model) on algae events predicted with rf

#---------------- visualize results -------------------
df_only_algae_rf_pred<-fread(file = "/home/rstudio/marine_data/results/pred_algae_species_mlr_from_pred_rf_algae_2nd_batch.csv",check.names = F)

df_only_algae_rf_pred<-fread(file = "/home/rstudio/marine_data/results/pred_algae_species_rf_from_pred_rf_algae_2nd_batch.csv",check.names = F)

# UMAP to visualize algae
df_only_algae_rf_pred_only_trained_species<-df_only_algae_rf_pred[df_only_algae_rf_pred$species_predicted!="Species_not_trained",]

umap_df<-umap_analysis(df_only_algae_rf_pred_only_trained_species,downs = 0.001,n_neigh = 15,min_dist = 0.5,n_epch = 200,
                       spread_val = 6,metric_select = "cosine") # best tuning parameters

get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "Species",size_p = 0.5,
                 legend_title = "Species",show_legend = T,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25)


get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "species_predicted",size_p = 0.5,
                 legend_title = "Species",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25)


get_scatter_plot(df = umap_df,x_var = "UMAP1",y_var = "UMAP2", color_var = "class",size_p = 0.5,
                 legend_title = "Classes",show_legend = F,x_lab = "UMAP1",y_lab = "UMAP2",
                 size_title_x = 30,size_title_y = 30,size_axis_text = 25)


#------- bar plots F1 scores
# mlr
species_names<-str_remove(names(all_scores),"Class: ")
f1_scores<-as.vector(all_scores)
df_scores<-as.data.frame(cbind(species_names,f1_scores),stringsAsFactors = F)
df_scores<-df_scores[order(df_scores$f1_scores,decreasing = T),]
df_scores$f1_scores<-as.numeric(df_scores$f1_scores)
df_scores$species_names<-factor(df_scores$species_names,levels = df_scores$species_names)
species_names_order_mlr<-df_scores$species_names

# rf
species_names<-str_remove(names(all_scores),"Class: ")
f1_scores<-as.vector(all_scores)
df_scores<-as.data.frame(cbind(species_names,f1_scores),stringsAsFactors = F)
df_scores<-df_scores[order(df_scores$f1_scores,decreasing = T),]
df_scores$f1_scores<-as.numeric(df_scores$f1_scores)
df_scores$species_names<-factor(df_scores$species_names,levels = species_names_order_mlr)
df_scores$species_names

get_bar_plot(df = df_scores,x_var = "species_names",y_var = "f1_scores",size_axis_text=18,size_title_x=20,size_title_y=20,
             x_lab = "Species",y_lab = "F1 score")
