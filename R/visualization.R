
get_scatter_plot<-function(df,x_var="PC1",y_var="PC2",color_var="Species",size_p=1,size_axis_text=18,
                           size_title_x=20,size_title_y=20,legend_title="Species",show_legend=F,
                           x_lab="PC1",y_lab="PC2",shape_var=NULL,col_single_algae_species=NULL,color_as_prob=F){
  if(is.null(col_single_algae_species)==F){
    inds_select_species<-which(df$Species %in% col_single_algae_species & df$type == "algae")
    df$Species[-inds_select_species]<-"Other species"
  }
  if(is.null(shape_var)==T){
    gg_plot<-ggplot(data = df, aes(x=unlist(df[,x_var,with=F]), y=unlist(df[,y_var,with=F]))) + geom_point(aes(color=unlist(df[,color_var,with=F])),size=size_p)
    
  }else{
    gg_plot<-ggplot(data = df, aes(x=unlist(df[,x_var,with=F]), y=unlist(df[,y_var,with=F]))) + geom_point(aes(color=unlist(df[,color_var,with=F]),shape = unlist(df[,shape_var,with=F])),size=size_p)
  }
  
  gg_plot<- gg_plot + theme(axis.text=element_text(size=size_axis_text),axis.title.x = element_text(size = size_title_x,face="bold"),
                                    axis.title.y = element_text(size=size_title_y,face="bold"))
  
  gg_plot <- gg_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  gg_plot<- gg_plot +  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  gg_plot<-gg_plot + theme(legend.key.size = unit(1, 'cm'),legend.title = element_text(size=20),legend.text = element_text(size=20))
  gg_plot<-gg_plot + guides(color = guide_legend(override.aes = list(size = 10))) + labs(color=legend_title,shape=legend_title)
  if(color_var=="measure"){
    df$measure<-str_remove(df$measure,"na_")
    gg_plot<-gg_plot + scale_color_manual(values= c("m1"="#009E73","m2"="#0072B2","m3"="#E69F00"))
  }else if(color_var=="type" || color_var=="type_predicted"){
    gg_plot<-gg_plot + scale_color_manual(values= c("algae"="#E69F00","unknown"="#0072B2"))
    
  }else if(color_var=="class"){
    gg_plot<-gg_plot + scale_color_manual(values= c("Chlorophyceae"="#355517","Eustigmatophyceae"="#FF7F00",
                                                    "Cryptophyceae"="#FDBF6F","Trebouxiophyceae"="#156804","Prymnesiophyceae"="#1F78B4",
                                                    "Chlorodendrophyceae"="#B2DF8A","Dinophyceae"="#CAB2D6","Cyanophyceae"="#A6CEE3",
                                                    "Ulvophyceae"="#00CC00","Euglenophyceae"="#FF00FF","Chromeridophyceae"="#7E3F1C",
                                                    "Olisthodiscophyceae"="#FB9A99","Bacillariophyceae"="#FFFF00","Mamiellophyceae"="#6A3D9A",
                                                    "Mediophyceae"="#E31A1C","Dictyocophyceae"="#00CC99"))
  }
  
  if(show_legend==F){
    gg_plot<- gg_plot + theme(legend.position="none")
  }
  gg_plot<-gg_plot + xlab(x_lab) + ylab(y_lab)
  if(color_as_prob==T){
    gg_plot<-gg_plot + scale_color_gradientn(colors = c("lavender", "orchid", "mediumorchid", "purple", "purple4"))  
      
  }
  return(gg_plot)
  
}

get_line_plot<-function(df,x_var="PC1",y_var="PC2",color_var="Species",size_axis_text=18,size_title_x=20,size_title_y=20,legend_title="Species",show_legend=F,
                           x_lab="PC1",y_lab="PC2",col_single_species=NULL){
  

  gg_plot<-ggplot(df, aes(x=unlist(df[,x_var,with=F]),y=unlist(df[,y_var,with=F]))) + geom_line(aes(color=unlist(df[,color_var,with=F]),group=unlist(df[,color_var,with=F]),
                                                                                                    linetype = unlist(df[,color_var,with=F]),linewidth = unlist(df[,color_var,with=F])))
  
  gg_plot<- gg_plot + theme(axis.text=element_text(size=size_axis_text),axis.title.x = element_text(size = size_title_x,face="bold"),
                            axis.title.y = element_text(size=size_title_y,face="bold"))
  
  gg_plot <- gg_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  
  if(show_legend==F){
    gg_plot<- gg_plot + theme(legend.position="none")
  }
  
  gg_plot<-gg_plot + xlab(x_lab) + ylab(y_lab)
  if(is.null(col_single_species)==F){
    # set color for selected line
    vec_color<-rep("#0072B2",length(unique(df$Species)))
    df_colors<-as.data.frame(cbind(unique(df$Species),vec_color),stringsAsFactors=F)
    colnames(df_colors)<-c("Species","color")
    ind<-grep(col_single_species,df_colors$Species)
    df_colors$color[ind]<-"#E69F00"
    color_vector <- setNames(df_colors$color, df_colors$Species)
    print(color_vector)
    gg_plot<- gg_plot + scale_color_manual(values=color_vector)
    
    # set linetype for selected line
    
    vec_linetype<-rep("solid",length(unique(df$Species)))
    df_linetypes<-as.data.frame(cbind(unique(df$Species),vec_linetype),stringsAsFactors=F)
    colnames(df_linetypes)<-c("Species","linetype")
    ind<-grep(col_single_species,df_linetypes$Species)
    df_linetypes$linetype[ind]<-"dashed"
    linetype_vector <- setNames(df_linetypes$linetype, df_linetypes$Species)
    gg_plot<- gg_plot + scale_linetype_manual(values=linetype_vector)
    
    # set linewidth for selected line
    
    vec_linewidth<-rep(0.5,length(unique(df$Species)))
    df_linewidths<-as.data.frame(cbind(unique(df$Species),vec_linewidth),stringsAsFactors=F)
      
    colnames(df_linewidths)<-c("Species","linewidth")
    df_linewidths$linewidth<-as.numeric(df_linewidths$linewidth)
    ind<-grep(col_single_species,df_linewidths$Species)
    df_linewidths$linewidth[ind]<-2
    linewidth_vector <- setNames(df_linewidths$linewidth, df_linewidths$Species)
    print(linewidth_vector)
    gg_plot<- gg_plot + scale_linewidth_manual(values=linewidth_vector)
    
    
  }else{
    vec_linetype<-rep("solid",length(unique(df$Species)))
    df_linetypes<-as.data.frame(cbind(unique(df$Species),vec_linetype),stringsAsFactors=F)
    colnames(df_linetypes)<-c("Species","linetype")
    linetype_vector <- setNames(df_linetypes$linetype, df_linetypes$Species)
    gg_plot<- gg_plot + scale_linetype_manual(values=linetype_vector)
  }
  gg_plot<- gg_plot + scale_x_discrete(expand = c(0.02, 0.02))
  return(gg_plot)
  
}

# function to confusiom matrix of prediction vs origina
get_confusion_matrix<-function(y_pred,y_true,size_axis_text=18,
                               size_title_x=20,size_title_y=20,show_legend=F){
  
  conf_matrix <- table(Predicted = y_pred, Actual = y_true)
  conf_df <- as.data.frame(conf_matrix)
  conf_df$Proportion <- conf_df$Freq / sum(conf_df$Freq)
  conf_df$Label <- sprintf("%d\n(%.1f%%)", conf_df$Freq, conf_df$Proportion * 100)
  gg_plot <-ggplot(data = conf_df, aes(x = Actual, y = Predicted, fill = Proportion)) + geom_tile(color = "white") 
  gg_plot<- gg_plot + scale_fill_gradient(low = "#0072B2", high = "#E69F00") + geom_text(aes(label = Label), color = "white", size = 6) 
  gg_plot<- gg_plot + labs(title = NULL, x = "Reference", y = "Predicted") + theme_minimal() 
  gg_plot<- gg_plot + theme(axis.text=element_text(size=size_axis_text),axis.title.x = element_text(size = size_title_x,face="bold"),
                            axis.title.y = element_text(size=size_title_y,face="bold"), axis.ticks.y.right = element_blank(), 
                            axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank(), axis.text.y.right = element_blank(),
                            panel.grid.major = element_blank())
  gg_plot<-gg_plot + theme(legend.key.size = unit(1, 'cm'),legend.title = element_text(size=20),legend.text = element_text(size=20))
  
  if(show_legend==F){
    gg_plot<- gg_plot + theme(legend.position="none")
  }
  
  
  return(gg_plot)
}

# function to plot line plot results cross val
get_line_plot_v2<-function(df,x_var="PC1",y_var="PC2",color_var="Species",size_axis_text=18,size_title_x=20,size_title_y=20,legend_title="Species",show_legend=F,
                        x_lab="PC1",y_lab="PC2"){
  
  
  gg_plot<-ggplot(df, aes(x=df[,x_var],y=df[,y_var])) + geom_line(aes(color=df[,color_var],group=df[,color_var],
                                                                      linetype = df[,color_var]),linewidth=2)
  
  gg_plot<- gg_plot + theme(axis.text=element_text(size=size_axis_text),axis.title.x = element_text(size = size_title_x,face="bold"),
                            axis.title.y = element_text(size=size_title_y,face="bold"))
  
  gg_plot <- gg_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  
  if(show_legend==F){
    gg_plot<- gg_plot + theme(legend.position="none")
  }
  
  gg_plot<-gg_plot + xlab(x_lab) + ylab(y_lab) + labs(color=legend_title,linetype=legend_title)

  gg_plot<-gg_plot + theme(legend.key.size = unit(2, 'cm'),legend.title = element_text(size=20),legend.text = element_text(size=20))
  
  return(gg_plot)
  
}

# function to plot bar plot
get_bar_plot<-function(df,x_var="PC1",y_var="PC2",size_axis_text=18,size_title_x=20,size_title_y=20,legend_title="Species",show_legend=F,
                       x_lab="PC1",y_lab="PC2"){
  
  df<-as.data.frame(df)
  gg_plot<-ggplot(df, aes(x=df[,x_var],y=df[,y_var])) + geom_bar(stat="identity",fill="#E69F00")
  
  gg_plot<- gg_plot + theme(axis.text=element_text(size=size_axis_text),axis.title.x = element_text(size = size_title_x,face="bold"),
                            axis.title.y = element_text(size=size_title_y,face="bold"))
  
  gg_plot <- gg_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  
  if(show_legend==F){
    gg_plot<- gg_plot + theme(legend.position="none")
  }
  
  gg_plot<-gg_plot + xlab(x_lab) + ylab(y_lab) + labs(color=legend_title,linetype=legend_title)
  
  gg_plot<-gg_plot + theme(legend.key.size = unit(2, 'cm'),legend.title = element_text(size=20),legend.text = element_text(size=20))
  gg_plot<- gg_plot + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  return(gg_plot)
  
}
