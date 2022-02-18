library("strap")
library("stringr")
library("ips")
library("ggtree")
require(treeio)
require(ggplot2)
require(ggtree)
library("reshape2")
library("ggstance")
library("diversitree")
library("phytools")
library("ggimage")
library("ggpubr")
library("dplyr")

###############################################
#set different downsampling schemes
###############################################
percentage <- c(0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,
                0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95)

dt <- data.frame()

for(k in 1:length(percentage)){
  
  loop <- 0
  accuracy <- c()
  rate <- c()
  while(loop <= 9){
    ##implementation
    #node_count <- 0
    tree<- read.nexus("meannode.tree")
    tree$root.time<-max(dist.nodes(tree))
    groups <- read.table("table.txt", sep='',
                         col.names = c('ID','group'),
                         header = FALSE, stringsAsFactors = FALSE)
    gdata <- as.data.frame(groups)
    
    gui_id <- grep("Guinea",groups$group)
    sie_id <- grep("SierraLeone",groups$group)
    
    tree1 <- full_join(as_tibble(tree),gdata,by = c('label'='ID'))
    tree3 <- as.phylo(tree1)
    
    
    ran_gen <- function(gui_per,sie_per,tree1,tree3){
      
      n_gui <- floor(length(gui_id)*gui_per)
      n_sie <- floor(length(sie_id)*sie_per)
      
      ran_gui <- sample.int(length(gui_id),n_gui)
      ran_sie <- sample.int(length(sie_id),n_sie)
      ran_gui_sort <- sort(ran_gui)
      ran_sie_sort <- sort(ran_sie)
      
      drop_id <- c()
      if(length(ran_gui_sort)>0){
        for(i in 1:length(ran_gui_sort)){
          drop_id[i] = gui_id[ran_gui_sort[i]]
        }
      }
      
      if(length(ran_sie_sort)>0){
        for(i in 1:length(ran_sie_sort)){
          drop_id[i] = sie_id[ran_sie_sort[i]]
        }
      }
      tip <- gdata$group
      
      if(length(drop_id)>=1){
        d <- c()
        rm_chr <- c()
        for(i in 1:length(drop_id)){
          d[i] <- tree3$tip.label[drop_id[i]]
          rm_chr[i] <- tree1$group[drop_id[i]]
        }
        tip <- tip[-drop_id]}
      
      tree5 <- drop.tip(tree3, d,trim.internal = TRUE)
      
      for(i in 1:length(tree5$edge.length)){
        if(tree5$edge.length[i] < 0){
          tree5$edge.length[i] = 0
        }
      }
      results <- list("tree"=tree5, "tip"=tip,"drop"=d)
      return(results)
    }
    ###############################################
    #change here from SIE to GUI
    gui_per <- 0.0
    sie_per <- percentage[i]
    print(sie_per)
    
    ###############################################
    results <- ran_gen(gui_per=gui_per,sie_per=sie_per,tree1=tree1,tree3=tree3)
    tree5 <- results$tree
    tip <- results$tip
    drop <- results$drop
    
    ans <- ace(tip,tree5, type = 'd')
    
    rfp <- ans$rates
    
    ID <- tree5$tip.label
    group <- tip
    df <- data.frame(ID,group)
    
    tree10 <- full_join(as_tibble(tree5),df,by = c('label'='ID'))
    
    j <- apply(ans$lik.anc, 1, which.max)
    len <- length(j)
    num_tip <- len+1
    for(i in num_tip+1:(num_tip+len)){
      if(i-num_tip < num_tip){
        #print(i-num_tip)
        if(j[i-num_tip]==1){
          tree10$group[i] = "Guinea"
        }
        else if(j[i-num_tip]==2){
          tree10$group[i] = "Liberia"
        }
        else if(j[i-num_tip]==3){
          tree10$group[i] = "Mali"
        }
        else if(j[i-num_tip]==4){
          tree10$group[i] = "SierraLeone"
        }
        else{
          tree10$group[i] = "Unknown"
        }
      }
    }
    
    tree_cool <- as.treedata(tree10)
    
    down_random_sie<-ggtree(tree_cool,right=TRUE,mrsd = "2015-01-31",options(ignore.negative.edge=TRUE))+
      geom_point(aes(color=group))+theme_tree2()+ggtitle("Downsample Tree, Downsample Sierra Leone")+ 
      theme(text = element_text(size = 10))
    
    ############################
    
    #ace for non-downsampling
    
    treefull <- read.nexus("meannode.tree")
    
    #set negative branch to 0
    for(i in 1:522){
      if(treefull$edge.length[i] < 0){
        treefull$edge.length[i] = 0
      }
    }
    
    groups <- read.table("table.txt", sep='',
                         col.names = c('ID','group'),
                         header = FALSE, stringsAsFactors = FALSE)
    gdata <- as.data.frame(groups)
    
    
    tree1full <- full_join(as_tibble(treefull),gdata,by = c('label'='ID'))
    tree2full <- as.treedata(tree1full)
    
    tree3full <- as.phylo(tree1full)
    
    #try to make node names
    try_make_node <- tree3full
    
    try_tree <- makeNodeLabel(tree3full,method = "number")
    
    
    #######
    t <- tree3full$tip.label
    
    t <- tree1full$group[1:262]
    
    
    ansfull <- ace(t,treefull, type = 'd')
    jfull <- apply(ansfull$lik.anc, 1, which.max)
    for(i in 263:523){
      if(jfull[i-262]==1){
        tree1full$group[i] = "Guinea"
        try_tree$node.label[i-262] = "Guinea"
      }
      else if(jfull[i-262]==2){
        tree1full$group[i] = "Liberia"
        try_tree$node.label[i-262] = "Liberia"
      }
      else if(jfull[i-262]==3){
        tree1full$group[i] = "Mali"
        try_tree$node.label[i-262] = "Mali"
      }
      else if(jfull[i-262]==4){
        tree1full$group[i] = "SierraLeone"
        try_tree$node.label[i-262] = "SierraLeone"
      }
      else{
        tree1full$group[i] = "Unknown"
        try_tree$node.label[i-262] = "Unknown"
      }
    }
    
    tree4full <- as.treedata(tree1full)
    
    true_tree<-ggtree(tree4full,right=TRUE,mrsd = "2015-01-31",options(ignore.negative.edge=TRUE))+
       geom_point(aes(color=group))+theme_tree2()+ggtitle("Initial Tree")
    
    treephylofull <- as.phylo(tree1full)
    treedropfull <- drop.tip(treephylofull, drop,trim.internal=TRUE,subtree=TRUE)
    
    try_drop <- drop.tip(try_tree,drop,trim.internal=TRUE)
    t_drop <-as_tibble(try_drop)
    
    
    drop_data <- gdata
    for(i in 1:length(drop)){
      drop_data <- drop_data[-c((grep(drop[i],drop_data$ID,fixed = TRUE))),]
    }
    
    ID <- t_drop$label
    group <- c(tip,try_drop$node.label)
    
    df_drop <- data.frame(ID,group)
    t_drop <- full_join(as_tibble(try_drop),df_drop,by=c('label'='ID'))
    try_cool <- as.treedata(try_drop)
    
    true_down_random<-ggtree(try_cool,right=TRUE,mrsd = "2015-01-31",options(ignore.negative.edge=TRUE))+
      geom_point(aes(color=group))+theme_tree2()+ggtitle("True Tree")
    
    
    #ggarrange(p3,p1,p2,nrow = 1 ,common.legend = TRUE,legend = "right")
    
    rate_for_plot <- rfp
    acc <- 0
    for(i in 1:length(group)){
      if(tree_cool[i]$group == group[i]){
        acc <- acc+1
      }
    }
    ac <- acc/length(group)
    
    
    
    loop <- loop+1
    accuracy[loop] <- ac
    rate[loop] <- rate_for_plot
  }
  percent <- rep(as.character(percentage[k]),10)
  df <- data.frame(accuracy,rate,percent,node_location)
  dt <- rbind(dt,df)
}

###############################################
#generating bar plots
###############################################

p1 <- ggplot(dt,aes(x=percent,y=accuracy)) + geom_boxplot()+ggtitle("Downsample SierraLeone")

p2 <- ggplot(dt,aes(x=percent,y=rate)) + geom_boxplot()+ggtitle("Downsample SierraLeone")

ggarrange(p1,p2,nrow = 1)


