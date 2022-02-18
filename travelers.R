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
#read the tree file 
#insert group names (states) into the tree
#preparation before applying ace function
###############################################


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


###############################################
#apply ace function 
###############################################

t <- tree3full$tip.label
t <- tree1full$group[1:262]
ansfull <- ace(t,treefull, type = 'd')

###############################################
#use lik.anc data to find the internal states
###############################################

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

###############################################
#plot the full tree
###############################################
tree4full <- as.treedata(tree1full)

true_tree<-ggtree(tree4full,right=TRUE,mrsd = "2015-01-31",options(ignore.negative.edge=TRUE))+
  geom_point(aes(color=group))+theme_tree2()+ggtitle("Full Tree, Full Reconstruction") + theme(text = element_text(size = 10))

################################################
#use true_tree to determine which tip where there was a recent state change (parent or parent's parent)

#extract data from the true tree
###############################################

dt <- true_tree$data

tip <- c()

#for loop
#parent & parent's parent
#for(i in 1:262){
#  state <- dt$group[i]
#  if(dt$group[dt$parent[i]]!=state|dt$group[dt$parent[dt$parent[i]]]!=state){
#    tip <- append(tip,i)
#  }
#}

#parent 
#for(i in 1:262){
#  state <- dt$group[i]
#  if(dt$group[dt$parent[i]]!=state){
#    tip <- append(tip,i)
#  }
#}

#3 generations
for(i in 1:262){
  state <- dt$group[i]
  if(dt$group[dt$parent[i]]!=state|dt$group[dt$parent[dt$parent[i]]]!=state|dt$group[dt$parent[dt$parent[dt$parent[i]]]]!=state){
    tip <- append(tip,i)
  }
}

for(i in 1:length(tip)){
  print(i)
  print(dt$group[tip[i]])
}
#drop unknown and Mali cases
#not_consider <- c(6,7,23,24)
#not_consider <- c(2,3,12,13)
not_consider <- c(8,12,31,32)
tip <- tip[-not_consider]


tree<- read.nexus("meannode.tree")
tree$root.time<-max(dist.nodes(tree))
groups <- read.table("table.txt", sep='',
                     col.names = c('ID','group'),
                     header = FALSE, stringsAsFactors = FALSE)
gdata <- as.data.frame(groups)

tree1 <- full_join(as_tibble(tree),gdata,by = c('label'='ID'))
tree3 <- as.phylo(tree1)


drop_traverlers <- function(tree1,tree3,tip){
  d <- c()
  rm_chr <- c()
  for(i in 1:length(tip)){
    d[i] <- tree3$tip.label[tip[i]]
    rm_chr[i] <- tree1$group[tip[i]]
  }
  name <- tip
  tip <- gdata$group
  tip <- tip[-name]
  
  tree5 <- drop.tip(tree3, d,trim.internal = TRUE)
  
  for(i in 1:length(tree5$edge.length)){
    if(tree5$edge.length[i] < 0){
      tree5$edge.length[i] = 0
    }
  }
  results <- list("tree"=tree5, "tip"=tip,"drop"=d)
  return(results)
}

results <- drop_traverlers(tree1=tree1,tree3=tree3,tip=tip)
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

down_travelers<-ggtree(tree_cool,right=TRUE,mrsd = "2015-01-31",options(ignore.negative.edge=TRUE))+
  geom_point(aes(color=group))+theme_tree2()+ggtitle("Downsample Tree, Full Reconstruction")+ theme(text = element_text(size = 10))

#(true_tree + down_travelers)+plot_layout(guides = "collect")

#############################################################
treephylofull <- as.phylo(tree1full)
treedropfull <- drop.tip(treephylofull, drop,trim.internal=TRUE,subtree=TRUE)

try_drop <- drop.tip(try_tree,drop,trim.internal=TRUE)
t_drop <-as_tibble(try_drop)
#try_cool <- as.treedata(try_drop)


drop_data <- gdata
for(i in 1:length(drop)){
  drop_data <- drop_data[-c((grep(drop[i],drop_data$ID,fixed = TRUE))),]
}

ID <- t_drop$label
group <- c(tip,try_drop$node.label)

df_drop <- data.frame(ID,group)
t_drop <- full_join(as_tibble(try_drop),df_drop,by=c('label'='ID'))
try_123 <- as.treedata(try_drop)


true_down_traverlers<-ggtree(try_123,right=TRUE,mrsd = "2015-01-31",options(ignore.negative.edge=TRUE))+
  geom_point(aes(color=group))+theme_tree2()+ggtitle("Preferentially sampling traverlers")+ theme(text = element_text(size = 10))


(true_tree + down_travelers+true_down_traverlers)+plot_layout(guides = "collect")


