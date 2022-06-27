
support = read.delim("data/tree/concord.cf.stat", header=T, comment.char="#") # astral tree
astral.tree <- read.tree("data/tree/concord.cf.tree")

phy = astral.tree
cf.file="data/tree/concord.cf.stat"
value.check = T
legend=T
CF="gCF"
col.palette="RdYlBu"

library(RColorBrewer); 

col.pallette

  # set the plot layout
  if(legend==T){layout(matrix(1:2,ncol=2), width = c(4,1),height = c(1,1))}
  
  # read in the gCF file
  cf.data = read.delim(cf.file, header=T, comment.char="#")
  
  # rename the first column from 'ID' to 'child', 
  # the naming convention for branches in IQTREE is not the same for edges in APE
  # but it is the same as the number of the child node below an edge in APE, so we can use that
  names(cf.data)[1] = "child"
  # add root node, and apply an arbitrary value, this may just be funny behavior because of 
  cf.data[1,"gCF"] <- 100; # arbitrary CF value for 2nd node
  cf.data[1,"sCF"] <- 100;
  cf.data <- rbind(cf.data, c(min(cf.data$child)-1, 100, NA,NA,NA,NA,NA)) # arbitrary CF value for root node
  
  # View(rbind(cf.data, c(min(cf.data$child)-1, 100, NA,NA,NA,NA,NA)))
  
  # create a data frame of all the edges in the tree, listed by (parent, child) nodes
  all.edges <- data.frame(parent=phy$edge[,1], child=phy$edge[,2])
  
  # extract just the branches that lead to tips
  tip.edges <- subset(all.edges, all.edges$child <= Ntip(phy))
  tip.edges$gCF <- 100 # give these branches arbitrary values of 100
  tip.edges$sCF <- 100 # give these branches arbitrary values of 100
  
  
  # extract just the internal branches
  int.edges <- subset(all.edges, all.edges$child > Ntip(phy))
  
  # get the branch-appropriate concordance factors
  int.edges <- left_join(int.edges, cf.data)
  int.edges <- int.edges[,c("parent", "child", "gCF", "sCF")] # subset to just the data we care about
  
  # combine the tip and internal branch information, sort it by the original order!
  combo.edges <- rbind(tip.edges, int.edges)
  combo.edges <- combo.edges[order(match(combo.edges$child, all.edges$child)),]     
  
  # create a color ramp of your choice
  clz <- colorRampPalette(brewer.pal(6, col.palette)); point.colors <- clz(100)
  # if(CF=="gCF"){edge.cols <- point.colors[combo.edges$gCF]} # match the CF values to colors
  if(CF=="gCF"){edge.cols <- c(point.colors[combo.edges$gCF])} 
  if(CF=="sCF"){edge.cols <- point.colors[combo.edges$sCF]} # match the CF values to colors
  
  combo.edges$CF_col <- edge.cols # add it to the data frame
  
  # make the branches we don't have CFs for black or colored
  if (!is.null(terminal.color)) {
    combo.edges[which(combo.edges$parent == min(combo.edges$parent) & combo.edges$child > Ntip(phy)),]$CF_col <- terminal.color
    #combo.edges[which(combo.edges$parent == min(combo.edges$parent)+1 & combo.edges$child > Ntip(phy)),]$gCF_col <- "black"
    combo.edges[which(combo.edges$child <= Ntip(phy)),]$CF_col <- terminal.color
  }
  
  # if you want to change the line types (it doesn't look great so I'm ignoring this)
  #combo.edges$gCF_lty <- "solid"
  #combo.edges[which(combo.edges$child <= Ntip(phy)),]$gCF_lty <- "dotted"
  
  # provide a check to make sure the values all add up
  if (value.check==TRUE) {tip.list <- combo.edges[which(combo.edges$child <= Ntip(phy)),]
  node.list <- table(tip.list$parent); node.list <- which(node.list > 1)
  target.node <- as.numeric(names(node.list)[sample(1:length(node.list),1)])
  chosen.tips <- getDescendants(phy, target.node)
  #t1 <- phy$tip.label[chosen.tips[1]]; t2 <- phy$tip.label[chosen.tips[2]]; tips <- c(t1,t2)
  if(CF=="gCF"){
    true.gCF <- cf.data[which(cf.data$child==target.node),"gCF"]
    target.gCF <- combo.edges[which(combo.edges$child==target.node) ,"gCF"]
  }
  if(CF=="sCF"){
    true.gCF <- cf.data[which(cf.data$child==target.node),"sCF"]
    target.gCF <- combo.edges[which(combo.edges$child==target.node) ,"sCF"]
  }
  #print(paste("node", target.node, ", the parent of tips", t1, "&", t2, "should have a gCF value of", target.gCF))
  if(!true.gCF==target.gCF){print("something funny about your tree, this plot ain't right")}
  else if(true.gCF==target.gCF){print("all good in the hood")}
  }
  
  
  library(stringr)
  subset_tips <- phy$tip.label[stringr::str_detect(phy$tip.label, pattern = "Anilios")]
  subset_tips2 <- phy$tip.label[stringr::str_detect(phy$tip.label, pattern = "Acutotyphlops")]
  subset_tips3 <- phy$tip.label[stringr::str_detect(phy$tip.label, pattern = "Ramphotyphlops")]
  subset_tips4 <- phy$tip.label[stringr::str_detect(phy$tip.label, pattern = "Indotyphlops")]
  phy_subset <- keep.tip(phy, tip = NULL)
  
  
  dev.off()
  # finally, plot it all.
  
  plot.phylo(phy, edge.color=combo.edges$CF_col, edge.width=3); # can add 'edge.lty=combo.edges$gCF_lty' 
  axisPhylo()
  
  
  
  if(legend==T){color.bar(colorRampPalette(brewer.pal(6, col.palette))(100))}