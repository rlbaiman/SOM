library(maps)
library(mapdata)
library(maps)
library(mapdata)
library(ggplot2)
library(raster) 
library(rgdal) 
library(ncdf4)
library(kohonen)
library(RColorBrewer)
library(RCurl)
library(lattice)
library(MASS)
library("ape")

###Read in the csv file for matrix from SOM_make_matrix.R. Each row should be a timestep of data 
file<-read.csv('SOM_Resources/latadjust_Anom_SOM_matrix.csv',stringsAsFactors = FALSE)
matrix2<-as.matrix(file)
matrix<-matrix2[1:nrow(matrix2),2:ncol(matrix2)] # remove extra column
rm(file)
rm(matrix2)

##################################################
#Step by step run of SOM function
#############################################################################
#SOM_function(matrix,1,4,4,5,500,.05,0)
dimx=4
dimy=4 
radius=4
iterations=500
alpha_start=.05
alpha_end=0

som_grid<-somgrid(xdim=dimx, ydim=dimy, topo='rectangular', neighbourhood.fct = 'bubble')

som_model<- som(matrix, grid=som_grid, radius=radius ,rlen=iterations, alpha=c(alpha_start,alpha_end), dist.fcts='euclidean')

nodes<-getCodes(som_model)

write.csv(nodes,'out_path')

node_assignment<-som_model$unit.classif 

write.csv(node_assignment,'Ice sheets & Climate Dropbox/Becca Baiman/Fall_2021/SOM_testing_again/node_assignment')

predictions<-predict(som_model, newdata=matrix)
predictions2<-predict(som_model)

plot(som_model, type='dist.neighbours')
plot(som_model, 'mapping')
#################################################################################

###Function that runs SOM and outputs everything we want###
SOM_function<-function(matrix, test_num, dimx, dimy, radius, iterations,alpha_start, alpha_end){
  
  timestart<-Sys.time()
  #Create the SOM grid dimensions that you want
  som_grid<-somgrid(xdim=dimx, ydim=dimy, topo='rectangular', neighbourhood.fct = 'gaussian')
  
  #Run the SOM
  som_model<- som(matrix, grid=som_grid, radius=radius ,rlen=iterations, alpha=c(alpha_start,alpha_end), dist.fcts='euclidean')
  
  name=paste0('Nodes_test',test_num,'.csv')
  nodes<-getCodes(som_model)
  write.csv(nodes,paste('Ice sheets & Climate Dropbox/Becca Baiman/Fall2021/SOM_testing_again/',name))
  
  ## RMSD
  RMSD<-mean(sqrt(som_model6$distances/(11773)))
  
  ##MAD
  node_assignment<-som_model$unit.classif 
  count=0
  for (i in 1:length(node_assignment)){
    MAD_timestep<-sum(abs(matrix[i,]-nodes[node_assignment[i],]))/11773
    count<-count+MAD_timestep
  }
  MAD<-count/4006
  
  #mean distance to closest nodes. Can be used to see how many rlen's were necessary
  par(mar=c(5,6,4,1)+.1)
  plot(som_model, type='changes')
  
  
  #Sammon plot
  lattice.options(default.theme = standard.theme(color = FALSE))
  
  codes<-getCodes(som_model, idx=1)
  hc<-cutree(hclust(dist(codes)),k=9)
  #add a map of thick lines representing clusters
  add.cluster.boundaries(som_model, hc)
  
  dis<-data.frame(dist(nodes))
  sammon<-sammon(dis)
  plot(sammon$points, type='n')
  text(sammon$points, labels=as.character(1:nrow(codes)))
  
  #returns
  print(paste('Root Mean Squared:',round(RMSD,2)))
  print(paste('Mean Absolute Difference',round(MAD,2)))
  runtime<-round(difftime(Sys.time(),timestart),2)
  print(runtime)
  
  
  #rm(som_model)
  rm(name)
  #rm(nodes)
  rm(runtime)
  rm(som_grid)
}



#####Example running function###
SOM_function(matrix,1,4,4,5,500,.05,0)



#############################################################
#Examples of different things you can do with SOM results


#Create the SOM grid dimensions that you want
som_grid<-somgrid(xdim=5, ydim=5, topo='rectangular', neighbourhood.fct = "gaussian")

#Run the SOM
som_model2<- som(matrix, grid=som_grid, radius=5 ,rlen=500, alpha=c(0.05,0), dist.fcts="euclidean")
runtime<-difftime(Sys.time(),timestart)


#Save the map nodes as nodes (each of these will have the same length
#as the matrix data rows)
nodes<-getCodes(som_model2, idx=1)

#Save the nodes as a csv file to map in python
write.csv(nodes,'Ice sheets & Climate Dropbox/Becca Baiman/Summer_2021/Testing_files/SOM_nodes_551.csv')


#list of which node each input data belongs to
node_assignment<-som_model2$unit.classif
write.csv(node_assignment,'Ice sheets & Climate Dropbox/Becca Baiman/Summer_2021/Testing_files/SOM_node_assignments_551.csv')

#list of distances between objects and their winning node
distances<-som_model2$distances

#mean distance to closest nodes. Can be used to see how many rlen's were necessary
par(mar=c(5,6,4,1)+.1)
plot(som_model2, type='changes')

#Shows how many timesteps are classified in each node
plot(som_model, type='count', main='Node Counts')


#U-matrix plot. Shows cluster of nodes that are like each other by coloring 
plot(som_model, type='dist.neighbours')

plot(som_model, type='mapping')

plot(som_model, type='quality')

#Error measures
## quantization error:
qe<-mean(som_model$distances)


######################################################################
#Make a Sammon Mapping of a result

sammon_plot<-function(matrix,xdim, ydim, radius, rlen, alpha_start, alpha_end){
  lattice.options(default.theme = standard.theme(color = FALSE))
  som_grid<-somgrid(xdim=xdim, ydim=ydim, topo='rectangular')
  som<- som(matrix, grid=som_grid, radius=radius,rlen=rlen, alpha=c(alpha_start,alpha_end), )
  
  #hierarchical clustering 
  codes<-getCodes(som, idx=1)
  hc<-cutree(hclust(dist(codes)),k=9)
  #add a map of thick lines representing clusters
  add.cluster.boundaries(som, hc)
  
  dis<-dist(codes)
  sammon<-sammon(dis)
  plot(sammon$points, type='n')
  segment(points[1])
  
  text(sammon$points, pos=1,labels=as.character(1:nrow(codes)))
  #lines from 1
  segments(points[1,1],points[1,2],points[2,1],points[2,2])
  segments(points[1,1],points[1,2],points[4,1],points[4,2])
  
  #lines from 2
  segments(points[2,1],points[2,2],points[3,1],points[3,2])
  segments(points[2,1],points[2,2],points[5,1],points[5,2])
  
  #lines from 3
  segments(points[3,1],points[3,2],points[6,1],points[6,2])
  
  #lines from 4
  segments(points[4,1],points[4,2],points[5,1],points[5,2])
  segments(points[4,1],points[4,2],points[7,1],points[7,2])
  
  #lines from 5
  segments(points[5,1],points[5,2],points[6,1],points[6,2])
  segments(points[5,1],points[5,2],points[8,1],points[8,2])
  
  #lines from 6
  segments(points[6,1],points[6,2],points[9,1],points[9,2])
  
  #lines from 7
  segments(points[7,1],points[7,2],points[8,1],points[8,2])
  
  #lines from 8
  segments(points[8,1],points[8,2],points[9,1],points[9,2])
  
}


sammon_plot(matrix, 4,4,2,200,.01,.001)







lattice.options(default.theme = standard.theme(color = FALSE))
som_grid<-somgrid(xdim=3, ydim=3, topo='rectangular')

som<- som(matrix, grid=som_grid, radius=2,rlen=200, alpha=c(.01,.001), )

#hierarchical clustering 
codes<-getCodes(som_model, idx=1)
hc<-cutree(hclust(dist(codes)),k=9)
#add a map of thick lines representing clusters
add.cluster.boundaries(som_model, hc)

codes<-getCodes(som_model,idx=1)
dis<-dist(codes)
sammon<-sammon(dis)
plot(sammon$points, type='n')


text(sammon$points, labels=as.character(1:nrow(codes)))








hc<-cutree(hclust(dist(codes)),k=9)
#add a map of thick lines representing clusters
add.cluster.boundaries(som, hc)

codes<-getCodes(som_model2,idx=1)
dis<-dist(codes)
sammon<-sammon(dis)
points<-sammon$points


plot(sammon$points,pch=17)
text(sammon$points, pos=1,labels=as.character(1:nrow(codes)))

#lines from 1
segments(points[1,1],points[1,2],points[2,1],points[2,2])
segments(points[1,1],points[1,2],points[4,1],points[4,2])

#lines from 2
segments(points[2,1],points[2,2],points[3,1],points[3,2])
segments(points[2,1],points[2,2],points[5,1],points[5,2])

#lines from 3
segments(points[3,1],points[3,2],points[6,1],points[6,2])

#lines from 4
segments(points[4,1],points[4,2],points[5,1],points[5,2])
segments(points[4,1],points[4,2],points[7,1],points[7,2])

#lines from 5
segments(points[5,1],points[5,2],points[6,1],points[6,2])
segments(points[5,1],points[5,2],points[8,1],points[8,2])

#lines from 6
segments(points[6,1],points[6,2],points[9,1],points[9,2])

#lines from 7
segments(points[7,1],points[7,2],points[8,1],points[8,2])

#lines from 8
segments(points[8,1],points[8,2],points[9,1],points[9,2])






#########
#Test clustering nodes
#########
nodes<-read.table('Ice sheets & Climate Dropbox/Becca Baiman/Fall_2021/Good_nodes/Data/nodes.csv',header=TRUE, sep=',')

d<- dist(nodes)
cah <- hclust(d,method="ward.D")
plot(cah)
groups <- cutree(cah,k=6)
plot(som_model,bgcol = rainbow(6)[groups],
     main="Hierarchical Clusters")
add.cluster.boundaries(SOM,clustering=groups)


plot(as.phylo(cah), type = "unrooted", cex = 0.6,
     no.margin = TRUE)


plot(as.phylo(cah), type = "fan")

