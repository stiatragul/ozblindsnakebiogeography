library(raster)
library(sp)
library(rgeos)
library(maps)
library(mapproj)

#Each realm's shape was obtained by intersecting the realm shapes with the shape containing all varanid ranges.
#That way realms are not artificially close.

afr <- shapefile("afrotropical.shx")
aus <- shapefile("australian.shx")
pap <- shapefile("papuan.shx")
ori <- shapefile("oriental.shx")
pal <- shapefile("palearctic.shx")
sah <- shapefile("saharo_arabian.shx")
sin <- shapefile("sino_japanese.shx")
phi <- shapefile("philippines.shx")
mel <- shapefile("melanesian.shx")

dist.mat <- matrix(NA,nrow=9,ncol=9)
rownames(dist.mat) <- c("afr","aus","pap","ori","pal","sah","sin","phi","mel")
colnames(dist.mat) <- c("afr","aus","pap","ori","pal","sah","sin","phi","mel")
shp.list <- list(afr,aus,pap,ori,pal,sah,sin,phi,mel)
names(shp.list) <- c("afr","aus","pap","ori","pal","sah","sin","phi","mel")

shp.list.t <- list()
for (i in 1:length(shp.list)){
  shp.list.t[[i]] <- spTransform(shp.list[[i]], CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"))
}
names(shp.list.t) <- names(shp.list)

par(mfrow=c(3,3))
for (i in 1:length(shp.list)){
  plot(shp.list[[i]],col="red",main=names(shp.list)[i])
}
dev.off()

par(mfrow=c(3,3))
for (i in 1:length(shp.list.t)){
  plot(shp.list.t[[i]],col="red",main=names(shp.list.t)[i])
}
dev.off()

cuentas <- 0
for(i in 1:9){
  for (j in 1:9){
    dist.mat[i,j] <- gDistance(shp.list.t[[i]],shp.list.t[[j]])
    cuentas <- cuentas+1
    print(paste(cuentas," of 81 done",sep=""))
  }
}

dist.mat

dist.mat.t <- dist.mat/unique(sort(dist.mat))[2]
dist.mat.t[which(dist.mat.t==0)] <- 1
dist.mat.t <- round(dist.mat.t,2)
diag(dist.mat.t) <- 0

dist.mat.df <- as.data.frame(dist.mat.t)
dist.mat.df <- dist.mat.df[c(4,7,2,3,5,6,1,8,9),c(4,7,2,3,5,6,1,8,9)]

dist.mat.df[1,] <- apply(dist.mat.df[1:2,],2,FUN=min)
dist.mat.df <- dist.mat.df[-2,]

dist.mat.df[4,] <- apply(dist.mat.df[4:5,],2,FUN=min)
dist.mat.df <- dist.mat.df[-5,]

dist.mat.df[,1] <- apply(dist.mat.df[,1:2],1,FUN=min)
dist.mat.df <- dist.mat.df[,-2]

dist.mat.df[,4] <- apply(dist.mat.df[,4:5],1,FUN=min)
dist.mat.df <- dist.mat.df[,-5]

rownames(dist.mat.df) <- LETTERS[1:nrow(dist.mat.df)]
colnames(dist.mat.df) <- LETTERS[1:nrow(dist.mat.df)]

write.table(dist.mat.df,"realm_dist.txt",sep = "\t",quote = F,row.names = F)
