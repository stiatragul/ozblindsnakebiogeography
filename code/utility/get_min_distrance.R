library(raster)
library(sp)
library(sf)
library(maps)
library(mapproj)

# Each .shp file was extracted using QGIS 
# CRS == EPSG:4326 - WGS 84 - Geographic
# Australia is from WWF Biomes

# From Natural Earth -- Admin 0-Countries
# PNG (ex Bismarck) is manually cut out in QGIS
# Bismarck + Solomon Island is manually cut out in QGIS
# Lesser Sunda manually cut out from Carlos

reg.aus <- shapefile("data/bioregions_map/map/aus_biomes.shx")
reg.sun <- shapefile("data/bioregions_map/map/LesserSundas.shx")
reg.bis <- shapefile("data/bioregions_map/map/Bismarck.shx")
reg.spg <- shapefile("data/bioregions_map/map/Papua.shx")

shp.list <- list(subset(reg.aus,reg.aus$BIOME==unique(reg.aus$BIOME)[1]), # Tropical Grassland (Savannah)
                 subset(reg.aus,reg.aus$BIOME==unique(reg.aus$BIOME)[9]), # Temperate Forest
                 subset(reg.aus,reg.aus$BIOME==unique(reg.aus$BIOME)[2]), # Arid
                 subset(reg.aus,reg.aus$BIOME==unique(reg.aus$BIOME)[5]), # Mediterranean
                 subset(reg.aus,reg.aus$BIOME==unique(reg.aus$BIOME)[3]), # Tropical Forest
                 subset(reg.aus,reg.aus$BIOME==unique(reg.aus$BIOME)[8]), # Temperate Grasslands
                 reg.sun, # Lesser Sunda
                 reg.spg, # Papua
                 reg.bis) # Bismarck Archipelago

names(shp.list) <- c("TropicalGrasslands",
                     "TemperateForest",
                     "Arid",
                     "Mediterranean",
                     "TropicalForest",
                     "TemperateGrasslands",
                     "LesserSundaIslands",
                     "Papua",
                     "Bismarck")

# Distance matrix set up and naming
dist.mat <- matrix(NA,nrow=9,ncol=9)

# Check if matrix is the way we like
rownames(dist.mat) <- names(shp.list)
colnames(dist.mat) <- names(shp.list)
dist.mat

shp.list.t <- list()
for (i in 1:length(shp.list)){
  shp.list.t[[i]] <- spTransform(shp.list[[i]], CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"))
}
names(shp.list.t) <- names(shp.list)

# Plot the 9 bioregions in our study
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

# Calculate distance using rgeos::gDistance
cuentas <- 0
for(i in 1:9){
  for (j in 1:9){
    dist.mat[i,j] <- rgeos::gDistance(shp.list.t[[i]],shp.list.t[[j]])
    cuentas <- cuentas+1
    print(paste(cuentas," of 81 done",sep=""))
  }
}

# Rescale so lowest distance is = 1
lowest_distance <- unique(sort(dist.mat))[2]
dist.mat.rescale <- dist.mat/lowest_distance

# replace adjacent 
dist.mat.rescale[which(dist.mat.rescale==0)] <- 1
dist.mat.rescale <- round(dist.mat.rescale,2)
# Make distance to biome itself still 0
diag(dist.mat.rescale) <- 0

# Prepare dataframe
dist.mat.df <- as.data.frame(dist.mat.rescale)

# Rename it with letters
rownames(dist.mat.df) <- LETTERS[1:nrow(dist.mat.df)]
colnames(dist.mat.df) <- LETTERS[1:nrow(dist.mat.df)]

# Save as txt file for biogeoBEARS
write.table(dist.mat.df, file = "data/bears_txt/biome_distance.txt", sep = "\t", quote = F, row.names = F)
