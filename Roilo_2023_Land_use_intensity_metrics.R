###################################################
#
# Title: Roilo_2023_Land_use_intensity_metrics.R
# Purpose: calculate different land use intensity metrics commonly used in biodiversity models, using different spatial aggregation units (grid types)
# over the Mulde River Basin, Germany.
# Reference: Roilo et al. "A multidimensional approach is needed to better quantify land-use intensity in biodiversity models".
# Author: Stephanie Roilo, Technische Universität Dresden
# Date: last updated on March 13th, 2023.
#
###################################################

# set the language to EN
Sys.setenv(LANGUAGE="en")
# load packages
library(rgdal)
library(raster)
library(sf)
library(dplyr)
library(units)
library(terra)
library(tmap)
library(ggplot2)
library(ggridges)
library(tidyr)
library(moments)
library(mapview)

# function for selecting the polygon centroids, or, if they fall outside of the polygon, a random point within it (source: https://stackoverflow.com/questions/52522872/r-sf-package-centroid-within-polygon ):
st_centroid_within_poly <- function (poly) {
  # check if centroid is in polygon
  ctrd <- st_centroid(poly, of_largest_polygon = TRUE)
  in_poly <- diag(st_within(ctrd, poly, sparse = F))
  # replace geometries that are not within polygon with st_point_on_surface()
  st_geometry(ctrd[!in_poly,]) <- st_geometry(st_point_on_surface(poly[!in_poly,]))
  ctrd
}  

###### GRID PREPARATION --------------------------
# set working directory
setwd("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/")
# load the shapefile of the boundary of the Mulde River Basin
mulde = st_read("Mulde/Mulde_EPSG3035.shp") 
# create a square grid covering the study area
sqgrid = st_make_grid(mulde, cellsize = c(1000,1000), what = "polygons", square=T)  
sqgrid = st_as_sf(sqgrid)
sqgrid = sqgrid[which(st_intersects(sqgrid, mulde, sparse= F) == T),]
# add unique identifier for each square
sqgrid$FID = c(1:nrow(sqgrid))

# create a hexagonal grid covering the same area
hexgrid = st_make_grid(sqgrid, cellsize = c(1000,1000), what = "polygons", square=F) # square=F to make hexagonal grids
hexgrid = st_as_sf(hexgrid)
hexgrid = hexgrid[which(st_intersects(hexgrid, mulde, sparse= F) == T),]
# plot and check results spatially
mapview(hexgrid) + sqgrid
# add unique identifier for each hexagon
hexgrid$FID = c(1:nrow(hexgrid))

# create voronoi grid
# the land cover map by Preidl et al. (2020). Introducing APiC for regionalised land cover mapping on the national scale using Sentinel-2A imagery. Remote Sensing of Environment, 240, 111673 was
# cropped to the study area extent and converted to polygons in QGIS. Now let's load it in R and keep only the polygons within our study region
luc = st_read("Land cover/APiC_Agricultural-Land-Cover-Germany_RSE-2020/u2018_clc2018_v2020_20u1_fgdb.shp")
luc = st_transform(luc, 3035) %>%  st_intersection(mulde)
luc = luc[,c("DN", "geometry")]
luc$area_m2 = set_units(st_area(luc), NULL)
# how many polygons are these? (Square grid has 6158; hexagonal one has 7070)
nrow(luc)  # 265567
mean(luc$area_m2)  # 21890.59 
median(luc$area_m2)  # 400
# eliminate polygons with area smaller than 12.5 ha, to reduce number of polygons; this corresponds to 1/8 of the area of a square cell of 1km2
lucf = luc[luc$area_m2>= 125000,]
mapview(lucf, zcol="DN")
# remove unclassified land-cover polygons
lucf = lucf[lucf$DN>0,]
## now calculate the centroid (or, if falling outside the polygon, a random inner point) of each polygon 
cents = st_centroid_within_poly(lucf)
# compute voronoi polygons for each centroid
vor = st_voronoi(st_combine(st_geometry(cents)))
vor = st_collection_extract(vor) 
vor = vor[unlist(st_intersects(cents, vor))]
vor = st_set_geometry(cents, vor)
mapview(vor, zcol="DN") + cents
## voronoi polygons are now computed for a much larger area (the entire bounding box of the Mulde region); cut them to the Mulde extent
vor = st_intersection(vor, st_buffer(mulde, 100))
vor = vor[,c("DN", "geometry")]
rm(cents, lucf, st_centroid_within_poly)

#save grids to file
st_write(sqgrid, "Mulde/Mulde_1km_grid.shp")
st_write(hexgrid, "Mulde/Mulde_hexagonal_grid.shp")
st_write(vor, "Mulde/Voronoi_PreidlMap.shp")


###### COMPUTATION OF LAND USE METRICS  --------------------------
# now, compute the land use intensity metrics for each of the three grids
sqgrid = st_read("Mulde/Mulde_1km_grid.shp")
hexgrid = st_read("Mulde/Mulde_hexagonal_grid.shp")
vor = st_read("Mulde/Voronoi_PreidlMap.shp")

# load Integrated Administration and Control System (IACS) data for the year 2019
in9 = st_read("INVEKOS_data/edited/InVeKoS_2019_Sc.shp")
# load file with crop-codes and corresponding crop names
codes <- read.table("INVEKOS_data/NC_LIST_SN_2020_sc.txt", header = T, sep="\t", stringsAsFactors = FALSE)
codes[1,1] <- "050"
# match crop name and group to the IACS shapefile
in9$crop <- codes$NC_BEZ[match(in9$HKCODE, codes$NC)]
in9$BNK <- codes$BNK[match(in9$HKCODE, codes$NC)]

# create a function that calculates proportion of area within a grid that is covered by a second layer, 
# and that adds it as an attribute column of name layer_name to the grid
CalcProp = function(grid, layer, layer_name) {
  grid$FID = c(1:nrow(grid))
  int = st_intersection(grid, layer)
  int$area = set_units(st_area(int), NULL)
  arabdf = aggregate(area ~ FID, FUN=sum, data=int)   
  grid[, layer_name] <- arabdf$area[match(grid$FID, arabdf$FID)]/set_units(st_area(grid), NULL)
  grid[is.na(grid[, layer_name]), layer_name] <- 0
  return(grid)
}

# 1. calculate proportion of utilized agricultural area (UAA)
sqgrid = CalcProp(sqgrid, in9, "UAA")
hexgrid = CalcProp(hexgrid, in9, "UAA")
vor = CalcProp(vor, in9, "UAA")

# 2. calculate % arable land 
arab = in9[in9$BNK=="AL", ]  # filter the IACS data to retain only field parcels that are arable land
sqgrid = CalcProp(sqgrid, arab, "Arable")
hexgrid = CalcProp(hexgrid, arab, "Arable")
vor = CalcProp(vor, arab, "Arable")

# 3. calculate proportion of arable land on total UAA
sqgrid$Arab_UAA = ifelse(sqgrid$UAA==0, 0, sqgrid$Arable / sqgrid$UAA)
hexgrid$Arab_UAA = ifelse(hexgrid$UAA==0, 0, hexgrid$Arable / hexgrid$UAA)
vor$Arab_UAA = ifelse(vor$UAA==0, 0, vor$Arable / vor$UAA)

# 4. calculate mean field size
# calculate area of each field parcel and then convert fields to centroids
in9$area = st_area(in9)
fieldc = st_centroid(in9)
# in the square grid
test = st_join(fieldc, sqgrid[,c("FID", "geometry"),])
f_df = aggregate(area ~ FID, FUN=mean, data=test)   
sqgrid[, "FieldSize"] <- f_df$area[match(sqgrid$FID, f_df$FID)]/10000 #convert the area in hectares
sqgrid[is.na(sqgrid[, "FieldSize"]), "FieldSize"] <- 0
# in the hexagonal grid
test = st_join(fieldc, hexgrid[,c("FID", "geometry"),])
f_df = aggregate(area ~ FID, FUN=mean, data=test)   
hexgrid[, "FieldSize"] <- f_df$area[match(hexgrid$FID, f_df$FID)]/10000 #convert the area in hectares
hexgrid[is.na(hexgrid[, "FieldSize"]), "FieldSize"] <- 0
# in the voronoi grid
test = st_join(fieldc, vor[,c("FID", "geometry"),])
f_df = aggregate(area ~ FID, FUN=mean, data=test)   
vor[, "FieldSize"] <- f_df$area[match(vor$FID, f_df$FID)]/10000 #convert the area in hectares
vor[is.na(vor[, "FieldSize"]), "FieldSize"] <- 0

# 5. compute mean nitrogen input index, based on the land-use intensity map by Temme & Verburg (2011). Mapping and modelling of changes in agricultural intensity in Europe. Agriculture, Ecosystems & Environment, 140(1-2), 46-56.
LUI = rast("C:/Users/sroilo/Desktop/BESTMAP documents/Data/Europe/Land Use Intensity map 2000/inthier2000_2.asc")
crs(LUI) <- "+proj=aea +lat_1=32.5 +lat_2=54.5 +lat_0=51.40 +lon_0=22.65 +x_0=0 +y_0=0 +ellps=WGS72 +units=m +no_defs"
LUI = project(LUI, "epsg:3035", method="near")
# compute a general LUI indicator encompassing both cropland and grassland, based on the kg of N input per ha
lui_ar = classify(LUI, rcl= matrix(c(0,1,2,3,4,5,0,1,2,3,1,2), ncol=2), othersNA=F)  # NEW CLASSES: 1 = < 50 kg N/ha, class 2 = 50-150 kg N/ha, class 3 = > 150 kg N/ha
sqgrid$LUI = terra::extract(x=lui_ar, y=vect(sqgrid), na.rm=T, method="simple", weights=T, fun="mean", touches=T)[,2]
hexgrid$LUI = terra::extract(x=lui_ar, y=vect(hexgrid), na.rm=T, method="simple", weights=T, fun="mean", touches=T)[,2]
vor$LUI = terra::extract(x=lui_ar, y=vect(vor), na.rm=T, method="simple", weights=T, fun="mean", touches=T)[,2]
mapview(vor, zcol="LUI")

# 6. calculate proportion of intensively managed land, i.e. UAA without Agri-Environmental Schemes (AES) or Ecological Focus Areas (EFA)
intens = st_make_valid(in9[is.na(in9$AUM_INFO1), ])
vec_int = vect(intens)
# now remove also linear EFA areas
efa9 = vect("INVEKOS_data/edited/EFA_2019.shp")
diff = vec_int - efa9   
diff2 = st_as_sf(diff)  #reconvert in sf object and plot to check everything went well
mapview(diff2)
# remove all field-wide EFA
diff2 = diff2[is.na(diff2$EFA_TYP),]
sqgrid = CalcProp(sqgrid, diff2, "Intens_farm")
hexgrid = CalcProp(hexgrid, diff2, "Intens_farm")
vor = CalcProp(vor, diff2, "Intens_farm")

# 7. calculate proportion of conventionally managed land, i.e. NOT organic
NOorganic = in9[in9$MERKM_OKO=="F", ]
sqgrid = CalcProp(sqgrid, NOorganic, "Conv_farm")
hexgrid = CalcProp(hexgrid, NOorganic, "Conv_farm")
vor = CalcProp(vor, NOorganic, "Conv_farm")

# 8. calculate agricultural land use diversity as a Shannon diversity index of different crops and grassland types
# in the square grid
int = st_intersection(sqgrid, in9)
int$area = set_units(st_area(int), NULL)
for ( i in unique(int$FID)) {
  subset = int[int$FID == i, ]  # subset of intersected fields with one grid cell
  cropdf = aggregate(area ~ HKCODE, data= subset, FUN = sum)
  cropdf$prop = cropdf$area/1000000
  sqgrid$ALU_div[sqgrid$FID==i] = -sum(sapply(X = cropdf$prop, FUN = function(x) {x * log(x)}, simplify = "vector"))
}
sqgrid$ALU_div[is.na(sqgrid$ALU_div)] <-0
mapview(sqgrid, zcol="ALU_div")
# in the hexagonal grid
int = st_intersection(hexgrid, in9)
int$area = set_units(st_area(int), NULL)
hex_area = set_units(st_area(hexgrid[1,]), NULL)
for ( i in unique(int$FID)) {
  subset = int[int$FID == i, ]  # subset of intersected fields with one grid cell
  cropdf = aggregate(area ~ HKCODE, data= subset, FUN = sum)
  cropdf$prop = cropdf$area/hex_area
  hexgrid$ALU_div[hexgrid$FID==i] = -sum(sapply(X = cropdf$prop, FUN = function(x) {x * log(x)}, simplify = "vector"))
}
hexgrid$ALU_div[is.na(hexgrid$ALU_div)] <-0
mapview(hexgrid, zcol="ALU_div")
# in the voronoi grid
int = st_intersection(vor, in9)
int$area = set_units(st_area(int), NULL)
for ( i in unique(int$FID)) {
  subset = int[int$FID == i, ]  # subset of intersected fields with one grid cell
  cropdf = aggregate(area ~ HKCODE, data= subset, FUN = sum)
  cropdf$prop = cropdf$area/set_units(st_area(vor[vor$FID==i,]), NULL)
  vor$ALU_div[vor$FID==i] = -sum(sapply(X = cropdf$prop, FUN = function(x) {x * log(x)}, simplify = "vector"))
}
vor$ALU_div[is.na(vor$ALU_div)] <-0
mapview(vor, zcol="ALU_div")

# 9. calculate land-use/land cover heterogeneity from the land cover map by Preidl et al. (2020)
# load vectorised Preidl map
lcm = st_read("Land cover/APiC_Agricultural-Land-Cover-Germany_RSE-2020/u2018_clc2018_v2020_20u1_fgdb.shp") %>% st_transform(3035)
int = st_intersection(sqgrid, lcm)
int$area = set_units(st_area(int), NULL)
for ( i in unique(int$FID)) {
  subset = int[int$FID == i, ]  
  cropdf = aggregate(area ~ DN, data= subset, FUN = sum)
  cropdf$prop = cropdf$area/1000000
  sqgrid$LULC_div[sqgrid$FID==i] = -sum(sapply(X = cropdf$prop, FUN = function(x) {x * log(x)}, simplify = "vector"))
}
sqgrid$LULC_div[is.na(sqgrid$LULC_div)] <-0
# in the hexagonal grid
int = st_intersection(hexgrid, lcm)
int$area = set_units(st_area(int), NULL)
hex_area = set_units(st_area(hexgrid[1,]), NULL)
for ( i in unique(int$FID)) {
  subset = int[int$FID == i, ]  
  cropdf = aggregate(area ~ DN, data= subset, FUN = sum)
  cropdf$prop = cropdf$area/hex_area
  hexgrid$LULC_div[hexgrid$FID==i] = -sum(sapply(X = cropdf$prop, FUN = function(x) {x * log(x)}, simplify = "vector"))
}
hexgrid$LULC_div[is.na(hexgrid$LULC_div)] <-0
# in the voronoi grid
vor = vor[,-1]  # remove old DN information
int = st_intersection(vor, lcm)
int$area = set_units(st_area(int), NULL)
for ( i in unique(int$FID)) {
  subset = int[int$FID == i, ]  
  cropdf = aggregate(area ~ DN, data= subset, FUN = sum)
  cropdf$prop = cropdf$area/set_units(st_area(vor[vor$FID==i,]), NULL)
  vor$LULC_div[vor$FID==i] = -sum(sapply(X = cropdf$prop, FUN = function(x) {x * log(x)}, simplify = "vector"))
}
vor$LULC_div[is.na(vor$LULC_div)] <-0

# 10. calculate land-cover heterogeneity index based on the S2GLC map by Malinowski et al. "Automated production of a land cover/use map of Europe based on Sentinel-2 imagery." Remote Sensing 12.21 (2020): 3523.
map = st_read("Land cover/S2GLC_2017/S2GLC_polygonized.shp")
int = st_intersection(sqgrid, map)
int$area = set_units(st_area(int), NULL)
for ( i in unique(int$FID)) {
  subset = int[int$FID == i, ]  
  cropdf = aggregate(area ~ DN, data= subset, FUN = sum)
  cropdf$prop = cropdf$area/1000000
  sqgrid$LC_div[sqgrid$FID==i] = -sum(sapply(X = cropdf$prop, FUN = function(x) {x * log(x)}, simplify = "vector"))
}
sqgrid$LC_div[is.na(sqgrid$LC_div)] <-0
# in the hexagonal grid
int = st_intersection(hexgrid, map)
int$area = set_units(st_area(int), NULL)
hex_area = set_units(st_area(hexgrid[1,]), NULL)
for ( i in unique(int$FID)) {
  subset = int[int$FID == i, ]  
  cropdf = aggregate(area ~ DN, data= subset, FUN = sum)
  cropdf$prop = cropdf$area/hex_area
  hexgrid$LC_div[hexgrid$FID==i] = -sum(sapply(X = cropdf$prop, FUN = function(x) {x * log(x)}, simplify = "vector"))
}
hexgrid$LC_div[is.na(hexgrid$LC_div)] <-0
# in the voronoi grid
int = st_intersection(vor, lcm)
int$area = set_units(st_area(int), NULL)
for ( i in unique(int$FID)) {
  subset = int[int$FID == i, ]  
  cropdf = aggregate(area ~ DN, data= subset, FUN = sum)
  cropdf$prop = cropdf$area/set_units(st_area(vor[vor$FID==i,]), NULL)
  vor$LC_div[vor$FID==i] = -sum(sapply(X = cropdf$prop, FUN = function(x) {x * log(x)}, simplify = "vector"))
}
vor$LC_div[is.na(vor$LC_div)] <-0

# rescale the variables to that they all range between 0 and 1, and calculate the mathematical opposite of the diversity/heterogeneity variables
# drop the geometry to ease computations
squ = st_drop_geometry(sqgrid); hex = st_drop_geometry(hexgrid); voro = st_drop_geometry(vor)
for ( i in c(2:ncol(squ))) {
  max_n = max(squ[,i])
  squ[,i] <- squ[,i]/max_n
  max_n = max(hex[,i])
  hex[,i] <- hex[,i]/max_n
  max_n = max(voro[,i])
  voro[,i] <- voro[,i]/max_n
}
# now, compute the opposite of the diversity/heterogeneity variables, to get a homogeneity index
squ$ALU_div = 1 - squ$ALU_div; hex$ALU_div = 1 - hex$ALU_div; voro$ALU_div = 1 - voro$ALU_div;
squ$LCLU_div = 1 - squ$LCLU_div; hex$LCLU_div = 1 - hex$LCLU_div; voro$LCLU_div = 1 - voro$LCLU_div;
squ$LC_div = 1 - squ$LC_div; hex$LC_div = 1 - hex$LC_div; voro$LC_div = 1 - voro$LC_div;
# update names of metrics
names(squ) = c("FID", "UAA","Arable","Arable_UAA", "Field_size", "N_input","Intens_farm", "Conv_farm" , "ALU_homog" , "LULC_homog", "LC_homog")
names(hex) = c("FID", "UAA","Arable","Arable_UAA", "Field_size", "N_input","Intens_farm", "Conv_farm" , "ALU_homog" , "LULC_homog", "LC_homog")
names(voro) = c("FID", "UAA","Arable","Arable_UAA", "Field_size", "N_input","Intens_farm", "Conv_farm" , "ALU_homog" , "LULC_homog", "LC_homog")

# rejoin with the geometries and save to file
sqgrid = dplyr::bind_cols(sqgrid[,"geometry"], squ)
hexgrid = dplyr::bind_cols(hexgrid[,"geometry"], hex)
vor = dplyr::bind_cols(vor[,"geometry"], voro)
st_write(sqgrid, "Land use intensity/metrics/Square_metric_20220908.shp")
st_write(hexgrid, "Land use intensity/metrics/Hexagonal_metric_20220908.shp")
st_write(vor, "Land use intensity/metrics/Voronoi_metric_20220908.shp")


### MEAN LUI AND SD ACROSS METRICS  -----------------
# now, calculate mean and SD of metrics within the same grid
# use the dataframes without geometry
squ$mean = rowMeans(squ[,c(2:11)])
squ$SD = apply(squ[,c(2:11)],1, sd, na.rm = TRUE)
# join with the geometry and filter out cells with no UAA
squaa = bind_cols(sqgrid, squ[,12:15])
squaa = squaa[squaa$UAA>0,]
# plot and save maps to file
m1 = tm_shape(squaa) + tm_polygons("mean", palette="Reds", title = "Mean", border.alpha = 0.1) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
tmap_save(m1, width= 10, height=7, filename="Land use intensity/images/maps/Square_mean.jpeg")
m2 = tm_shape(squaa) + tm_polygons("SD", palette="Blues", title = "Standard deviation", border.alpha = 0.1) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
tmap_save(m2, width= 10, height=7, filename="Land use intensity/images/maps/Square_SD.jpeg")

## do the same with the hexagonal and voronoi grid
hex$mean = rowMeans(hex[,c(2:11)])
hex$SD = apply(hex[,c(2:11)],1, sd, na.rm = TRUE)
# merge with geometry and get rid of cells with no UAA
hexuaa = bind_cols(hexgrid, hex[,12:15])
hexuaa = hexuaa[hexuaa$UAA>0,]
# plot and save maps to file
m1 = tm_shape(hexuaa) + tm_polygons("mean", palette="Reds", title = "Mean", border.alpha = 0.1) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
tmap_save(m1, width= 10, height=7, filename="Land use intensity/images/maps/Hexagon_mean.jpeg")
m2 = tm_shape(hexuaa) + tm_polygons("SD", palette="Blues", title = "Standard deviation", border.alpha = 0.1) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
tmap_save(m2, width= 10, height=7, filename="Land use intensity/images/maps/Hexagon_SD.jpeg")

# voronoi grid
voro$mean = rowMeans(voro[,c(2:11)])
voro$SD = apply(voro[,c(2:11)],1, sd, na.rm = TRUE)
# mrege with geometry and get rid of cells with no UAA
vorouaa = bind_cols(vor, voro[,12:15])
vorouaa = vorouaa[vorouaa$UAA>0,]
# plot and save maps to file
m1 = tm_shape(vorouaa) + tm_polygons("mean", palette="Reds", title = "Mean", border.alpha = 0.1) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
tmap_save(m1, width= 10, height=7, filename="Land use intensity/images/maps/Voronoi_mean.jpeg")
m2 = tm_shape(vorouaa) + tm_polygons("SD", palette="Blues", title = "Standard deviation", border.alpha = 0.1) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
tmap_save(m2, width= 10, height=7, filename="Land use intensity/images/maps/Voronoi_SD.jpeg")

# save newly prepared grids to file
st_write(squaa, "Land use intensity/metrics/SqUAA.shp")
st_write(hexuaa, "Land use intensity/metrics/HexUAA.shp")
st_write(vorouaa, "Land use intensity/metrics/VoroUAA.shp")

#### MAKE PLOTS of the different grids & metrics  -----------------------------------------
library(tmap)

squ1 = tm_shape(sqgrid) + tm_polygons("Arable", palette="viridis", title = expression("% arable land")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
squ2 = tm_shape(sqgrid) + tm_polygons("UAA", palette="viridis", title = expression("% Utilized Agricultural Area")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
squ3 = tm_shape(sqgrid) + tm_polygons("Arb_UAA", palette="viridis", title = expression("% arable land on total UAA")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
squ4 = tm_shape(sqgrid) + tm_polygons("Fild_sz", palette="viridis", title = expression("Mean field size (ha)")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
squ5 = tm_shape(sqgrid) + tm_polygons("N_input", palette="viridis", title = expression("Nitrogen input (index)")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
squ6 = tm_shape(sqgrid) + tm_polygons("Intns_f", palette="viridis", title = expression("% intensive farming")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
squ7 = tm_shape(sqgrid) + tm_polygons("Cnv_frm", palette="viridis", title = expression("% conventional farming")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
squ8 = tm_shape(sqgrid) + tm_polygons("ALU_hmg", palette="viridis", title = expression("ALU homogeneity")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
squ9 = tm_shape(sqgrid) + tm_polygons("LCLU_hm", palette="viridis", title = expression("LULC homogeneity")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
squ10 = tm_shape(sqgrid) + tm_polygons("LC_homg", palette="viridis", title = expression("LC homogeneity")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
# combine all maps and save to file
all = tmap_arrange(squ1, squ2, squ3, squ4, squ5, squ6, squ7, squ8, squ9, squ10)
tmap_save(all, width= 20, height=20, filename="Land use intensity/Square_grid_202210.jpeg")

## hexagonal grid
hex1 = tm_shape(hexgrid) + tm_polygons("Arable", palette="viridis", title = expression("% arable land")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
hex2 = tm_shape(hexgrid) + tm_polygons("UAA", palette="viridis", title = expression("% Utilized Agricultural Area")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
hex3 = tm_shape(hexgrid) + tm_polygons("Arb_UAA", palette="viridis",  title = expression("% arable land on total UAA")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
hex4 = tm_shape(hexgrid) + tm_polygons("Fild_sz", palette="viridis", title = expression("Mean field size (ha)")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
hex5 = tm_shape(hexgrid) + tm_polygons("N_input", palette="viridis", title = expression("Nitrogen input (index)")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
hex6 = tm_shape(hexgrid) + tm_polygons("Intns_f", palette="viridis", title = expression("% intensive farming")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
hex7 = tm_shape(hexgrid) + tm_polygons("Cnv_frm", palette="viridis", title = expression("% conventional farming")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
hex8 = tm_shape(hexgrid) + tm_polygons("ALU_hmg", palette="viridis", title = expression("ALU homogeneity")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
hex9 = tm_shape(hexgrid) + tm_polygons("LCLU_hm", palette="viridis", title = expression("LULC homogeneity")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
hex10 = tm_shape(hexgrid) + tm_polygons("LC_homg", palette="viridis", title = expression("LC homogeneityy")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
# combine all maps and save to file
all = tmap_arrange(hex1, hex2, hex3, hex4, hex5, hex6, hex7, hex8, hex9, hex10)
tmap_save(all, width= 20, height=20, filename="Land use intensity/Hexagonal_grid_202210.jpeg")

## voronoi grid
vor1 = tm_shape(vor) + tm_polygons("Arable", palette="viridis", title = expression("% arable land")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
vor2 = tm_shape(vor) + tm_polygons("UAA", palette="viridis", title = expression("% Utilized Agricultural Area")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
vor3 = tm_shape(vor) + tm_polygons("Arb_UAA", palette="viridis", title = expression("% arable land on total UAA")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
vor4 = tm_shape(vor) + tm_polygons("Fild_sz", palette="viridis", title = expression("Mean field size (ha)")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
vor5 = tm_shape(vor) + tm_polygons("N_input", palette="viridis", title = expression("Nitrogen input (index)")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
vor6 = tm_shape(vor) + tm_polygons("Intns_f", palette="viridis", title = expression("% intensive farming")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
vor7 = tm_shape(vor) + tm_polygons("Cnv_frm", palette="viridis", title = expression("% conventional farming")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
vor8 = tm_shape(vor) + tm_polygons("ALU_hmg", palette="viridis", title = expression("ALU homogeneity")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
vor9 = tm_shape(vor) + tm_polygons("LCLU_hm", palette="viridis", title = expression("LULC homogeneity")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
vor10 = tm_shape(vor) + tm_polygons("LC_homg", palette="viridis", title = expression("LC homogeneity")) + 
  tm_compass(type = "8star", position = c("center", "top")) +
  tm_scale_bar(breaks = c(0,15,30), text.size = 1, position = c("right", "bottom")) + 
  tm_layout(legend.position = c("right", "top"), legend.text.size = 1, legend.title.size=1.5)
# combine all maps and save to file
all = tmap_arrange(vor1, vor2, vor3, vor4, vor5, vor6, vor7, vor8, vor9, vor10)
tmap_save(all, width= 20, height=20, filename="Land use intensity/Voronoi_grid_202210.jpeg")


rm(list=ls())
setwd("C:/Users/sroilo/Documents")
