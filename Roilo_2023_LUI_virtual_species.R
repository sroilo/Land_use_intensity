###################################################
#
# Title: Roilo_2023_LUI_virtual_species.R
# Purpose: create three virtual species with known land-use intensity (LUI)-species relationships to test how using different subsets of LUI metrics 
#          affects the outcomes of biodiversity models (we here used Generalised Additive Models to model the virtual species occurrence).
# Reference: Roilo et al. "Quantifying agricultural land-use intensity for spatial biodiversity modelling: effects of different metrics and spatial aggregation methods".
# Author: Stephanie Roilo, Technische Universität Dresden
# Date: last updated on September 15th, 2023.
#
###################################################

# set the language to EN
Sys.setenv(LANGUAGE="en")
# load packages
library(virtualspecies)
library(sp)
library(sf)
library(raster)
library(dplyr)
library(terra)
library(psych)   # to plot nice pairs.plots
library(MuMIn)   # for automatic model selection
library(DHARMa)  # to simulate residuals and check model fit
library(ncf)  # to check for spatial autocorrelation in the model residuals
library(mgcv)   # for fitting GAM models
library(ggplot2)   # for plotting 
library(ggridges)  # for plotting ridgeplots
library(tidyr)    # for handling data

# check the virtualspecies tutorial here: http://borisleroy.com/files/virtualspecies-tutorial.html 
# load the previously computed LUI metrics in the square grid, and prepare the raster layers to be used to set up the virtual species
sqgrid = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/metrics/Square_metric_20220908.shp")
r = rast(crs="epsg:3035",extent=ext(sqgrid), resolution=c(1000,1000))

# We design three virtual species: a grassland species, an arable land species, and a wetland species
### VIRTUAL GRASSLAND SPECIES --------------------------------------------------------------------
# convert the LUI metrics from the square grid to raster
Intens_f = rasterize(vect(sqgrid), r, field="Intns_f")
LULC_homog = rasterize(vect(sqgrid), r, field="LCLU_hm")
Arable_UAA = rasterize(vect(sqgrid), r, field="Arb_UAA")

# combine everything in a rasterstack (of type RasterStack, raster package)
env = c(Arable_UAA, Intens_f, LULC_homog)  
names(env) = c("Arable_UAA", "Intens_f", "LULC_homog")
env =  env %>% raster::stack()
plot(env)
# define the response curves of the virtual species to the LUI metrics 
vgrass_parameters = formatFunctions(Arable_UAA = c(fun = 'logisticFun', beta = 0.6, alpha = 0.1),
                                    LULC_homog = c(fun = 'linearFun', a=-0.95, b=1),
                                    Intens_f = c(fun = 'linearFun', a = -0.8, b = 1)) 
# create the virtual species
vgrass <- generateSpFromFun(raster.stack = env, 
                            parameters = vgrass_parameters, 
                            species.type="multiplicative",
                            rescale = TRUE, rescale.each.response = F, plot = TRUE)
# plot species response plots
plotResponse(vgrass)
# save suitability raster to file
writeRaster(vgrass$suitab.raster, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Grassland_suitability_20230305.tif")

#convert the suitability raster into a binary presence-absence raster using a logistic curve
set.seed(13)
vgrass_PA <- convertToPA(vgrass, PA.method = "probability", prob.method = "logistic", beta = 0.3, alpha = -0.05, plot=T)
# save rasters to file
writeRaster(vgrass_PA$probability.of.occurrence, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Grassland_ProbOcc_20230305.tif")
writeRaster(vgrass_PA$pa.raster, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Grassland_PA_20230305.tif")

# sample virtual species occurrence at 500 random points, only in grid cells which have UAA>0
uaa = as(sqgrid[sqgrid$UAA>0,], "Spatial")
set.seed(13)
vgrass_points <- sampleOccurrences(vgrass_PA, n = 500, type = "presence-absence", sampling.area = uaa, 
                                 detection.probability = 1, extract.probability=TRUE, plot=T)
# convert to sf object
vgrass_points <- vgrass_points$sample.points %>% st_as_sf(coords = c("x", "y"), crs=3035)
sum(vgrass_points$Observed)  # 99 presence points out of 500
# extract the LUI metrics at the 500 randomly sampled points
vgrass_points = st_join(vgrass_points, sqgrid, left=T)
names(vgrass_points)[4:14] <- c("FID_SQ", "UAA_SQ","Arable_SQ", "Arable_UAA_SQ", "Field_size_SQ", "N_input_SQ", "Intens_f_SQ", "Conv_farm_SQ", 
                              "ALU_homog_SQ", "LULC_homog_SQ", "LC_homog_SQ")
names(vgrass_points)[2] <- "presence"
# check if points fall in NA areas of the grid
sum(is.na(vgrass_points$FID_SQ))  # 0
# add spatial coordinates to the dataframe
vgrass_points = bind_cols(st_coordinates(vgrass_points), vgrass_points)
# save to file
st_write(vgrass_points, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Grassland_sampled_points_202030305.gpkg")

### VIRTUAL ARABLE SPECIES -------------------------------------------------
# convert the LUI metrics from the square grid to raster 
Arable = rasterize(vect(sqgrid), r, field="Arable")
Field_size = rasterize(vect(sqgrid), r, field="Fild_sz") # see Gayer et al 2019 "but not skylarks" publication
N_input = rasterize(vect(sqgrid), r, field="N_input")   

# combine everything in a rasterstack (of type RasterStack, raster package)
env = c(Arable, N_input, Field_size)  
names(env) = c("Arable", "N_input", "Field_size")
env =  env %>% raster::stack()
plot(env)
# define the response curves of the virtual species to the LUI metrics 
varab_parameters = formatFunctions(Arable = c(fun = 'quadraticFun', a = -1, b = 2, c = 0),
                                   N_input = c(fun = 'linearFun', a=-0.8, b=1),
                                   Field_size = c(fun = 'logisticFun', beta = 0.1, alpha = -0.05))
# create the virtual species
varab <- generateSpFromFun(raster.stack = env, 
                            parameters = varab_parameters, 
                            species.type="multiplicative",
                            rescale = TRUE, rescale.each.response = F, plot = TRUE)
# plot species response plots
plotResponse(varab)
# save suitability raster to file
writeRaster(varab$suitab.raster, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Arable_suitability_20230908.tif")

#convert the suitability raster into a binary presence-absence raster using a logistic curve
set.seed(13)
varab_PA <- convertToPA(varab, PA.method = "probability", prob.method = "logistic", beta = 0.1, alpha = -0.05, plot=T)
# save rasters to file
writeRaster(varab_PA$probability.of.occurrence, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Arable_ProbOcc_20230908.tif")
writeRaster(varab_PA$pa.raster, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Arable_PA_20230908.tif")

# sample virtual species occurrence at 500 random points, only in grid cells which have UAA>0
uaa = as(sqgrid[sqgrid$UAA>0,], "Spatial")
set.seed(50)
varab_points <- sampleOccurrences(varab_PA, n = 500, type = "presence-absence", sampling.area = uaa, 
                                   detection.probability = 1, extract.probability=TRUE, plot=T)
# convert to sf object
varab_points <- varab_points$sample.points %>% st_as_sf(coords = c("x", "y"), crs=3035)
sum(varab_points$Observed)  # 215 presence points out of 500

# extract the LUI metrics at the 500 randomly sampled points
varab_points = st_join(varab_points, sqgrid, left=T)
names(varab_points)[4:14] <- c("FID_SQ", "UAA_SQ","Arable_SQ", "Arable_UAA_SQ", "Field_size_SQ", "N_input_SQ", "Intens_f_SQ", "Conv_farm_SQ", 
                                "ALU_homog_SQ", "LULC_homog_SQ", "LC_homog_SQ")
names(varab_points)[2] <- "presence"
# check if points fall in NA areas of the grid
sum(is.na(varab_points$FID_SQ))  # 0
# add spatial coordinates to the dataframe
varab_points = bind_cols(st_coordinates(varab_points), varab_points)
# save to file
st_write(varab_points, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Arable_sampled_points_202030908.gpkg")

### VIRTUAL WETLAND SPECIES ------------
# convert the LUI metrics from the square grid to raster 
Arable_UAA = rasterize(vect(sqgrid), r, field="Arb_UAA")  #likes (wet) grasslands along rivers, etc.
Conv_farm = rasterize(vect(sqgrid), r, field="Cnv_frm")   # likes organic farming - less pollution in the water and herbicides in the grass
LC_homog = rasterize(vect(sqgrid), r, field="LC_homg")    # like diversity of land covers (must have water and also tussocky vegetation and grassland)

# combine everything in a rasterstack (of type RasterStack, raster package)
env = c(Arable_UAA, Conv_farm, LC_homog)  
names(env) = c("Arable_UAA", "Conv_farm", "LC_homog")
env =  env %>% raster::stack()
plot(env)
# define the response curves of the virtual species to the LUI metrics 
vwet_parameters = formatFunctions(Arable_UAA = c(fun = 'linearFun', a=-1, b=1),  
                                  Conv_farm = c(fun = 'quadraticFun', a = -1, b = 0, c = 1),
                                  LC_homog = c(fun = 'logisticFun', beta = 0.6, alpha = 0.1) )
# create the virtual species
vwet <- generateSpFromFun(raster.stack = env, 
                          parameters = vwet_parameters, 
                          species.type="multiplicative",
                          rescale = TRUE, rescale.each.response = F, plot = TRUE)
# plot species response plots
plotResponse(vwet)
# save suitability raster to file
writeRaster(vwet$suitab.raster, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Wetland_suitability_20230908.tif")

#convert the suitability raster into a binary presence-absence raster using a logistic curve
set.seed(13)
vwet_PA <- convertToPA(vwet, PA.method = "probability", prob.method = "logistic", beta = 0.2, alpha = -0.05, plot=T)
# save rasters to file
writeRaster(vwet_PA$probability.of.occurrence, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Wetland_ProbOcc_20230908.tif")
writeRaster(vwet_PA$pa.raster, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Wetland_PA_20230908.tif")

# sample virtual species occurrence at 500 random points, only in grid cells which have UAA>0
uaa = as(sqgrid[sqgrid$UAA>0,], "Spatial")
set.seed(100)
vwet_points <- sampleOccurrences(vwet_PA, n = 500, type = "presence-absence", sampling.area = uaa, 
                                   detection.probability = 1, extract.probability=TRUE, plot=T)
# convert to sf object
vwet_points <- vwet_points$sample.points %>% st_as_sf(coords = c("x", "y"), crs=3035)
sum(vwet_points$Observed)  # 204 presence points out of 500

# extract the LUI metrics at the 500 randomly sampled points
vwet_points = st_join(vwet_points, sqgrid, left=T)
names(vwet_points)[4:14] <- c("FID_SQ", "UAA_SQ","Arable_SQ", "Arable_UAA_SQ", "Field_size_SQ", "N_input_SQ", "Intens_f_SQ", "Conv_farm_SQ", 
                                "ALU_homog_SQ", "LULC_homog_SQ", "LC_homog_SQ")
names(vwet_points)[2] <- "presence"
# check if points fall in NA areas of the grid
sum(is.na(vwet_points$FID_SQ))  # 0
# add spatial coordinates to the dataframe
vwet_points = bind_cols(st_coordinates(vwet_points), vwet_points)
# save to file
st_write(vwet_points, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Wetland_sampled_points_202030908.gpkg")


### MODELLING the grassland species -----------------------
# load sampled presence/absence points for the species
vgrass_sq = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Grassland_sampled_points_202030305.gpkg") %>% st_drop_geometry()
sum(vgrass_sq$presence) # 500 points, 99 presences
# check correlation among variables
plotd = vgrass_sq[order(vgrass_sq$presence, decreasing=F),]
names(plotd)[7:16] = c( "UAA","Arable","Arable_UAA", "Field_size", "N_input","Intens_f","Conv_farm","ALU_homog","LULC_homog", "LC_homog")
pairs.panels(plotd[, c("presence", "UAA", "Arable", "Arable_UAA", "Conv_farm", "Intens_f", "N_input", 
                       "ALU_homog", "LULC_homog", "LC_homog", "Field_size")], bg=c("blue","yellow")[as.factor(plotd$presence)], pch=21,
             method = "spearman", hist.col = "darkorange", density = TRUE, ellipses = FALSE, main= "Spearman correlation coefficients - sampled points for the grassland species in the square grid" )
# prepare data for GAM modelling
dat = vgrass_sq[,c(1,2,4,7:16)]
# create a correlation matrix for the dataset
cmat <- abs(cor(dat[,4:13], method="spearman")) <= 0.7
cmat[!lower.tri(cmat, diag=F)] <- NA 
# rename variables in the correlation matrix to match the names in the full model to be dredged
row.names(cmat) = paste0("s(", row.names(cmat), ", k = 5)")
colnames(cmat) = paste0("s(", colnames(cmat), ", k = 5)")

# fit full model and dredge
fits = gam(presence ~ s(UAA_SQ,k=5) + s(Arable_SQ,k=5) + s(Arable_UAA_SQ,k=5) + s(Field_size_SQ,k=5) + s(N_input_SQ,k=5) + s(Intens_f_SQ,k=5) + s(Conv_farm_SQ,k=5) + 
             s(ALU_homog_SQ,k=5) + s(LULC_homog_SQ,k=5) + s(LC_homog_SQ,k=5), data= dat, family="binomial", na.action="na.fail")
dr1 = dredge(fits, subset=cmat, extra = c("R^2", "adjR^2"))  
write.table(dr1, sep=";", dec=",", row.names=F, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/Virtual_Grassland_dredge_20230909_GAMcor.csv")

# change names of dat columns for the plots
names(dat)[4:13] = c("UAA","Arable","Arable_UAA","Field_size", "N_input","Intens_f","Conv_farm","ALU_homog","LULC_homog","LC_homog") 
# fit overall best model
fits = gam(presence ~ s(Arable_UAA, k=5) + s(Intens_f, k=5) + s(LULC_homog, k=5),
           data= dat, family="binomial", na.action="na.fail")
summary(fits)  # all highly significant
# diagnostic plots
simulateResiduals(fits, plot=T)  # OK
max.dist = max(dist(dat[,c("X","Y")])) * 2/3
scor = ncf::spline.correlog(x = dat$X, y = dat$Y, z = fits$residuals, resamp=100, xmax = max.dist)
plot(scor)   # NO PROBLEM !
# plot partial effect plots
par(mfrow=c(1,4))
plot(fits,residuals=FALSE,all.terms=TRUE,shade=TRUE,shade.col=2, seWithMean=FALSE, ylim=c(-18,13)) 

# best model with only land use and management metrics: 
sbest1 = gam(presence ~ s(Arable_UAA,k=5) +s(UAA, k=5), data= dat, family="binomial", na.action="na.fail")
summary(sbest1)   # both very significant!
par(mfrow=c(1,4))
plot(sbest1,residuals=FALSE,all.terms=TRUE,shade=TRUE,shade.col=2, seWithMean=FALSE, ylim=c(-18,13)) 
simulateResiduals(sbest1, plot=T)  # OK
scor = ncf::spline.correlog(x = dat$X, y = dat$Y, z = sbest1$residuals, resamp=100, xmax = max.dist)
plot(scor)   # NO PROBLEM !

# best model with only landscape structure and management: 
sbest2 = gam(presence ~ s(Intens_f,k=5) + s(N_input,k=5) + s(LULC_homog,k=5) + s(Field_size,k=5), data= dat, family="binomial", na.action="na.fail")
summary(sbest2)   # all significant
par(mfrow=c(1,4))
plot(sbest2,residuals=FALSE,all.terms=TRUE,shade=TRUE,shade.col=2, seWithMean=FALSE, ylim=c(-18,13)) 
simulateResiduals(sbest2, plot=T) # OK
scor = ncf::spline.correlog(x = dat$X, y = dat$Y, z = sbest2$residuals, resamp=100, xmax = max.dist)
plot(scor)   # NO PROBLEM !

# best model with only land use landscape structure: 
sbest3 = gam(presence ~  s(UAA,k=5) + s(Arable_UAA,k=5) + s(LULC_homog,k=5), data= dat, family="binomial", na.action="na.fail")
summary(sbest3)   # all significant
par(mfrow=c(1,4))
plot(sbest3,residuals=FALSE,all.terms=TRUE,shade=TRUE,shade.col=2, seWithMean=FALSE, ylim=c(-18,13)) 
simulateResiduals(sbest3, plot=T) # OK
scor = ncf::spline.correlog(x = dat$X, y = dat$Y, z = sbest3$residuals, resamp=100, xmax = max.dist)
plot(scor)   # NO PROBLEM !

### MODELLING the arable species -----------------------
# load sampled presence/absence points for the species
varab_sq = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Arable_sampled_points_202030908.gpkg") %>% st_drop_geometry()
sum(varab_sq$presence) # 500 points, 215 presences
# check correlation among variables
plotd = varab_sq[order(varab_sq$presence, decreasing=F),]
names(plotd)[7:16] = c( "UAA","Arable","Arable_UAA", "Field_size", "N_input","Intens_f","Conv_farm","ALU_homog","LULC_homog", "LC_homog")
pairs.panels(plotd[, c("presence", "UAA", "Arable", "Arable_UAA", "Conv_farm", "Intens_f", "N_input", 
                       "ALU_homog", "LULC_homog", "LC_homog", "Field_size")], bg=c("blue","yellow")[as.factor(plotd$presence)], pch=21,
             method = "spearman", hist.col = "darkorange", density = TRUE, ellipses = FALSE, main= "Spearman correlation coefficients - sampled points for the arable land species in the square grid" )
# prepare data for GAM modelling
dat = varab_sq[,c(1,2,4,7:16)]
# create a correlation matrix for the dataset
cmat <- abs(cor(dat[,4:13], method="spearman")) <= 0.7
cmat[!lower.tri(cmat, diag=F)] <- NA 
# rename variables in the correlation matrix to match the names in the full model to be dredged
row.names(cmat) = paste0("s(", row.names(cmat), ", k = 5)")
colnames(cmat) = paste0("s(", colnames(cmat), ", k = 5)")

# fit full model and dredge
fits = gam(presence ~ s(UAA_SQ,k=5) + s(Arable_SQ,k=5) + s(Arable_UAA_SQ,k=5) + s(Field_size_SQ,k=5) + s(N_input_SQ,k=5) + s(Intens_f_SQ,k=5) + s(Conv_farm_SQ,k=5) + 
             s(ALU_homog_SQ,k=5) + s(LULC_homog_SQ,k=5) + s(LC_homog_SQ,k=5), data= dat, family="binomial", na.action="na.fail")
dr1 = dredge(fits, subset=cmat, extra = c("R^2", "adjR^2")) 
write.table(dr1, sep=";", dec=",", row.names=F, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/Virtual_Arable_dredge_20230908_GAMcor.csv")

# change names of dat columns for the plots
names(dat)[4:13] = c("UAA","Arable","Arable_UAA","Field_size", "N_input","Intens_f","Conv_farm","ALU_homog","LULC_homog","LC_homog") 
# fit overall best model
fits = gam(presence ~  s(Arable, k=5) + s(N_input, k=5) + s(Field_size, k=5),
           data= dat, family="binomial", na.action="na.fail")
summary(fits)  
# diagnostic plots
simulateResiduals(fits, plot=T)  # OK
max.dist = max(dist(dat[,c("X","Y")])) * 2/3
scor = ncf::spline.correlog(x = dat$X, y = dat$Y, z = fits$residuals, resamp=100, xmax = max.dist)
plot(scor)   # NO PROBLEM !
# plot partial effect plots
par(mfrow=c(1,4))
plot(fits,residuals=FALSE,all.terms=TRUE,shade=TRUE,shade.col=2, seWithMean=FALSE, ylim=c(-5,10)) 

# best model with only land use and management metrics: 
sbest1 = gam(presence ~ s(Arable,k=5) + s(N_input,k=5), data= dat, family="binomial", na.action="na.fail")
summary(sbest1)   
par(mfrow=c(1,4))
plot(sbest1,residuals=FALSE,all.terms=TRUE,shade=TRUE,shade.col=2, seWithMean=FALSE, ylim=c(-5,10)) 
simulateResiduals(sbest1, plot=T)  # OK
scor = ncf::spline.correlog(x = dat$X, y = dat$Y, z = sbest1$residuals, resamp=100, xmax = max.dist)
plot(scor)   # NO PROBLEM !

# best model with only landscape structure and management: 
sbest2 = gam(presence ~ s(Intens_f,k=5) + s(N_input,k=5) + s(LULC_homog,k=5) + s(Field_size, k=5), data= dat, family="binomial", na.action="na.fail")
summary(sbest2)   # all significant
par(mfrow=c(1,4))
plot(sbest2,residuals=FALSE,all.terms=TRUE,shade=TRUE,shade.col=2, seWithMean=FALSE, ylim=c(-5,10)) 
simulateResiduals(sbest2, plot=T) # OK
scor = ncf::spline.correlog(x = dat$X, y = dat$Y, z = sbest2$residuals, resamp=100, xmax = max.dist)
plot(scor)   # NO PROBLEM !

# best model with only land use landscape structure: 
sbest3 = gam(presence ~  s(Arable,k=5) + s(Field_size, k=5), data= dat, family="binomial", na.action="na.fail")
summary(sbest3)   # all significant
par(mfrow=c(1,4))
plot(sbest3,residuals=FALSE,all.terms=TRUE,shade=TRUE,shade.col=2, seWithMean=FALSE, ylim=c(-5,10)) 
# diagnostic plots with DHARMa package
simulateResiduals(sbest3, plot=T) # OK
scor = ncf::spline.correlog(x = dat$X, y = dat$Y, z = sbest3$residuals, resamp=100, xmax = max.dist)
plot(scor)   # NO PROBLEM !

### MODELLING the wetland species -----------------------
# load sampled presence/absence points for the species
vwet_sq = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Wetland_sampled_points_202030908.gpkg") %>% st_drop_geometry()
sum(vwet_sq$presence) # 500 points, 204 presences
# check correlation among variables
plotd = vwet_sq[order(vwet_sq$presence, decreasing=F),]
names(plotd)[7:16] = c( "UAA","Arable","Arable_UAA", "Field_size", "N_input","Intens_f","Conv_farm","ALU_homog","LULC_homog", "LC_homog")
pairs.panels(plotd[, c("presence", "UAA", "Arable", "Arable_UAA", "Conv_farm", "Intens_f", "N_input", 
                       "ALU_homog", "LULC_homog", "LC_homog", "Field_size")], bg=c("blue","yellow")[as.factor(plotd$presence)], pch=21,
             method = "spearman", hist.col = "darkorange", density = TRUE, ellipses = FALSE, main= "Spearman correlation coefficients - sampled points for the wetland species in the square grid" )
# prepare data for GAM modelling
dat = vwet_sq[,c(1,2,4,7:16)]
# create a correlation matrix for the dataset
cmat <- abs(cor(dat[,4:13], method="spearman")) <= 0.7
cmat[!lower.tri(cmat, diag=F)] <- NA 
# rename variables in the correlation matrix to match the names in the full model to be dredged
row.names(cmat) = paste0("s(", row.names(cmat), ", k = 5)")
colnames(cmat) = paste0("s(", colnames(cmat), ", k = 5)")

# fit full model and dredge
fits = gam(presence ~ s(UAA_SQ,k=5) + s(Arable_SQ,k=5) + s(Arable_UAA_SQ,k=5) + s(Field_size_SQ,k=5) + s(N_input_SQ,k=5) + s(Intens_f_SQ,k=5) + s(Conv_farm_SQ,k=5) + 
             s(ALU_homog_SQ,k=5) + s(LULC_homog_SQ,k=5) + s(LC_homog_SQ,k=5), data= dat, family="binomial", na.action="na.fail")
dr1 = dredge(fits, subset=cmat, extra = c("R^2", "adjR^2")) 
write.table(dr1, sep=";", dec=",", row.names=F, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/Virtual_Wetland_dredge_20230908_GAMcor.csv")
# repeat without outliers (with z-score >3) in Field_size
#abs(scale(dat$Field_size))  # data points with Field_size higher than 0.2 are outliers
#dat = dat[-which(dat$Field_size_SQ>0.2),]
#write.table(dr1, sep=";", dec=",", row.names=F, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/Virtual_Wetland_dredge_20230908_noOutliers.csv")

# change names of dat columns for the plots
names(dat)[4:13] = c("UAA","Arable","Arable_UAA","Field_size", "N_input","Intens_f","Conv_farm","ALU_homog","LULC_homog","LC_homog") 
# fit overall best model
fits = gam(presence ~ s(Arable_UAA, k=5) + s(Conv_farm, k=5) + s(LC_homog, k=5) + s(Field_size, k=5),
           data= dat, family="binomial", na.action="na.fail")
summary(fits)  
simulateResiduals(fits, plot=T)  # two outliers detected, check if outlier test is significant
testOutliers(fits) # outliers detected, but test not significant (n.s.)
max.dist = max(dist(dat[,c("X","Y")])) * 2/3
scor = ncf::spline.correlog(x = dat$X, y = dat$Y, z = fits$residuals, resamp=100, xmax = max.dist)
plot(scor)   # NO PROBLEM !
# plot partial effect plots
par(mfrow=c(1,5))
plot(fits,residuals=FALSE,all.terms=TRUE,shade=TRUE,shade.col=2, seWithMean=FALSE, ylim=c(-25,10)) 

# best model with only land use and management metrics: 
sbest1 = gam(presence ~ s(UAA,k=5) + s(Arable_UAA,k=5) + s(N_input,k=5), data= dat, family="binomial", na.action="na.fail")
summary(sbest1)   # Arable and Intens_f NOT significant!
par(mfrow=c(1,5))
plot(sbest1,residuals=FALSE,all.terms=TRUE,shade=TRUE,shade.col=2, seWithMean=FALSE, ylim=c(-25,10)) 
simulateResiduals(sbest1, plot=T)  # OK
scor = ncf::spline.correlog(x = dat$X, y = dat$Y, z = sbest1$residuals, resamp=100, xmax = max.dist)
plot(scor)   # NO PROBLEM !

# best model with only landscape structure and management: 
sbest2 = gam(presence ~ s(Conv_farm,k=5) + s(N_input,k=5) + s(LC_homog,k=5) + s(LULC_homog,k=5) + s(Field_size, k=5), data= dat, family="binomial", na.action="na.fail")
summary(sbest2)   # all significant
par(mfrow=c(1,5))
plot(sbest2,residuals=FALSE,all.terms=TRUE,shade=TRUE,shade.col=2, seWithMean=FALSE, ylim=c(-25,10)) 
simulateResiduals(sbest2, plot=T) # OK
scor = ncf::spline.correlog(x = dat$X, y = dat$Y, z = sbest2$residuals, resamp=100, xmax = max.dist)
plot(scor)   # NO PROBLEM !

# best model with only land use landscape structure: 
sbest3 = gam(presence ~ s(UAA,k=5) + s(Arable_UAA,k=5) + s(LC_homog,k=5) + s(Field_size,k=5), data= dat, family="binomial", na.action="na.fail")
summary(sbest3)   # all significant
par(mfrow=c(1,5))
plot(sbest3,residuals=FALSE,all.terms=TRUE,shade=TRUE,shade.col=2, seWithMean=FALSE, ylim=c(-25,10)) 
simulateResiduals(sbest3, plot=T) # outliers detected
testOutliers(sbest3)  # significant
scor = ncf::spline.correlog(x = dat$X, y = dat$Y, z = sbest3$residuals, resamp=100, xmax = max.dist)
plot(scor)   # NO PROBLEM !


## RIDGELINE PLOT GRASSLAND SPECIES ----------------
# load the sampled points
vspec = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Grassland_sampled_points_202030305.gpkg")
# join with LUI metrics from the hexagonal and voronoi grids (square grid info is already in the shapefile)
hexgrid = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/metrics/Hexagonal_metric_20220908.shp")
vor = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/metrics/Voronoi_metric_20220908.shp")
vspec = st_join(vspec, hexgrid, left=T)
vspec = st_join(vspec, vor, left=T)
names(vspec)[17:38] <- c("FID_HX","UAA_HX","Arable_HX", "Arable_UAA_HX", "Field_size_HX", "N_input_HX", "Intens_f_HX", "Conv_farm_HX", "ALU_homog_HX", "LULC_homog_HX", "LC_homog_HX",
                        "FID_VO","UAA_VO","Arable_VO", "Arable_UAA_VO", "Field_size_VO", "N_input_VO", "Intens_f_VO", "Conv_farm_VO", "ALU_homog_VO", "LULC_homog_VO", "LC_homog_VO")
# save to file
vsave = vspec[,c(1,2,4,7:16,18:27,29:38)]
st_write(vsave, "C:/Users/sroilo/Desktop/BESTMAP documents/Papers/Land use intensity/data_upload/VirtualSp_Grassland_sampled_points_20240120.gpkg")
# filter to retain only the presence points of the virtual species
vspec = vspec[vspec$presence==1,]  # NOTE: rows are more than 99 now, as some fall in between two hexagons, or two voronoi polygons -> some points are joined to multiple geometries
# drop geometry and reorder the variables
vspec = st_drop_geometry(vspec[,c(1,2,4,7:16,18:27,29:38)])
# convert from wide to long format
allpts2 = gather(vspec, key = Metric, Intensity, UAA_SQ:LC_homog_VO)
allpts2$Grid = "Square grid"
allpts2$Grid[grepl(pattern="_HX", x=allpts2$Metric)] <- "Hexagonal grid"
allpts2$Grid[grepl(pattern="_VO", x=allpts2$Metric)] <- "Voronoi grid"
allpts2 <- allpts2 %>%  mutate(Grid=factor(Grid,levels=unique(allpts2$Grid)) )
allpts2$Metric = gsub(pattern="_SQ", replacement="", 
                      gsub(pattern="_HX", replacement="",     
                           gsub(pattern="_VO", replacement="", allpts2$Metric)))
allpts2 <- allpts2 %>%  mutate(Metric=factor(Metric,levels= rev(c("UAA", "Arable", "Arable_UAA", "Conv_farm", "Intens_f", "N_input", "ALU_homog",
                                                                  "LC_homog", "LULC_homog", "Field_size"))) )
# plot ridgeplots for all grid types and all metrics (exported dimension: 700x550)
# ridgeplot tutorials: https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html 
# and https://jtr13.github.io/cc19/ridgeline-plots.html 
ggplot(allpts2, aes(x=Intensity, y=Metric, fill = factor(stat(quantile)))) +    
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,#scale=0.9    
    quantiles = 5, quantile_lines = TRUE, jittered_points = TRUE,aes(point_alpha=0.4, point_shape=20, point_stroke = 0.5, point_color="dark_red")) +   
  scale_fill_viridis_d(name = "Quantiles", alpha=0.5) +labs(x="Intensity",y='') +facet_wrap( ~Grid, labeller = "label_both")
ggsave(filename="C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/images/Fig2.jpg", width= 2400, height= 1600, units="px", plot = last_plot())  # change name of output file according to metric

## RIDGELINE PLOT ARABLE LAND SPECIES ----------------
# load the sampled points
vspec = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Arable_sampled_points_202030908.gpkg")
# join with LUI metrics from the hexagonal and voronoi grids (square grid info is already in the shapefile)
hexgrid = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/metrics/Hexagonal_metric_20220908.shp")
vor = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/metrics/Voronoi_metric_20220908.shp")
vspec = st_join(vspec, hexgrid, left=T)
vspec = st_join(vspec, vor, left=T)
names(vspec)[17:38] <- c("FID_HX","UAA_HX","Arable_HX", "Arable_UAA_HX", "Field_size_HX", "N_input_HX", "Intens_f_HX", "Conv_farm_HX", "ALU_homog_HX", "LULC_homog_HX", "LC_homog_HX",
                         "FID_VO","UAA_VO","Arable_VO", "Arable_UAA_VO", "Field_size_VO", "N_input_VO", "Intens_f_VO", "Conv_farm_VO", "ALU_homog_VO", "LULC_homog_VO", "LC_homog_VO")
# save to file
vsave = vspec[,c(1,2,4,7:16,18:27,29:38)]
st_write(vsave, "C:/Users/sroilo/Desktop/BESTMAP documents/Papers/Land use intensity/data_upload/VirtualSp_Arable_sampled_points_20240120.gpkg")
# filter to retain only the presence points of the virtual species
vspec = vspec[vspec$presence==1,]  # NOTE: rows are more than 215 now, as some fall in between two hexagons, or two voronoi polygons -> some points are joined to multiple geometries
# drop geometry and reorder the variables
vspec = st_drop_geometry(vspec[,c(1,2,4,7:16,18:27,29:38)])
# convert from wide to long format
allpts2 = gather(vspec, key = Metric, Intensity, UAA_SQ:LC_homog_VO)
allpts2$Grid = "Square grid"
allpts2$Grid[grepl(pattern="_HX", x=allpts2$Metric)] <- "Hexagonal grid"
allpts2$Grid[grepl(pattern="_VO", x=allpts2$Metric)] <- "Voronoi grid"
allpts2 <- allpts2 %>%  mutate(Grid=factor(Grid,levels=unique(allpts2$Grid)) )
allpts2$Metric = gsub(pattern="_SQ", replacement="", 
                      gsub(pattern="_HX", replacement="",     
                           gsub(pattern="_VO", replacement="", allpts2$Metric)))
allpts2 <- allpts2 %>%  mutate(Metric=factor(Metric,levels= rev(c("UAA", "Arable", "Arable_UAA", "Conv_farm", "Intens_f", "N_input", "ALU_homog",
                                                                 "LC_homog", "LULC_homog", "Field_size"))) )
# plot ridgeplots  (exported dimension: 700x550)
ggplot(allpts2, aes(x=Intensity, y=Metric, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,#scale=0.9    
    quantiles = 5, quantile_lines = TRUE, jittered_points = TRUE,aes(point_alpha=0.4, point_shape=20, point_stroke = 0.5, point_color="dark_red")) +   
    scale_fill_viridis_d(name = "Quantiles", alpha=0.5) +labs(x="Intensity",y='') +facet_wrap( ~Grid, labeller = "label_both")

## RIDGELINE PLOT WETLAND SPECIES ----------------
# load the sampled points
vspec = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Wetland_sampled_points_202030908.gpkg")
# join with LUI metrics from the hexagonal and voronoi grids (square grid info is already in the shapefile)
hexgrid = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/metrics/Hexagonal_metric_20220908.shp")
vor = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/metrics/Voronoi_metric_20220908.shp")
vspec = st_join(vspec, hexgrid, left=T)
vspec = st_join(vspec, vor, left=T)
names(vspec)[17:38] <- c("FID_HX","UAA_HX","Arable_HX", "Arable_UAA_HX", "Field_size_HX", "N_input_HX", "Intens_f_HX", "Conv_farm_HX", "ALU_homog_HX", "LULC_homog_HX", "LC_homog_HX",
                         "FID_VO","UAA_VO","Arable_VO", "Arable_UAA_VO", "Field_size_VO", "N_input_VO", "Intens_f_VO", "Conv_farm_VO", "ALU_homog_VO", "LULC_homog_VO", "LC_homog_VO")
# save to file
vsave = vspec[,c(1,2,4,7:16,18:27,29:38)]
st_write(vsave, "C:/Users/sroilo/Desktop/BESTMAP documents/Papers/Land use intensity/data_upload/VirtualSp_Wetland_sampled_points_20240120.gpkg")
# filter to retain only the presence points of the virtual species
vspec = vspec[vspec$presence==1,]  # NOTE: rows are more than 215 now, as some fall in between two hexagons, or two voronoi polygons -> some points are joined to multiple geometries
# drop geometry and reorder the variables
vspec = st_drop_geometry(vspec[,c(1,2,4,7:16,18:27,29:38)])
# convert from wide to long format
allpts2 = gather(vspec, key = Metric, Intensity, UAA_SQ:LC_homog_VO)
allpts2$Grid = "Square grid"
allpts2$Grid[grepl(pattern="_HX", x=allpts2$Metric)] <- "Hexagonal grid"
allpts2$Grid[grepl(pattern="_VO", x=allpts2$Metric)] <- "Voronoi grid"
allpts2 <- allpts2 %>%  mutate(Grid=factor(Grid,levels=unique(allpts2$Grid)) )
allpts2$Metric = gsub(pattern="_SQ", replacement="", 
                      gsub(pattern="_HX", replacement="",     
                           gsub(pattern="_VO", replacement="", allpts2$Metric)))
allpts2 <- allpts2 %>%  mutate(Metric=factor(Metric,levels= rev(c("UAA", "Arable", "Arable_UAA", "Conv_farm", "Intens_f", "N_input", "ALU_homog",
                                                                  "LC_homog", "LULC_homog", "Field_size"))) )
# plot ridgeplots  (exported dimension: 700x550)
ggplot(allpts2, aes(x=Intensity, y=Metric, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,#scale=0.9    
    quantiles = 5, quantile_lines = TRUE, jittered_points = TRUE,aes(point_alpha=0.4, point_shape=20, point_stroke = 0.5, point_color="dark_red")) +   
  scale_fill_viridis_d(name = "Quantiles", alpha=0.5) +labs(x="Intensity",y='') +facet_wrap( ~Grid, labeller = "label_both")


##### DISTANCE measures across the density distributions ------------
# compute the symmetrical Kullback-Leibler distance for each pair of memtrics, within and across grids
library(philentropy)
# function to compute symmetrical kullback_leibler_distance for each pair of distributions
sym_KLD = function(distr1, distr2, unit, epsilon) {
  kld1 = kullback_leibler_distance(distr1, distr2, testNA=F, unit=unit, epsilon = epsilon)
  kld2 = kullback_leibler_distance(distr2, distr1, testNA=F, unit=unit, epsilon = epsilon)
  skld = 0.5*sum(kld1,kld2)
  return(skld)
}

# retain only metrics in the dataframe, get rid of FID and other variables, and retain only presence points
allmets = vspec[vspec$presence==1,c(4:33)]
#remove potential NAs, and reorder columns to match the order of the metrics as presented in table 1 of the manuscript
allmets = na.omit(allmets[,c("UAA_SQ","Arable_SQ", "Arable_UAA_SQ", "Conv_farm_SQ","Intens_f_SQ","N_input_SQ","ALU_homog_SQ", "LULC_homog_SQ", "LC_homog_SQ", "Field_size_SQ",
                             "UAA_HX","Arable_HX", "Arable_UAA_HX", "Conv_farm_HX", "Intens_f_HX", "N_input_HX", "ALU_homog_HX", "LULC_homog_HX", "LC_homog_HX", "Field_size_HX",
                             "UAA_VO","Arable_VO", "Arable_UAA_VO", "Conv_farm_VO", "Intens_f_VO", "N_input_VO","ALU_homog_VO", "LULC_homog_VO", "LC_homog_VO", "Field_size_VO" )])
#  calculate KLD for each pair of metrics
kld_df = data.frame(row.names = names(allmets))
for ( i in c(1:ncol(allmets))) {
  met_i = names(allmets)[i]
  kld_df[,met_i] <- NA
  for (j in c(ifelse(i==30, 30, i+1):ncol(allmets))) {
    kld_df[j,i] <- sym_KLD(allmets[,i], allmets[,j], unit="log", epsilon = 0.00001)
  }
}
# write results to table
write.table(kld_df, sep=";", dec=",", row.names = T, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/Kullback_leibler_GrasslandSp_20230308.csv")
#write.table(kld_df, sep=";", dec=",", row.names = T, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/Kullback_leibler_WetlandSp_20230911.csv")
#write.table(kld_df, sep=";", dec=",", row.names = T, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/Kullback_leibler_ArableSp_20230911.csv")

### CORRELATION STRUCTURES ------------------------
# create a cross-grid cross-metric correlation matrix
library(corrplot)
cmat <- cor(allmets, method="spearman")
testRes = cor.mtest(cmat, conf.level = 0.95)
# Spearman correlation values
corrplot(cmat, order = 'original', method="color", addCoef.col = 'black', tl.pos = 'd',  tl.srt=45,
         type="upper",tl.cex=0.6, number.cex=0.6, col = COL2('PiYG'),
         title= "Spearman correlation matrix", mar=c(2,2,2,2)) 
# significance levels
corrplot(cmat, order = 'original',  method="color", p.mat = testRes$p, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.8, insig = 'label_sig', pch.col = 'black',
         tl.pos = 'a', type="upper",tl.cex=0.8,  col = COL2('PiYG'),title= "Spearman correlation matrix", mar=c(2,2,2,2)) 

# make paired plots for the sampled grid cells - hexagonal grid
vspec = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Grassland_sampled_points_202030305.gpkg")
hexgrid = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/metrics/Hexagonal_metric_20220908.shp")
vspec = st_join(vspec[,c("X","Y","presence")], hexgrid, left=T)
ploth = vspec[order(vspec$presence, decreasing=F),]
names(ploth)[5:14] = c( "UAA","Arable","Arable_UAA", "Field_size", "N_input","Intens_f","Conv_farm","ALU_homog","LULC_homog", "LC_homog")
ploth = st_drop_geometry(ploth[-which(duplicated(ploth[,c("X","Y")])),])
pairs.panels(ploth[, c("presence", "UAA", "Arable", "Arable_UAA", "Conv_farm", "Intens_f", "N_input", 
                       "ALU_homog", "LULC_homog", "LC_homog", "Field_size")], bg=c("blue","yellow")[as.factor(ploth$presence)], pch=21,
             method = "spearman", hist.col = "darkorange", density = TRUE, ellipses = FALSE, 
             main= "Spearman correlation coefficient - sampled points for the grassland species in the hexagonal grid" )

# make paired plots for the sampled grid cells - voronoi grid
vspec = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_Grassland_sampled_points_202030305.gpkg")
vor = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/metrics/Voronoi_metric_20220908.shp")
vspec = st_join(vspec[,c("X","Y","presence")], vor, left=T)
ploth = vspec[order(vspec$presence, decreasing=F),]
names(ploth)[5:14] = c( "UAA","Arable","Arable_UAA", "Field_size", "N_input","Intens_f","Conv_farm","ALU_homog","LULC_homog", "LC_homog")
ploth = st_drop_geometry(ploth)  # no point duplicates here
pairs.panels(ploth[, c("presence", "UAA", "Arable", "Arable_UAA", "Conv_farm", "Intens_f", "N_input", 
                       "ALU_homog", "LULC_homog", "LC_homog", "Field_size")], bg=c("blue","yellow")[as.factor(ploth$presence)], pch=21,
             method = "spearman", hist.col = "darkorange", density = TRUE, ellipses = FALSE, 
             main= "Spearman correlation coefficient - sampled points for the grassland species in the voronoi grid" )

# re-run the section of code above for each of the three virtual species, after loading the respective dataset.

rm(list=ls())
