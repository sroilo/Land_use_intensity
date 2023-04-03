###################################################
#
# Title: Roilo_2023_LUI_virtual_species.R
# Purpose: create a virtual species with known land-use intensity (LUI)-species relationships to test how using different subsets of LUI metrics affects the outcomes of 
# biodiversity models (we here used Generalised Additive Models to model the virtual species occurrence).
# Reference: Roilo et al. "A multidimensional approach is needed to better quantify land-use intensity in biodiversity models".
# Author: Stephanie Roilo, Technische Universität Dresden
# Date: last updated on March 13th, 2023.
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
library(mapview)
library(terra)
library(psych)   # to plot nice pairs.plots
library(MuMIn)
library(DHARMa)  # to simulate residuals and check model fit
library(ncf)  # to check for spatial autocorrelation in the model residuals
library(mgcv)
library(ggplot2)   # for plotting ridgeplots
library(ggridges)  # for plotting ridgeplots
library(tidyr)    # for handling data

# CREATE A VIRTUAL SPECIES -------------------------
# check the virtualspecies tutorial here: http://borisleroy.com/files/virtualspecies-tutorial.html 
# load the previously computed LUI metrics in the square grid, and prepare the raster layers to be used to set up the virtual species
sqgrid = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/metrics/Square_metric_20220908.shp")
r = rast(crs="epsg:3035",extent=ext(sqgrid), resolution=c(1000,1000))
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
vsax_parameters = formatFunctions(Arable_UAA = c(fun = 'logisticFun', beta = 0.6, alpha = 0.1),
                                  LULC_homog = c(fun = 'linearFun', a=-0.95, b=1),
                                  Intens_f = c(fun = 'linearFun', a = -0.8, b = 1)) 
# create the virtual species
vsax <- generateSpFromFun(raster.stack = env, 
                          parameters = vsax_parameters, 
                          species.type="multiplicative",
                          rescale = TRUE, rescale.each.response = F, plot = TRUE)
# plot species response plots
plotResponse(vsax)
# save suitability raster to file
writeRaster(vsax$suitab.raster, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_suitability_20230305F.tif")

#convert the suitability raster into a binary presence-absence raster using a logistic curve
set.seed(13)
vsax_PA <- convertToPA(vsax, PA.method = "probability", prob.method = "logistic", beta = 0.3, alpha = -0.05, plot=T)
# save rasters to file
writeRaster(vsax_PA$probability.of.occurrence, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_ProbOcc_20230305F.tif")
writeRaster(vsax_PA$pa.raster, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_PA_20230305F.tif")

# sample virtual species occurrence at 500 random points, only in grid cells which have UAA>0
uaa = as(sqgrid[sqgrid$UAA>0,], "Spatial")
set.seed(13)
vsax_points <- sampleOccurrences(vsax_PA, n = 500, type = "presence-absence", sampling.area = uaa, 
                                 detection.probability = 1, extract.probability=TRUE, plot=T)
# convert to sf object
vsax_points <- vsax_points$sample.points %>% st_as_sf(coords = c("x", "y"), crs=3035)
sum(vsax_points$Observed)  # 99 presence points out of 500

#### MODELLING VIRTUAL SPECIES OCCURRENCE --------------------
# extract the LUI metrics at the 500 randomly sampled points
vsax_points = st_join(vsax_points, sqgrid, left=T)
names(vsax_points)[4:14] <- c("FID_SQ", "UAA_SQ","Arable_SQ", "Arable_UAA_SQ", "Field_size_SQ", "N_input_SQ", "Intens_f_SQ", "Conv_farm_SQ", 
                              "ALU_homog_SQ", "LULC_homog_SQ", "LC_homog_SQ")
names(vsax_points)[2] <- "presence"
# check if points fall in NA areas of the grid
sum(is.na(vsax_points$FID_SQ))  # 0
# add spatial coordinates to the dataframe
vsax_points = bind_cols(st_coordinates(vsax_points), vsax_points)
# save to file
#st_write(vsax_points, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_sampled_points_SQ_202030305F.gpkg")

### MODELLING -----------------------
sax_sq = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_sampled_points_SQ_202030305F.gpkg") %>% st_drop_geometry()
sum(sax_sq$presence) # 500 points, 99 presences
# check correlation among variables
plotd = sax_sq[order(sax_sq$presence, decreasing=F),]
names(plotd)[7:16] = c( "UAA","Arable","Arable_UAA", "Field_size", "N_input","Intens_f","Conv_farm","ALU_homog","LULC_homog", "LC_homog")
pairs.panels(plotd[, c("presence", "UAA", "Arable", "Arable_UAA", "Conv_farm", "Intens_f", "N_input", 
                       "ALU_homog", "LULC_homog", "LC_homog", "Field_size")], bg=c("blue","yellow")[as.factor(plotd$presence)], pch=21,
             method = "spearman", hist.col = "darkorange", density = TRUE, ellipses = FALSE, main= "Spearman correlation coefficients - sampled points in the square grid" )
# prepare data for GAM modelling
dat = sax_sq[,c(1,2,4,7:16)]
cmat <- abs(cor(dat[,4:13], method="spearman")) <= 0.7
cmat[!lower.tri(cmat, diag=F)] <- NA  

# fit full model and dredge
fits = gam(presence ~ s(UAA_SQ,k=5) + s(Arable_SQ,k=5) + s(Arable_UAA_SQ,k=5) + s(Field_size_SQ,k=5) + s(N_input_SQ,k=5) + s(Intens_f_SQ,k=5) + s(Conv_farm_SQ,k=5) + 
             s(ALU_homog_SQ,k=5) + s(LULC_homog_SQ,k=5) + s(LC_homog_SQ,k=5), data= dat, family="binomial", na.action="na.fail")
dr1 = dredge(fits, subset=cmat, extra = c("R^2", "adjR^2"))  
write.table(dr1, sep=";", dec=",", row.names=F, "C:/Users/sroilo/Desktop/BESTMAP documents/Papers/Land use intensity/revision/Virtual_square_dredge_20230306_GAM.csv")

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
plot(fits,residuals=FALSE,all.terms=TRUE,shade=TRUE,shade.col=2, seWithMean=FALSE) 

# best model with only land use and management metrics: 
sbest1 = gam(presence ~ s(UAA, k=5) + s(Arable,k=5) + s(Arable_UAA,k=5)+ s(Intens_f,k=5), data= dat, family="binomial", na.action="na.fail")
summary(sbest1)   # Arable and Intens_f NOT significant!
par(mfrow=c(1,4))
plot(sbest1,residuals=FALSE,all.terms=TRUE,shade=TRUE,shade.col=2, seWithMean=FALSE) 
simulateResiduals(sbest1, plot=T)  # OK
scor = ncf::spline.correlog(x = dat$X, y = dat$Y, z = sbest1$residuals, resamp=100, xmax = max.dist)
plot(scor)   # NO PROBLEM !

# best model with only landscape structure and management: 
sbest2 = gam(presence ~ s(Intens_f,k=5) + s(N_input,k=5) + s(LULC_homog,k=5) + s(Field_size,k=5), data= dat, family="binomial", na.action="na.fail")
summary(sbest2)   # all significant
par(mfrow=c(1,4))
plot(sbest2,residuals=FALSE,all.terms=TRUE,shade=TRUE,shade.col=2, seWithMean=FALSE) 
# diagnostic plots with DHARMa package
simulateResiduals(sbest2, plot=T) # OK
scor = ncf::spline.correlog(x = dat$X, y = dat$Y, z = sbest2$residuals, resamp=100, xmax = max.dist)
plot(scor)   # NO PROBLEM !

# best model with only land use landscape structure: 
sbest3 = gam(presence ~  s(UAA,k=5) + s(Arable_UAA,k=5) + s(LULC_homog,k=5), data= dat, family="binomial", na.action="na.fail")
summary(sbest3)   # all significant
par(mfrow=c(1,4))
plot(sbest3,residuals=FALSE,all.terms=TRUE,shade=TRUE,shade.col=2, seWithMean=FALSE) 
# diagnostic plots with DHARMa package
simulateResiduals(sbest3, plot=T) # OK
scor = ncf::spline.correlog(x = dat$X, y = dat$Y, z = sbest3$residuals, resamp=100, xmax = max.dist)
plot(scor)   # NO PROBLEM !


### RIDGEPLOTS of LUI values at the presence points of the virtual species -------
# load the sampled points
vspec = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/virtualSpecies/VirtualSp_sampled_points_SQ_202030305F.gpkg")
# join with LUI metrics from the hexagonal and voronoi grids (square grid info is already in the shapefile)
hexgrid = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/metrics/Hexagonal_metric_20220908.shp")
vor = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/metrics/Voronoi_metric_20220908.shp")
vspec = st_join(vspec, hexgrid, left=T)
vspec = st_join(vspec, vor, left=T)
names(vspec)[17:38] <- c("FID_HX","UAA_HX","Arable_HX", "Arable_UAA_HX", "Field_size_HX", "N_input_HX", "Intens_f_HX", "Conv_farm_HX", "ALU_homog_HX", "LULC_homog_HX", "LC_homog_HX",
                        "FID_VO","UAA_VO","Arable_VO", "Arable_UAA_VO", "Field_size_VO", "N_input_VO", "Intens_f_VO", "Conv_farm_VO", "ALU_homog_VO", "LULC_homog_VO", "LC_homog_VO")
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
                                                              "LULC_homog", "LC_homog", "Field_size"))) )
# plot ridgeplots for all grid types and all metrics
# ridgeplot tutorials: https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html 
# and https://jtr13.github.io/cc19/ridgeline-plots.html 
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
sym_KLD(vspec$UAA_SQ, vspec$UAA_VO, unit="log", epsilon = 0.00001) # works!
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
write.table(kld_df, sep=";", dec=",", row.names = T, "C:/Users/sroilo/Desktop/BESTMAP documents/Papers/Land use intensity/revision/Kullback_leibler_divergence_20230308.csv")

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
vspec = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Papers/Land use intensity/revision/VirtualSp_sampled_points_SQ_202030305F.gpkg")
hexgrid = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/metrics/Hexagonal_metric_20220908.shp")
vspec = st_join(vspec[,c("X","Y","presence")], hexgrid, left=T)
# plot paired plots (square grid one was already done previously)
ploth = vspec[order(vspec$presence, decreasing=F),]
names(ploth)[5:14] = c( "UAA","Arable","Arable_UAA", "Field_size", "N_input","Intens_f","Conv_farm","ALU_homog","LULC_homog", "LC_homog")
ploth = st_drop_geometry(ploth[-which(duplicated(ploth[,c("X","Y")])),])
pairs.panels(ploth[, c("presence", "UAA", "Arable", "Arable_UAA", "Conv_farm", "Intens_f", "N_input", 
                       "ALU_homog", "LULC_homog", "LC_homog", "Field_size")], bg=c("blue","yellow")[as.factor(ploth$presence)], pch=21,
             method = "spearman", hist.col = "darkorange", density = TRUE, ellipses = FALSE, 
             main= "Spearman correlation coefficient - sampled points in the hexagonal grid" )

# make paired plots for the sampled grid cells - voronoi grid
vspec = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Papers/Land use intensity/revision/VirtualSp_sampled_points_SQ_202030305F.gpkg")
vor = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land use intensity/metrics/Voronoi_metric_20220908.shp")
vspec = st_join(vspec[,c("X","Y","presence")], vor, left=T)
# plot paired plots (square grid one was already done previously)
ploth = vspec[order(vspec$presence, decreasing=F),]
names(ploth)[5:14] = c( "UAA","Arable","Arable_UAA", "Field_size", "N_input","Intens_f","Conv_farm","ALU_homog","LULC_homog", "LC_homog")
ploth = st_drop_geometry(ploth)  # no point duplicates here
pairs.panels(ploth[, c("presence", "UAA", "Arable", "Arable_UAA", "Conv_farm", "Intens_f", "N_input", 
                       "ALU_homog", "LULC_homog", "LC_homog", "Field_size")], bg=c("blue","yellow")[as.factor(ploth$presence)], pch=21,
             method = "spearman", hist.col = "darkorange", density = TRUE, ellipses = FALSE, 
             main= "Spearman correlation coefficient - sampled points in the voronoi grid" )
