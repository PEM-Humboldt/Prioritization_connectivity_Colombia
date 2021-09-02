rm(list = ls())
gc()


#library(devtools)
#library(remotes)
library(rgeos)
library(Makurhini)
library(raster)
library(rgdal)
library(xfun)
library(sf)
library(spatialEco)
library(geogrid)
library(rgeos)

#Set working directory
setwd("")

#load ecorregions
ecorregions_COL<-readOGR("./Regiones_bioticas/regiones_bioticos/Regiones_Bióticas.shp")

#load human footprint
human.footprint<-raster("./human_footprint/IHEH_2018.tif")

#create ecorregion list
eco_list<-split(ecorregions_COL,ecorregions_COL$Region)

#protected areas shapefile
PNN<- readOGR("./RUNAP/RUNAP_ajust_proj.shp")
PNN <- gBuffer(PNN, byid=TRUE, width=0)
PNN<-st_as_sf(PNN)
#calculate pa area
PNN$area<-as.numeric(st_area(PNN))
#select columns
PNN<-PNN[c("id_pnn","area")]

#Colombia forest
forest_col<-raster("forest/Bosque_col_300m.tif")
#water index
water_col<-raster("Ecosystem_services/water_index_300m.tif")
#carbon index
carbon_col<-raster("Ecosystem_services/carbon_index_300m.tif")
#species richness raster
species_richness<-raster("rasters_aves/aves-shp/mapas_riqueza/riqueza_aves_elev_eco_prj.tif")
#corridor raster (generated in linkage mapper)
corridors<-raster("connectivity/corridors/corridors_col_200k.tif")
#land rent raster
cost<-raster("costo_Colombia/total_rent_300m.tif")

#create list of rasters
r<-list(forest_col,water_col,carbon_col,species_richness,corridors,cost,
        human.footprint)

#make sure all rasters have the same extent
extend_all =
  function(rasters){
    extent(Reduce(extend,rasters))
  }

extent_r<-extend_all(r)

re = lapply(r, function(r){extend(r,extent_r)})


#calculate metrics for all pus within each ecorregions
results_conn<-list()

for (i in 1:length(eco_list)
     ){
  timestamp()
  #crop data to ecorregion extent
  #25000: buffer around ecorregion
  #crop hf to ecorregion extent
  hf_eco<-crop(human.footprint,extent(eco_list[[i]]) + 25000,snap = "near")
  #crop all rasters to ecorregion extent
  allrast<-crop(stack(re),extent(eco_list[[i]]) + 25000,snap = "near")
  #convert PNN to raster and crop to ecorregion extent
  PNN_eco<-crop(as(PNN,"Spatial"),extent(eco_list[[i]]) + 25000,snap = "near")
  PNN.raster<- rasterize(PNN_eco,hf_eco, mask=TRUE) #crops to polygon edge & converts to raster
  PNN.raster[which(!is.na(PNN.raster[]))]<-1
  #plot(PNN.raster)
  names(PNN.raster)<-"PNN"
  PNN.raster[which(is.na(PNN.raster[]))]<-0
  PNN.raster<-extend(PNN.raster,allrast)
  allrast<-stack(allrast,PNN.raster)
  
  eco_list[[i]] <- gBuffer(eco_list[[i]], byid=TRUE, width=0)
  
  #create grid PUs per ecorregion (hexagonal here)
  HexPts <-spsample(eco_list[[i]], type="hexagonal", cellsize=1200)
  HexPols <- HexPoints2SpatialPolygons(HexPts)
  #add ID to each PU
  HexPols$ID<-1:length(HexPols@polygons)
  
  #mean removing NAs
  mean_narm <- function(x, na.rm = T) { 
    mean(x)
  } 
  
  #Zonal statistics per PU
  system.time(ex <- spatialEco::zonal.stats(st_as_sf(HexPols),
                                            allrast, stats = "mean_narm"))
  
  #add colnames
  colnames(ex)<-names(allrast)
  #add data to grid shapefile
  HexPols@data<-cbind(HexPols@data,ex)
  
  #create data frame
  s_po_df<-as.data.frame(HexPols)
  
  colnames(s_po_df)<-c("ID","Bosque_col_300m",          
                       "water_index_300m","carbon_index_300m","biod",
                       "corridors_col_200k","total_rent_300m","IHEH_2018",                
                       "PNN" )
  
  s_po_df$Region<-eco_list[[i]]$Region
  
  #save dataframe
  write.csv(s_po_df,paste0("./Regiones_bioticas/","grid_",eco_list[[i]]$Region,"_calculations.csv"))
  
  #save PU grid shapefile
  writeOGR(HexPols,"./Regiones_bioticas/grids_eco",
           paste0("grid_eco","_", eco_list[[i]]$Region,"_","comp_PNN"),
           driver="ESRI Shapefile")
  
  print(i)
  gc()
  
}
  