#____________________________________#
# Code written by Felipe Suarez    #
# Version :  01-07-2020            #
#_____________________________________#

#clear workspace
rm(list=ls())

library(SpaDES)
library(raster)
#library(spatial.tools)
library(xfun)
library(rgdal)
library(ReIns)


setwd("C:/Data/50Reefs/Felipe/iucn/Cambio_climatico/refugios_2030/2021-2040_ssp585_10")

# load sp rasters

myfiles = list.files(pattern="*.tif") 

list_raster<-list()
for (i in 1:length(myfiles)){
  print(file_ext(myfiles[[i]]))
  if(file_ext(myfiles[[i]]) == "xml"|file_ext(myfiles[[i]]) == "dbf"
     |file_ext(myfiles[[i]]) == "ovr"|file_ext(myfiles[[i]]) == "cpg"
     |file_ext(myfiles[[i]]) == "tfw"|file_ext(myfiles[[i]]) == "lock"){
    list_raster[[i]]<-NULL
  }else{
    list_raster[[i]]<-myfiles[[i]]
    
  }
}

list_raster<-plyr::compact(list_raster)

#READ RUNAP data in raster format
setwd("C:/Data/50Reefs/Felipe/iucn/RUNAP/RUNAP_rasters")
PNN<-raster("./RUNAP_2018_wgs1km.tif")

#get raster resolution to calculate areas:

#pixel size from raster(when using projected coordinated system: 
res_r<-res(PNN)[[1]]^2
#or specify it manually:
res_r<-1 #km2



#this loop generates the statistics for all the species
sp_stats<-list()#list to save stats per species
  
for (j in 1:length(list_raster)){#list_pol is the list with all the focus_sps
    setwd("C:/Data/50Reefs/Felipe/iucn/Cambio_climatico/refugios_2030/2021-2040_ssp585_10")
    focus_sp<-raster(list_raster[[j]])#read each sp raster
    if (is.null(intersect(extent(PNN), focus_sp))){#check if rasters intersect
      dfsp<-as.data.frame(matrix(nrow = 1,ncol = 3))
      colnames(dfsp)<-c("PNN","freq")
      dfsp$species<-names(focus_sp)
      dfsp$year_PNN<-strsplit(names(PNN),"_")[[1]][[2]]
    }else{
      #convert raster to df
      cells<-rasterToPoints(focus_sp)
      cells<-as.data.frame(cells)
      cells<-subset(cells,cells[3]==1)
      #extract PNN values that intersect the sp distribution area
      PNN_data<-raster::extract(PNN,cells[1:2])
      cells$PNN<-PNN_data
      if(nrow(cells)>0){
        hola<-plyr::count(cells[4])
        hola$species<-names(focus_sp)
        #hola$year_PNN<-strsplit(names(PNN),"_")[[1]][[2]]
      }
    }
    #assign values depending on whether pixels are inside or outside the PNN
    hola[which(is.na(hola$PNN)),"in_PNN"]<-0
    hola[which(!is.na(hola$PNN)),"in_PNN"]<-1
    
    #area inside and outside PNN
    count_pnn<-aggregate(hola$freq,list(PNN=hola$in_PNN),"sum")
    total_area<-sum(count_pnn$x)
    
    total_protected<-count_pnn[which(count_pnn$PNN == "1"),2]
    if(length(total_protected)>0){
      hola$per_protected<-(total_protected/total_area)*100
    }else{
      hola$per_protected<-0
    }
    
    #area in km2
    hola$total_area<-total_area
    hola$total_area_km2<-hola$total_area*res_r
    
    hola<-hola[!duplicated(hola[5:7]),]
    hola<-hola[c("species","per_protected","total_area","total_area_km2")]
    
    print(paste(j))
    sp_stats[[j]]<-hola
    gc()
  }

#df with the following columns:
stats_allsp<-do.call(rbind,sp_stats)

#species: sp file name
#per_protected: percentage of the sp distribution area in Colombia that is protected
#total_area_km2: area of  sp distribution in km2

#Calculate percentage needed to protect#########

#prop_achieved: proportion of distribution area to be protected that is actually protected
#achieved: 0 if the prop_achieved is <0.9, 1 if it is >0.9

df<-stats_allsp

#per_to_protect: percentage of the sp distribution area to be to protect
df$per_to_protect <- (1 -(pexp(df$total_area/1e5,rate=0.5)))*100

df$prot_gap<-(df$per_to_protect - df$per_protected)/100

df[which(df$prot_gap<0),"prot_gap"]<-0

write.csv(stats_allsp,"statistics_sp_2021_2040_ssp126_10.csv")

#Create raster to map species based on conservation gap###

df_con<-subset(df,df$prot_gap > 0)
head(df_con)
list_raster<-df_con$species


focus_sp<-list_raster[[1]]
rs1<-raster(focus_sp)
df_sp<-subset(df_con,df_con$species == focus_sp)
rs1[which(rs1[] > 0)]<-df_sp$prot_gap

for(i in 2:length(list_raster)){
  focus_sp<-list_raster[[i]]
  df_sp<-subset(df_con,df_con$species == focus_sp)
  rast<-raster(focus_sp)
  rast<-extend(rast,rs1)
  rs1<-extend(rs1,rast)
  rast[which(rast[] > 0)]<-df_sp$prot_gap
  init_rast<-stack(rs1,rast)
  rs1 <- calc(init_rast, sum, na.rm = T)
  print(i)
  gc()
}

