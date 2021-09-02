rm(list = ls())
gc()

library(devtools)
library(remotes)
library(rgeos)
#install_github("connectscape/Makurhini", dependencies = TRUE, upgrade = "never")
library(Makurhini)
library(raster)
library(rgdal)
library(xfun)
library(sf)
library(geogrid)
library(rgeos)
library(plyr)

setwd("")

#Read all data####

#Ecorregiones

ecorregions_COL<-readOGR("./Regiones_bioticas/regiones_bioticos/Regiones_Bióticas.shp")
eco_list<-split(ecorregions_COL,ecorregions_COL$Region)


#human footprint

human.footprint<-raster(paste0("./human_footprint","/","IHEH_2018.tif"))

#Protected area network

PNN<- readOGR("./RUNAP/RUNAP_WDPA_prj.shp")
PNN <- gBuffer(PNN, byid=TRUE, width=0)
PNN<-st_as_sf(PNN)
#calculate area of each PA
PNN$area<-st_area(PNN)
PNN$area<-as.numeric(PNN$area)
#add ID to PA
PNN$id_pnn<-1:nrow(PNN)
PNN<-PNN[c("id_pnn","area")]

#prepare grid calculation data
grid_cals<-list.files("./Regiones_bioticas/grid_calculations",full.names = T)
grid_cals<-lapply(grid_cals,read.csv)
alldata<-do.call(rbind.fill,grid_cals)
head(alldata)


#define scenario combinations####
#create dataframe showing objective and indices combinations 
#names must be the same as the ones used in the prioritization analysis
#this is done to call all the prioritization files

scenarios<-data.frame(objective_1= c("carbon_ratio","carbon_ratio","carbon_ratio",
                                     "water_ratio","water_ratio","corridors_rent"),
                      objective_2= c("water_ratio","corridors_rent","biod_ratio",
                                     "corridors_rent",
                                     "biod_ratio","biod_ratio"),
                      index_1 = c("carbon_index_300m","carbon_index_300m","carbon_index_300m",
                                  "water_index_300m","water_index_300m","corridors"),
                      index_2 = c("water_index_300m","corridors","biod",
                                  "corridors","biod","biod"))

#Calculate connectivity indices

for(i in 1:length(eco_list)){
  
  #select data ecorregion of interest
  s_po_df<-subset(alldata,alldata$Region == eco_list[[i]]$Region)
  
  #get PUs outside PNN
  s_po_df<-subset(s_po_df,s_po_df$PNN == 0)

  #get PA data inside ecorregion of interest
  
  pnn_eco<-readOGR(paste0("./Regiones_bioticas/pnn_per_eco/eco_",
                          i,"_PNN.shp"))
  
  #get PU shapefile (created in file "calculate_metrics_per_pu)
  
  grid_eco<-readOGR(paste0("./Regiones_bioticas/grids_eco/grid_eco_",
                           i, "_comp_PNN.shp"))
  
  #Get PUs outside PAs
  grid_eco<-subset(grid_eco,grid_eco$PNN == 0)
  
  
  conn_results<-list()
  for (j in 1:nrow(scenarios)) {
    
    sc_obj1<-read.csv(paste0("calculations_ranking_analysis/50k/target_30/",
                             "eco_",i,"_",scenarios[j,1],"_",scenarios[j,2],
                             "_weight_1_30.csv"))
    sc_obj2<-read.csv(paste0("calculations_ranking_analysis/50k/target_30/",
                             "eco_",i,"_",scenarios[j,1],"_",scenarios[j,2],
                             "_weight_last.csv"))
    
    obj1_PU<-subset(grid_eco,grid_eco$ID %in% sc_obj1$ID)
    obj1_PU$diss<-1
    obj1_PU$sce<-scenarios[j,1]
    obj2_PU<-subset(grid_eco,grid_eco$ID %in% sc_obj2$ID)
    obj2_PU$diss<-1
    obj2_PU$sce<-scenarios[j,2]
    
    
    #classify according to HF categories
    obj1_PU[which(obj1_PU$IHEH_20 <= 15),"cat"]<-"Natural"
    obj1_PU[which(obj1_PU$IHEH_20 > 15 & obj1_PU$IHEH_20 <= 40),"cat"]<-"Low"
    obj1_PU[which(obj1_PU$IHEH_20 > 40 & obj1_PU$IHEH_20 <= 60),"cat"]<-"Mid"
    obj1_PU[which(obj1_PU$IHEH_20 > 60),"cat"]<-"High"
    
    obj2_PU[which(obj2_PU$IHEH_20 <= 15),"cat"]<-"Natural"
    obj2_PU[which(obj2_PU$IHEH_20 > 15 & obj2_PU$IHEH_20 <= 40),"cat"]<-"Low"
    obj2_PU[which(obj2_PU$IHEH_20 > 40 & obj2_PU$IHEH_20 <= 60),"cat"]<-"Mid"
    obj2_PU[which(obj2_PU$IHEH_20 > 60),"cat"]<-"High"
    
    weight_list<-list(obj1_PU,obj2_PU)
    df<-list()
    for(s in 1:length(weight_list)){
      eco_newareas<-gUnaryUnion(weight_list[[s]],id = weight_list[[s]]@data$diss)
      crs(eco_newareas)<-crs(pnn_eco)
      eco_newareas<-st_as_sf(eco_newareas)
      #st_crs(eco_newareas)<-st_crs(pnn_eco)
      #eco_newareas<-st_as_sf(eco_newareas)
      eco_newareas<-st_cast(eco_newareas,"POLYGON")
      eco_newareas$area<-as.numeric(st_area(eco_newareas))
      
      eco_newareas$id_pnn<-paste0("new_", 1:nrow(eco_newareas))
      
      #new_PNN<-eco_newareas
      new_PNN<-rbind(st_as_sf(pnn_eco["id_pnn"]),eco_newareas["id_pnn"])
      
      eco_focus<-st_as_sf(eco_list[[i]])
      
      
      #plot(eco_focus,add = T)
      #plot(eco_newareas[1],add = T, col = "red")
      #plot(new_PNN[1],add = T, col = "green")
      
      #calculate connectivity for complementary areas####
      #define distance thresholds
      dist.thr<-c(1000,5000,10000,50000)
      system.time(conn_comp<-try(Makurhini::MK_ProtConn(new_PNN,
                                                    eco_focus,
                                                    thintersect = 0,
                                                    distance = list(type = "centroid", 
                                                                    resistance = human.footprint),
                                                    #attribute = "area",
                                                    #area_unit= "m2",
                                                    keep = 0.5,
                                                    probability = 0.5,
                                                    transboundary = 50000,
                                                    distance_thresholds= dist.thr)))
      results_conn<-list()
      if(isTRUE(class(conn_comp)=="try-error")| isTRUE(class(conn_comp)=="try-error")){
        results_conn[[i]]<-as.data.frame(matrix(nrow = 1,ncol = 21,data = "NA"))
        colnames(results_conn[[i]])<-colnames(results_conn[[2]])
        results_conn[[i]]$ECO<-eco_list[[i]]$FIRST_Dist
      }else{
        for(d in 1:length(dist.thr)){
          test2<-as.data.frame(conn_comp[[d]])
          colnames(test2)<-c("Index","Value","ProtConn indicator","Percentage")
          complement<-as.data.frame(t(test2[c(3,4)]))
          colnames(complement)<-test2$`ProtConn indicator`
          complement<-complement[2,]
          complement$ECO<-eco_list[[i]]$FIRST_Dist
          complement$year<-strsplit(names(human.footprint),"_")[[1]][[2]]
          complement$EC_PC<-test2$Value[[1]]
          complement$PC<-test2$Value[[2]]
          complement$sce<-weight_list[[s]]$sce[[1]]
          complement$distance<-dist.thr[d]
          results_conn[[d]]<-complement
          gc()
        }
      }
      df[[s]]<-do.call(rbind,results_conn)
    }
    
    conn_results[[j]]<-do.call(rbind,df)
    
    #MUCH FASTER WHEN THE POLYGON IS DISSOLVED
    #SIMPLIFY POLYGON
    #HERE PROTCONN_WITHIN AND PROTCONN_CONTIG 
    #GIVE DIFFERENTE RESULTS IF THE POLYGON IS DISSOLVED OR NOT
    #EVERYTHING ELSE IS THE SAME
    
    gc()
    timestamp()
    #print(j)
    systime2 <- Sys.time() 
  }
  
  #calculate connectivity without complementary areas####
  
  system.time(conn_cur<-try(Makurhini::MK_ProtConn(pnn_eco,
                                                    eco_focus,
                                                    thintersect = 0,
                                                    distance = list(type = "centroid", 
                                                                    resistance = human.footprint),
                                                    #attribute = "area",
                                                    #area_unit= "m2",
                                                    keep = 0.5,
                                                    probability = 0.5,
                                                    transboundary = 50000,
                                                    distance_thresholds= dist.thr)))
  results_conn_cur<-list()
  if(isTRUE(class(conn_comp)=="try-error")| isTRUE(class(conn_comp)=="try-error")){
    results_conn_cur[[i]]<-as.data.frame(matrix(nrow = 1,ncol = 21,data = "NA"))
    colnames(results_conn_cur[[i]])<-colnames(results_conn_cur[[2]])
    results_conn_cur[[i]]$ECO<-eco_list[[i]]$FIRST_Dist
  }else{
    for(d in 1:length(dist.thr)){
      test2<-as.data.frame(conn_cur[[d]])
      colnames(test2)<-c("Index","Value","ProtConn indicator","Percentage")
      complement<-as.data.frame(t(test2[c(3,4)]))
      colnames(complement)<-test2$`ProtConn indicator`
      complement<-complement[2,]
      complement$ECO<-eco_list[[i]]$FIRST_Dist
      complement$year<-strsplit(names(human.footprint),"_")[[1]][[2]]
      complement$EC_PC<-test2$Value[[1]]
      complement$PC<-test2$Value[[2]]
      complement$sce<-"PNN"
      complement$distance<-dist.thr[d]
      results_conn_cur[[d]]<-complement
      gc()
    }
  }
  
  #merge connectivity results
  results_conn<-do.call(rbind,conn_results)
  conn_init<-do.call(rbind,results_conn_cur)

  allresults<-rbind(results_conn,conn_init)  
  
  write.csv(allresults,paste0("results_conn_t30_eco_",i,".csv"))
  }
  