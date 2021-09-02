rm(list = ls())#remove all the elements (clear workspace)
gc()

library(tidyverse)
library(magrittr)
library(rio) ## Packege to read large tables
library(plyr) ## to manupulate tables
library(foreign) 
library(data.table)
library(sf)
library(raster)
library(gdistance)
library(rgeos)
library(rgdal)

setwd("")


#load ecorregions

ecorregions_COL<-readOGR("./Regiones_bioticas/regiones_bioticos/Regiones_Bióticas.shp")

#list of ecorregions
eco_list<-split(ecorregions_COL,ecorregions_COL$Region)

#human footprint

human.footprint<-raster(paste0("./human_footprint","/","IHEH_2018.tif"))

#PNN
PNN<- readOGR("./RUNAP/RUNAP_ajust_proj.shp")
PNN <- gBuffer(PNN, byid=TRUE, width=0)
PNN<-st_as_sf(PNN)
PNN$area<-st_area(PNN)
PNN$area<-as.numeric(PNN$area)
PNN<-PNN[c("id_pnn","area")]

#prepare grid calculation data
grid_cals<-list.files("./Regiones_bioticas/grid_calculations",full.names = T)
grid_cals<-lapply(grid_cals,read.csv)
alldata<-do.call(rbind.fill,grid_cals)
head(alldata)

#define objective function
# Scenario-- theta can be only from 0-1 as it is a weight... 'para' = weight  
para <- seq(0,1,0.1)

#objective function:
#remember that data needs to be normalized
obj_fun<-function(para,#weigths
                  obj_1,#objective 1
                  obj_2 #objective 2
){
  obj_para <-(para*obj_1)+ ((1-para)*obj_2)
  return(obj_para)
}

#generate objective combinations####
#the dataframe will allow to run scenarios
#by calling the objectives of interest
#column names represent the objective values that will be generated
# in the loop for each ecorregion
#
scenarios<-data.frame(objective_1= c("carbon_ratio","carbon_ratio","carbon_ratio",
                                     "water_ratio","water_ratio","corridors_rent"),
                      objective_2= c("water_ratio","corridors_rent","biod_ratio",
                                     "corridors_rent",
                                     "biod_ratio","biod_ratio"),
                      index_1 = c("carbon_index_300m","carbon_index_300m","carbon_index_300m",
                                  "water_index_300m","water_index_300m","corridors"),
                      index_2 = c("water_index_300m","corridors","biod",
                                  "corridors","biod","biod"))


# % to protect
target = 30 

for(i in 2:length(eco_list)){
  
  #get human footprint per ecorregion####
  hf_eco<-crop(human.footprint,extent(eco_list[[i]]),snap = "near")
  #protected areas inside ecorregions
  PNN_eco<-crop(as(PNN,"Spatial"),extent(eco_list[[i]]),snap = "near")
  PNN.polygon<- rasterize(PNN_eco,hf_eco, mask=TRUE) #crops to polygon edge & converts to raster
  PNN.polygon[which(!is.na(PNN.polygon[]))]<-1
  
  #create raster ecorregion######
  
  ecorregion.raster<- rasterize(eco_list[[i]],hf_eco, mask=TRUE)
  ecorregion.raster[which(is.na(ecorregion.raster[]))]<--1
  eco_raster<-ecorregion.raster
  eco_raster[which(ecorregion.raster[]>-1)]<-1
  eco_raster[which(PNN.polygon[]>0)]<-2
  eco_raster[which(ecorregion.raster[]==-1)]<--1
  
  #extract values
  ecorregion_df<-as.data.frame(rasterToPoints(eco_raster))
  
  #calculate proportion area protected per ecorregion####
  
  eco_protected<-nrow(subset(ecorregion_df,ecorregion_df[3] == 2))
  eco_total<-nrow(subset(ecorregion_df,ecorregion_df[3] > 0))
  prop_protected<-(eco_protected/eco_total)*100
  
  #read data of planning units####
  
  s_po_df<-subset(alldata,alldata$Region == eco_list[[i]]$Region)
  
  #remove 0s from lant rent values so we do not have to divide by 0
  hello<-subset(s_po_df, s_po_df$total_rent_300m > 0 & s_po_df$total_rent_300m < 1)
  
  if(nrow(hello)> 0){
    s_po_df[which(s_po_df$total_rent_300m == 0),"total_rent_300m"]<-min(hello$total_rent_300m)
  }else{
    s_po_df[which(s_po_df$total_rent_300m == 0),"total_rent_300m"]<-0.1
  }
  
  # function to standardize data####
  range01 <- function(x, na.rm = T){(x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))}
  
  #standardize corridor values####
  
  #corridor values correspond to the distance in m from the best corridor
  #the larger the value in the original table, the farthest it is from the
  #potential corridors. We convert NAs and cells > 50 km
  #to the maximum distance value (200 km) to ensure
  #that these cells have the lowest chance of being selected in 
  #the prioritization analysis
  
  s_po_df[which(is.na(s_po_df$corridors_col_200k)),"corridors_col_200k"]<-2e5
  s_po_df[which(s_po_df$corridors_col_200k > 5e4),"corridors_col_200k"]<-2e5
  
  #transform corridor values
  #transform to meters
  s_po_df$corridors_col_200k<-s_po_df$corridors_col_200k/1000
  s_po_df$corridors<- range01(s_po_df$corridors_col_200k)
  s_po_df$corridors<-1-s_po_df$corridors
  
  #only consider forest corridors
  #s_po_df$corridors<-s_po_df$Bosque_col_300m*s_po_df$corridors
  
  #divide by rent
  s_po_df$corridors_rent<-s_po_df$corridors/s_po_df$total_rent_300m
  s_po_df[which(s_po_df$corridors_col_200k == 2e5/1000),"corridors_rent"]<-min(s_po_df$corridors_rent, na.rm = T)
  
  
  #standardize ES and biodiversity values####
  s_po_df$carbon_ratio<-s_po_df$carbon_index_300m/s_po_df$total_rent_300m
  s_po_df$water_ratio<-s_po_df$water_index_300m/s_po_df$total_rent_300m
  #biod values
  s_po_df$biod_ratio<- range01(s_po_df$biod)
  s_po_df$biod_ratio<-s_po_df$biod_ratio/s_po_df$total_rent_300m
  
  #change NAs to 0
  s_po_df[which(is.na(s_po_df$carbon_ratio)),"carbon_ratio"]<-0
  s_po_df[which(is.na(s_po_df$water_ratio)),"water_ratio"]<-0
  s_po_df[which(is.na(s_po_df$biod_ratio)),"biod_ratio"]<-0
  
  #-------------------- Set up benefit-cost ratio targeting model -----------------
  #budget here represents the total area that needs to be protected to reach target
  budget <- (eco_total*((target - prop_protected)/100))*0.3*0.3
  print(budget)
  print(prop_protected)
  
  #READ GRIDS and PNN around ecorregion###
  grid_eco<-readOGR(paste0("./Regiones_bioticas/grids_eco/grid_eco_",
                           eco_list[[i]]$Region, "_comp_PNN.shp"))
  
  #subset planning units outside PNN
  grid_eco<-subset(grid_eco,grid_eco$PNN == 0)
  
  #select PUs outside PAs
  s_po_df<-subset(s_po_df,s_po_df$PNN == 0)

  
  #pb <- txtProgressBar(min = 0, max = length(budget), style = 3) #To check the progress
  
  #Apply objective function to all the scenario combinations
  systime1 <- Sys.time()
  conn_results<-list()
  for (j in 1:nrow(scenarios)) {
    #apply function to all dataframe rows using all the parameter values
    vector_t<- apply(s_po_df, 1, function(x)obj_fun(para,x[scenarios[j,1]],x[scenarios[j,2]])) #apply obj function
    Scenario_ES<-(as.data.frame(vector_t))
    Scenario_ES<-as.data.frame(t(Scenario_ES))
    colnames(Scenario_ES)<-para
    Scenario_ES$ID<-s_po_df$ID
    Scenario_ES<-merge(Scenario_ES,s_po_df, by = "ID") #
    
    #we are not using a present net value; use the area of each PU
    #as the "budget spent"
    Scenario_ES$overall_npv<-1.2
    
    #setTxtProgressBar(pb, b)
    df_para<-list() 
    df_para_subscenario<-list()
    #To save results for pareto (obj_1, obj_2, theta and budget)
    pareto_results<-data.frame(matrix(0,0,4)) 
    colnames(pareto_results)<-c(paste0(colnames(s_po_df[scenarios[j,1]]),"_sum"), paste0(colnames(s_po_df[scenarios[j,2]]),"_sum"), "para", "budget")
    
    #rank PUs based on the results of the objective function for each weight
    for (p in 2:(length(para) + 1)){
      Scenario_ES<-Scenario_ES[order(-Scenario_ES[p]),]
      bud_spent<-Scenario_ES[1,"overall_npv"]
      sum_budget<-1
      
      while (bud_spent <= budget){
        bud_spent<-bud_spent + Scenario_ES[sum_budget + 1,"overall_npv"] 
        sum_budget<-sum_budget + 1
      }
      
      #select PUs until the target value (e.g. 30 %) is reached
      df_selected<-Scenario_ES[1:sum_budget,c(1,p,length(Scenario_ES))]
      #record area to reach target
      df_selected$budget<-budget
      #record weigth value
      df_selected$para<-para[p -1]
      colnames(df_selected)[2]<-"ES_scenario"
      #save data frame as object of list
      df_para[[p -1]]<-df_selected
      
      #select PUs ranked highest for obj1
      obj1_pu_df_selected <- merge(df_selected, s_po_df[c("ID",colnames(s_po_df[scenarios[j,3]]))], by= "ID")
      obj1_sum <- sum(obj1_pu_df_selected[length(obj1_pu_df_selected)])
      if(is.na(obj1_sum))obj1_sum<-0
      
      #select PUs ranked highest for obj2
      obj2_pu_df_selected <- merge(df_selected, s_po_df[c("ID",colnames(s_po_df[scenarios[j,4]]))], by="ID")
      obj2_sum <- sum(obj2_pu_df_selected[length(obj2_pu_df_selected)]) 
      if(is.na(obj2_sum))obj2_sum<-0
      # results to create pareto graph
      pareto_results[p-1,]<-c(obj1_sum, obj2_sum, para[p-1], budget)
    }
    
    #save dataframes with selected units for each objective
    #df_para[[length(df_para)]] means that the last parameter was used
    #in our case is 1 and it means that more weight was given to
    #objective 1 over objective 2
    #"_weight_last.csv" means that the second objective in the file name
    #has the highest weight, whereas "_weight_1_30.csv" means that the 
    #first scenario has the highest weight
    
    write.csv(pareto_results, paste0("calculations_ranking_analysis/50k/target_30/",
                                     "eco_",eco_list[[i]]$Region,"_",colnames(s_po_df[scenarios[j,1]]),
                                     "_",colnames(s_po_df[scenarios[j,2]]),
                                     "_para_30.csv"))
    
    write.csv(df_para[[length(df_para)]], paste0("calculations_ranking_analysis/50k/target_30/",
                                                 "eco_",eco_list[[i]]$Region,"_",colnames(s_po_df[scenarios[j,1]]),"_",colnames(s_po_df[scenarios[j,2]]),
                                                 "_weight_1_30.csv"))
    write.csv(df_para[[1]], paste0("calculations_ranking_analysis/50k/target_30/",
                                   "eco_",eco_list[[i]]$Region,"_",colnames(s_po_df[scenarios[j,1]]),"_",colnames(s_po_df[scenarios[j,2]]),
                                   "_weight_last.csv"))
  
    
  }
  
  #select PUs from each scenario and save them as shp
  
  g_po<-subset(grid_eco,grid_eco$ID %in% df_para[[length(df_para)]]$ID)
  g_po_2<-subset(grid_eco,grid_eco$ID %in% df_para[[1]]$ID)
  
  writeOGR(g_po_2,"./calculations_ranking_analysis/50k/target_30/corridors_shp_30",
           paste0("eco_",i,"_",colnames(s_po_df[scenarios[j,2]]),"_",
                  colnames(s_po_df[scenarios[j,1]]),"_30"),
           driver="ESRI Shapefile")
  writeOGR(g_po,"./calculations_ranking_analysis/50k/target_30/corridors_shp_30",
           paste0("eco_",i,"_",colnames(s_po_df[scenarios[j,1]]),"_",
                  colnames(s_po_df[scenarios[j,2]]),"_30"),
           driver="ESRI Shapefile")
  print(i)
}
