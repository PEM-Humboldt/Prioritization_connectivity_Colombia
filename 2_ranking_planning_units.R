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

setwd("C:/Data/50Reefs/Felipe/iucn")


#load ecorregions

#ecorregions_COL<-readOGR("./Regiones_bioticas/regiones_bioticos/Regiones_Bióticas.shp")

#prepare grid calculation data
grid_cals_f<-list.files("./NATGEO/grid_calculations_NATGEO",full.names = T)
grid_cals<-lapply(grid_cals_f,read.csv)
#alldata<-do.call(rbind.fill,grid_cals)
#head(alldata)

reg<-strsplit(grid_cals_f,"grid_")
reg<-lapply(reg,function(x)strsplit(x[[3]],"_calculations"))
reg<-lapply(reg,function(x)return(x[[1]][[1]]))

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
scenarios<-data.frame(objective_1= c("sp_126_ratio","sp_585_ratio",
                                     "sp_present_ratio"),      
                      objective_2= c("corridors_rent"),
                      valuesobj1 = c("sp_126","sp_585",
                                     "sp_present"),
                      valuesobj2 = c("corridors"))

# % to protect
target = 30 

for(i in 2:length(grid_cals)){
  
  ecorregion_df<-grid_cals[[i]]
  #calculate proportion area protected per ecorregion####
  
  eco_protected<-nrow(subset(ecorregion_df,!is.na(ecorregion_df$PNN)))
  eco_total<-nrow(ecorregion_df)
  prop_protected<-(eco_protected/eco_total)*100
  
  #read data of planning units####
  
  ecorregion_df$Region<-reg[[i]]
  
  s_po_df<-ecorregion_df
  
  s_po_df<-rename(s_po_df,"total_rent_300m" = "MEAN")
  
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
  
  s_po_df[which(is.na(s_po_df$corridors)),"corridors"]<-2e5
  s_po_df[which(s_po_df$corridors > 5e4),"corridors"]<-2e5
  
  #transform corridor values
  #transform to meters
  s_po_df$corridors<-s_po_df$corridors/1000
  s_po_df$corridors<- range01(s_po_df$corridors)
  s_po_df$corridors<-1-s_po_df$corridors
  
  #only consider forest corridors
  #s_po_df$corridors<-s_po_df$Bosque_col_300m*s_po_df$corridors
  
  #divide by rent
  s_po_df$corridors_rent<-s_po_df$corridors/s_po_df$total_rent_300m
  s_po_df[which(s_po_df$corridors == 2e5/1000),"corridors_rent"]<-min(s_po_df$corridors_rent, na.rm = T)
  
  
  #standardize ES and biodiversity values####
  #biod values
  s_po_df$sp_126_ratio<- range01(s_po_df$sp_126)
  s_po_df$sp_126_ratio<-s_po_df$sp_126/s_po_df$total_rent_300m
  
  s_po_df$sp_585_ratio<- range01(s_po_df$sp_585)
  s_po_df$sp_585_ratio<-s_po_df$sp_585/s_po_df$total_rent_300m
  
  s_po_df$sp_present_ratio<- range01(s_po_df$sp_present)
  s_po_df$sp_present_ratio<-s_po_df$sp_present/s_po_df$total_rent_300m
  

  #change NAs to 0
  s_po_df[which(is.na(s_po_df$corridors)),"corridors"]<-0
  s_po_df[which(is.na(s_po_df$sp_126)),"sp_126"]<-0
  s_po_df[which(is.na(s_po_df$sp_585)),"sp_585"]<-0
  s_po_df[which(is.na(s_po_df$sp_present)),"sp_present"]<-0
  
  #change NAs to 0
  s_po_df[which(is.na(s_po_df$corridors_rent)),"corridors_rent"]<-0
  s_po_df[which(is.na(s_po_df$sp_126_ratio)),"sp_126_ratio"]<-0
  s_po_df[which(is.na(s_po_df$sp_585_ratio)),"sp_585_ratio"]<-0
  s_po_df[which(is.na(s_po_df$sp_present_ratio)),"sp_present_ratio"]<-0
  
  #-------------------- Set up benefit-cost ratio targeting model -----------------
  #budget here represents the total area that needs to be protected to reach target
  budget <- (eco_total*((target - prop_protected)/100))#*0.3*0.3
  #print(budget)
  #print(prop_protected)
  
  s_po_df<-subset(s_po_df,is.na(s_po_df$PNN))
  
  #####
  
  grid_eco<-readOGR(paste0("NATGEO/grids_eco_NATGEO/grid_eco_",
                           s_po_df$Region[[1]],"_calc.shp"))

  #Apply objective function to all the scenario combinations
  systime1 <- Sys.time()
  conn_results<-list()
  for (j in 1:nrow(scenarios)) {
    #apply function to all dataframe rows using all the parameter values
    vector_t<- apply(s_po_df[c(scenarios[j,1],scenarios[j,2])], 1, function(x)obj_fun(para,x[scenarios[j,1]],x[scenarios[j,2]])) #apply obj function
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
    
    write.csv(Scenario_ES,paste0("NATGEO/ranking_analysis/",
                                 "eco_",s_po_df$Region[[1]],"_",colnames(s_po_df[scenarios[j,1]]),
                                 "_",colnames(s_po_df[scenarios[j,2]]),
                                 "_allpara.csv"))
    
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
    
    write.csv(pareto_results, paste0("NATGEO/ranking_analysis/",
                                     "eco_",s_po_df$Region[[1]],"_",colnames(s_po_df[scenarios[j,1]]),
                                     "_",colnames(s_po_df[scenarios[j,2]]),
                                     "_para_30.csv"))
    
    write.csv(df_para[[length(df_para)]], paste0("NATGEO/ranking_analysis/",
                                                 "eco_",s_po_df$Region[[1]],"_",colnames(s_po_df[scenarios[j,1]]),"_",colnames(s_po_df[scenarios[j,2]]),
                                                 "_weight_1_30.csv"))
    write.csv(df_para[[1]], paste0("NATGEO/ranking_analysis/",
                                   "eco_",s_po_df$Region[[1]],"_",colnames(s_po_df[scenarios[j,1]]),"_",colnames(s_po_df[scenarios[j,2]]),
                                   "_weight_last.csv"))
    
    grid_eco<-readOGR(paste0("NATGEO/grids_eco_NATGEO/grid_eco_",
                             s_po_df$Region[[1]],"_calc.shp"))
    
    
    g_po<-subset(grid_eco,grid_eco$ID %in% df_para[[length(df_para)]]$ID)
    g_po_2<-subset(grid_eco,grid_eco$ID %in% df_para[[1]]$ID)
    
    writeOGR(g_po_2,"NATGEO/corridors_shp_30_1",
             paste0("eco_",s_po_df$Region[[1]],"_",colnames(s_po_df[scenarios[j,2]]),"_",
                    colnames(s_po_df[scenarios[j,1]]),"_30"),
             driver="ESRI Shapefile")
    writeOGR(g_po,"NATGEO/corridors_shp_30_1",
             paste0("eco_",s_po_df$Region[[1]],"_",colnames(s_po_df[scenarios[j,1]]),"_",
                    colnames(s_po_df[scenarios[j,2]]),"_30"),
             driver="ESRI Shapefile")
    
  }
  
  #select PUs from each scenario and save them as shp
  
  print(i)
}
