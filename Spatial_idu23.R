library(sf)
library(spdep)
library(whitebox)
library(car)
ibrary(spatialreg)
library(gstat)

###############################
###Section 1 - Data preparation and variable calculation
###############################
#Upload and prepare data file
Dados_SpatialCorkoak<-read_excel("C:/Users/caven/Google Drive (paulofirmino@e-isa.ulisboa.pt)/2023-Doutoramento/SpatialCorkoak/Dados_SpatialCorkOak.xlsx")
Dados_SpatialCorkoak<-st_as_sf(Dados_SpatialCorkoak, coords=c("coordsX","coordsY"))
parcela<-Dados_SpatialCorkOak[which(Dados_SpatialCorkOak$Site == "A"),] ## Select on of the 4 properties (A, B, C, D) - A and B are Santarém, C and D are Castelo Branco 
parcela<-st_as_sf(parcela, coords=c("coordsX","coordsY"),crs = 3763) #transform to sf object
parcela<-parcela[which(parcela$mais_velha == 0),]  ##Remove clearly older trees present in the plantation area
parcela_vivas<-parcela[which(parcela$Morta == 0),]  ##Remove dead trees
parcela_vivas<-parcela_vivas[-which(parcela_vivas$du_annual_growth == 0),] #do not run for A, only site no du_annual_growth == 0
#"parcela_vivas" will be used when considering individual trees

#1.1.1 Geographical position indices calculation (TPI and TRI)
#f <- matrix(1, nrow=5, ncol=5)
#TPI <- focal(Altimetria, w=f, fun=function(x, ...) x[5] - mean(x[-5]), pad=TRUE, padValue=NA);plot(TPI)
#TRI <- focal(Altimetria, w=f, fun=function(x) sum(abs(x[-5]-x[5]))/8, pad=TRUE, padValue=NA);plot(TRI)

#1.1.2 Geographical position indices calculation (TWI)
# TWI calculation based on the methodology of https://vt-hydroinformatics.github.io/rgeoraster.html
# Prepare DEM for Hydrology Analyses
#  wbt_breach_depressions_least_cost(dem = "/path/parcela_Altimetria.tif", output = "/path/parcela_DTM_breach.tif", dist = 5,  fill = TRUE)
#  wbt_fill_depressions_wang_and_liu(dem = "/path/parcela_DTM_breach.tif",output = "/path/parcela_DTM_breach_fill.tif")
# Visualize and correct filled sinks and breached depression
#  filled_breached <- raster("/path/parcela_DTM_breach_fill.tif")
#  plot(filled_breached)
#    difference <- Altimetria - filled_breached
#  difference[difference == 0] <- NA
#D infinity flow accumulation (alternative flow accumulation may be calculated from D infinity method D8 Flow Accumulation)
#  wbt_d_inf_flow_accumulation("/path/parcela_DTM_breach_fill.tif","/path/Infinit_FlowAccum.tif")
#  dinf <- raster("/path/Infinit_FlowAccum.tif");  plot(dinf)
#Calculate Specific Contributing Area
#  wbt_d_inf_flow_accumulation(input = "/path/parcela_DTM_breach_fill.tif",output = "/path/FlowAccum2.tif", out_type = "Specific Contributing Area")
#Calculate slope or use slope mapped tif from data
#  wbt_slope(dem = "/path/parcela_DTM_breach_fill.tif",output = "/path/parcela_slope_degrees.tif",units = "degrees")
#Calculate topographic wetness index
#  wbt_wetness_index(sca = "/path/FlowAccum2.tif",slope = "/path/parcela_slope_degrees.tif", output = "/path/TWI.tif")              
#  twi <- raster("/path/TWI.tif");  plot(twi)

########
#1.2 Dependent variables - Calculate diameter related variables and attribute explanatory values to points
# Dependent variable is considered in the ways: 1) as the individual diameter annual growth; 2) as the basal area annual growth of a group of closely located trees
#Annual growth is considered so that all plots are normalized by age
########
#1.2.1 Individual diameter annual growth
#Already calculated in dataframe as column "du_annual_growth"
#Attribute exploratory variable value according to the pixel where tree is positioned
#parcela$Cea_1m_1px<-extract(Cea_1m,parcela) #Variable Example, already on the dataframe

###############################
##2. Spatial autocorrelation analysis
###############################
#Section 2 - Spatial autocorrelation analysis
#Spatial autocorrelation analysis is applied to each plot separately, with the following steps: 
#1) Spatial matrix definition;
#2) Spatial weights definition; 
#3) Global Moran's I analysis
#4) Correlogram (Moran's I statistics) analysis
#5) Local Moran's I analysis

########
##2.1 Spatial matrix definition
########             
#Although plantations are regular, trees are not at exact distances to make preferable the using of polygons for analysis.
#Its opted using points for the analysis. Output of 2.1 is nb object.
#Three methods are tested to choice the proper spatial matrix # Distance based neighbours; k neighbours; Tree area of influence sobreposition
#2.1.1 Distance base neighbours listw objects creation #To find the most adequate distances for the lag, we calculate the listw object from 2m to 15m 
#Simultaneously we test the most adequate weights 1) Row-normalized binary (rnorm_bin); Binary (bin); Row-normalized Inverted distance weights (rnorm_IDW); Inverted distance weights (IDW)
# empty vector for storage of Moran's I statistic values
moran_I_rnorm_bin<- c()
moran_I_bin <- c()
moran_I_rnorm_idw<-c()
moran_I_idw <- c()

#Cycle to run all Moran's I tests for each combination of lag/weight - slow 
for (d in seq(2, 15, 1)) {
  rnorm_bin <- nb2listw(dnearneigh(parcela_vivas, d1 = 0, d2 = d), style = "W", zero.policy = TRUE) #Row-normalized binary weights (style ="w")
  bin <- nb2listw(dnearneigh(parcela_vivas, d1 = 0, d2 = d), style = "B", zero.policy = TRUE)   #Binary weights (style ="B")
  #To use that idw weights, it requires nb2listwdist
  rnorm_idw <-nb2listwdist(dnearneigh(parcela_vivas, d1=0, d2=d),parcela_vivas,type="idw",style="W", zero.policy = TRUE) #Row-normalized idw (style ="w")
  idw <-nb2listwdist(nb2listw(dnearneigh(parcela_vivas, d1=0, d2=d),zero.policy = TRUE)$n,parcela_vivas,type="idw",zero.policy = TRUE) #idw (style ="Raw")
  #Calculate Moran's I statistics
  moran1 <- moran.mc(parcela_vivas$du_annual_growth, rnorm_bin, nsim = 999, zero.policy = TRUE)
  moran2 <- moran.mc(parcela_vivas$du_annual_growth, bin, nsim = 999, zero.policy = TRUE)
  moran3 <- moran.mc(parcela_vivas$du_annual_growth, rnorm_idw, nsim = 999, zero.policy = TRUE)
  moran4 <- moran.mc(parcela_vivas$du_annual_growth, idw, nsim = 999, zero.policy = TRUE)
  #Add the value to the empty vector
  moran_I_rnorm <- c(moran_I_rnorm, moran1$statistic)
  moran_I_bin <- c(moran_I_bin, moran2$statistic)
  moran_I_idw <- c(moran_I_idw, moran3$statistic)
  moran_I_rnorm_idw <- c(moran_I_rnorm_idw, moran4$statistic)
}
#Finalize the dataframe for plotting
moran_I_Adist <-  data.frame(moran1 =  moran_I_rnorm,moran2 =  moran_I_bin,moran3 =  moran_I_rnorm_idw,moran4 = moran_I_idw, distance = seq(2, 15, 1))
#Lock the dataframe to into a dataframe for each site, to not have to run again
#moran_I_Adist_4w <- data.frame(moran1 =  moran_I_rnorm,moran2 =  moran_I_bin,moran3 =  moran_I_rnorm_idw,moran4 = moran_I_idw, distance = seq(2, 15, 1))
#moran_I_Bdist_4w <- data.frame(moran1 =  moran_I_rnorm,moran2 =  moran_I_bin,moran3 =  moran_I_rnorm_idw,moran4 = moran_I_idw, distance = seq(2, 15, 1))
#moran_I_Cdist_4w <- data.frame(moran1 =  moran_I_rnorm,moran2 =  moran_I_bin,moran3 =  moran_I_rnorm_idw,moran4 = moran_I_idw, distance = seq(2, 15, 1))
#moran_I_Ddist_4w <- data.frame(moran1 =  moran_I_rnorm,moran2 =  moran_I_bin,moran3 =  moran_I_rnorm_idw,moran4 = moran_I_idw, distance = seq(2, 15, 1))
#moran_I_abcddist_4w <- data.frame(moran1 =  moran_I_rnorm,moran2 =  moran_I_bin,moran3 =  moran_I_rnorm_idw,moran4 = moran_I_idw, distance = seq(2, 15, 1))
#Plot and compare lags/weights for and betwee each site
plot<-(ggplot() + 
  geom_point() +
    geom_line(data = moran_I_abcddist_4w, aes(x = distance, y = moran1, color = "red")) +
    geom_line(data = moran_I_abcddist_4w, aes(x = distance, y = moran2, color = "black")) +
    geom_line(data = moran_I_abcddist_4w, aes(x = distance, y = moran3, color = "blue")) +
    geom_line(data = moran_I_abcddist_4w, aes(x = distance, y = moran4, color = "green")) +
    scale_color_identity(name = "Weights",
                         breaks = c("red", "black", "blue","green"),
                         labels = c("Row.norm Binary","Binary","Row.norm IDW","IDW"),
                         guide = "legend")+
  ylab("Moran's I statistic")+
  ggtitle("All")+theme(plot.title = element_text(hjust = 0.5),legend.justification = "left")+
  scale_y_continuous(
    limits = c(0, 0.60),    # Set the y-axis range
    breaks = seq(0, 0.60, 0.05))+    # Set custom breaks
    scale_x_continuous(
      limits = c(0, 15),    # Set the y-axis range
      breaks = seq(0, 15, 2),    # Set custom breaks
 )
 )
#lock plots
#plotd_a<-plot_a
#plotd_b<-plot_b
#plotd_c<-plot_c
#plotd_d<-plot_d
#plotd_abcd<-plot_abcd
grid.arrange(plot_a, plot_b, plot_c, plot_d, nrow = 2, ncol=2 )#Plot the four sites 
#######
#Evaluation:
#Best distance-base lags is between 5-7m for any site
#Row-normalized weights have lower Moran I values, but ensure spatial weights statistics assumptions to correctly run models
#Row-normalized binary or row-normalized weights are most adequate depending on the site

#2.1.2 k-neighbours base listw objects creation #To find the most adequate number of neighbourslag, we calculate the listw object for k=4, k=8 (and k=12 to see the evolution)
#Run similar cycle for k-neighbours now
# empty vector for storage of Moran's I values
moran_I_rnorm_bin<- c()
moran_I_bin <- c()
moran_I_rnorm_idw<-c()
moran_I_idw <- c()
for (k in seq(4, 12, 4)) {
  rnorm_bin <- nb2listw(knn2nb(knearneigh(parcela_vivas, k=k)), style = "W", zero.policy = TRUE)
  bin <- nb2listw(knn2nb(knearneigh(parcela_vivas, k=k)), style = "B", zero.policy = TRUE)
  rnorm_IDW <-nb2listwdist(knn2nb(knearneigh(parcela_vivas, k=k)),parcela_vivas,type="idw",style="W",zero.policy = TRUE)
  IDW <-nb2listwdist(knn2nb(knearneigh(parcela_vivas, k=k)),parcela_vivas,type="idw",zero.policy = TRUE)
  moran1 <- moran.mc(parcela_vivas$du_annual_growth, rnorm_bin, nsim = 999, zero.policy = TRUE)
  moran2 <- moran.mc(parcela_vivas$du_annual_growth, bin, nsim = 999, zero.policy = TRUE)
  moran3 <- moran.mc(parcela_vivas$du_annual_growth, rnorm_idw, nsim = 999, zero.policy = TRUE)
  moran4 <- moran.mc(parcela_vivas$du_annual_growth, idw, nsim = 999, zero.policy = TRUE)
  moran_I_rnorm_bin <- c(rnorm_bin, moran1$statistic)
  moran_I_bin <- c(moran_I_bin, moran2$statistic)
  moran_I_rnorm_idw <- c(moran_I_rnorm_idw, moran3$statistic)
  moran_I_idw <- c(moran_I_idw, moran4$statistic)
  }
#Lock the values
#moran_I_Akneigh_4w <- data.frame(moran1=moran_I_rnorm,moran2=moran_I_bin,moran3=moran_I_rnorm_idw,moran4=moran_I_idw,neighbours=seq(4, 12, 4))
#moran_I_Bkneigh_4w <- data.frame(moran1=moran_I_rnorm,moran2=moran_I_bin,moran3=moran_I_rnorm_idw,moran4=moran_I_idw,neighbours=seq(4, 12, 4))
#moran_I_Ckneigh_4w <- data.frame(moran1=moran_I_rnorm,moran2=moran_I_bin,moran3=moran_I_rnorm_idw,moran4=moran_I_idw,neighbours=seq(4, 12, 4))
#moran_I_Dkneigh_4w <- data.frame(moran1=moran_I_rnorm,moran2=moran_I_bin,moran3=moran_I_rnorm_idw,moran4=moran_I_idw,neighbours=seq(4, 12, 4))
#moran_I_abcdkneigh_4w <- data.frame(moran1=moran_I_rnorm,moran2=moran_I_bin,moran3=moran_I_rnorm_idw,moran4=moran_I_idw,neighbours=seq(4, 12, 4))
plot<-(ggplot() + 
           geom_point() +
           geom_line(data = moran_I_abcdkneigh_4w, aes(x = neighbours, y = moran1, color = "red")) +
           geom_line(data = moran_I_abcdkneigh_4w, aes(x = neighbours, y = moran2, color = "black")) +
           geom_line(data = moran_I_abcdkneigh_4w, aes(x = neighbours, y = moran3, color = "blue")) +
           geom_line(data = moran_I_abcdkneigh_4w, aes(x = neighbours, y = moran4, color = "green")) +
           scale_color_identity(name = "Weights",
                                breaks = c("red", "black", "blue","green"),
                                labels = c("Row.norm Binary","Binary","Row.norm IDW","IDW"),
                                guide = "legend")+
           ylab("Moran's I statistic")+
           ggtitle("All")+theme(plot.title = element_text(hjust = 0.5),legend.justification = "left")+
           scale_y_continuous(
             limits = c(0, 0.60),    # Set the y-axis range
             breaks = seq(0, 0.60, 0.02),    # Set custom breaks
           )
)
#plotk_a<-plot
#plotk_b<-plot
#plotk_c<-plot
#plotk_d<-plot
#plotd_abcd<-plot
grid.arrange(plot_a, plot_b, plot_c, plot_d, nrow = 2, ncol=2 )#Plot the four sites 
#######
#Evaluation:
#K-neighbours-based lags is more adequate for k=4 for any site
#Values can be close to the distance-base lags. A table has to me create for a direct comparison.













#dnearneigh(parcela_vivas, d1=0, d2=8)
#dnearneigh(parcela_vivas, d1=0, d2=10)
#dnearneigh(parcela_vivas, d1=0, d2=15)
#2.1.2 k neighbours - Selected 4 and 8 neighbours.
#knn2nb(knearneigh(parcela_vivas, k=4))
#knn2nb(knearneigh(parcela_vivas, k=8))
#2.1.3  Tree area of influence sobreposition
#Requires 3 steps: 
#1)Calculating individual tree area of influence;
#2)Creating of neighbour/distance table; 
#3)Identifying the trees with sobreposed areas of influence and providing a value for them as neighbours

#2.1.3.1  Calculating individual tree area of influence, according Paulo et al. 2016
parcela_vivas$du_dug<-""
parcela_vivas$du_dug<-parcela_vivas$du/(sqrt(sum((parcela_vivas$du)^2)/length(parcela_vivas$du)))
parcela_vivas$du_dug<-as.numeric(parcela_vivas$du_dug)
parcela_vivas$influence_area_est<-2.5*(28.502*exp((-66.436+4.201*(parcela_vivas$du_dug))/(19.817+parcela_vivas$du))/2) 
plot(parcela_vivas["influence_area_est"],pch=16)

#2.1.3.2 Creating of neighbour/distance table; 
parcela_vivas$continId<-c(1:dim(parcela_vivas)[1])#remake the continId for parcela_vivas, after removing dead elements
xy<-st_coordinates(parcela_vivas)
matrix<-as.matrix(dist(xy, "euclidean"), labels=TRUE) #Distance matrix
colnames(matrix) <- rownames(matrix) <- parcela_vivas$continId #Name rows and columns
melt <- melt(matrix) #unmaked the matrix to columns
M101<-data.frame(parcela_vivas)
melt_0<-melt[which(melt$value!=0),]  #removes the distances=zero
submelt01<-melt[melt_0$Var1 %in% M101$continId & melt_0$Var2 %in% M101$continId,] #creates a table with Var1, Var2 and the distance value between then
submelt01<-submelt01[which(submelt01$value<33),] #Removes any lines with distace > 33m, which is the highest distance between two neighbouring trees. This condition fastens the process.
submelt01<-submelt01[which(submelt01$value!=0),] #Removes distances =0
mt<-merge(submelt01, M101[,c("continId","influence_area_est" )], by.x="Var1", by.y="continId") #Creates a diameter column for Var1, depending on tree number
mtt<-merge(mt, M101[,c("continId","influence_area_est" )], by.x="Var2", by.y="continId") #Creates a diameter column for Var2, depending on tree number
names(mtt)<-c(names(mtt)[1:3], "influence_area_est2",  "influence_area_est1") #change the names of the two new columns
head(mtt)

#2.1.3.3 Identifying the trees with sobreposed areas of influence and providing a value for them as neighbours
matrix<-matrix(nrow=dim(parcela_vivas)[1], ncol=dim(parcela_vivas)[1]) #starts a blank matrix for the cycle
for (k in c(1:nrow(mtt))){ 
  valor1<-mtt[k,1]
  valor2<-mtt[k,2]
  linha<- mtt[k,]
  matrix[valor1,valor2]<- ifelse((linha$influence_area_est1 + linha$influence_area_est2 > linha$value),1,0)
}
#asks if the sum of radius of area of influence of two trees is higher than their distance. 
#If true, they are neighbours. The value attributed if neighbours can work as a weight. In this case is 1, similar to binary case.

matrix[is.na(matrix)] = 0 #turns NA to zeros
colnames(matrix) <- rownames(matrix) <- parcela_vivas$continId #give names to new matrix
#View(matrix) #check if result is correct
W<- mat2listw(matrix) #creates a listw object from this matrix, ready to be applied, with style=matrix itself.
WW<-mat2listw(matrix,style = "W")

###Repeat the cycle about attribute the distance between instead of 1. 
matrix<-matrix(nrow=dim(parcela_vivas)[1], ncol=dim(parcela_vivas)[1]) #starts a blank matrix for the cycle
for (k in c(1:nrow(mtt))){ 
  valor1<-mtt[k,1]
  valor2<-mtt[k,2]
  linha<- mtt[k,]
  matrix[valor1,valor2]<- ifelse((linha$influence_area_est1 + linha$influence_area_est2 > linha$value),linha$value,0)
}
matrix[is.na(matrix)] = 0 #turns NA to zeros
colnames(matrix) <- rownames(matrix) <- parcela_vivas$continId #give names to new matrix
#View(matrix) #check if result is correct
Wdist<- mat2listw(matrix) #creates a listw object from this matrix, this time with dist as weights
WdistW<- mat2listw(matrix,style="W")

moran.test(parcela_vivas$du_annual_growth, W, zero.policy = TRUE)
moran.test(parcela_vivas$du_annual_growth, WW, zero.policy = TRUE)
moran.test(parcela_vivas$du_annual_growth, Wdist, zero.policy = TRUE)
moran.test(parcela_vivas$du_annual_growth, WdistW, zero.policy = TRUE)
####2.1.4  Plot and compare results
par(mfrow=c(1,1));par(cex=0.6, pch=16)
plot(dnearneigh(parcela_vivas, d1=0, d2=5), coord=(parcela_vivas$geometry), col="red");title(main=paste("Distance based neighbours - 5m"),cex.main=2)
plot(dnearneigh(parcela_vivas, d1=0, d2=6), coord=(parcela_vivas$geometry), col="red");title(main=paste("Distance based neighbours - 6m"),cex.main=2)
plot(dnearneigh(parcela_vivas, d1=0, d2=7), coord=(parcela_vivas$geometry), col="red");title(main=paste("Distance based neighbours - 7m"),cex.main=2)
plot(knn2nb(knearneigh(parcela_vivas, k=4)), coord=(parcela_vivas$geometry), col="red");title(main=paste("k neighbours - 4"),cex.main=2)
plot(knn2nb(knearneigh(parcela_vivas, k=8)), coord=(parcela_vivas$geometry), col="red");title(main=paste("k neighbours - 8"),cex.main=2)
plot(W,coord=(parcela_vivas$geometry),col="red"); title(main=paste("Area of Influence approach"),cex.main=2)

######
#Evaluation: 
#8m distance and k=4 neighbours seem to underestimate the amount of links per tree
#15m distance seems to overestimate the amount of links per tree    
#10 distance and k=8 seem to have a most realistic number and distribution of neighbours
#Area of Influence approach seems to work better on site A and B, where tree size is more regular
#Area of Influence approach may overestimate the neighbours in larger trees
#Area of Influence approach underestimate spatial autocorrelation in very small trees (very few or no neighbours) in C and D
########

##2.1.2 Spatial weights matrix definition
########    
##Three methods are tested to choose weights: 1) row-normalized (W); 2) Binary (B); Inverse of distance (idw)
#The three methods generate distinct spatial weight matrices, an object class listw
#Example of a previously defined nb object
#nb2listw(dnearneigh(parcela_vivas, d1=0, d2=8),zero.policy = TRUE, style="W") #row-normalized
#nb2listw(dnearneigh(parcela_vivas, d1=0, d2=8),zero.policy = TRUE,style="B") #binary
#nb2listwdist(nb2listw(dnearneigh(parcela_vivas, d1=0, d2=8),zero.policy = TRUE)) #idw        

#To select the most suitable weights, they are compared on the following spatial analysis functions

########
##2.2 Global Moran'I test
########  
# A table can be created fro each plot to compare the multiple methods outputs: Morans'I statistic and respective p-value.
#Example of the output table, plot A:
#                         Weights    Row-normalized          Inverse Distance                    Binary        
#                          Method Moran I statistic p-value Moran I statistic p-value Moran I statistic p-value
#  Distance based neighbours - 8m              0.16       0             0.187       0             0.145       0
# Distance based neighbours - 10m             0.129       0             0.154       0              0.12       0
# Distance based neighbours - 15m             0.129       0             0.136       0             0.112       0
#                k neighbours - 4             0.158       0             0.193       0             0.158       0
#                k neighbours - 8             0.127       0             0.156       0             0.127       0
#      Area of Influence approach              0.07   4e-05             0.102       0             0.103       0

tabela<-data.frame()
tabela[1,1]<- "Method"
tabela[2,1]<- "Distance based neighbours - 5m"
tabela[3,1]<- "Distance based neighbours - 6m"
tabela[4,1]<- "Distance based neighbours - 7m"
tabela[5,1]<- " k neighbours - 4"
tabela[6,1]<- " k neighbours - 8"
tabela[7,1]<- "Area of Influence approach"
tabela[1,2]<-"Moran I statistic";tabela[1,3]<-"p-value";tabela[1,4]<-"Moran I statistic";tabela[1,5]<-"p-value";tabela[1,6]<-"Moran I statistic";tabela[1,7]<-"p-value";
#distance base tests - row-normalized weights; zero.policy=TRUE must be present since there is always points with no (living) neighbours
tabela[2,2]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(dnearneigh(parcela_vivas, d1=0, d2=5),zero.policy = TRUE),zero.policy=TRUE)$estimate[1],3)
tabela[2,3]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(dnearneigh(parcela_vivas, d1=0, d2=5),zero.policy = TRUE),zero.policy=TRUE)$p.value,5)
tabela[3,2]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(dnearneigh(parcela_vivas, d1=0, d2=6),zero.policy = TRUE),zero.policy=TRUE)$estimate[1],3)
tabela[3,3]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(dnearneigh(parcela_vivas, d1=0, d2=6),zero.policy = TRUE),zero.policy=TRUE)$p.value,5)
tabela[4,2]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(dnearneigh(parcela_vivas, d1=0, d2=7),zero.policy = TRUE),zero.policy=TRUE)$estimate[1],3)
tabela[4,3]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(dnearneigh(parcela_vivas, d1=0, d2=7),zero.policy = TRUE),zero.policy=TRUE)$p.value,5)
#distance base tests - inverse of distance weights; can be used nb2listwdist function
tabela[2,4]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listwdist(nb2listw(dnearneigh(parcela_vivas, d1=0, d2=5),zero.policy = TRUE)$n,parcela_vivas,type="idw",zero.policy = TRUE),zero.policy=TRUE)$estimate[1],3)
tabela[2,5]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listwdist(nb2listw(dnearneigh(parcela_vivas, d1=0, d2=5),zero.policy = TRUE)$n,parcela_vivas,type="idw",zero.policy = TRUE),zero.policy=TRUE)$p.value,5)
tabela[3,4]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listwdist(nb2listw(dnearneigh(parcela_vivas, d1=0, d2=6),zero.policy = TRUE)$n,parcela_vivas,type="idw",zero.policy = TRUE),zero.policy=TRUE)$estimate[1],3)
tabela[3,5]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listwdist(nb2listw(dnearneigh(parcela_vivas, d1=0, d2=6),zero.policy = TRUE)$n,parcela_vivas,type="idw",zero.policy = TRUE),zero.policy=TRUE)$p.value,5)
tabela[4,4]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listwdist(nb2listw(dnearneigh(parcela_vivas, d1=0, d2=7),zero.policy = TRUE)$n,parcela_vivas,type="idw",zero.policy = TRUE),zero.policy=TRUE)$estimate[1],3)
tabela[4,5]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listwdist(nb2listw(dnearneigh(parcela_vivas, d1=0, d2=7),zero.policy = TRUE)$n,parcela_vivas,type="idw",zero.policy = TRUE),zero.policy=TRUE)$p.value,5)
#distance base tests - binary weights
tabela[2,6]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(dnearneigh(parcela_vivas, d1=0, d2=5),zero.policy = TRUE,style="B"),zero.policy=TRUE)$estimate[1],3)
tabela[2,7]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(dnearneigh(parcela_vivas, d1=0, d2=5),zero.policy = TRUE,style="B"),zero.policy=TRUE)$p.value,5)
tabela[3,6]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(dnearneigh(parcela_vivas, d1=0, d2=6),zero.policy = TRUE,style="B"),zero.policy=TRUE)$estimate[1],3)
tabela[3,7]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(dnearneigh(parcela_vivas, d1=0, d2=6),zero.policy = TRUE,style="B"),zero.policy=TRUE)$p.value,5)
tabela[4,6]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(dnearneigh(parcela_vivas, d1=0, d2=7),zero.policy = TRUE,style="B"),zero.policy=TRUE)$estimate[1],3)
tabela[4,7]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(dnearneigh(parcela_vivas, d1=0, d2=7),zero.policy = TRUE,style="B"),zero.policy=TRUE)$p.value,5)
#k neighbours - row-normalized weights
tabela[5,2]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(knn2nb(knearneigh(parcela_vivas, k=4)),style="W"),zero.policy=TRUE)$estimate[1],3)
tabela[5,3]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(knn2nb(knearneigh(parcela_vivas, k=4)),style="W"),zero.policy=TRUE)$p.value,5)
tabela[6,2]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(knn2nb(knearneigh(parcela_vivas, k=8)),style="W"),zero.policy=TRUE)$estimate[1],3)
tabela[6,3]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(knn2nb(knearneigh(parcela_vivas, k=8)),style="W"),zero.policy=TRUE)$p.value,5)
#k neighbours - inverse of distance weights 
tabela[5,4]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listwdist(nb2listw(knn2nb(knearneigh(parcela_vivas, k=4)))$n,parcela_vivas,type="idw"),zero.policy=TRUE)$estimate[1],3)
tabela[5,5]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listwdist(nb2listw(knn2nb(knearneigh(parcela_vivas, k=4)))$n,parcela_vivas,type="idw"),zero.policy=TRUE)$p.value,5)
tabela[6,4]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listwdist(nb2listw(knn2nb(knearneigh(parcela_vivas, k=8)))$n,parcela_vivas,type="idw"),zero.policy=TRUE)$estimate[1],3)
tabela[6,5]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listwdist(nb2listw(knn2nb(knearneigh(parcela_vivas, k=8)))$n,parcela_vivas,type="idw"),zero.policy=TRUE)$p.value,5)
#k neighbours - binary weights
tabela[5,6]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(knn2nb(knearneigh(parcela_vivas, k=4)),style="B"),zero.policy=TRUE)$estimate[1],3)
tabela[5,7]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(knn2nb(knearneigh(parcela_vivas, k=4)),style="B"),zero.policy=TRUE)$p.value,5)
tabela[6,6]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(knn2nb(knearneigh(parcela_vivas, k=8)),style="B"),zero.policy=TRUE)$estimate[1],3)
tabela[6,7]<-round(moran.test(parcela_vivas$du_annual_growth,nb2listw(knn2nb(knearneigh(parcela_vivas, k=8)),style="B"),zero.policy=TRUE)$p.value,5)
#Area of influence neighbours - three weights (WW, Wdist, W)
tabela[7,2]<-round(moran.test(parcela_vivas$du_annual_growth,WW, zero.policy=TRUE)$estimate[1],3)
tabela[7,3]<-round(moran.test(parcela_vivas$du_annual_growth,WW, zero.policy=TRUE)$p.value,5)
tabela[7,4]<-round(moran.test(parcela_vivas$du_annual_growth,Wdist, zero.policy=TRUE)$estimate[1],3)
tabela[7,5]<-round(moran.test(parcela_vivas$du_annual_growth,Wdist, zero.policy=TRUE)$p.value,5)
tabela[7,6]<-round(moran.test(parcela_vivas$du_annual_growth,W, zero.policy=TRUE)$estimate[1],3)
tabela[7,7]<-round(moran.test(parcela_vivas$du_annual_growth,W, zero.policy=TRUE)$p.value,5)
nomes<-c("Weights","Row-normalized","","Inverse Distance","", "Binary","");colnames(tabela)<-nomes
#tabelaA<-tabela #save the table
tabelaA

listw_A<-nb2listw(dnearneigh(parcela_vivas, d1=0, d2=6),zero.policy = TRUE,style="W")
listw_B<-nb2listwdist(nb2listw(dnearneigh(parcela_vivas, d1=0, d2=7),zero.policy = TRUE)$n,parcela_vivas,type="idw",style="W",zero.policy = TRUE)
listw_C<-nb2listwdist(nb2listw(knn2nb(knearneigh(parcela_vivas, k=4)))$n,parcela_vivas,type="idw",style="W")
listw_D<-nb2listwdist(nb2listw(knn2nb(knearneigh(parcela_vivas, k=4)))$n,parcela_vivas,type="idw",style="W") 
listw_all<-nb2listwdist(nb2listw(knn2nb(knearneigh(parcela_vivas, k=4)))$n,parcela_vivas,type="idw",style="W") 
  
moran.test(parcela_vivas$du_annual_growth,listw_A, zero.policy=TRUE)
moran.test(parcela_vivas$du_annual_growth,listw_B, zero.policy=TRUE)
moran.test(parcela_vivas$du_annual_growth,listw_C, zero.policy=TRUE)
moran.test(parcela_vivas$du_annual_growth,listw_D, zero.policy=TRUE)
moran.test(parcela_vivas$du_annual_growth,listw_all, zero.policy=TRUE)

teste<-dnearneigh(parcela_vivas, d1=0, d2=7)
summary(teste)
teste<-knearneigh(parcela_vivas, k=4)
teste$nn

histA_class<-c(0,1,2,3,4,5)
histA_freqabs<-c(64,251,458,47,9,2)
histB_class<-c(0,1,2,3,4,5)
histB_freqabs<-c(108,356,554,93,17,3)
#histC_class<-c(0,1,2,3,4,5,6,7,8,9,10,11,12)
#histC_freqabs<-c(6,25,93,134,157,159,157,151,111,80,53,25,7)

bardplotD<-barplot(histD_freqabs, xlab="neighbours", names.arg = histD_class)
bardplotA<-barplot(histA_freqabs, xlab="neighbours", names.arg = histA_class)
bardplotB<-barplot(histB_freqabs, xlab="neighbours", names.arg = histB_class)
bardplotC<-barplot(histC_freqabs, xlab="neighbours", names.arg = histC_class)
par(mfrow = c(2, 2))
bardplotA<-barplot(histA_freqabs, xlab="number of neighbours", names.arg = histD_class, main="A")
bardplotB<-barplot(histB_freqabs, xlab="number of neighbours", names.arg = histD_class, main="B")
bardplotC<-barplot(histC_freqabs, xlab="number of neighbours", names.arg = histC_class, main="C")
bardplotD<-barplot(histD_freqabs, xlab="number of neighbours", names.arg = histD_class, main="D")

teste<-as.numeric(teste)
summary(teste)
hist(teste )
#######
#Evaluation:
#Generally similar results on Moran I statistics on distance based and k-neighbours spatial matrices. 
#Area of influence method with distinct results, and depends on the plot
#Row-normalized and Binary weights produce same results on k-neighbours methods, as expected
#idw weights show consistent slightly higher Moran's I statistics, except for site C (?)
#8m distance base and k4 neighbours show slightly higher Moran's I statistics for any weights


########
##2.3 Correlogram
########
par(mfrow=c(2,3))
plot(sp.correlogram(knn2nb(knearneigh(parcela_vivas, k=4)), parcela_vivas$du_annual_growth, method="I" ,order=9,zero.policy=TRUE),ylim=c(0,0.35), main="k4 neighbours")
plot(sp.correlogram(knn2nb(knearneigh(parcela_vivas, k=8)), parcela_vivas$du_annual_growth, method="I" ,order=9,zero.policy=TRUE),ylim=c(0,0.35), main="k8 neighbours")
plot(sp.correlogram(dnearneigh(parcela_vivas, d1=0, d2=8), parcela_vivas$du_annual_growth, method="I" ,order=9,zero.policy=TRUE),ylim=c(0,0.35), main="dist 8m")
plot(sp.correlogram(dnearneigh(parcela_vivas, d1=0, d2=10), parcela_vivas$du_annual_growth, method="I" ,order=9,zero.policy=TRUE),ylim=c(0,0.35), main="dist 10m")
plot(sp.correlogram(dnearneigh(parcela_vivas, d1=0, d2=15), parcela_vivas$du_annual_growth, method="I" ,order=9,zero.policy=TRUE),ylim=c(0,0.35) ,main="dist 15m")
plot(sp.correlogram(W$neighbours, parcela_vivas$du_annual_growth, method="I" ,order=9,zero.policy=TRUE),ylim=c(0,0.35) ,main="Area of Influence")
par(mfrow=c(1,1))
correl_a<-plot(sp.correlogram(dnearneigh(parcela_vivas, d1=0, d2=6), parcela_vivas$du_annual_growth, method="I" ,order=9,zero.policy=TRUE),ylim=c(0,0.40) ,main="A")
correl_b<-plot(sp.correlogram(dnearneigh(parcela_vivas, d1=0, d2=7), parcela_vivas$du_annual_growth, method="I" ,order=9,zero.policy=TRUE),ylim=c(0,0.40) ,main="B")
correl_C<-plot(sp.correlogram(dnearneigh(parcela_vivas, d1=0, d2=7), parcela_vivas$du_annual_growth, method="I" ,order=30,zero.policy=TRUE),ylim=c(0,0.40) ,main="C")
correl_d<-plot(sp.correlogram(dnearneigh(parcela_vivas, d1=0, d2=5), parcela_vivas$du_annual_growth, method="I" ,order=9,zero.policy=TRUE),ylim=c(0,0.40) ,main="D")

#######
#Evaluation:
#Lower k neighbours or distance show highest error bars, values tend to stay similar in lags
#K-neighbours and distance approach tend to show similar results
#Area of influence approach produce distinct results. Its more difficult to explain the lag system with this nb object
#Checking Area of influence approach correlogram $cardnos argument show very variable number of neighbours
par(mfrow=c(1,1),cex=0.8)
teste<-(moran.plot(parcela_vivas$du_annual_growth,listw=nb2listw(dnearneigh(parcela_vivas, d1=0, d2=8),zero.policy = TRUE), main="dist 8m"))
teste<-moran.plot(parcela_vivas$du_annual_growth,listw=nb2listw(dnearneigh(parcela_vivas, d1=0, d2=10),zero.policy = TRUE),  main="dist 10m")
teste<-moran.plot(parcela_vivas$du_annual_growth,listw=nb2listw(dnearneigh(parcela_vivas, d1=0, d2=15),zero.policy = TRUE), main="dist 15m")
teste<-moran.plot(parcela_vivas$du_annual_growth,listw=nb2listw(knn2nb(knearneigh(parcela_vivas, k=4)),zero.policy = TRUE), main="k4 neighbours")
teste<-moran.plot(parcela_vivas$du_annual_growth,listw=nb2listw(knn2nb(knearneigh(parcela_vivas, k=8)),zero.policy = TRUE), main="k8 neighbours")
moran.plot(parcela_vivas$du_annual_growth,listw=WW,zero.policy = TRUE, main="Area of influence approach")
mtext("Site D", side = 3, line = - 1.5, outer = TRUE)

print (inf.rows <- which(rowSums(teste$is.inf) == TRUE))

##########################################################################
# empty vector for storage of Moran's I values
moran_I_rnorm <- c()
moran_I_bin <- c()
moran_I_idw <- c()
moran_I_rnorm_idw<-c()
# loop d through a sequence ranging from 50 to 2000
for (d in seq(2, 15, 1)) {
  Rownorm <- nb2listw(dnearneigh(parcela_vivas, d1 = 0, d2 = d), style = "W", zero.policy = TRUE)
  Binary <- nb2listw(dnearneigh(parcela_vivas, d1 = 0, d2 = d), style = "B", zero.policy = TRUE)
  IDW <-nb2listwdist(nb2listw(dnearneigh(parcela_vivas, d1=0, d2=d),zero.policy = TRUE)$n,parcela_vivas,type="idw",zero.policy = TRUE)
  rnorm_IDW <-nb2listwdist(nb2listw(dnearneigh(parcela_vivas, d1=0, d2=d),zero.policy = TRUE)$n,parcela_vivas,type="idw",style="W", zero.policy = TRUE)
   moran1 <- moran.mc(parcela_vivas$du_annual_growth, Rownorm, nsim = 999, zero.policy = TRUE)
  moran2 <- moran.mc(parcela_vivas$du_annual_growth, Binary, nsim = 999, zero.policy = TRUE)
  moran3 <- moran.mc(parcela_vivas$du_annual_growth, IDW, nsim = 999, zero.policy = TRUE)
  moran4 <- moran.mc(parcela_vivas$du_annual_growth, rnorm_IDW, nsim = 999, zero.policy = TRUE)
  moran_I_rnorm <- c(moran_I_rnorm, moran1$statistic)
  moran_I_bin <- c(moran_I_bin, moran2$statistic)
  moran_I_idw <- c(moran_I_idw, moran3$statistic)
  moran_I_rnorm_idw <- c(moran_I_rnorm_idw, moran4$statistic)
}

#moran_I_Adist <- data.frame(moran1 =  moran_I_rnorm,moran2 =  moran_I_bin,moran3 =  moran_I_idw, distance = seq(2, 15, 1))
#moran_I_Adist_4w <- data.frame(moran1 =  moran_I_rnorm,moran2 =  moran_I_bin,moran3 =  moran_I_idw,moran4 = moran_I_rnorm_idw, distance = seq(2, 15, 1))
#moran_I_Bdist_4w <- data.frame(moran1 =  moran_I_rnorm,moran2 =  moran_I_bin,moran3 =  moran_I_idw,moran4 = moran_I_rnorm_idw, distance = seq(2, 15, 1))
#moran_I_Cdist_4w <- data.frame(moran1 =  moran_I_rnorm,moran2 =  moran_I_bin,moran3 =  moran_I_idw,moran4 = moran_I_rnorm_idw, distance = seq(2, 15, 1))
#moran_I_Ddist_4w <- data.frame(moran1 =  moran_I_rnorm,moran2 =  moran_I_bin,moran3 =  moran_I_idw,moran4 = moran_I_rnorm_idw, distance = seq(2, 15, 1))
#moran_I_alldist_4w <- data.frame(moran1 =  moran_I_rnorm,moran2 =  moran_I_bin,moran3 =  moran_I_idw,moran4 = moran_I_rnorm_idw, distance = seq(2, 15, 1))


#moran_I_Bdist <- data.frame(moran1 =  moran_I_rnorm,moran2 =  moran_I_bin,moran3 =  moran_I_idw, distance = seq(2, 15, 1))
#moran_I_Cdist <- data.frame(moran1 =  moran_I_rnorm,moran2 =  moran_I_bin,moran3 =  moran_I_idw, distance = seq(2, 15, 1))
#moran_I_Ddist <- data.frame(moran1 =  moran_I_rnorm,moran2 =  moran_I_bin,moran3 =  moran_I_idw, distance = seq(2, 15, 1))

plot_all<-(ggplot() + 
  geom_point() +
    geom_line(data = moran_I_alldist_4w, aes(x = distance, y = moran2, color = "black")) +
    geom_line(data = moran_I_alldist_4w, aes(x = distance, y = moran3, color = "blue")) +
    geom_line(data = moran_I_alldist_4w, aes(x = distance, y = moran1, color = "red")) +
    geom_line(data = moran_I_alldist_4w, aes(x = distance, y = moran4, color = "green")) +
    scale_color_identity(name = "Weights",
                         breaks = c("black", "blue", "red","green"),
                         labels = c("Binary", "IDW","Row.norm Binary","Row.norm IDW"),
                         guide = "legend")+
  ylab("Moran's I statistic")+
  ggtitle("All")+theme(plot.title = element_text(hjust = 0.5),legend.justification = "left")+
  scale_y_continuous(
    limits = c(0, 0.60),    # Set the y-axis range
    breaks = seq(0, 0.60, 0.05))+    # Set custom breaks
    scale_x_continuous(
      limits = c(0, 15),    # Set the y-axis range
      breaks = seq(0, 15, 2),    # Set custom breaks
  )
 )
plot_a
grid.arrange(plot_a, plot_b, plot_c, plot_d, nrow = 2, ncol=2 )



# empty vector for storage of Moran's I values
moran_I <- c()
moran_I_rnorm<-c()
moran_I_bin<-c()
moran_I_idw<-c()
moran_I_rnorm_idw<-c()
# loop d through a sequence ranging from 50 to 2000
for (k in seq(4, 12, 4)) {
  Rownorm <- nb2listw(knn2nb(knearneigh(parcela_vivas, k=k)), style = "W", zero.policy = TRUE)
  Binary <- nb2listw(knn2nb(knearneigh(parcela_vivas, k=k)), style = "B", zero.policy = TRUE)
  IDW <-nb2listwdist(nb2listw(knn2nb(knearneigh(parcela_vivas, k=k)),zero.policy = TRUE)$n,parcela_vivas,type="idw",zero.policy = TRUE)
  rnorm_IDW <-nb2listwdist(nb2listw(knn2nb(knearneigh(parcela_vivas, k=k)),zero.policy = TRUE)$n,parcela_vivas,type="idw",style="W",zero.policy = TRUE)
  moran1 <- moran.mc(parcela_vivas$du_annual_growth, Rownorm, nsim = 999, zero.policy = TRUE)
  moran2 <- moran.mc(parcela_vivas$du_annual_growth, Binary, nsim = 999, zero.policy = TRUE)
  moran3 <- moran.mc(parcela_vivas$du_annual_growth, IDW, nsim = 999, zero.policy = TRUE)
  moran4 <- moran.mc(parcela_vivas$du_annual_growth, rnorm_IDW, nsim = 999, zero.policy = TRUE)
  moran_I_rnorm <- c(moran_I_rnorm, moran1$statistic)
  moran_I_bin <- c(moran_I_bin, moran2$statistic)
  moran_I_idw <- c(moran_I_idw, moran3$statistic)
  moran_I_rnorm_idw <- c(moran_I_rnorm_idw, moran4$statistic)
 }
#moran_I_Akneigh_4w <- data.frame(moran1=moran_I_rnorm,moran2=moran_I_bin,moran3=moran_I_idw,moran4=moran_I_rnorm_idw,neighbours=seq(4, 12, 4))
#moran_I_Bkneigh_4w <- data.frame(moran1=moran_I_rnorm,moran2=moran_I_bin,moran3=moran_I_idw,moran4=moran_I_rnorm_idw,neighbours=seq(4, 12, 4))
#moran_I_Ckneigh_4w <- data.frame(moran1=moran_I_rnorm,moran2=moran_I_bin,moran3=moran_I_idw,moran4=moran_I_rnorm_idw,neighbours=seq(4, 12, 4))
#moran_I_Dkneigh_4w <- data.frame(moran1=moran_I_rnorm,moran2=moran_I_bin,moran3=moran_I_idw,moran4=moran_I_rnorm_idw,neighbours=seq(4, 12, 4))
#moran_I_allkneigh_4w <- data.frame(moran1=moran_I_rnorm,moran2=moran_I_bin,moran3=moran_I_idw,moran4=moran_I_rnorm_idw,neighbours=seq(4, 12, 4))


#moran_I_Bneigh_teste <- data.frame(moran1=moran_I_rnorm,moran2=moran_I_bin,moran3=moran_I_idw,moran4=moran_I_rnorm_idw,neighbours=seq(4, 12, 4))

plot_k<-(ggplot() + 
           geom_point() +
           geom_line(data = moran_I_allkneigh_4w, aes(x = neighbours, y = moran2, color = "black")) +
           geom_line(data = moran_I_allkneigh_4w, aes(x = neighbours, y = moran3, color = "blue")) +
           geom_line(data = moran_I_allkneigh_4w, aes(x = neighbours, y = moran1, color = "red")) +
           geom_line(data = moran_I_allkneigh_4w, aes(x = neighbours, y = moran4, color = "green")) +
           scale_color_identity(name = "Weights",
                                breaks = c("black", "blue", "red","green"),
                                labels = c("Binary", "IDW","Row.norm Binary","Row.norm IDW"),
                                guide = "legend")+
           ylab("Moran's I statistic")+
           ggtitle("All")+theme(plot.title = element_text(hjust = 0.5),legend.justification = "left")+
           scale_y_continuous(
             limits = c(0, 0.60),    # Set the y-axis range
             breaks = seq(0, 0.60, 0.02),    # Set custom breaks
           )
)
plot_k

#plot_k<-(ggplot(moran_I, aes(x = neighbours, y = moran)) + 
#  geom_point() +
#    geom_line())+
#  ylab("Moran's I statistic")+
#  ggtitle("Contiguity lag")+theme(plot.title = element_text(hjust = 0.5))+
#  scale_y_continuous(
#    limits = c(0, 0.20),    # Set the y-axis range
#    breaks = seq(0, 0.25, 0.02),    # Set custom breaks
#  )
grid.arrange(plot_d, plot_k, nrow = 1, top = "Site A" )


##############################################################3
########
##2.4 Local Moran's I plot
########
#A plot can be created to check 1) the precision of the higher and lower clusters; 2) if they clusters are in accordance to field observations
#A plot is created for the six neighbours approaches used before

#localm<-localmoran(parcela_vivas$du_annual_growth, listw=nb2listw(dnearneigh(parcela_vivas, d1=0, d2=8),zero.policy = TRUE,style="W"),alternative = "greater",zero.policy = TRUE)
#localm<-localmoran(parcela_vivas$du_annual_growth, listw=nb2listw(dnearneigh(parcela_vivas, d1=0, d2=10),zero.policy = TRUE,style="W"),alternative = "greater",zero.policy = TRUE)
#localm<-localmoran(parcela_vivas$du_annual_growth, listw=nb2listw(dnearneigh(parcela_vivas, d1=0, d2=15),zero.policy = TRUE,style="W"),alternative = "greater",zero.policy = TRUE)
#localm<-localmoran(parcela_vivas$du_annual_growth, listw=nb2listw(knn2nb(knearneigh(parcela_vivas, k=4)),zero.policy = TRUE,style="W"),alternative = "greater",zero.policy = TRUE)
#localm<-localmoran(parcela_vivas$du_annual_growth, listw=nb2listw(knn2nb(knearneigh(parcela_vivas, k=8)),zero.policy = TRUE,style="W"),alternative = "greater",zero.policy = TRUE)
#localm<-localmoran(parcela_vivas$du_annual_growth, listw=WW,alternative = "greater")

localm<-localmoran(parcela_vivas$du_annual_growth, listw=nb2listwdist(nb2listw(dnearneigh(parcela_vivas, d1=0, d2=6),zero.policy = TRUE)$n,parcela_vivas,type="idw",zero.policy = TRUE),alternative = "greater")

#Prepare an object for the plot
local_plot<-parcela_vivas  
quadrant <- vector(mode = "numeric", length = nrow(localm)) #prepares a vector to receive a categorical value according to each quadrant
m.qualification <- local_plot$du_annual_growth - mean(local_plot$du_annual_growth) # centers the variable of interest around its mean
m.local <- localm[, "Ii"] - mean(na.omit(localm[, "Ii"])) # centers the local Moran's around the mean
signif <- 0.1# significance threshold - important for minimizing the outliers low-low and high-high  

# builds a data quadrant
quadrant[m.qualification > 0 & m.local > 0] <- 4 #high-high values
quadrant[m.qualification < 0 & m.local < 0] <- 1 #low-low values
quadrant[m.qualification < 0 & m.local > 0] <- 2 #low-high values
quadrant[m.qualification > 0 & m.local < 0] <- 3 #high-low values
quadrant[localm[, "Pr(z > E(Ii))"] > signif] <- 0 #other values

# plot in r
brks <- c(0, 1, 2, 3, 4)
colors <- c("grey", "blue", "orange", "red", "darkgreen")
plot(local_plot[1],  col = colors[findInterval(quadrant, brks, all.inside = FALSE)],  pch = c(1, 18, 16, 16, 16)[as.factor(quadrant)])
legend("bottomleft", legend = c("insignificant","low-low","low-high","high-low","high-high"),fill=colors,bty="n")

shape_lisa <- shape %>% 
  mutate(nb = st_neighbors(geometry),
         wts = st_weights(nb),
         lag_SID79 = st_lag(SID79, nb, wts),
         lisa = categorize_lisa(SID79, lag_SID79))

# report results
ggplot(data = shape_lisa) +
  geom_sf(aes(fill = lisa))

localm[,"Ii"]
summary(localm)
teste<-localm[,"Ii"] <- localm[,"Ii"] - mean(localm[,"Ii"], na.rm = TRUE) 

teste2 <- parcela_vivas$du_annual_growth - mean(parcela_vivas$du_annual_growth, na.rm = TRUE) 

#######
#Evaluation:
#k-4 and 8m dist produces a plot with low focus in comparing to others. This is most relevant for plot A, where there is less spatial autocorrelation
#k-8, 10m, and 15m seem to produce intended clusters, with a tradeoff in precision/clusters size
#Area of influence approach fails at plotting realistic clusters

##################################################################################################################################
#### Spatial weights matrix selection according to results from neighbours network, correlogram, global moran's I and local moran's I.
#nb neighbours network --  k4, dist 8m, dist 15m and AIA were considered less adequate
#correlograms -- k4, dist 8m and AIA were considered less adequate        
#Global Moran's I test -- AIA was considered less adequate   
#Local Moran's I test -- k4, dist 8m and AIA were considered less adequate     

#Due to potentially being more accurate/balanced in describing the data, two spatial matrices will continue to the modelling phase: k-8 and dist 10m.


###############################
##3. Modelling
###############################
#Section 3 - Non-spatial and spatial modelling
#Non-spatial and spatial models will be fitted to each plot separately, and to all data, with the following steps: 
#1) Non-spatial linear modelling - Individual tree and tree group methods
#2) Spatial linear modelling; 

########
##3.1 Non-spatial modelling 
########
##3.3.1 Individual tree linear modelling
########                
X<-parcela_vivas
st_geometry(X)=NULL
f<-(du_annual_growth~Cea_0.5m_1px+Cea_1m_1px+Altimetria_1px+slope_1px+cos_aspect_1px+TRI_1px+TPI_tree_1px+TWI_tree_1px)

#Automatic stepwise selection of variables
mod<-(step(lm(f, data=X)))
summary(mod)
vif(mod) #check multicollinearity
plot(mod) #plot residuals

#Manual selection of variables - add biological interpretation to variable selection
#summary(lm(du_annual_growth~Cea_0.5m_1px+Cea_1m_1px+TPI_tree_1px+TWI_tree_1px, data=X))  #A
#summary(lm(du_annual_growth~Cea_1m_1px+slope_1px+TRI_1px, data=X))  #B
#summary(lm(du_annual_growth~Cea_0.5m_1px+Altimetria_1px, data=X))  #C
#summary(lm(du_annual_growth~Cea_0.5m_1px+Altimetria_1px+slope_1px+TWI_tree_1px, data=X))  #D
#summary(lm(du_annual_growth~Cea_0.5m_1px+Altimetria_1px+slope_1px+TPI_tree_1px, data=X)) #all

vif(lm(du_annual_growth~Cea_0.5m_1px+Altimetria_1px+slope_1px+TPI_tree_1px, data=X))
plot(lm(du_annual_growth~Cea_0.5m_1px+Altimetria_1px+slope_1px+TPI_tree_1px, data=X))
AIC(lm(du_annual_growth~Cea_0.5m_1px+Altimetria_1px+slope_1px+TPI_tree_1px, data=X))
boxplot(X$du_annual_growth)
########
##3.3.2  Group linear modelling
######## 
###Two methods were tested for creat tree group units: 1) groups of k-8 neighbours (groups of 9 trees); 2) groups by intersection of tree in fixed area grids.
X<-arvores
f<-(annuak_BA~Cea_0.5m+Cea_1m+Elev+slope+cos_aspect+TRI+TPI+TWI)
mod2<-(step(lm(f, data=X)))
summary(mod2)
vif(mod2) #check multicollinearity
plot(mod2) #plot residuals





