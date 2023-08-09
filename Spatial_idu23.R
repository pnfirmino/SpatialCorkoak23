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
df<-read_excel("C:/Users/caven/Google Drive (paulofirmino@e-isa.ulisboa.pt)/2023-Doutoramento/SpatialCorkoak/Dados_SpatialCorkOak.xlsx")
dfxy<-st_as_sf(df, coords=c("coordsX","coordsY")) #sf object
dfxy<-dfxy[which(dfxy$mais_velha == 0),]  ##Remove clearly older trees present in the plantation area
dfxy<-dfxy[which(dfxy$Morta == 0),]  ##Remove dead trees
dfxy<-dfxy[-which(dfxy$du_annual_growth == 0),] #remove trees with no diameter

## Need to calculate dug to obtain estimated influece area for each tree
dfxy$du_dug<-""; dfxy$du_dug<-dfxy$du/(sqrt(sum((dfxy$du)^2)/length(dfxy$du)))
dfxy$influence_area_est<-2.5*(28.502*exp((-66.436+4.201*(dfxy$du_dug))/(19.817+dfxy$du))/2)  Calculating individual tree area of influence, according Paulo et al. 2016

#Fix each area as an individual dataset and to the fifth s the total data
A<-dfxy[which(dfxy$Site == "A"),]; A$continId<-c(1:dim(A)[1]) #Requires remaking the Continuous Id since observations were removed (dead, older,...)
B<-dfxy[which(dfxy$Site == "B"),]; B$continId<-c(1:dim(B)[1])
C<-dfxy[which(dfxy$Site == "C"),]; C$continId<-c(1:dim(C)[1])
D<-dfxy[which(dfxy$Site == "D"),] ;D$continId<-c(1:dim(D)[1])
ABCD<-dfxy; ABCD$continId<-c(1:dim(ABCD)[1])

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
#parcela$Cea_1m_1px<-extract(Cea_1m,dfxy) #Variable Example, already on the dataframe

###############################
##2. Spatial autocorrelation analysis
###############################
#Section 2 - Spatial autocorrelation analysis
#Spatial autocorrelation analysis is applied to each plot separately, with the following steps: 
#1) Spatial matrix definition;
#2) Spatial weights definition; 
#3)  Moran's I analysis
#4) Correlogram (Moran's I statistics) analysis

########
##2.1 Spatial weights matrix definition
########             
#Although plantations are regular, trees are not at exact distances, not being preferable using polygons for analysis.
#Its opted to use points for the analysis. Output of 2.1 is nb object.
#Three methods are tested to choice the proper spatial matrix # Distance based neighbours; k neighbours; Tree area of influence sobreposition

## The section 2 is put as cycle, with 5 loops for datasets A, B, C, D and ABCD
#Prepare objects for cycle
p<-list(A,B,C,D,ABCD)
code<-list("A","B","C","D","ABCD")
tabela<-list("tabelaA","tabelaB","tabelaC","tabelaD","tabelaABCD")
list<-c()
kneigh<-c()

for(i in c(1:length(p))){ 
 
  
#2.1.1 Distance base neighbours listw objects creation #To find the most adequate distances for the lag, we calculate the listw object from 2m to 15m 
#Simultaneously we test the most adequate weights 1) Row-normalized binary (rnorm_bin); Binary (bin); Row-normalized Inverted distance weights (rnorm_IDW); Inverted distance weights (IDW)

#Cycle to run all Moran's I tests for each combination of lag/weight - slow 
# empty vector for storage of Moran's I statistic values
moran_I_rnorm_bin<- c()
moran_I_bin <- c()
moran_I_rnorm_idw<-c()
moran_I_idw <- c()
  
for (d in seq(2, 15, 1)) {
Rownorm <- nb2listw(dnearneigh( p[[i]], d1 = 0, d2 = d), style = "W", zero.policy = TRUE) #Row-normalized binary weights (style ="w")
Binary <- nb2listw(dnearneigh( p[[i]], d1 = 0, d2 = d), style = "B", zero.policy = TRUE) #Binary weights (style ="B")
#To use that idw weights, it requires nb2listwdist
IDW <-nb2listwdist(dnearneigh( p[[i]], d1=0, d2=d),p[[i]],type="idw",zero.policy = TRUE) #Row-normalized idw (style ="w")
rnorm_IDW <-nb2listwdist(dnearneigh(p[[i]], d1=0, d2=d),p[[i]],type="idw",style="W", zero.policy = TRUE) #idw (style ="Raw")
#Calculate Moran's I statistics
moran1 <- moran.mc(p[[i]]$du_annual_growth, Rownorm, nsim = 999, zero.policy = TRUE)
moran2 <- moran.mc(p[[i]]$du_annual_growth, Binary, nsim = 999, zero.policy = TRUE)
moran3 <- moran.mc(p[[i]]$du_annual_growth, IDW, nsim = 999, zero.policy = TRUE)
moran4 <- moran.mc(p[[i]]$du_annual_growth, rnorm_IDW, nsim = 999, zero.policy = TRUE) #Add the value to the empty vector
moran_I_rnorm_bin <- c(moran_I_rnorm_bin, moran1$statistic)
moran_I_bin <- c(moran_I_bin, moran2$statistic)
moran_I_idw <- c(moran_I_idw, moran3$statistic)
moran_I_rnorm_idw <- c(moran_I_rnorm_idw, moran4$statistic)
}
#Finalize the dataframe for plotting
moran_I_alldist_4w <- data.frame(moran1 =  moran_I_rnorm,moran2 =  moran_I_bin,moran3 =  moran_I_idw,moran4 = moran_I_rnorm_idw, distance = seq(2, 15, 1))
  #Plot and compare lags/weights for and betwee each site
 dist[[i]]<-(ggplot() + 
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
              ggtitle(paste0(code[[i]]))+theme(plot.title = element_text(hjust = 0.5),legend.justification = "left")+
              scale_y_continuous(
                limits = c(0, 0.60),    # Set the y-axis range
                breaks = seq(0, 0.60, 0.05))+    # Set custom breaks
              scale_x_continuous(
                limits = c(0, 15),    # Set the y-axis range
                breaks = seq(0, 15, 2),    # Set custom breaks
              )
)

#2.1.2 k-neighbours base listw objects creation #To find the most adequate number of neighbourslag, we calculate the listw object for k=4, k=8 (and k=12 to see the evolution)
#Run similar cycle for k-neighbours now
# empty vector for storage of Moran's I values
moran_I <- c()
moran_I_rnorm<-c()
moran_I_bin<-c()
moran_I_idw<-c()
moran_I_rnorm_idw<-c()
for (k in seq(4, 12, 4)) {
  Rownorm <- nb2listw(knn2nb(knearneigh(p[[i]], k=k)), style = "W", zero.policy = TRUE)
  Binary <- nb2listw(knn2nb(knearneigh(p[[i]], k=k)), style = "B", zero.policy = TRUE)
  IDW <-nb2listwdist(knn2nb(knearneigh(p[[i]], k=k)),p[[i]],type="idw",zero.policy = TRUE)
  rnorm_IDW <-nb2listwdist(knn2nb(knearneigh(p[[i]], k=k)),p[[i]],type="idw",style="W",zero.policy = TRUE)
  moran1 <- moran.mc(p[[i]]$du_annual_growth, Rownorm, nsim = 999, zero.policy = TRUE)
  moran2 <- moran.mc(p[[i]]$du_annual_growth, Binary, nsim = 999, zero.policy = TRUE)
  moran3 <- moran.mc(p[[i]]$du_annual_growth, IDW, nsim = 999, zero.policy = TRUE)
  moran4 <- moran.mc(p[[i]]$du_annual_growth, rnorm_IDW, nsim = 999, zero.policy = TRUE)
  moran_I_rnorm <- c(moran_I_rnorm, moran1$statistic)
  moran_I_bin <- c(moran_I_bin, moran2$statistic)
  moran_I_idw <- c(moran_I_idw, moran3$statistic)
  moran_I_rnorm_idw <- c(moran_I_rnorm_idw, moran4$statistic)
}
moran_I_allkneigh_4w <- data.frame(moran1 =  moran_I_rnorm,moran2 =  moran_I_bin,moran3 =  moran_I_idw,moran4 = moran_I_rnorm_idw, neighbours = seq(4, 12, 4))

kneigh[[i]]<-(ggplot() + 
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
           ggtitle(paste0(code[[i]]))+theme(plot.title = element_text(hjust = 0.5),legend.justification = "left")+
           scale_y_continuous(
             limits = c(0, 0.60),    # Set the y-axis range
             breaks = seq(0, 0.60, 0.02),    # Set custom breaks
           )
)
  
#2.1.3 Identifying the trees with sobreposed areas of influence and providing a value for them as neighbours
#Creating of neighbour/distance table; 
  matrix<-as.matrix(dist(st_coordinates(p[[i]]), "euclidean"), labels=TRUE) #Distance matrix
  colnames(matrix) <- rownames(matrix) <- p[[i]]$continId #Name rows and columns
  melt <- melt(matrix) #unmaked the matrix to columns
  melt_0<-melt[which(melt$value!=0),]  #removes the distances=zero
  M101<-data.frame(p[[i]]) #change to a different object to avoid errors from p[[i]]
  submelt01<-melt[melt_0$Var1 %in%  M101$continId & melt_0$Var2 %in%  M101$continId,] #creates a table with Var1, Var2 and the distance value between then
  submelt01<-submelt01[which(submelt01$value<33),] #Removes any lines with distance > 33m, which is the highest distance between two neighbouring trees. This condition fastens the process.
  submelt01<-submelt01[which(submelt01$value!=0),] #Removes distances =0
  mt<-merge(submelt01, M101[,c("continId","influence_area_est" )], by.x="Var1", by.y="continId") #Creates a diameter column for Var1, depending on tree number
  mtt<-merge(mt, M101[,c("continId","influence_area_est" )], by.x="Var2", by.y="continId") #Creates a diameter column for Var2, depending on tree number
  names(mtt)<-c(names(mtt)[1:3], "influence_area_est2",  "influence_area_est1") #change the names of the two new columns
  head(mtt)
  
  # Identifying the trees with sobreposed areas of influence and providing a value for them as neighbours
  matrix<-matrix(nrow=dim(M101)[1], ncol=dim(M101)[1]) #starts a blank matrix for the cycle
  for (k in c(1:nrow(mtt))){ 
    valor1<-mtt[k,1]
    valor2<-mtt[k,2]
    linha<- mtt[k,]
    matrix[valor1,valor2]<- ifelse((linha$influence_area_est1 + linha$influence_area_est2 > linha$value),1,0)
  }
  #asks if the sum of radius of area of influence of two trees is higher than their distance. 
  #If true, they are neighbours. The value attributed if neighbours can work as a weight. In this case is 1, similar to binary case.
  
  matrix[is.na(matrix)] = 0 #turns NA to zeros
  colnames(matrix) <- rownames(matrix) <- p[[i]]$continId #give names to new matrix
  #View(matrix) #check if result is correct
  W<- mat2listw(matrix) #creates a listw object from this matrix, ready to be applied, with style=matrix itself.
  WW<-mat2listw(matrix,style = "W") #creates a listw object from this matrix, but row-normalized
  ###Repeat the cycle about attribute the distance between instead of 1. 
  matrix<-matrix(nrow=dim(p[[i]])[1], ncol=dim(p[[i]])[1]) #starts a blank matrix for the cycle
  for (k in c(1:nrow(mtt))){ 
    valor1<-mtt[k,1]
    valor2<-mtt[k,2]
    linha<- mtt[k,]
    matrix[valor1,valor2]<- ifelse((linha$influence_area_est1 + linha$influence_area_est2 > linha$value),linha$value,0)
  }
  matrix[is.na(matrix)] = 0 #turns NA to zeros
  colnames(matrix) <- rownames(matrix) <- p[[i]]$continId #give names to new matrix
  #View(matrix) #check if result is correct
  Wdist<- mat2listw(matrix) #creates a listw object from this matrix, this time with dist as weights
  WdistW<- mat2listw(matrix,style="W") #row normalized idw weights

 ####2.1.4  Plot and compare results /Global Global Moran'I test

# A table can be created fro each plot to compare the multiple methods outputs: Morans'I statistic and respective p-value.
#Example of the output table, plot A:
#                         Weights    Row-normalized Binary          Binary                  Row-normalized idw              idw                 
#                          Method Moran I statistic p-value Moran I statistic p-value Moran I statistic p-value  Moran I statistic p-value
#  Distance based neighbours - 5m             XXX       0             XXX       0             XXX         0            XXX         0
# Distance based neighbours -  6m             XXX       0              XXX       0             XXX        0            XXX         0
# Distance based neighbours -  7m             XXX       0             XXX       0             XXX         0            XXX         0
#                k neighbours - 4             XXX       0             XXX       0             XXX         0            XXX         0
#                k neighbours - 8             XXX       0             XXX       0             XXX         0            XXX         0
#      Area of Influence approach             XXX   4e-05             XXX       0             XXX         0            XXX         0

  
tab<-data.frame()
tab[1,1]<- "Method"
tab[2,1]<- "Distance based neighbours - 5m"
tab[3,1]<- "Distance based neighbours - 6m"
tab[4,1]<- "Distance based neighbours - 7m"
tab[5,1]<- " k neighbours - 4"
tab[6,1]<- " k neighbours - 8"
tab[7,1]<- "Area of Influence approach"
tab[1,2]<-"Moran I statistic";tab[1,3]<-"p-value";tab[1,4]<-"Moran I statistic";tab[1,5]<-"p-value";tab[1,6]<-"Moran I statistic";tab[1,7]<-"p-value";
#distance base tests - row-normalized weights; zero.policy=TRUE must be present since there is always points with no (living) neighbours
tab[2,2]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(dnearneigh(p[[i]], d1=0, d2=5),zero.policy = TRUE),zero.policy=TRUE)$estimate[1],3)
tab[2,3]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(dnearneigh(p[[i]], d1=0, d2=5),zero.policy = TRUE),zero.policy=TRUE)$p.value,5)
tab[3,2]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(dnearneigh(p[[i]], d1=0, d2=6),zero.policy = TRUE),zero.policy=TRUE)$estimate[1],3)
tab[3,3]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(dnearneigh(p[[i]], d1=0, d2=6),zero.policy = TRUE),zero.policy=TRUE)$p.value,5)
tab[4,2]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(dnearneigh(p[[i]], d1=0, d2=7),zero.policy = TRUE),zero.policy=TRUE)$estimate[1],3)
tab[4,3]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(dnearneigh(p[[i]], d1=0, d2=7),zero.policy = TRUE),zero.policy=TRUE)$p.value,5)
#distance base tests - inverse of distance weights; can be used nb2listwdist function
tab[2,4]<-round(moran.test(p[[i]]$du_annual_growth,nb2listwdist(dnearneigh(p[[i]], d1=0, d2=5),p[[i]],type="idw",zero.policy = TRUE),zero.policy=TRUE)$estimate[1],3)
tab[2,5]<-round(moran.test(p[[i]]$du_annual_growth,nb2listwdist(dnearneigh(p[[i]], d1=0, d2=5),p[[i]],type="idw",zero.policy = TRUE),zero.policy=TRUE)$p.value,5)
tab[3,4]<-round(moran.test(p[[i]]$du_annual_growth,nb2listwdist(dnearneigh(p[[i]], d1=0, d2=6),p[[i]],type="idw",zero.policy = TRUE),zero.policy=TRUE)$estimate[1],3)
tab[3,5]<-round(moran.test(p[[i]]$du_annual_growth,nb2listwdist(dnearneigh(p[[i]], d1=0, d2=6),p[[i]],type="idw",zero.policy = TRUE),zero.policy=TRUE)$p.value,5)
tab[4,4]<-round(moran.test(p[[i]]$du_annual_growth,nb2listwdist(dnearneigh(p[[i]], d1=0, d2=7),p[[i]],type="idw",zero.policy = TRUE),zero.policy=TRUE)$estimate[1],3)
tab[4,5]<-round(moran.test(p[[i]]$du_annual_growth,nb2listwdist(dnearneigh(p[[i]], d1=0, d2=7),p[[i]],type="idw",zero.policy = TRUE),zero.policy=TRUE)$p.value,5)
#distance base tests - binary weights
tab[2,6]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(dnearneigh(p[[i]], d1=0, d2=5),zero.policy = TRUE,style="B"),zero.policy=TRUE)$estimate[1],3)
tab[2,7]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(dnearneigh(p[[i]], d1=0, d2=5),zero.policy = TRUE,style="B"),zero.policy=TRUE)$p.value,5)
tab[3,6]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(dnearneigh(p[[i]], d1=0, d2=6),zero.policy = TRUE,style="B"),zero.policy=TRUE)$estimate[1],3)
tab[3,7]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(dnearneigh(p[[i]], d1=0, d2=6),zero.policy = TRUE,style="B"),zero.policy=TRUE)$p.value,5)
tab[4,6]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(dnearneigh(p[[i]], d1=0, d2=7),zero.policy = TRUE,style="B"),zero.policy=TRUE)$estimate[1],3)
tab[4,7]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(dnearneigh(p[[i]], d1=0, d2=7),zero.policy = TRUE,style="B"),zero.policy=TRUE)$p.value,5)
#k neighbours - row-normalized weights
tab[5,2]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(knn2nb(knearneigh(p[[i]], k=4)),style="W"),zero.policy=TRUE)$estimate[1],3)
tab[5,3]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(knn2nb(knearneigh(p[[i]], k=4)),style="W"),zero.policy=TRUE)$p.value,5)
tab[6,2]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(knn2nb(knearneigh(p[[i]], k=8)),style="W"),zero.policy=TRUE)$estimate[1],3)
tab[6,3]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(knn2nb(knearneigh(p[[i]], k=8)),style="W"),zero.policy=TRUE)$p.value,5)
#k neighbours - inverse of distance weights 
tab[5,4]<-round(moran.test(p[[i]]$du_annual_growth,nb2listwdist(knn2nb(knearneigh(p[[i]], k=4)),p[[i]],type="idw"),zero.policy=TRUE)$estimate[1],3)
tab[5,5]<-round(moran.test(p[[i]]$du_annual_growth,nb2listwdist(knn2nb(knearneigh(p[[i]], k=4)),p[[i]],type="idw"),zero.policy=TRUE)$p.value,5)
tab[6,4]<-round(moran.test(p[[i]]$du_annual_growth,nb2listwdist(knn2nb(knearneigh(p[[i]], k=8)),p[[i]],type="idw"),zero.policy=TRUE)$estimate[1],3)
tab[6,5]<-round(moran.test(p[[i]]$du_annual_growth,nb2listwdist(knn2nb(knearneigh(p[[i]], k=8)),p[[i]],type="idw"),zero.policy=TRUE)$p.value,5)
#k neighbours - binary weights
tab[5,6]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(knn2nb(knearneigh(p[[i]], k=4)),style="B"),zero.policy=TRUE)$estimate[1],3)
tab[5,7]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(knn2nb(knearneigh(p[[i]], k=4)),style="B"),zero.policy=TRUE)$p.value,5)
tab[6,6]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(knn2nb(knearneigh(p[[i]], k=8)),style="B"),zero.policy=TRUE)$estimate[1],3)
tab[6,7]<-round(moran.test(p[[i]]$du_annual_growth,nb2listw(knn2nb(knearneigh(p[[i]], k=8)),style="B"),zero.policy=TRUE)$p.value,5)
#Area of influence neighbours - three weights (WW, Wdist, W)
tab[7,2]<-round(moran.test(p[[i]]$du_annual_growth,WW, zero.policy=TRUE)$estimate[1],3)
tab[7,3]<-round(moran.test(p[[i]]$du_annual_growth,WW, zero.policy=TRUE)$p.value,5)
tab[7,4]<-round(moran.test(p[[i]]$du_annual_growth,Wdist, zero.policy=TRUE)$estimate[1],3)
tab[7,5]<-round(moran.test(p[[i]]$du_annual_growth,Wdist, zero.policy=TRUE)$p.value,5)
tab[7,6]<-round(moran.test(p[[i]]$du_annual_growth,W, zero.policy=TRUE)$estimate[1],3)
tab[7,7]<-round(moran.test(p[[i]]$du_annual_growth,W, zero.policy=TRUE)$p.value,5)
nomes<-c("Weights","Row-normalized","","Inverse Distance","", "Binary","");colnames(tab)<-nomes
tabela[[i]]<-tab #save the table

####2.1.5  Plot and compare results
par(mfrow=c(2,3));par(cex=0.6, pch=16)
plot(dnearneigh(p[[i]], d1=0, d2=5), coord=(p[[i]]$geometry), col="red");title(main=paste("Distance based neighbours - 5m"),cex.main=2)
plot(dnearneigh(p[[i]], d1=0, d2=6), coord=(p[[i]]$geometry), col="red", lwd = 2);title(main=paste("Distance based neighbours - 6m"),cex.main=2)
plot(dnearneigh(p[[i]], d1=0, d2=7), coord=(p[[i]]$geometry), col="red", lwd = 3);title(main=paste("Distance based neighbours - 7m"),cex.main=2)
plot(knn2nb(knearneigh(p[[i]], k=4)), coord=(p[[i]]$geometry), col="red");title(main=paste("k neighbours - 4"),cex.main=2)
plot(knn2nb(knearneigh(p[[i]], k=8)), coord=(p[[i]]$geometry), col="red");title(main=paste("k neighbours - 8"),cex.main=2)
plot(WW,coord=(p[[i]]$geometry),col="red"); title(main=paste("Area of Influence approach"),cex.main=2)

####2.1.6 Correlogram
par(mfrow=c(2,3))
plot(sp.correlogram(knn2nb(knearneigh(p[[i]], k=4)), p[[i]]$du_annual_growth, method="I" ,order=9,zero.policy=TRUE),ylim=c(0,0.35), main="k4 neighbours")
plot(sp.correlogram(knn2nb(knearneigh(p[[i]], k=8)), p[[i]]$du_annual_growth, method="I" ,order=9,zero.policy=TRUE),ylim=c(0,0.35), main="k8 neighbours")
plot(sp.correlogram(dnearneigh(p[[i]], d1=0, d2=5), p[[i]]$du_annual_growth, method="I" ,order=9,zero.policy=TRUE),ylim=c(0,0.35), main="dist 5m")
plot(sp.correlogram(dnearneigh(p[[i]], d1=0, d2=6), p[[i]]$du_annual_growth, method="I" ,order=9,zero.policy=TRUE),ylim=c(0,0.35), main="dist 6m")
plot(sp.correlogram(dnearneigh(p[[i]], d1=0, d2=7), p[[i]]$du_annual_growth, method="I" ,order=9,zero.policy=TRUE),ylim=c(0,0.35) ,main="dist 7m")
plot(sp.correlogram(WW$neighbours, p[[i]]$du_annual_growth, method="I" ,order=9,zero.policy=TRUE),ylim=c(0,0.35) ,main="Area of Influence")
} #  End of Spatial weights matrix cycle

############
####2.1.7 Results
############
## Output objects
#Plots of Moran's I statistic values from distance lags (2 to 15m) for individual plots A,B,C,D
grid.arrange(dist[[1]], dist[[2]], dist[[3]], dist[[4]], ncol=2 )
#Plots of Moran's I statistic values from distance lags (2 to 15m) for ABCD
dist[[5]]
#Plots of Moran's I statistic values from  k neighbours lags (4 to 8m) for individual plots A,B,C,D
grid.arrange(kneigh[[1]], kneigh[[2]], kneigh[[3]], kneigh[[4]], ncol=2 )
#Plots of Moran's I statistic values from  k neighbours lags (4 to 8m) for ABCD
kneigh[[5]]
#Tables of Moran's I statistics for lags and weights
tabela[[1]];tabela[[2]]; tabela[[3]]; tabela[[4]]; tabela[[5]]

#Evaluation:
#Best distance-base lags is between 5-7m for any site
#Row-normalized weights have lower Moran I values, but ensure spatial weights statistics assumptions to correctly run models
#Row-normalized binary or row-normalized weights are most adequate depending on the site
#K-neighbours-based lags is more adequate for k=4 for any site
#Values can be close to the distance-base lags. A table has be created for a direct comparison.

#Store the choosen most adequate spatial weights matrices for modelling phase
listw_A<-nb2listw(dnearneigh(A, d1=0, d2=6),zero.policy = TRUE,style="W")
listw_B<-nb2listwdist(dnearneigh(B, d1=0, d2=7),B,type="idw",style="W",zero.policy = TRUE)
listw_C<-nb2listwdist(nb2listw(knn2nb(knearneigh(C, k=4)))$n,C,type="idw",style="W")
listw_D<-nb2listwdist(nb2listw(knn2nb(knearneigh(D, k=4)))$n,D,type="idw",style="W") 
listw_ABCD<-nb2listwdist(nb2listw(knn2nb(knearneigh(ABCD, k=4)))$n,ABCD,type="idw",style="W")



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
X<-dfxy
st_geometry(X)=NULL
f<-(du_annual_growth~Cea_0.5m_1px+Cea_1m_1px+Altimetria_1px+slope_1px+cos_aspect_1px+TRI_1px+TPI_tree_1px+TWI_tree_1px)

#Automatic stepwise selection of variables
mod<-(step(lm(f, data=X)))
summary(mod)
vif(mod) #check multicollinearity
plot(mod) #plot residuals

#Manual selection of variables - add biological interpretation to variable selection
#List of full model and most adequate submodel, to use on spatial modelling phase
modelA=lm(du_annual_growth~Cea_0.5m_1px+Cea_1m_1px+Altimetria_1px+slope_1px+cos_aspect_1px+TRI_1px+TPI_tree_1px+TWI_tree_1px,data=X)
submodelA=lm(du_annual_growth~Cea_0.5m_1px+Cea_1m_1px+TPI_tree_1px+TWI_tree_1px,data=X)
modelB=lm(du_annual_growth~Cea_0.5m_1px+Cea_1m_1px+Altimetria_1px+slope_1px+cos_aspect_1px+TRI_1px+TPI_tree_1px+TWI_tree_1px,data=X)
submodelB=lm(du_annual_growth~slope_1px+TRI_1px,data=X)
modelC=lm(du_annual_growth~Cea_0.5m_1px+Cea_1m_1px+Altimetria_1px+slope_1px+cos_aspect_1px+TRI_1px+TPI_tree_1px+TWI_tree_1px,data=X)
submodelC=lm(du_annual_growth~Cea_0.5m_1px+Altimetria_1px,data=X)
modelD=lm(du_annual_growth~Cea_0.5m_1px+Cea_1m_1px+Altimetria_1px+slope_1px+cos_aspect_1px+TRI_1px+TPI_tree_1px+TWI_tree_1px,data=X)
submodelD=lm(du_annual_growth~Cea_1m_1px+Altimetria_1px+slope_1px+TWI_tree_1px,data=X)
modelABCD=lm(du_annual_growth~Cea_0.5m_1px+Cea_1m_1px+Altimetria_1px+slope_1px+cos_aspect_1px+TRI_1px+TPI_tree_1px+TWI_tree_1px,data=X)
submodelABCD=lm(du_annual_growth~Cea_0.5m_1px+Cea_1m_1px+Altimetria_1px+slope_1px+TPI_tree_1px,data=X)

#Model assumptions verification
#vif(lm(du_annual_growth~Cea_0.5m_1px+Altimetria_1px+slope_1px+TPI_tree_1px, data=X)) #check multicollinearity
#plot(lm(du_annual_growth~Cea_0.5m_1px+Altimetria_1px+slope_1px+TPI_tree_1px, data=X)) #check model assumptions
#moran.test(model$residuals,listw_A ,zero.policy = TRUE) #Check if there is autocorrelation on the residuals
#AIC(lm(du_annual_growth~Cea_0.5m_1px+Altimetria_1px+slope_1px+TPI_tree_1px, data=X)) #check 

#####################################################################

########
##3.2 Spatial modelling
######## 
#Choosing the most suitable model between OLS, Spatial Lag Model (SLM), Spatial Error Model (SEM), and Simultaneous Autocorrelated Model (SAR)


submodel<-list(submodelA, submodelB, submodelC, submodelD, submodelABCD)
spatialwm<-list(listw_A,listw_B,listw_C,listw_D,listw_ABCD )
splag<-list("splagA","splagB","splagC","splagD","splagABCD")
sperror<-list("sperrorA","sperrorB","sperrorC","sperrorD","sperrorABCD")
spSAC<-list("spSACA","spSACB","spSACC","spSACD","spSACABCD")

for(i in c(1:length(p))){ 
  
  lm.LMtests(submodel[[i]],spatialwm[[i]], zero.policy = TRUE, test="all") #Check the scores from Lagrange Modifier tests and SARMA

  #Run models so a likelihood tests may be applied. Slow process, but is once required once.
  splag[[i]]<-lagsarlm(submodel[[i]],data=p[[i]],listw=spatialwm[[i]],zero.policy = TRUE)
  sperror[[i]]<-errorsarlm(submodel[[i]],data=p[[i]],listw=spatialwm[[i]],zero.policy = TRUE)
  spSAC[[i]]<-sacsarlm(submodel[[i]],data=p[[i]],listw=spatialwm[[i]],zero.policy = TRUE)

 LR.Sarlm(spatiallag,spatialSAC) #SLM vs SAR
 LR.Sarlm(spatialerr,spatialSAC) #SEM vs SAR 
  }
#By examining LMtests scores, some models are similar so a likelihood test should also be tried (this emplied that all spatial models were run previously)
#######
#Evaluation:
#All OLS models have spatial autocorrelations in residuals
#SEM is as good as SAR for site A, so SEM is selected
#SAR is more adequate for remaining sites

#Spatial Models Applied
#A
summary(sperror[[1]])
AIC(sperror[[1]])
round(1-(sum((p[[1]]$du_annual_growth-(sperror[[1]]$fitted.values))^2)/sum((p[[1]]$du_annual_growth-mean(p[[1]]$du_annual_growth))^2)),4) #R2 
Ypred<-fitted(sperror[[1]]) ; plot(Ypred~p[[1]]$du_annual_growth)
moran.test(sperror[[1]]$residuals,listw_A ,zero.policy = TRUE) #check spatial autocorrelation from residuals
#B
summary(spSAC[[2]])
round(1-(sum((p[[2]]$du_annual_growth-(spSAC[[2]]$fitted.values))^2)/sum((p[[2]]$du_annual_growth-mean(p[[2]]$du_annual_growth))^2)),4)
Ypred<-fitted(spSAC[[2]]) ; plot(Ypred~p[[2]]$du_annual_growth)
moran.test(spSAC[[2]]$residuals,listw_B ,zero.policy = TRUE)                
#C
spatialSACC<-sacsarlm(lm(du_annual_growth~Cea_0.5m_1px+Altimetria_1px,data=p[[3]],listw=listw_C,zero.policy = TRUE)
summary(spatialSACC)
round(1-(sum((p[[3]]$du_annual_growth-(spatialSACC$fitted.values))^2)/sum((p[[3]]$du_annual_growth-mean(p[[3]]$du_annual_growth))^2)),4)
Ypred<-fitted(spatialSACC) ; plot(Ypred~p[[3]]$du_annual_growth)
moran.test(spatialSACC$residuals,listw_C ,zero.policy = TRUE)
#D
spatialSACD<-sacsarlm(model,data=p[[4]],listw=listw_D,zero.policy = TRUE)
summary(spatialSACD)
round(1-(sum((p[[4]]$du_annual_growth-(spatialSACD$fitted.values))^2)/sum((p[[4]]$du_annual_growth-mean(p[[4]]$du_annual_growth))^2)),4)
Ypred<-fitted(spatialSACD) ; plot(Ypred~p[[4]]$du_annual_growth)
moran.test(spatialSACD$residuals,listw_D ,zero.policy = TRUE)
#ABCD
spatialSACABCD<-sacsarlm( lm(du_annual_growth~Cea_0.5m_1px+Altimetria_1px+slope_1px,data=parcela_vivas),
  data=p[[5]],listw=listw_ABCD,zero.policy = TRUE)
summary(spatialSACABCD)
Ypred<-fitted(spatialSACABCD)
plot(Ypred~p[[5]]$du_annual_growth)
round(1-(sum((p[[5]]$du_annual_growth-(spatialSACABCD$fitted.values))^2)/sum((p[[5]]$du_annual_growth-mean(p[[5]]$du_annual_growth))^2)),4)
moran.test(spatialSACABCD$residuals,listw_ABCD ,zero.policy = TRUE)


















EXTRA - Removed from paper
##############################################################3
########
##2.3 Local Moran's I plot
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


