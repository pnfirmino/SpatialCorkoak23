library(sf)
library(spdep)
library(whitebox)
library(car)
ibrary(spatialreg)
library(gstat)


###Section 1 - Data preparation and variable calculation
#Upload and prepare data file
Dados_SpatialCorkOak_9_maio_23 <- read_excel("/path/Dados_SpatialCorkOak_9_maio_23.xlsx")
parcela<-Dados_SpatialCorkOak1[which(Dados_SpatialCorkOak1$Site == "A"),] ## Select on of the 4 properties (A, B, C, D) - A and B are Santarém, C and D are Castelo Branco 
parcela<-parcela[which(parcela$mais_velha == 0),]  ##Remove clearly older trees present in the plantation area

########
### Variable calculation - these variable are already calculated on the sf object. Code displayed in this section is only for informative purposes
########

#1.1 Geographical position indices calculation (TPI and TRI)
f <- matrix(1, nrow=5, ncol=5)
TPI <- focal(Altimetria, w=f, fun=function(x, ...) x[5] - mean(x[-5]), pad=TRUE, padValue=NA)
plot(TPI)
TRI <- focal(Altimetria, w=f, fun=function(x) sum(abs(x[-5]-x[5]))/8, pad=TRUE, padValue=NA)
plot(TRI)

#1.2 Geographical position indices calculation (TWI)
# TWI calculation based on the methodology of https://vt-hydroinformatics.github.io/rgeoraster.html
  # Prepare DEM for Hydrology Analyses
  wbt_breach_depressions_least_cost(
  dem = "/path/parcela_Altimetria.tif",
  output = "/path/parcela_DTM_breach.tif",
  dist = 5,
  fill = TRUE)
  # wbt_fill_depressions_wang_and_liu(
  dem = "/path/parcela_DTM_breach.tif",
  output = "/path/parcela_DTM_breach_fill.tif",
  )
  # Visualize and correct filled sinks and breached depression
  filled_breached <- raster(/path/"parcela_DTM_breach_fill.tif")
  plot(filled_breached)
    difference <- Altimetria - filled_breached
  difference[difference == 0] <- NA
  #D infinity flow accumulation (alternative flow accumulation may be calculated from D infinity method D8 Flow Accumulation)
  wbt_d_inf_flow_accumulation("/path/parcela_DTM_breach_fill.tif",
                              "/path/Infinit_FlowAccum.tif")
  dinf <- raster("/path/Infinit_FlowAccum.tif")
  plot(dinf)
  #Calculate Specific Contributing Area
  wbt_d_inf_flow_accumulation(input = "/path/parcela_DTM_breach_fill.tif",
                            output = "/path/FlowAccum2.tif",
                            out_type = "Specific Contributing Area")
  #Calculate slope or use slope mapped tif from data
  wbt_slope(dem = "/path/parcela_DTM_breach_fill.tif",
          output = "/path/parcela_slope_degrees.tif",
          units = "degrees")
  #Calculate topographic wetness index
  wbt_wetness_index(sca = "/path/FlowAccum2.tif",
                  slope = "/path/parcela_slope_degrees.tif",
                  output = "/path/TWI.tif")
  twi <- raster("/path/TWI.tif")
  plot(twi)

########
#1.3 Calculate diameter related variables and attribute explanatory values to points
# Response variable is considered in the ways: 1) as the individual tree diameter; 2) as the basal area of a group of trees
########
#1.3.1 - Attribute explanatory values to points (individual tree positions)
#Attribute exploratory variable value according to the pixel where tree is positioned
parcela$Cea_1m_1px<-extract(Cea_1m,parcela) #Variable Example
########
#1.3.2


########
#Section 2 - Spatial autocorrelation analysis
#Spatial autocorrelation analysis is applied to each plot separately, with the following steps: 
#1) Spatial matrix definition;
#2) Spatial weights definition; 
#3) Global Moran's I analysis
#4) Correlogram (Moran's I statistics) analysis
#5) Local Moran's I analysis

#2.1 Spatial matrix definition
#Although plantations are regular, trees are not at exact distances to make preferable the using of polygons for analysis.
#Its opted using points for the analysis. Output of 2.1 is nb object.
#Three methods are tested to choice the proper spatial matrix 
#1) Distance based neighbours; k neighbours; Tree area of influence sobreposition
#2.1.1 Distance base neighbours nb objects creation - Selected 8, 10 and 15m. The idea is capture links according to plantation spacing
dnearneigh(parcela_vivas, d1=0, d2=8)
dnearneigh(parcela_vivas, d1=0, d2=10)
dnearneigh(parcela_vivas, d1=0, d2=15)
#2.1.2 k neighbours - Selected 4 and 8 neighbours.
knn2nb(knearneigh(parcela_vivas, k=4))
knn2nb(knearneigh(parcela_vivas, k=8))
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
matrix<-matrix(nrow=length(parcela_vivas), ncol=length(parcela_vivas)) #starts a blank matrix for the cycle
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
View(matrix) #check if result is correct
W<- mat2listw(matrix) #creates a listw object from this matrix, ready to be applied.


####2.1.4  Plot and compare results
par(mfrow=c(2,3));par(cex=0.4, pch=16)
plot(dnearneigh(parcela_vivas, d1=0, d2=8), coord=(parcela_vivas$geometry), col="red");title(main=paste("Distance based neighbours - 8m"),cex.main=2)
plot(dnearneigh(parcela_vivas, d1=0, d2=10), coord=(parcela_vivas$geometry), col="red");title(main=paste("Distance based neighbours -10m"),cex.main=2)
plot(dnearneigh(parcela_vivas, d1=0, d2=15), coord=(parcela_vivas$geometry), col="red");title(main=paste("Distance based neighbours - 15m"),cex.main=2)
plot(knn2nb(knearneigh(parcela_vivas, k=4)), coord=(parcela_vivas$geometry), col="red");title(main=paste("k neighbours - 4"),cex.main=2)
plot(knn2nb(knearneigh(parcela_vivas, k=8)), coord=(parcela_vivas$geometry), col="red");title(main=paste("k neighbours - 8"),cex.main=2)
plot(W,coord=(parcela_vivas$geometry),col="red"); title(main=paste("Area of Influence approach"),cex.main=2)












1.3.2
#######################################
#Considering du_mean response variable

##### Creating variables considering a moving window of the subject tree and the 8 closer neighbours (queen's case continuity)

#Remove lines with NA's in TRI/TPI, respective to points very close to the limits of the area and the 5x5 matrix from TRI is not completetly available
A <- A[!is.na(A$TRI_1px), ]  

#Calculation of du_mean and attribution of values to tree points

#Creating a dataframe for each subject tree and respective closest 8 neighbours (queen case)
valor_du_id <- data.frame(matrix(NA,  ncol = 39))
valor_du_1 <- data.frame(matrix(NA,  ncol = 39))
valor_du_2 <- data.frame(matrix(NA,  ncol = 39))
valor_du_3 <- data.frame(matrix(NA,  ncol = 39))
valor_du_4 <- data.frame(matrix(NA,  ncol = 39))
valor_du_5 <- data.frame(matrix(NA,  ncol = 39))
valor_du_6 <- data.frame(matrix(NA,  ncol = 39))
valor_du_7 <- data.frame(matrix(NA,  ncol = 39))
valor_du_8 <- data.frame(matrix(NA,  ncol = 39))

#Find each subject tree and respective closest 8 neighbours (queen case)
A$continId<-c(1:dim(A)[1])
for(i in c(1:nrow(A))) {
  valor_du_id[i,]<-A[i,]
  valor_du_1[i,]<-(A[which(A$continId == ((knearneigh(A, k=8)$nn[i]))),])
  valor_du_2[i,]<-(A[which(A$continId == ((knearneigh(A, k=8)$nn[i,2]))),])
  valor_du_3[i,]<-(A[which(A$continId == ((knearneigh(A, k=8)$nn[i,3]))),])
  valor_du_4[i,]<-(A[which(A$continId == ((knearneigh(A, k=8)$nn[i,4]))),])
  valor_du_5[i,]<-(A[which(A$continId == ((knearneigh(A, k=8)$nn[i,5]))),])
  valor_du_6[i,]<-(A[which(A$continId == ((knearneigh(A, k=8)$nn[i,6]))),])
  valor_du_7[i,]<-(A[which(A$continId == ((knearneigh(A, k=8)$nn[i,7]))),])
  valor_du_8[i,]<-(A[which(A$continId == ((knearneigh(A, k=8)$nn[i,8]))),])
}
head(valor_du_mixid)
colnames(valor_du_id)<-names(A)
colnames(valor_du_1)<-names(A)
colnames(valor_du_2)<-names(A)
colnames(valor_du_3)<-names(A)
colnames(valor_du_4)<-names(A)
colnames(valor_du_5)<-names(A)
colnames(valor_du_6)<-names(A)
colnames(valor_du_7)<-names(A)
colnames(valor_du_8)<-names(A)

#Calculate mean_du
names(A)
A$mean_du_mix<-((valor_du_id[,"du"]+valor_du_1[,"du"]+valor_du_2[,"du"]+valor_du_3[,"du"]+
                           valor_du_4[,"du"]+valor_du_5[,"du"]+valor_du_6[,"du"]+valor_du_7[,"du"]+
                           valor_du_8[,"du"])/9)
plot(A["mean_du_mix"], pch=16)

#Calculate values of terrain variables attributes to points
A$Cea_0.5_mix<-((valor_du_id[,"Cea_0.5m_1px"]+valor_du_1[,"Cea_0.5m_1px"]+valor_du_2[,"Cea_0.5m_1px"]+valor_du_3[,"Cea_0.5m_1px"]+  #Example of Cea_0.5
                           valor_du_4[,"Cea_0.5m_1px"]+valor_du_5[,"Cea_0.5m_1px"]+valor_du_6[,"Cea_0.5m_1px"]+valor_du_7[,"Cea_0.5m_1px"]+
                           valor_du_8[,"Cea_0.5m_1px"])/9)

#######################################
#From now on is required to consider two datasets, the original with full data; A2 where dead trees are removed
######################################
#Spatial autocorrelation analysis

A2<-A[which(A$Morta == 0),] ## Select on of the 4 properties (A, B, C, D)

#Definition of Spatial Matrix
#Selection of neighbours according to overlapping areas of influence

A2$du_dug<-""
A2$du_dug<-A2$du/(sqrt(sum((A2$du)^2)/length(A2$du)))
A2$du_dug<-as.numeric(A2$du_dug)
A2$influence_area_est<-2.5*(28.502*exp((-66.436+4.201*(A2$du_dug))/(19.817+A2$du))/2) 
plot(A2["influence_area_est"],pch=16)

xyA=st_coordinates(A2)
matrixA<-dist(xyA, "euclidean") #Calculas as distancias em classe dist
matrixA=as.matrix(matrixA, labels=TRUE) #Passa o objecto de classe dist para matriz
colnames(matrixA) <- rownames(matrixA) <- A2$continId
meltA <- melt(matrixA) #desfaz a matrix em colunas apenas
summary(meltA)

melt_0<-melt[which(meltA$value!=0),]  
submelt01<-meltA[melt_0$Var1 %in% M101$continId & melt_0$Var2 %in% A2$continId,]
submelt01<-submelt01[which(submelt01$value<33),]   #minimize the possible neighbours to run faster, where 33 is the max possible distance between tree overlapping areas of influence
submelt01<-submelt01[which(submelt01$value!=0),]
mt<-merge(submelt01, A2[,c("continId","influence_area_est" )], by.x="Var1", by.y="continId") #Criar uma coluna de diametros para a Var1, dependendo do numero da arvore
mtt<-merge(mt, M101[,c("continId","influence_area_est" )], by.x="Var2", by.y="continId")  #Criar uma coluna de diametros para a Var2, dependendo do numero da arvore
names(mtt)<-c(names(mtt)[1:3], "influence_area_est2",  "influence_area_est1") #mudar os nomes das colunas da novas da dataframe
head(mtt)
matrixA<-matrix(nrow=nrow(A2), ncol=nrow(A2))
for (k in c(1:nrow(mtt))){
  valor1<-mtt[k,1]
  valor2<-mtt[k,2]
  linha<- mtt[k,]
  #  matrixA[valor1,valor2]<- ifelse((linha$influence_area_est1 + linha$influence_area_est2 > linha$value),valor2,0)
  matrixA[valor1,valor2]<- ifelse((linha$influence_area_est1 + linha$influence_area_est2 > linha$value),1,0)
  }
matrixA[is.na(matrixA)] = 0
colnames(matrixA) <- rownames(matrixA) <- A2$continId
View(matrixA)
WA<- mat2listw(matrixA)
par(cex=0.5, pch=16)
plot(WA,coord=st_coordinates(A2), col="red")

#Global Moran's I test
moran.test(A2$du, listw=WA,zero.policy=TRUE)

#Correlogram construction and plotting 
 sp.correlogram(Wa, var=A2$du, method="I", order=8)

#Local Moran's I Analysis

  nlist <- knn2nb(knearneigh(A2,k=8))   #Neighbours are defined as the 8 closest neighbours, as it make more sense for local moran tests
  Wk8 <- nb2listw(A2,style="W")
  
  localm<-localmoran(A2$du, listw=Wk8,alternative = "greater")
  A2$localmoran2<-localm[,"Pr(z > 0)"]
  quadrant <- vector(mode="numeric",length=nrow(localm))
  m.qualification <- A2$du - mean(A2$du)  # centers the variable of interest around its mean
  m.local <- localm[,1] - mean(localm[,1])    # centers the local Moran's around the mean
  signif <- 0.1 # significance threshold
  quadrant[m.qualification >0 & m.local>0] <- 4  # builds a data quadrant
  quadrant[m.qualification <0 & m.local<0] <- 1      
  quadrant[m.qualification <0 & m.local>0] <- 2
  quadrant[m.qualification >0 & m.local<0] <- 3
  quadrant[A2$localmoran2>signif] <- 0   
  # plot 
  brks <- c(0,1,2,3,4)
  colors <- c("grey","black","orange","red","darkgreen")
  plot(A2[1],col=colors[findInterval(quadrant,brks,all.inside=FALSE)], pch=c(1,18,16,16,16)[as.factor(quadrant)],main="")
  legend("bottomleft", legend = c("insignificant","low-low","low-high","high-low","high-high"),
        fill=colors,bty="n")
       
#######################################
# Non-spatial modelling

X<-A
st_geometry(X)=NULL
summary(lm(f , data=X))
vif(lm(funcao, data=X))
plot(lm(funcao, data=X))
AIC(lm(funcao , data=X))
  
 f<-as.formula("du~  Cea_1m_1px+  Cea_0.5m_1px+  slope_1px+  Elev_1px+  TRI_1px+  TPI_1px+  cos_aspect_1px+  TWI_1px")        #function for linear models on individual tree diameter for singular sites
 f<-as.formula("du_annual_growth~ Cea_1_mean+  Cea_0.5_mean+  slope_mean+  Elev_mean+  TRI_mean+ cos_aspect_mean+ TWI_mean")  #function for linear models on individual tree diameter for all data
 f<-as.formula("du_mean~  Cea_1m_mean+  Cea_0.5m_1px+  slope_1px+  Elev_1px+  TRI_1px+  TPI_1px+  cos_aspect_1px+  TWI_1px")  #function for linear models on mean tree diameter for singular sites
 f<-as.formula("mean_du_annual_growth~ Cea_1_mean+  Cea_0.5_mean+  slope_mean+  Elev_mean+  TRI_mean+ cos_aspect_mean+ TWI_mean") #function for linear models on mean tree diameter for all data
 
#######################################
### Spatial modelling

duplicated(A$geometry)
 W <- nb2listw(WA$neighbours, style="W",zero.policy = TRUE)
 f <- as.formula("du~Cea_1m_1px+  Cea_0.5m_1px+  slope_1px+  Elev_1px+  TRI_1px+  TPI_1px+  cos_aspect_1px+  TWI_1px")

 
 mod.lag <- lagsarlm(f,data=A2,listw=Wa,zero.policy = TRUE)
 summary(mod.lag)
 
 mod.err <- errorsarlm(f,data=A2,listw=Wa,listw=W,zero.policy = TRUE)
 summary(mod.err)
 
 LR.sarlm(m,mod.lag)
 LR.sarlm(m,mod.err)
 AIC(mod.lag,mod.err)
 
 mod.sac <- sacsarlm(f,data=A2,listw=Wa,zero.policy = TRUE)
 LR.sarlm(mod.sac,mod.lag)
 LR.sarlm(mod.sac,mod.err)
 AIC(mod.lag,mod.err,mod.sac,n)
 summary(mod.sac)
 
 mod.sac2 <- sacsarlm(du ~ Cea_0.5m_1px+slope_1px+TRI_1px+TPI_1px+TWI_1px,data=A2,listw=W,zero.policy = TRUE) ##stepwise removal of variables by the user
 summary(mod.sac2)
 AIC(mod.lag,mod.err,mod.sac, mod.sac2)
 hist(mod.sac2$residuals)
 qqnorm(mod.sac2$residuals)
 qqline(mod.sac2$residuals)
 plot(mod.sac2$residuals ~ fitted(mod.sac2), data=A2)

 

