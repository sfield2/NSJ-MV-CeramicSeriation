
############## NORTHERN SAN JUAN - MESA VERDE CERAMIC SERIATION ################



#### ***NOTE ****
#### If you don't want to build calibration seriation, skip STEP #1.1 & #1.2 
#### and instead import seriated calibration count data (cal) in STEP #1.3



#### ENVIRONMENTS & DATA #######################################################
packages <-c('tidyverse','sp','ggplot2','rgdal','forecast','colorspace','FactoMineR','factoextra','dplyr','ggpubr','lsa','tabula','ca','raster','reshape2')
for(p in packages) if(p %in% rownames(installed.packages()) == F) { install.packages(p) }
for(p in packages) suppressPackageStartupMessages(library(p,quietly=T,character.only=T))

setwd('/Users/seanp/Documents/R/DATING_COMPARISON')
theme_set(theme_bw())

#### data
# calib <- read.csv("http.csv")
# cal <- 
# df <- 

calib <- read.csv("./DATA/VEP_Calibration.csv",header=T,fileEncoding = 'UTF-8-BOM')

df<- read.csv("./DATA/sample_counts.csv",header=T,fileEncoding = 'UTF-8-BOM')

######### STEP #1.1: RUN CORRESPONDENCE ANALYSIS ON CALIBRATION DATA ###########
#### Limit calibration data to NSJ/MV diagnostic types 
calib_diag <- subset(calib, select=c(Chapin.Gray,Moccasin.Gray,Mancos.Gray,
                                     Mancos.Corrugated,Dolores.Corrugated,Mesa.Verde.Corrugated,
                                     Chapin.B.W,Piedra.B.W,Cortez.B.W,
                                     Mancos.B.W,McElmo.B.W,Mesa.Verde.B.W,Abajo.R.O,
                                     Bluff.B.R, Deadmans.B.R))
calib_diag[is.na(calib_diag)] <- 0

#### Run correspondence analysis on limited data
rownames(calib_diag) <- calib$Site.Name
res.ca <- CA(calib_diag, graph = FALSE)
p <- fviz_ca_biplot(res.ca, repel = TRUE,max.overlaps = 15)

#### Run K-means Cluster analysis to determine maximum number
#### of clusters with at least three sites per cluster
row <- get_ca_row(res.ca)
data <- as.data.frame(row$coord)
df <- data[,1:2]
km <- kmeans(df, centers = 10, nstart = 25)
p2<-fviz_cluster(km, data = df)
plot(p2)

#### Group sites within their clusters based on K-Means
data <- p2$data
data$cluster <- as.vector(data$cluster)


#### Example of re-assigning clusters to chronological order
#### Lower right-hand group is earlier in time
calib_groups<- data%>%
  mutate(reassign_cluster = ifelse(cluster %in% 10,1,
                            ifelse(cluster %in% 1, 2,
                            ifelse(cluster %in% 6, 3,
                            ifelse(cluster %in% 7, 4,
                            ifelse(cluster %in% 4, 5,
                            ifelse(cluster %in% 8, 6,
                            ifelse(cluster %in% 9, 7,
                            ifelse(cluster %in% 3, 8,
                            ifelse(cluster %in% 2, 9,
                            ifelse(cluster %in% 5, 10,NA)))))))))))

#### Merge group number with original calibration data
calib_diag$calibration_groups <- calib_groups$reassign_cluster

######### STEP #1.2: CALCULATE TOTAL NUMBER OF SHERDS IN CALIBRATION GROUPS ####
#### Create dataframe for calibration groups
cal<- as.data.frame(matrix(NA,10,15))
colnames(cal) <- colnames(calib_diag)[1:15]

#### Execute loop that sums total number of diagnostic type counts for each group
for(j in 1:10){
  for (i in 1:15){
    cal[j,i]<-sum(calib_diag[,i][calib_diag$calibration_groups == j])
  }
}

######### STEP #1.3: CALCULATE TYPE PROPORTIONS IN CALIBRATION GROUPS ##########
#### Can import seriated calibration count data (cal) here if you want

# cal <- read.csv("http.csv")
cal <- read.csv("./DATA/calibration_counts.csv",header=T,fileEncoding = 'UTF-8-BOM')
cal_allprop <-prop.table(as.matrix(cal),1)*100

######### STEP #2.0: CALCULATE TOTAL NUMBER OFSHERDS IN SAMPLE DATA ############
df[is.na(df)] <- 0
#### Limit survey data to diagnostic types (Wetherill considered as Mancos BW)
df_diag <- subset(df, select=c(Site.Number,Chapin.Gray,Moccasin.Gray,Mancos.Gray,
                                   Mancos.Corrugated,Dolores.Corrugated,Mesa.Verde.Corrugated,
                                   Chapin.B.W,Piedra.B.W,Cortez.B.W,
                                   McElmo.B.W,Mesa.Verde.B.W,Abajo.R.O,
                                   Bluff.B.R,Deadman.s.B.R))%>%
  mutate(Mancos.B.W = df$Mancos.B.W+df$Wetherill.B.W)%>%
  replace(is.na(.),0)%>%
  group_by(Site.Number) %>% 
  summarise_each(funs(sum))


## reorder so types are ordered chronologically by ware
df_diag <- df_diag[, c(1,2,3,4,5,6,7,8,9,10,16,11,12,13,14,15)]

######### STEP #2.1: CALCULATE PROPORTION OF SHERDS IN SAMPLE DATA #############
rownames(df_diag) <-NULL
df_diag <- column_to_rownames(df_diag,"Site.Number")
df_allprop <-prop.table(as.matrix(df_diag),1)*100

######### STEP #3.0: MEASURE DIFFERENCE BETWEEN SURVEY DATA AND CALIBRATION #### 
#### Create Brainerd-Robinson function (Peeples 2011)
BR <- function(x) {
  rd <- dim(x)[1]
  results <- matrix(0,rd,rd)
  for (s1 in 1:rd) {
    for (s2 in 1:rd) {
      x1Temp <- as.numeric(x[s1, ])
      x2Temp <- as.numeric(x[s2, ])
      br.temp <- 0
      results[s1,s2] <- 200 - (sum(abs(x1Temp - x2Temp)))}}
  row.names(results) <- row.names(x)
  colnames(results) <- row.names(x)
  return(results)}

#### Merge calibration and sample data
test <-rbind(cal_allprop,df_allprop)

#### Run BR and extract relevant data
test_br_out <- as.data.frame(BR(test))
test_br_out <- test_br_out[11:nrow(test),1:10]


### normalize
test_br_zscore<- as.data.frame(matrix(NA,nrow(df_diag),nrow(cal)))
colnames(test_br_zscore) <- 1:10
rownames(test_br_zscore) <- rownames(df_diag)

n <- nrow(test_br_zscore)
n2 <- ncol(test_br_zscore)

for(i in 1:n){
  for(j in 1:n2){
    test_br_zscore[i,j] <- (test_br_out[i,j]-rowMeans(test_br_out[i,]))/sd(test_br_out[i,])
  }
}


######### STEP #4.0: TURN INTO TIME SERIES  ####################################
ranges<- c("0700-0799","0800-0879","0880-0919","0920-0979","0980-1019", "1020-1059",
           "1060-1099","1100-1179","1180-1249","1250-1280")
period <- 1:10
ranges <- data.frame(ranges,period)

test_br_zscore <- rownames_to_column(test_br_zscore, "Site.Number")

test_br_zscore_ts<- test_br_zscore%>%
  gather(key="Timeperiod", value = "Value",`1`:`10`)%>%
  mutate(range = ranges$ranges[match(Timeperiod,ranges$period)])%>%
  mutate(start=substr(range,1L,4L))%>%
  mutate(end=substr(range,6L,9L))%>%
  rowwise() %>%
  do(data.frame(Site = .$Site.Number, phase = .$Timeperiod, range = .$range, start =.$start, end= .$end,
                Value=.$Value,Year=seq(.$start,.$end,by=1)))

######### STEP #4.1: IMPORT TREE RINGS  ########################################
tr <- read.csv("./DATA/Treerings.csv",header=T,fileEncoding = 'UTF-8-BOM')


######### STEP #4.2: PLOT
png('/Users/seanp/Documents/R/DATING_COMPARISON/FIGURES/SERIATION_REDEUX/comps2.png',height=1250,width=1800)
ggplot()+
  geom_hline(yintercept = 0,linetype="dashed",alpha=0.35,col="salmon",size=1.25)+
  geom_vline(data=tr,aes(xintercept=Year),linetype="dashed",col="#B3B5FF",size=1.75)+
  geom_line(data=test_br_zscore_ts,aes(Year, Value),col="#4381BB",size=1.75)+
  geom_smooth(data=test_br_zscore_ts,aes(Year,Value),col="#97C0B7B3",size=1.5,se=F)+
  facet_wrap(~Site,scales="free_x")+
    theme(legend.position="none",
          axis.text=element_text(size=25,color="black"),
          axis.title.y=element_text(size=28,color="black"),
          axis.title.x=element_blank(),
          strip.text=element_text(size=32,face="bold",color="black"))
dev.off()



####  CITATIONS

# Peeples, Matthew A.
# 2011 R Script for Calculating the Brainerd-Robinson Coefficient of Similarity 
# and Assessing Sampling Error. Electronic document
# http://www.mattpeeples.net/br.html, accessed March 17, 2023.
