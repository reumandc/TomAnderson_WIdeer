###Climate Data Cleaning and Management Steps

#load location/gps files
noaa.coord<-read.csv("Data/noaa_latlongUpd.csv")
noaa.coord<-noaa.coord[,c(4,5,6,7,8,22:24)]
usda<-read.csv("Data/WI_USDA_Sections.csv")

#rename counties
levels(noaa.coord$COUNTY_N_1)[levels(noaa.coord$COUNTY_N_1)=="Saint Croix"] <- "St. Croix"
levels(noaa.coord$COUNTY_N_1)[levels(noaa.coord$COUNTY_N_1)=="Fond du Lac"] <- "Fond Du Lac"
levels(usda$County)[levels(usda$County)=="Saint Croix"] <- "St. Croix"
levels(usda$County)[levels(usda$County)=="Fond du Lac"] <- "Fond Du Lac"

#Load temp data
temp1<-read.csv("Data/Cty_NOAA_Temp1981_2015_set1.csv",na.strings=c("","unknown","-9999"))
temp2<-read.csv("Data/Cty_NOAA_Temp1981_2015_set2.csv",na.strings=c("","unknown","-9999"))
temp3<-read.csv("Data/Cty_NOAA_Temp1981_2015_set3.csv",na.strings=c("","unknown","-9999"))
temp4<-read.csv("Data/Cty_NOAA_Temp1981_2015_set4.csv",na.strings=c("","unknown","-9999"))
temp5<-read.csv("Data/Cty_NOAA_Temp1981_2015_set5.csv",na.strings=c("","unknown","-9999"))
temp6<-read.csv("Data/Cty_NOAA_Temp1981_2015_set6.csv",na.strings=c("","unknown","-9999"))
temp7<-read.csv("Data/Cty_NOAA_Temp2016_2017.csv",na.strings=c("","unknown","-9999"))

wi.temp<-rbind(temp1,temp2,temp3,temp4,temp5,temp6)
wi.temp<-wi.temp[,c(1:7,12,17,22)]
wi.temp<-rbind(wi.temp,temp7)
wi.temp<-wi.temp[!duplicated(wi.temp),]

#Convert dates, and generate Julian date column
wi.temp$DATE<-gsub("(^\\d{4})(\\d{2})(\\d{2}$)", "\\1-\\2-\\3", wi.temp$DATE)
wi.temp$Jdate<-as.numeric(format(as.Date(as.character(wi.temp$DATE),"%Y-%m-%d"),"%j"))
wi.temp$Year<-as.numeric((format(as.Date(wi.temp$DATE,format="%Y-%m-%d"), "%Y")))
wi.temp$Month<-as.numeric((format(as.Date(wi.temp$DATE,format="%Y-%m-%d"), "%m")))

#Load precip data
pcp81_85<-read.csv("Data/Cty_NOAA_Precip1981_1985.csv",na.strings=c("","unknown","-9999"))#in standard units
pcp81_85$PRCP<-pcp81_85$PRCP*25.4
pcp81_85$SNOW<-pcp81_85$SNOW*25.4
pcp81_85$SNWD<-pcp81_85$SNWD*25.4
pcp86_91<-read.csv("Data/Cty_NOAA_Precip1986_1991.csv",na.strings=c("","unknown","-9999"))
pcp92_97<-read.csv("Data/Cty_NOAA_Precip1992_1997.csv",na.strings=c("","unknown","-9999"))
pcp98_03<-read.csv("Data/Cty_NOAA_Precip1998_2003.csv",na.strings=c("","unknown","-9999"))
pcp04_09<-read.csv("Data/Cty_NOAA_Precip2004_2009.csv",na.strings=c("","unknown","-9999"))
pcp10_15<-read.csv("Data/Cty_NOAA_Precip2010_2015.csv",na.strings=c("","unknown","-9999"))
pcp16_17<-read.csv("Data/Cty_NOAA_Precip2016_2017.csv",na.strings=c("","unknown","-9999"))

wi.pcp<-rbind(pcp81_85,pcp86_91,pcp92_97,pcp98_03,pcp04_09,pcp10_15)  #combine all files into one
wi.pcp<-wi.pcp[,c(1:6,27,32,37)]  #filter to relevant columns,ignores 
wi.pcp<-rbind(wi.pcp,pcp16_17)

#Convert dates, and generate Julian date column
wi.pcp$DATE<-gsub("(^\\d{4})(\\d{2})(\\d{2}$)", "\\1-\\2-\\3", wi.pcp$DATE)
wi.pcp$Jdate<-as.numeric(format(as.Date(as.character(wi.pcp$DATE),"%Y-%m-%d"),"%j"))
wi.pcp$Year<-as.numeric((format(as.Date(wi.pcp$DATE,format="%Y-%m-%d"), "%Y")))
wi.pcp$Month<-as.numeric((format(as.Date(wi.pcp$DATE,format="%Y-%m-%d"), "%m")))

#merge precip and temp
wi.clim<-merge(wi.pcp,wi.temp,by=c("STATION","STATION_NAME","DATE","Jdate","Year","Month"),all=T)
wi.clim<-merge(wi.clim,noaa.coord[,c(1,8)],by=c("STATION"),all.x=T)
names(wi.clim)[20]<-"COUNTY_NAM"
sum(is.na(wi.clim$COUNTY_NAM))

noaa.latlong<-wi.clim[!duplicated(wi.clim$STATION),] #get unique lat longs for each station

#determine data sample sizes
if(samplesize.flag=="y"){
  source('Code/testSampleSize.R')
}

#assign counties to stations on the borders/in Lake MI
wi.clim$COUNTY_NAM[wi.clim$STATION=="GHCND:USC00110203"]<-"Kenosha"
wi.clim$COUNTY_NAM[wi.clim$STATION=="GHCND:USC00477464"]<-"Polk"
wi.clim$COUNTY_NAM[wi.clim$STATION=="GHCND:USC00473038"]<-"Vernon"
wi.clim$COUNTY_NAM[wi.clim$STATION=="GHCND:USC00476922"]<-"Racine"
wi.clim$COUNTY_NAM[wi.clim$STATION=="GHCND:USC00476764"]<-"Ozaukee"
wi.clim$COUNTY_NAM[wi.clim$STATION=="GHCND:USC00478672"]<-"Manitowoc"
wi.clim$COUNTY_NAM[wi.clim$STATION=="GHCND:USC00478589"]<-"Trempealeau"
wi.clim$COUNTY_NAM[wi.clim$STATION=="GHCND:USC00214418"]<-"La Crosse"
wi.clim$COUNTY_NAM[wi.clim$STATION=="GHCND:USC00477121"]<-"Marathon"

#average over all data points per county by Julian day
cty.min<-aggregate(TMIN~COUNTY_NAM+Year+Month+Jdate,wi.clim,FUN="mean")
cty.max<-aggregate(TMAX~COUNTY_NAM+Year+Month+Jdate,wi.clim,FUN="mean")
cty.pcp<-aggregate(PRCP~COUNTY_NAM+Year+Month+Jdate,wi.clim,FUN="mean")
cty.snwd<-aggregate(SNWD~COUNTY_NAM+Year+Month+Jdate,wi.clim,FUN="mean")
cty.clim<-merge(cty.min,cty.max,by=c("COUNTY_NAM","Year","Month","Jdate"),all=T)
cty.clim<-merge(cty.clim,cty.pcp,by=c("COUNTY_NAM","Year","Month","Jdate"),all=T)
cty.clim<-merge(cty.clim,cty.snwd,by=c("COUNTY_NAM","Year","Month","Jdate"),all=T)

#drop Menomonie county- very little deer or climate data
cty.clim<-cty.clim[!(cty.clim$COUNTY_NAM=="Menominee"),]

#average data over months
cty.month<-ddply(cty.clim,c("COUNTY_NAM","Year","Month"),.fun = summarise,Tmin=mean(TMIN,na.rm=T),Tmax=mean(TMAX,na.rm=T),Pcp=mean(PRCP,na.rm=T),Snwd=mean(SNWD,na.rm=T),.drop=F)
cty.month<-cty.month[!(cty.month$Year==2017 & cty.month$Month>5),]

#drop Menomonie county- very little deer or climate data
cty.month<-cty.month[!(cty.month$COUNTY_NAM=="Menominee"),]

#Merge daily climate data by county with usda district data
usda.clim<-merge(cty.clim,usda,by.x="COUNTY_NAM",by.y="County",all=T)
#Merge monthly climate data by county with usda district data
usda.month<-merge(cty.month,usda,by.x="COUNTY_NAM",by.y="County",all=T)

# County Climate Data -----------------------------------------------------
if(scale.flag=="county"||scale.flag=="both"){

# openingday<-read.csv("Data/OpeningDay.csv")
# openingday<-openingday[order(openingday$Year),]
# openingday$Jdate<-as.numeric(format(as.Date(as.character(openingday$Date),"%m/%d/%Y"),"%j"))
# openday.dat<-openingday$Jdate[openingday$Year%in%c(minyear:maxyear)]
# 
# #Hunting Season Weather
# years<-unique(minyear:maxyear)
# hunt.clim<-list()
# mean.hunt<-matrix(NA,ncol=length(unique(cty.clim$COUNTY_NAM)),nrow=length(minyear:maxyear))
# ctynames<-sort(unique(cty.clim$COUNTY_NAM))
# respvars<-names(cty.clim)[5:8]
# for(resp in respvars){
#   for(ctyname in ctynames){
#     for(year in years){
#       tmp<-cty.clim[,resp][cty.clim$COUNTY_NAM==ctyname & 
#                              cty.clim$Year==year & 
#                              #as.numeric(cty.clim$Jdate)%in%c(openday.dat[year-1980]:(openday.dat[year-1980]+8))]
#                              as.numeric(cty.clim$Jdate)%in%c(openday.dat[year-min(years)+1]:(openday.dat[year-min(years)+1]+8))]
#       if(length(tmp)==0){
#         mean.hunt[years%in%year,ctynames%in%ctyname]<-NA
#       }
#       else{
#         mean.hunt[years%in%year,ctynames%in%ctyname]<-mean(tmp,na.rm=T)
#       }
#     }
#     hunt.clim[[resp]]<-mean.hunt
#   }
# }
# 
# mean.day<-round(mean(openday.dat),0)
# fixedtmax<-matrix(NA,ncol=length(ctynames),nrow=length(years))
# for(ctyname in ctynames){
#   for(year in years){
#     tmp<-cty.clim[,'TMAX'][cty.clim$COUNTY_NAM==ctyname & 
#                            cty.clim$Year==year & 
#                            as.numeric(cty.clim$Jdate)%in%c(mean.day:(mean.day+8))]
# 
#     if(length(tmp)==0){
#       fixedtmax[years%in%year,ctynames%in%ctyname]<-NA
#     }
#     else{
#       fixedtmax[years%in%year,ctynames%in%ctyname]<-mean(tmp,na.rm=T)
#     }
#   }
# }
# hunt.clim[[5]]<-fixedtmax

# hunt.clim.dt<-list()
# for(j in 1:length(hunt.clim)){
#   tmp.mat<-matrix(NA,nrow=length(minyear:maxyear),ncol=length(unique(cty.clim$COUNTY_NAM)))
#   for(i in 1:dim(hunt.clim[[1]])[2]){
#     if(sum(is.na(hunt.clim[[j]][,i]))>0.5*dim(hunt.clim[[j]])[1]){
#       tmp.mat[,i]<-NA
#     }
#     else{
#       tmp<-resid(lm(hunt.clim[[j]][,i]~years,na.action=na.exclude))
#       tmp<-tmp/sd(tmp,na.rm=T)
#       tmp.mat[,i]<-tmp
#     }
#   }
#   hunt.clim.dt[[j]]<-t(tmp.mat)
# }
# names(hunt.clim.dt)<-c("TMIN","TMAX","PRCP","SNWD","FixedTMAX")
# 
# #generate subsets of raw data for southern and northern counties
# hunt.clim<-lapply(hunt.clim,function(x){colnames(x)<-unique(cty.clim$COUNTY_NAM);x})
# hunt.clim<-lapply(hunt.clim,function(x){x[,order(colnames(x))]})
# hunt.climN<-lapply(hunt.clim,function(x){x[,northern.cty]})
# hunt.climS<-lapply(hunt.clim,function(x){x[,!(colnames(x)%in%northern.cty)]})
# 
# #generate subsets of detrended data for southern and northern counties
# hunt.clim.dt<-lapply(hunt.clim.dt,function(x){row.names(x)<-unique(cty.clim$COUNTY_NAM);x})
# hunt.clim.dt<-lapply(hunt.clim.dt,function(x){x[order(row.names(x)),]})
# hunt.clim.dtN<-lapply(hunt.clim.dt,function(x){x[northern.cty,]})
# hunt.clim.dtS<-lapply(hunt.clim.dt,function(x){x[southern.cty,]})

#NAs occur where the data doesn't exist
#NaNs are the data is missing for that county/year for the hunting season (i.e. mean of 10 NAs)

#Winter weather conditions
winter.tmp<-matrix(NA, ncol=length(unique(cty.month$COUNTY_NAM)),nrow=length(minyear:maxyear))
winter.clim<-list()
years<-minyear:maxyear
ctynames<-sort(unique(cty.month$COUNTY_NAM))
respvars<-names(cty.month)[4:7]

for(resp in respvars){
  for(ctyname in ctynames){      
    for(year in years){
      tmp<-cty.month[cty.month$COUNTY_NAM==ctyname & cty.month$Year==year & cty.month$Month<=3,] 
      tmp<-rbind(tmp, cty.month[cty.month$COUNTY_NAM==ctyname & cty.month$Year==(year-1) & cty.month$Month>=12,])
      
      if(length(tmp)==0){
        winter.tmp[years %in% year,ctynames %in% ctyname]<-NA
      }
      else{
        tmp1<-mean(tmp[,resp],na.rm=T)
        winter.tmp[years %in% year,ctynames %in% ctyname]<-tmp1
      }
    }
  }
  winter.clim[[resp]]<-winter.tmp
}

#calculate winter severity index (wsi) for each location
#wsi is the sum of number of days below 0 temp, and number of days with greater than 18" (45cm) snow
wsi.mat<-matrix(NA, ncol=length(unique(cty.month$COUNTY_NAM)),nrow=length(minyear:maxyear))
years<-minyear:maxyear
ctynames<-sort(unique(cty.month$COUNTY_NAM))

for(ctyname in ctynames){      
  for(year in years){
    tmp<-cty.clim[cty.clim$COUNTY_NAM==ctyname & cty.clim$Year==year & cty.clim$Month<=4,] 
    tmp<-rbind(tmp, cty.clim[cty.clim$COUNTY_NAM==ctyname & cty.clim$Year==(year-1) & cty.clim$Month>=12,])

    if(nrow(tmp)==0){
      wsi.mat[years %in% year,ctynames %in% ctyname]<-NA
    }
    else{
      if(sum(is.na(tmp$TMIN))>(0.10*dim(tmp)[1])||sum(is.na(tmp$SNWD))>(0.10*dim(tmp)[1])){
        wsi.mat[years %in% year,ctynames %in% ctyname]<-NA
      }
      else{
        air<-sum(tmp$TMIN<(-17.778),na.rm=T)
        snow<-sum(tmp$SNWD/10>45.72,na.rm=T)
        wsi<-air+snow
        wsi.mat[years %in% year,ctynames %in% ctyname]<-wsi
      }
    }
  }
}
colnames(wsi.mat)<-ctynames

#impute mean value if the number of NAs in a county is 2 or less
wsi.mat1<-matrix(wsi.mat,ncol=ncol(wsi.mat),nrow=nrow(wsi.mat))
for(i in 1:ncol(wsi.mat1)){
  if(sum(is.na(wsi.mat1[,i]))<3){
    wsi.mat1[is.na(wsi.mat1[,i]), i] <- mean(wsi.mat1[,i], na.rm = TRUE)
  }
}
winter.clim[[5]]<-wsi.mat1
names(winter.clim)<-c("Tmin","Tmax","Prcp","Snwd","WSI")
# #detrend WSI (uses values for imputed counties)
# wsi.mat.dt<-matrix(NA,nrow=nrow(wsi.mat1),ncol=ncol(wsi.mat1))
# for(i in 1:dim(wsi.mat1)[2]){
#   if(sum(is.na(wsi.mat1[,i]))>0.5*dim(wsi.mat1)[1]){
#     wsi.mat.dt[,i]<-NA
#   }
#   else{
#     tmp<-resid(lm(wsi.mat1[,i]~years,na.action=na.exclude))
#     tmp<-tmp/sd(tmp,na.rm=T)
#     wsi.mat.dt[,i]<-tmp
#   }
# }
# wsi.mat.dt<-t(wsi.mat.dt)
# winter.clim.dt[[5]]<-wsi.mat.dt #add to winter climate list

spring.tmp<-matrix(NA, ncol=length(unique(cty.month$COUNTY_NAM)),nrow=length(minyear:maxyear))
spring.clim<-list()
years<-minyear:maxyear
ctynames<-sort(unique(cty.month$COUNTY_NAM))
respvars<-names(cty.month)[4:7]

for(resp in respvars){
  for(ctyname in ctynames){      
    for(year in years){
      tmp<-cty.month[cty.month$COUNTY_NAM==ctyname & cty.month$Year==year & cty.month$Month>=3 & cty.month$Month<=5,] 

      if(length(tmp)==0){
        spring.tmp[years %in% year,ctynames %in% ctyname]<-NA
      }
      else{
        tmp1<-mean(tmp[,resp],na.rm=T)
        spring.tmp[years %in% year,ctynames %in% ctyname]<-tmp1
      }
    }
  }
  spring.clim[[resp]]<-spring.tmp
}
names(spring.clim)<-c("Tmin","Tmax","Prcp","Snwd")

if(detr.flag=="norm"){
  winter.clim.dt<-lapply(winter.clim,function(x){x<-reumannplatz::CleanData(t(x))$cleandat;x})
  spring.clim.dt<-lapply(spring.clim,function(x){x<-reumannplatz::CleanData(t(x))$cleandat;x})
  #hunt.clim.dt<-lapply(hunt.clim,function(x){x<-reumannplatz::CleanData(t(x))$cleandat;x})
  winter.clim.dt<-lapply(winter.clim.dt,function(x){rownames(x)<-colnames(wsi.mat);x})
  #hunt.clim.dt<-lapply(hunt.clim.dt,function(x){rownames(x)<-colnames(wsi.mat);x})
  names(winter.clim.dt)<-c("Tmin","Tmax","Prcp","Snwd","WSI")
  #names(hunt.clim.dt)<-c("Tmin","Tmax","Prcp","Snwd","FixedTmax")
}

if(detr.flag=="detr"){
  winter.clim.dt<-lapply(winter.clim,function(x){x<-reumannplatz::CleanData(t(x),normalize=F)$cleandat;x})
  spring.clim.dt<-lapply(spring.clim,function(x){x<-reumannplatz::CleanData(t(x),normalize=F)$cleandat;x})
  #hunt.clim.dt<-lapply(hunt.clim,function(x){x<-reumannplatz::CleanData(t(x),normalize=F)$cleandat;x})
  winter.clim.dt<-lapply(winter.clim.dt,function(x){rownames(x)<-colnames(wsi.mat);x})
  #hunt.clim.dt<-lapply(hunt.clim.dt,function(x){rownames(x)<-colnames(wsi.mat);x})
  names(winter.clim.dt)<-c("Tmin","Tmax","Prcp","Snwd","WSI")
  #names(hunt.clim.dt)<-c("Tmin","Tmax","Prcp","Snwd","FixedTmax")
}

if(detr.flag=="none"){
  winter.clim<-lapply(winter.clim,t)
  spring.clim<-lapply(spring.clim,t)
  #hunt.clim<-lapply(hunt.clim,t)
  winter.clim<-lapply(winter.clim,function(x){rownames(x)<-colnames(wsi.mat);x})
  #hunt.clim<-lapply(hunt.clim,function(x){rownames(x)<-colnames(wsi.mat);x})
  names(winter.clim)<-c("Tmin","Tmax","Prcp","Snwd","WSI")
  #names(hunt.clim)<-c("Tmin","Tmax","Prcp","Snwd","FixedTmax")
}

#detrend climate data
# winter.clim.dt<-list()
# for(j in 1:length(winter.clim)){
#   tmp.mat<-matrix(NA,nrow=length(minyear:maxyear),ncol=length(unique(cty.clim$COUNTY_NAM)))
#   for(i in 1:dim(winter.clim[[1]])[2]){
#     if(sum(is.na(winter.clim[[j]][,i]))>0.5*dim(winter.clim[[j]])[1]){
#       tmp.mat[,i]<-NA
#     }
#     else{
#       tmp<-resid(lm(winter.clim[[j]][,i]~years,na.action=na.exclude))
#       tmp<-tmp/sd(tmp,na.rm=T)
#       tmp.mat[,i]<-tmp
#     }
#   }
#   winter.clim.dt[[j]]<-t(tmp.mat)
# }
# names(winter.clim.dt)<-c("Tmin","Tmax","Prcp","Snwd","WSI")

#Generate subsets of winter weather for northern (N =18 counties) and southern WI (N =54 counties)
# winter.clim<-lapply(winter.clim,function(x){colnames(x)<-colnames(wsi.mat);x})
# winter.clim<-lapply(winter.clim,function(x){x[,order(colnames(x))]})
# winter.climN<-lapply(winter.clim,function(x){x[,northern.cty]})
# winter.climS<-lapply(winter.clim,function(x){x[,southern.cty]})
# 
# #Generate subsets of detrended winter weather for northern (N =18 counties) and southern WI (N =54 counties)
# winter.clim.dt<-lapply(winter.clim.dt,function(x){rownames(x)<-colnames(wsi.mat);x})
# winter.clim.dt<-lapply(winter.clim.dt,function(x){x[order(row.names(x)),]})
# winter.clim.dtN<-lapply(winter.clim.dt,function(x){x[northern.cty,]})
# winter.clim.dtS<-lapply(winter.clim.dt,function(x){x[southern.cty,]})

# apply(wsi.mat1,2,function(x){sum(is.na(x))})
# sum(!is.na(colMeans(wsi.mat1)))  #31 counties with < 10% missing, could add 11 counties
# sum(!is.na(rowMeans(wsi.mat)))  #31 counties with < 10% missing, could add 11 counties

}
#Notes
#counties with very incomplete data: st. croix. (2373- missing temp, not precip data),pepinX(3237),langlade. (9963),adamsX (4781)
#counites with some missingness: Iowa. (11601), Lafayette. (11553), KewauneeX (11736),marquette.(11870),eau claire. (11445)
#64 counties have complete time series

# USDA Climate Data -------------------------------------------------------
if(scale.flag=="usda"||scale.flag=="both"){
#Winter conditions for USDA districts
winter.tmp<-matrix(NA, ncol=length(unique(usda.month$Zone)),nrow=length(minyear:maxyear))
winter.clim.usda<-list()
years<-minyear:maxyear
usda.names<-sort(unique(usda.month$Zone))
respvars<-names(usda.month)[4:7]

for(resp in respvars){
  for(usdaname in usda.names){      
    for(year in years){
      tmp<-usda.month[usda.month$Zone==usdaname & usda.month$Year==year & usda.month$Month<=3,] 
      tmp<-rbind(tmp, usda.month[usda.month$Zone==usdaname & usda.month$Year==(year-1) & usda.month$Month>=12,])
      
      if(length(tmp)==0){
        winter.tmp[years %in% year,usda.names %in% usdaname]<-NA
      }
      else{
        tmp1<-mean(tmp[,resp],na.rm=T)
        winter.tmp[years %in% year,usda.names %in% usdaname]<-tmp1
      }
    }
  }
  winter.clim.usda[[resp]]<-winter.tmp
}

#calculate winter severity index (wsi) for each usda district
#wsi is the sum of number of days below 0 temp, and number of days with greater than 18" (45cm) snow
wsi.usda.mat<-matrix(NA, ncol=length(unique(usda.clim$Zone)),nrow=length(minyear:maxyear))
years<-minyear:maxyear
usda.names<-sort(unique(usda$Zone))

for(usdaname in usda.names){      
  for(year in years){
    tmp<-usda.clim[usda.clim$Zone==usdaname & usda.clim$Year==year & usda.clim$Month<=4,] 
    tmp<-rbind(tmp, usda.clim[usda.clim$Zone==usdaname & usda.clim$Year==(year-1) & usda.clim$Month>=12,])
    tmp<-aggregate(cbind(TMIN,TMAX,PRCP,SNWD)~Year+Month+Jdate,FUN="mean",data=tmp,na.rm=T) #average values across district
    if(nrow(tmp)==0){
      wsi.mat[years %in% year,usda.names %in% usdaname]<-NA
    }
    else{
      if(sum(is.na(tmp$TMIN))>(0.10*dim(tmp)[1])||sum(is.na(tmp$SNWD))>(0.10*dim(tmp)[1])){
        wsi.usda.mat[years %in% year,usda.names %in% usdaname]<-NA
      }
      else{
        air<-sum(tmp$TMIN<(-17.778),na.rm=T)
        snow<-sum(tmp$SNWD/10>45.72,na.rm=T)
        wsi<-air+snow
        wsi.usda.mat[years %in% year,usda.names %in% usdaname]<-wsi
      }
    }
  }
}
colnames(wsi.usda.mat)<-usda.names
winter.clim.usda[[5]]<-wsi.usda.mat
names(winter.clim.usda)<-c("Tmin","Tmax","Prcp","Snwd","WSI")

#detrend climate data
# wint.clim.usda.dt<-list()
# for(j in 1:length(winter.clim.usda)){
#   tmp.mat<-matrix(NA,nrow=length(minyear:maxyear),ncol=length(unique(usda.clim$Zone)))
#   for(i in 1:dim(winter.clim.usda[[1]])[2]){
#     if(sum(is.na(winter.clim.usda[[j]][,i]))>0.5*dim(winter.clim.usda[[j]])[1]){
#       tmp.mat[,i]<-NA
#     }
#     else{
#       tmp<-resid(lm(winter.clim.usda[[j]][,i]~years,na.action=na.exclude))
#       tmp<-tmp/sd(tmp,na.rm=T)
#       tmp.mat[,i]<-tmp
#     }
#   }
#   wint.clim.usda.dt[[j]]<-t(tmp.mat)
# }
# names(wint.clim.usda.dt)<-names(winter.clim.usda)

# #detrend WSI
# wsi.usda.mat.dt<-matrix(NA,nrow=nrow(wsi.usda.mat),ncol=ncol(wsi.usda.mat))
# for(i in 1:dim(wsi.usda.mat)[2]){
#   if(sum(is.na(wsi.usda.mat[,i]))>0.5*dim(wsi.usda.mat)[1]){
#     wsi.usda.mat.dt[,i]<-NA
#   }
#   else{
#     tmp<-resid(lm(wsi.usda.mat[,i]~years,na.action=na.exclude))
#     tmp<-tmp/sd(tmp,na.rm=T)
#     wsi.usda.mat.dt[,i]<-tmp
#   }
# }
# wsi.usda.mat.dt<-t(wsi.usda.mat.dt)
# wint.clim.usda.dt[[5]]<-wsi.usda.mat.dt
#names(wint.clim.usda.dt)<-c("Tmin","Tmax","Prcp","Snwd","WSI")

#Hunting conditions for USDA districts
# years<-unique(minyear:maxyear)
# hunt.clim.usda<-list()
# mean.hunt<-matrix(NA,ncol=length(unique(usda.clim$Zone)),nrow=length(minyear:maxyear))
# usda.names<-sort(unique(usda.clim$Zone))
# respvars<-names(usda.clim)[5:8]
# for(resp in respvars){
#   for(usdaname in usda.names){
#     for(year in years){
#       tmp<-usda.clim[,resp][usda.clim$Zone==usdaname & 
#                               usda.clim$Year==year & 
#                               as.numeric(usda.clim$Jdate)%in%c(openday.dat[year-min(years)+1]:(openday.dat[year-min(years)+1]+8))]
#       if(length(tmp)==0){
#         mean.hunt[years%in%year,usda.names%in%usdaname]<-NA
#       }
#       else{
#         mean.hunt[years%in%year,usda.names%in%usdaname]<-mean(tmp,na.rm=T)
#       }
#     }
#     hunt.clim.usda[[resp]]<-mean.hunt
#   }
# }
# 
# mean.day<-round(mean(openday.dat),0)
# fixedtmax<-matrix(NA,ncol=length(usda.names),nrow=length(years))
# for(usdaname in usda.names){
#   for(year in years){
#     tmp<-usda.clim[,'TMAX'][usda.clim$Zone==usdaname & 
#                               usda.clim$Year==year & 
#                              as.numeric(usda.clim$Jdate)%in%c(mean.day:(mean.day+8))]
#     
#     if(length(tmp)==0){
#       fixedtmax[years%in%year,usda.names%in%usdaname]<-NA
#     }
#     else{
#       fixedtmax[years%in%year,usda.names%in%usdaname]<-mean(tmp,na.rm=T)
#     }
#   }
# }
# hunt.clim.usda[[5]]<-fixedtmax


# hunt.clim.usda.dt<-list()
# for(j in 1:length(hunt.clim.usda)){
#   tmp.mat<-matrix(NA,nrow=length(minyear:maxyear),ncol=length(unique(usda.clim$Zone)))
#   for(i in 1:dim(hunt.clim.usda[[1]])[2]){
#     if(sum(is.na(hunt.clim.usda[[j]][,i]))>0.5*dim(hunt.clim.usda[[j]])[1]){
#       tmp.mat[,i]<-NA
#     }
#     else{
#       tmp<-resid(lm(hunt.clim.usda[[j]][,i]~years,na.action=na.exclude))
#       tmp<-tmp/sd(tmp,na.rm=T)
#       tmp.mat[,i]<-tmp
#     }
#   }
#   hunt.clim.usda.dt[[j]]<-t(tmp.mat)
# }
# names(hunt.clim.usda.dt)<-c("Tmin","Tmax","Prcp","Snwd","FixedTmax")

if(detr.flag=="norm"){
  winter.clim.usda<-lapply(winter.clim.usda,function(x){x<-reumannplatz::CleanData(t(x))$cleandat;x})
  #hunt.clim.usda<-lapply(hunt.clim.usda,function(x){x<-reumannplatz::CleanData(t(x))$cleandat;x})
  names(winter.clim.usda)<-c("Tmin","Tmax","Prcp","Snwd","WSI")
  #names(hunt.clim.usda)<-c("Tmin","Tmax","Prcp","Snwd","FixedTmax")
}

if(detr.flag=="detr"){
  winter.clim.usda<-lapply(winter.clim.usda,function(x){x<-reumannplatz::CleanData(t(x),normalize=F)$cleandat;x})
  #hunt.clim.usda<-lapply(hunt.clim.usda,function(x){x<-reumannplatz::CleanData(t(x),normalize=F)$cleandat;x})
  names(winter.clim.usda)<-c("Tmin","Tmax","Prcp","Snwd","WSI")
  #names(hunt.clim.usda)<-c("Tmin","Tmax","Prcp","Snwd","FixedTmax")
}

if(detr.flag=="none"){
  winter.clim.usda<-lapply(winter.clim.usda,function(x){x<-t(x);x})
  #hunt.clim.usda<-lapply(hunt.clim.usda,function(x){x<-t(x);x})
  names(winter.clim.usda)<-c("Tmin","Tmax","Prcp","Snwd","WSI")
  #names(hunt.clim.usda)<-c("Tmin","Tmax","Prcp","Snwd","FixedTmax")
}

}
# Large Scale Climate Patterns --------------------------------------------

#load large climate indice data
climate_indices <- read.csv("Data/climate_indices.csv")
mei<-read.csv("Data/mei.csv")
pdo<-read.csv("Data/pdo.csv")
npi<-read.csv("Data/npi.csv",na.strings = c("-999"))
enso <- read.csv("Data/ENSO.csv")

#compute yearly averages
mei$Mean<-rowMeans(mei[,2:dim(mei)[2]],na.rm=T)
npi$Mean<-rowMeans(npi[,2:dim(npi)[2]],na.rm=T)
npi$MeanZ<-scale(npi$Mean)#double check that this is right for z transformation
mean.clim<-aggregate(cbind(NAO,EA,WP,EP_NP,PNA,EA.WR)~Year,FUN="mean",data=climate_indices)
mean.pdo<-aggregate(Value~Year,pdo,FUN="mean")

#calculate winter means
## 1. Mean Winter NAO (Dec. through Mar.)
winter.nao<-rep(NA, length(minyear:maxyear))
years<-unique(climate_indices$Year[climate_indices$Year%in%c(minyear:maxyear)])
for(year in years){
  tmp<-climate_indices[climate_indices$Year==year & climate_indices$Month<=3,]
  tmp<-rbind(tmp, climate_indices[climate_indices$Year==(year-1) & climate_indices$Month>=12,])
  winter.nao[years %in% year]<-mean(tmp$NAO,na.rm=T)
}
winter.nao

## 2. Mean Winter PDO (Dec. through Mar.)
winter.pdo<-rep(NA, length(minyear:maxyear))
years<-unique(pdo$Year[pdo$Year%in%c(minyear:maxyear)])
for(year in years){
  tmp<-pdo[pdo$Year==year & pdo$Month<=3,]
  tmp<-rbind(tmp, pdo[pdo$Year==(year-1) & pdo$Month>=12,])
  winter.pdo[years %in% year]<-mean(tmp$Value,na.rm=T)
}
winter.pdo

## 3. Mean Winter MEI (Dec. through Mar.)
mei.long<-melt(mei[,1:13],id.vars=c("YEAR"),variable.name = "MONTHS",value.name = "MEI")
mei.long$Index<-rep(1:12,each=67)
winter.mei<-rep(NA, length(minyear:maxyear))
years<-unique(mei.long$YEAR[mei.long$YEAR%in%c(minyear:maxyear)])
for(year in years){
  tmp<-mei.long[mei.long$YEAR==year & mei.long$Index<=3,]
  tmp<-rbind(tmp, mei.long[mei.long$YEAR==(year-1) & mei.long$Index>=12,])
  winter.mei[years %in% year]<-mean(tmp$MEI,na.rm=T)
}

## 3. Mean Winter MEI (Dec. through Mar.)
enso.long<-melt(enso[,1:13],id.vars=c("Year"),variable.name = "Month",value.name = "enso")
enso.long$Index<-rep(1:12,each=147)
winter.enso<-rep(NA, length(minyear:maxyear))
years<-unique(enso.long$Year[enso.long$Year%in%c(minyear:maxyear)])
for(year in years){
  tmp<-enso.long[enso.long$Year==year & enso.long$Index<=3,]
  tmp<-rbind(tmp, enso.long[enso.long$Year==(year-1) & enso.long$Index>=12,])
  winter.enso[years %in% year]<-mean(tmp$enso,na.rm=T)
}

## 1. Mean summer NAO (Dec. through Mar.)
summer.nao<-rep(NA, length(minyear:maxyear))
years<-unique(climate_indices$Year[climate_indices$Year%in%c(minyear:maxyear)])
for(year in years){
  tmp<-climate_indices[climate_indices$Year==year & climate_indices$Month%in%c(4:11),]
  summer.nao[years %in% year]<-mean(tmp$NAO,na.rm=T)
}

## 2. Mean summer PDO (April through November)
summer.pdo<-rep(NA, length(minyear:maxyear))
years<-unique(pdo$Year[pdo$Year%in%c(minyear:maxyear)])
for(year in years){
  tmp<-pdo[pdo$Year==year & pdo$Month%in%c(4:11),]
  summer.pdo[years %in% year]<-mean(tmp$Value,na.rm=T)
}

## 3. Mean summer MEI (April through November)
mei.long<-melt(mei[,1:13],id.vars=c("YEAR"),variable.name = "MONTHS",value.name = "MEI")
mei.long$Index<-rep(1:12,each=67)
summer.mei<-rep(NA, length(minyear:maxyear))
years<-unique(mei.long$YEAR[mei.long$YEAR%in%c(minyear:maxyear)])
for(year in years){
  tmp<-mei.long[mei.long$YEAR==year & mei.long$Index%in%c(4:11),]
  summer.mei[years %in% year]<-mean(tmp$MEI,na.rm=T)
}

summer.enso<-rep(NA, length(minyear:maxyear))
years<-unique(enso.long$Year[enso.long$Year%in%c(minyear:maxyear)])
for(year in years){
  tmp<-enso.long[enso.long$Year==year & enso.long$Index%in%c(4:11),]
  summer.enso[years %in% year]<-mean(tmp$enso,na.rm=T)
}

#put into repeating matrices
if(scale.flag=="usda"||scale.flag=="both"){
  #nao.mat<-matrix(mean.clim$NAO[mean.clim$Year%in%c(minyear:maxyear)],ncol=length(minyear:maxyear),nrow=length(unique(usda.clim$Zone)),byrow=T)
  #pdo.mat<-matrix(mean.pdo$Value[mean.pdo$Year%in%c(minyear:maxyear)],ncol=length(minyear:maxyear),nrow=length(unique(usda.clim$Zone)),byrow=T)
  #mei.mat<-matrix(mei$Mean[mei$YEAR%in%c(minyear:maxyear)],ncol=length(minyear:maxyear),nrow=length(unique(usda.clim$Zone)),byrow=T)
  #winter indices
  win.nao.mat<-matrix(winter.nao[1:length(winter.nao)],ncol=length(minyear:maxyear),nrow=length(unique(usda.clim$Zone)),byrow=T)
  win.pdo.mat<-matrix(winter.pdo[1:length(winter.pdo)],ncol=length(minyear:maxyear),nrow=length(unique(usda.clim$Zone)),byrow=T)
  win.mei.mat<-matrix(winter.mei[1:length(winter.mei)],ncol=length(minyear:maxyear),nrow=length(unique(usda.clim$Zone)),byrow=T)
  win.enso.mat<-matrix(winter.enso[1:length(winter.enso)],ncol=length(minyear:maxyear),nrow=length(unique(usda.clim$Zone)),byrow=T)
  
  #summer indices
  sum.nao.mat<-matrix(summer.nao[1:length(summer.nao)],ncol=length(minyear:maxyear),nrow=length(unique(usda.clim$Zone)),byrow=T)
  sum.pdo.mat<-matrix(summer.pdo[1:length(summer.pdo)],ncol=length(minyear:maxyear),nrow=length(unique(usda.clim$Zone)),byrow=T)
  sum.mei.mat<-matrix(summer.mei[1:length(summer.mei)],ncol=length(minyear:maxyear),nrow=length(unique(usda.clim$Zone)),byrow=T)
  sum.enso.mat<-matrix(summer.enso[1:length(summer.enso)],ncol=length(minyear:maxyear),nrow=length(unique(usda.clim$Zone)),byrow=T)
  climindex.usda<-list(win.nao.mat=win.nao.mat,win.pdo.mat=win.pdo.mat,win.mei.mat=win.mei.mat,win.enso.mat=win.enso.mat,
                  sum.nao.mat=sum.nao.mat,sum.pdo.mat=sum.pdo.mat,sum.mei.mat=sum.mei.mat,sum.enso.mat=sum.enso.mat)#nao.mat=nao.mat,pdo.mat=pdo.mat,mei.mat=mei.mat,
  if(detr.flag=="detr"){
    climindex.usda<-lapply(climindex.usda,function(x){x<-reumannplatz::CleanData(x,normalize=F)$cleandat;x})
  }
  if(detr.flag=="norm"){
    climindex.usda<-lapply(climindex.usda,function(x){x<-reumannplatz::CleanData(x)$cleandat;x})
  }
  names(climindex.usda)<-c("WinterNAO","WinterPDO","WinterMEI","WinterENSO","SummerNAO","SummerPDO","SummerMEI","SummerENSO")
}
if(scale.flag=="county"||scale.flag=="both"){
  #nao.mat<-matrix(mean.clim$NAO[mean.clim$Year%in%c(minyear:maxyear)],ncol=length(minyear:maxyear),nrow=length(unique(cty.clim$COUNTY_NAM)),byrow=T)
  #pdo.mat<-matrix(mean.pdo$Value[mean.pdo$Year%in%c(minyear:maxyear)],ncol=length(minyear:maxyear),nrow=length(unique(cty.clim$COUNTY_NAM)),byrow=T)
  #mei.mat<-matrix(mei$Mean[mei$YEAR%in%c(minyear:maxyear)],ncol=length(minyear:maxyear),nrow=length(unique(cty.clim$COUNTY_NAM)),byrow=T)
  
  #winter indices
  win.nao.mat<-matrix(winter.nao[1:length(winter.nao)],ncol=length(minyear:maxyear),nrow=length(unique(cty.clim$COUNTY_NAM)),byrow=T)
  win.pdo.mat<-matrix(winter.pdo[1:length(winter.pdo)],ncol=length(minyear:maxyear),nrow=length(unique(cty.clim$COUNTY_NAM)),byrow=T) 
  win.mei.mat<-matrix(winter.mei[1:length(winter.mei)],ncol=length(minyear:maxyear),nrow=length(unique(cty.clim$COUNTY_NAM)),byrow=T) 
  win.enso.mat<-matrix(winter.enso[1:length(winter.enso)],ncol=length(minyear:maxyear),nrow=length(unique(cty.clim$COUNTY_NAM)),byrow=T)
  #summer indices
  sum.nao.mat<-matrix(summer.nao[1:length(summer.nao)],ncol=length(minyear:maxyear),nrow=length(unique(cty.clim$COUNTY_NAM)),byrow=T)
  sum.pdo.mat<-matrix(summer.pdo[1:length(summer.pdo)],ncol=length(minyear:maxyear),nrow=length(unique(cty.clim$COUNTY_NAM)),byrow=T)
  sum.mei.mat<-matrix(summer.mei[1:length(summer.mei)],ncol=length(minyear:maxyear),nrow=length(unique(cty.clim$COUNTY_NAM)),byrow=T)
  sum.enso.mat<-matrix(summer.enso[1:length(summer.enso)],ncol=length(minyear:maxyear),nrow=length(unique(cty.clim$COUNTY_NAM)),byrow=T)
  climindex<-list(win.nao.mat=win.nao.mat,win.pdo.mat=win.pdo.mat,win.mei.mat=win.mei.mat,win.enso.mat=win.enso.mat,
                  sum.nao.mat=sum.nao.mat,sum.pdo.mat=sum.pdo.mat,sum.mei.mat=sum.mei.mat,sum.enso.mat=sum.enso.mat)#nao.mat=nao.mat,pdo.mat=pdo.mat,mei.mat=mei.mat,
  climindex<-lapply(climindex,function(x){row.names(x)<-unique(cty.clim$COUNTY_NAM);x})

  if(detr.flag=="detr"){
    climindex<-lapply(climindex,function(x){x<-reumannplatz::CleanData(x,normalize=F)$cleandat;x})
  }
  if(detr.flag=="norm"){
    climindex<-lapply(climindex,function(x){x<-reumannplatz::CleanData(x)$cleandat;x})
  }
  names(climindex)<-c("WinterNAO","WinterPDO","WinterMEI","WinterENSO","SummerNAO","SummerPDO","SummerMEI","SummerENSO")
}



