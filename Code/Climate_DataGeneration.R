###Climate Data Cleaning and Management Steps

#load location/gps files
noaa.coord<-read.csv("Data/noaa_latlongUpd.csv")
noaa.coord<-noaa.coord[,c("STATION","STATION_NA","COUNTY_N_1")]
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
wi.temp<-wi.temp[,c("STATION","STATION_NAME","DATE","TMAX","TMIN")]
temp7<-temp7[,c("STATION","STATION_NAME","DATE","TMAX","TMIN")]
wi.temp<-rbind(wi.temp,temp7)
wi.temp<-wi.temp[!duplicated(wi.temp),]

#Convert dates, and generate Julian date column
wi.temp$DATE<-gsub("(^\\d{4})(\\d{2})(\\d{2}$)", "\\1-\\2-\\3", wi.temp$DATE)
wi.temp$Jdate<-as.numeric(format(as.Date(as.character(wi.temp$DATE),"%Y-%m-%d"),"%j"))
wi.temp$Year<-as.numeric((format(as.Date(wi.temp$DATE,format="%Y-%m-%d"), "%Y")))
wi.temp$Month<-as.numeric((format(as.Date(wi.temp$DATE,format="%Y-%m-%d"), "%m")))

#Load precip data
pcp81_85<-read.csv("Data/Cty_NOAA_Precip1981_1985.csv",na.strings=c("","unknown","-9999"))#in standard units
pcp81_85$PRCP<-pcp81_85$PRCP*25.4 #change the units
pcp81_85$SNOW<-pcp81_85$SNOW*25.4 #change the units
pcp81_85$SNWD<-pcp81_85$SNWD*25.4 #change the units
pcp86_91<-read.csv("Data/Cty_NOAA_Precip1986_1991.csv",na.strings=c("","unknown","-9999"))
pcp92_97<-read.csv("Data/Cty_NOAA_Precip1992_1997.csv",na.strings=c("","unknown","-9999"))
pcp98_03<-read.csv("Data/Cty_NOAA_Precip1998_2003.csv",na.strings=c("","unknown","-9999"))
pcp04_09<-read.csv("Data/Cty_NOAA_Precip2004_2009.csv",na.strings=c("","unknown","-9999"))
pcp10_15a<-read.csv("Data/Cty_NOAA_Precip2010_2015set1.csv",na.strings=c("","unknown","-9999"))
pcp10_15b<-read.csv("Data/Cty_NOAA_Precip2010_2015set2.csv",na.strings=c("","unknown","-9999"))
pcp16_17<-read.csv("Data/Cty_NOAA_Precip2016_2017.csv",na.strings=c("","unknown","-9999"))

wi.pcp<-rbind(pcp81_85,pcp86_91,pcp92_97,pcp98_03,pcp04_09,pcp10_15a,pcp10_15b)  #combine all files but '16-'17 into one
wi.pcp<-wi.pcp[,c("STATION","STATION_NAME","DATE","PRCP","SNWD")]  #filter to relevant columns 
pcp16_17<-pcp16_17[,c("STATION","STATION_NAME","DATE","PRCP","SNWD")]  #filter to relevant columns 
wi.pcp<-rbind(wi.pcp,pcp16_17)

#Convert dates, and generate Julian date column
wi.pcp$DATE<-gsub("(^\\d{4})(\\d{2})(\\d{2}$)", "\\1-\\2-\\3", wi.pcp$DATE)
wi.pcp$Jdate<-as.numeric(format(as.Date(as.character(wi.pcp$DATE),"%Y-%m-%d"),"%j"))
wi.pcp$Year<-as.numeric((format(as.Date(wi.pcp$DATE,format="%Y-%m-%d"), "%Y")))
wi.pcp$Month<-as.numeric((format(as.Date(wi.pcp$DATE,format="%Y-%m-%d"), "%m")))

#merge precip and temp
wi.clim<-merge(wi.pcp,wi.temp,by=c("STATION","STATION_NAME","DATE","Jdate","Year","Month"),all=T)
wi.clim<-merge(wi.clim,noaa.coord[,c("STATION","COUNTY_N_1")],by=c("STATION"),all.x=T)
names(wi.clim)[names(wi.clim)=="COUNTY_N_1"]<-"COUNTY_NAM"

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
cty.tmin <- by(cty.clim$TMIN, INDICES=list(cty.clim$COUNTY_NAM, cty.clim$Year,cty.clim$Month), FUN=mean,na.rm=T)
cty.tmin<-cbind(expand.grid(attributes(cty.tmin)$dimnames), as.vector(cty.tmin))
names(cty.tmin)<-c("COUNTY_NAM","Year","Month","Tmin")
cty.tmax <- by(cty.clim$TMAX, INDICES=list(cty.clim$COUNTY_NAM, cty.clim$Year,cty.clim$Month), FUN=mean,na.rm=T)
cty.tmax<-cbind(expand.grid(attributes(cty.tmax)$dimnames), as.vector(cty.tmax))
names(cty.tmax)<-c("COUNTY_NAM","Year","Month","Tmax")
cty.pcp <- by(cty.clim$PRCP, INDICES=list(cty.clim$COUNTY_NAM, cty.clim$Year,cty.clim$Month), FUN=mean,na.rm=T)
cty.pcp<-cbind(expand.grid(attributes(cty.pcp)$dimnames), as.vector(cty.pcp))
names(cty.pcp)<-c("COUNTY_NAM","Year","Month","Pcp")
cty.snwd <- by(cty.clim$SNWD, INDICES=list(cty.clim$COUNTY_NAM, cty.clim$Year,cty.clim$Month), FUN=mean,na.rm=T)
cty.snwd<-cbind(expand.grid(attributes(cty.snwd)$dimnames), as.vector(cty.snwd))
names(cty.snwd)<-c("COUNTY_NAM","Year","Month","Snwd")
cty.month<-merge(cty.tmin,cty.tmax,by=c("COUNTY_NAM","Year","Month"),all=T)
cty.month<-merge(cty.month,cty.pcp,by=c("COUNTY_NAM","Year","Month"),all=T)
cty.month<-merge(cty.month,cty.snwd,by=c("COUNTY_NAM","Year","Month"),all=T)
cty.month$Month<-as.numeric(cty.month$Month)

cty.month<-cty.month[!(cty.month$Year==2017 & cty.month$Month>5),]

#make sure Menomonie county is dropped
cty.month<-cty.month[!(cty.month$COUNTY_NAM=="Menominee"),]

#Merge daily climate data by county with usda district data
usda.clim<-merge(cty.clim,usda,by.x="COUNTY_NAM",by.y="County",all=T)
#Merge monthly climate data by county with usda district data
usda.month<-merge(cty.month,usda,by.x="COUNTY_NAM",by.y="County",all=T)

#remove large unneeded files
rm(temp1,temp2,temp3,temp4,temp5,temp6,temp7,wi.temp)
rm(pcp04_09,pcp10_15a,pcp10_15b,pcp16_17,pcp81_85,pcp86_91,pcp92_97,pcp98_03,wi.pcp)
rm(wi.clim,cty.min,cty.max,cty.snwd)

# County Climate Data -----------------------------------------------------
#Winter weather conditions
winter.tmp<-matrix(NA, ncol=length(unique(cty.month$COUNTY_NAM)),nrow=length(minyear:maxyear))
winter.clim<-list()
years<-minyear:maxyear
ctynames<-sort(unique(cty.month$COUNTY_NAM))
respvars<-c("Tmin","Tmax","Pcp","Snwd")

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
winter.clim<-lapply(winter.clim,t)
winter.clim<-lapply(winter.clim,function(x){rownames(x)<-colnames(wsi.mat);x})
names(winter.clim)<-c("Tmin","Tmax","Prcp","Snwd","WSI")
saveRDS(winter.clim,"Results/winter.clim.rds")

#Winter conditions for USDA districts
winter.tmp<-matrix(NA, ncol=length(unique(usda.month$Zone)),nrow=length(minyear:maxyear))
winter.clim.usda<-list()
years<-minyear:maxyear
usda.names<-sort(unique(usda.month$Zone))
respvars<-c("Tmin","Tmax","Pcp","Snwd")

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
winter.clim.usda<-lapply(winter.clim.usda,function(x){x<-t(x);x})
names(winter.clim.usda)<-c("Tmin","Tmax","Prcp","Snwd","WSI")
saveRDS(winter.clim.usda,"Results/winter.clim.usda.rds")

rm(cty.month,usda.month)

#load large climate indice data
climate_indices <- read.csv("Data/climate_indices.csv")
mei<-read.csv("Data/mei.csv")
pdo<-read.csv("Data/pdo.csv")

#compute yearly averages
mei$Mean<-rowMeans(mei[,2:dim(mei)[2]],na.rm=T)
mean.clim<-aggregate(cbind(NAO,EA,WP,EP_NP,PNA,EA.WR)~Year,FUN="mean",data=climate_indices)
mean.pdo<-aggregate(Value~Year,pdo,FUN="mean")

#calculate winter means
#Mean Winter NAO (Dec. through Mar.)
winter.nao<-rep(NA, length(minyear:maxyear))
years<-unique(climate_indices$Year[climate_indices$Year%in%c(minyear:maxyear)])
for(year in years){
  tmp<-climate_indices[climate_indices$Year==year & climate_indices$Month<=3,]
  tmp<-rbind(tmp, climate_indices[climate_indices$Year==(year-1) & climate_indices$Month>=12,])
  winter.nao[years %in% year]<-mean(tmp$NAO,na.rm=T)
}
winter.nao

#Mean Winter PDO (Dec. through Mar.)
winter.pdo<-rep(NA, length(minyear:maxyear))
years<-unique(pdo$Year[pdo$Year%in%c(minyear:maxyear)])
for(year in years){
  tmp<-pdo[pdo$Year==year & pdo$Month<=3,]
  tmp<-rbind(tmp, pdo[pdo$Year==(year-1) & pdo$Month>=12,])
  winter.pdo[years %in% year]<-mean(tmp$Value,na.rm=T)
}
winter.pdo

#Mean Winter MEI (Dec. through Mar.)
mei.long<-melt(mei[,1:13],id.vars=c("YEAR"),variable.name = "MONTHS",value.name = "MEI")
mei.long$Index<-rep(1:12,each=67)
winter.mei<-rep(NA, length(minyear:maxyear))
years<-unique(mei.long$YEAR[mei.long$YEAR%in%c(minyear:maxyear)])
for(year in years){
  tmp<-mei.long[mei.long$YEAR==year & mei.long$Index<=3,]
  tmp<-rbind(tmp, mei.long[mei.long$YEAR==(year-1) & mei.long$Index>=12,])
  winter.mei[years %in% year]<-mean(tmp$MEI,na.rm=T)
}

#Mean summer NAO (Dec. through Mar.)
summer.nao<-rep(NA, length(minyear:maxyear))
years<-unique(climate_indices$Year[climate_indices$Year%in%c(minyear:maxyear)])
for(year in years){
  tmp<-climate_indices[climate_indices$Year==year & climate_indices$Month%in%c(4:11),]
  summer.nao[years %in% year]<-mean(tmp$NAO,na.rm=T)
}

#Mean summer PDO (April through November)
summer.pdo<-rep(NA, length(minyear:maxyear))
years<-unique(pdo$Year[pdo$Year%in%c(minyear:maxyear)])
for(year in years){
  tmp<-pdo[pdo$Year==year & pdo$Month%in%c(4:11),]
  summer.pdo[years %in% year]<-mean(tmp$Value,na.rm=T)
}

#Mean summer MEI (April through November)
mei.long<-melt(mei[,1:13],id.vars=c("YEAR"),variable.name = "MONTHS",value.name = "MEI")
mei.long$Index<-rep(1:12,each=67)
summer.mei<-rep(NA, length(minyear:maxyear))
years<-unique(mei.long$YEAR[mei.long$YEAR%in%c(minyear:maxyear)])
for(year in years){
  tmp<-mei.long[mei.long$YEAR==year & mei.long$Index%in%c(4:11),]
  summer.mei[years %in% year]<-mean(tmp$MEI,na.rm=T)
}

#put into repeating matrices
win.nao.mat<-matrix(winter.nao[1:length(winter.nao)],ncol=length(minyear:maxyear),nrow=length(unique(usda.clim$Zone)),byrow=T)
win.pdo.mat<-matrix(winter.pdo[1:length(winter.pdo)],ncol=length(minyear:maxyear),nrow=length(unique(usda.clim$Zone)),byrow=T)
win.mei.mat<-matrix(winter.mei[1:length(winter.mei)],ncol=length(minyear:maxyear),nrow=length(unique(usda.clim$Zone)),byrow=T)

#summer indices
sum.nao.mat<-matrix(summer.nao[1:length(summer.nao)],ncol=length(minyear:maxyear),nrow=length(unique(usda.clim$Zone)),byrow=T)
sum.pdo.mat<-matrix(summer.pdo[1:length(summer.pdo)],ncol=length(minyear:maxyear),nrow=length(unique(usda.clim$Zone)),byrow=T)
sum.mei.mat<-matrix(summer.mei[1:length(summer.mei)],ncol=length(minyear:maxyear),nrow=length(unique(usda.clim$Zone)),byrow=T)
climindex.usda<-list(win.nao.mat=win.nao.mat,win.pdo.mat=win.pdo.mat,win.mei.mat=win.mei.mat,
                     sum.nao.mat=sum.nao.mat,sum.pdo.mat=sum.pdo.mat,sum.mei.mat=sum.mei.mat)#nao.mat=nao.mat,pdo.mat=pdo.mat,mei.mat=mei.mat,
names(climindex.usda)<-c("WinterNAO","WinterPDO","WinterMEI","SummerNAO","SummerPDO","SummerMEI")
saveRDS(climindex.usda,"Results/climindex.usda.rds")


#winter indices
win.nao.mat<-matrix(winter.nao[1:length(winter.nao)],ncol=length(minyear:maxyear),nrow=length(unique(cty.clim$COUNTY_NAM)),byrow=T)
win.pdo.mat<-matrix(winter.pdo[1:length(winter.pdo)],ncol=length(minyear:maxyear),nrow=length(unique(cty.clim$COUNTY_NAM)),byrow=T) 
win.mei.mat<-matrix(winter.mei[1:length(winter.mei)],ncol=length(minyear:maxyear),nrow=length(unique(cty.clim$COUNTY_NAM)),byrow=T) 
#summer indices
sum.nao.mat<-matrix(summer.nao[1:length(summer.nao)],ncol=length(minyear:maxyear),nrow=length(unique(cty.clim$COUNTY_NAM)),byrow=T)
sum.pdo.mat<-matrix(summer.pdo[1:length(summer.pdo)],ncol=length(minyear:maxyear),nrow=length(unique(cty.clim$COUNTY_NAM)),byrow=T)
sum.mei.mat<-matrix(summer.mei[1:length(summer.mei)],ncol=length(minyear:maxyear),nrow=length(unique(cty.clim$COUNTY_NAM)),byrow=T)
climindex<-list(win.nao.mat=win.nao.mat,win.pdo.mat=win.pdo.mat,win.mei.mat=win.mei.mat,
                sum.nao.mat=sum.nao.mat,sum.pdo.mat=sum.pdo.mat,sum.mei.mat=sum.mei.mat)
climindex<-lapply(climindex,function(x){row.names(x)<-unique(cty.clim$COUNTY_NAM);x})
names(climindex)<-c("WinterNAO","WinterPDO","WinterMEI","SummerNAO","SummerPDO","SummerMEI")
saveRDS(climindex,"Results/climindex.rds")

rm(cty.clim,usda.clim)
