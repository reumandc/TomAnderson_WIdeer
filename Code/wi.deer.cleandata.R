#Hunting, Deer and DVC Data Cleaning and Management Steps

#Load Data
wi.deer <- read.csv("Data/State_Kill.csv")
cty.deer <- read.csv("Data/County_KillTotal.csv")
cty.deer.buck<-read.csv("Data/County_Kill_BuckOnly.csv")
cty.deer.doe<-read.csv("Data/County_Kill_DoeOnly.csv")
cty.deer.gunbuck<-read.csv("Data/County_KillGunBuck.csv")
cty.deer.gundoe<-read.csv("Data/County_KillGunDoe.csv")
cty.deer.bowbuck<-read.csv("Data/County_KillBowBuck.csv")
cty.deer.bowdoe<-read.csv("Data/County_KillBowDoe.csv")
#gps <- read.csv("WI_Cty_Centroid.csv") #load if you want UTM coordinates rather than DD
gps <- read.csv("Data/WI_Cty_CentroidDD.csv")
crash <- read.csv("Data/DOT_Accidents1987_2016.csv")
adj.crash <- read.csv("Data/Adj_DOT_Accidents1987_2016.csv")
openingday<-read.csv("Data/OpeningDay.csv")
crop<-read.csv("Data/CornHarvest.csv")
usda<-read.csv("Data/WI_USDA_Sections.csv")
usda.centroid<-read.csv("Data/USDA_centroids.csv")
crop.cty<-merge(usda,crop,by="Zone")
traffic<-read.csv("Data/Cty_Traffic.csv")
abun.est<-read.csv("Data/DMU_DeerAbundance1981_2013.csv",na.strings = c(""))
dens.est<-read.csv("Data/DMU_DeerDensity1981_2013.csv",na.strings = c(""))
hunter.abun<-read.csv("Data/DMU_HunterAbundance1975_2013.csv",na.strings = c(""))
hunter.dens<-read.csv("Data/DMU_HunterDensity1975_2013.csv",na.strings = c(""))
hunter.dens1416<-read.csv("Data/hunterdens1416.csv",na.strings = c(""))
hunter.abun1416<-read.csv("Data/hunterabun1416.csv",na.strings = c(""))
cty.abun.1416 <- read.csv("Data/Cty_Abun_2014_2016.csv",na.strings = c(""))
permits <- read.csv("Data/permit_issued1970_2013.csv",na.strings = c("","."))

#CLEAN DATA
gps<-gps[order(gps$COUNTY_NAM),] #order counties
wi.deer$TotalKills<-rowSums(wi.deer[,c("TotalGun","TotalCross","TotalBow")],na.rm=T)#calculate state total kills across methods
wi.deer$TotalHunters<-rowSums(wi.deer[,c("GunHunters","BowHunters","CrossHunters")],na.rm=T)#calculate state total hunters across methods

#fix county names
levels(gps$COUNTY_NAM)[levels(gps$COUNTY_NAM)=="Saint Croix"] <- "St. Croix"
levels(gps$COUNTY_NAM)[levels(gps$COUNTY_NAM)=="Fond du Lac"] <- "Fond Du Lac"
levels(crash$County)[levels(crash$County)=="Sheboygon"] <- "Sheboygan"
levels(adj.crash$County)[levels(adj.crash$County)=="Sheboygon"] <- "Sheboygan"
levels(usda$County)[levels(usda$County)=="Saint Croix"] <- "St. Croix"
levels(usda$County)[levels(usda$County)=="Fond du Lac"] <- "Fond Du Lac"

#convert opening day of gun hunting to Julian
openingday$Jdate<-as.numeric(format(as.Date(as.character(openingday$Date),"%m/%d/%Y"),"%j"))

#eliminate dashes in DMU values
abun.est$Unit<-noquote(gsub("-", "", abun.est$Unit))
dens.est$Unit<-noquote(gsub("-", "", dens.est$Unit))

hunter.abun<-hunter.abun[,!names(hunter.abun)%in%c("X60_64","X65_69","X70_74","X1975","X1976","X1977","X1978","X1979","X1980","X")] #eliminate years prior to 1981, and notes column
hunter.dens<-hunter.dens[,!names(hunter.dens)%in%c("X60_64","X65_69","X70_74","X1975","X1976","X1977","X1978","X1979","X1980")] #eliminate years prior to 1981
hunter.dens1416$County<-as.factor(toTitleCase(tolower(as.character(hunter.dens1416$County)))) #make names match in capitalization
hunter.abun1416$County<-as.factor(toTitleCase(tolower(as.character(hunter.abun1416$County)))) #make names match in capitalization

cty.deer.gundoe[is.na(cty.deer.gundoe)]<-0 #replace NAs with 0 (2 in 1960)
cty.deer.gunbuck[is.na(cty.deer.gunbuck)]<-0 #replace NAs with 0 (1 in 1960)
cty.deer.bowdoe[is.na(cty.deer.bowdoe)]<-0 #replace NAs with 0 (1 in 1960)
cty.deer<-cty.deer[!cty.deer$County=="Total",] #drop totals row
cty.deer<-cty.deer[,!names(cty.deer)%in%c("Longitude","Latitude")] #drop gps columns

#Fill deer range estimates that were NA or wrong, based on deer range values from other files
dens.est$DeerRange[dens.est$Unit==5]<-226
dens.est$DeerRange[dens.est$Unit==35]<-409
dens.est$DeerRange[dens.est$Unit==38]<-388
dens.est$DeerRange[dens.est$Unit==39]<-411
dens.est$DeerRange[dens.est$Unit==70]<-273

#generate lists of spatial subsets of WI
northern.cty<-sort(c("Douglas","Bayfield","Ashland","Iron","Vilas","Forest","Florence","Marinette",
                     "Langlade","Lincoln","Price","Sawyer","Washburn","Burnett","Rusk","Taylor","Oneida",
                     "Barron","Oconto","Menominee","Shawano","Marathon","Clark","Chippewa","Polk","St. Croix",
                     "Dunn","Eau Claire"))
southern.cty<-unique(gps$COUNTY_NAM[!(gps$COUNTY_NAM%in%northern.cty)])
#cwd counties and zones
cwd<-c("Grant","Kenosha","Racine","Waukesha","Jefferson","Dane","Columbia",'Sauk',
       "Richland","Vernon","Crawford","Juneau","Iowa","Lafayette","Green","Rock","Walworth")#cwd counties
cwd.usda<-c("SE","SC","SW")#cwd zones

#convert DMU data to long format
abun.estL<-melt(abun.est,id.vars = c("Unit"),variable.name="Year",value.name="Abun")
dens.estL<-melt(dens.est,id.vars = c("Unit","DeerRange"),variable.name="Year",value.name="Dens")
hunter.abunL<-melt(hunter.abun,id.vars = c("Region","UnitGroup","Unit","DeerRangeOld","DeerRange2001"),variable.name="Year",value.name="HunterAbun")
hunter.densL<-melt(hunter.dens,id.vars = c("Region","UnitGroup","Unit","DeerRangeOld","DeerRange2001"),variable.name="Year",value.name="HunterDens")
abun.estL$Year<-as.numeric(sub('X', '', abun.estL$Year))
dens.estL$Year<-as.numeric(sub('X', '', dens.estL$Year))
hunter.abunL$Year<-as.numeric(sub('X', '', hunter.abunL$Year))
hunter.densL$Year<-as.numeric(sub('X', '', hunter.densL$Year))

#combine data with gps coordinates
cty.deer<-merge(gps[,names(gps)%in%c("COUNTY_NAM","Longitude","Latitude")],cty.deer,by.x=c("COUNTY_NAM"),by.y=c("County"),select.x=c("Latitude","Longitude"))
cty.deer.buck<-merge(gps[,names(gps)%in%c("COUNTY_NAM","Longitude","Latitude")],cty.deer.buck,by.x=c("COUNTY_NAM"),by.y=c("County"),select.x=c("Latitude","Longitude"))
cty.deer.doe<-merge(gps[,names(gps)%in%c("COUNTY_NAM","Longitude","Latitude")],cty.deer.doe,by.x=c("COUNTY_NAM"),by.y=c("County"),select.x=c("Latitude","Longitude"))
cty.deer.gunbuck<-merge(gps[,names(gps)%in%c("COUNTY_NAM","Longitude","Latitude")],cty.deer.gunbuck,by.x=c("COUNTY_NAM"),by.y=c("County"),select.x=c("Latitude","Longitude"))
cty.deer.gundoe<-merge(gps[,names(gps)%in%c("COUNTY_NAM","Longitude","Latitude")],cty.deer.gundoe,by.x=c("COUNTY_NAM"),by.y=c("County"),select.x=c("Latitude","Longitude"))
cty.deer.bowbuck<-merge(gps[,names(gps)%in%c("COUNTY_NAM","Longitude","Latitude")],cty.deer.bowbuck,by.x=c("COUNTY_NAM"),by.y=c("County"),select.x=c("Latitude","Longitude"))
cty.deer.bowdoe<-merge(gps[,names(gps)%in%c("COUNTY_NAM","Longitude","Latitude")],cty.deer.bowdoe,by.x=c("COUNTY_NAM"),by.y=c("County"),select.x=c("Latitude","Longitude"))
crash.gps<-merge(gps[,names(gps)%in%c("COUNTY_NAM","Longitude","Latitude")],crash,by.x=c("COUNTY_NAM"),by.y=c("County"))
adj.crash.gps<-merge(gps[,names(gps)%in%c("COUNTY_NAM","Longitude","Latitude")],adj.crash,by.x=c("COUNTY_NAM"),by.y=c("County"))
traffic.gps<-merge(gps[,names(gps)%in%c("COUNTY_NAM","Longitude","Latitude")],traffic,by.x=c("COUNTY_NAM"),by.y=c("COUNTY"))

#convert to long format
cty.deer1<-melt(cty.deer,id.vars = c("Longitude","Latitude","COUNTY_NAM"),variable.name="Year",value.name="Kills")
cty.deer.buck1<-melt(cty.deer.buck,id.vars = c("Longitude","Latitude","COUNTY_NAM"),variable.name="Year",value.name="Kills")
cty.deer.doe1<-melt(cty.deer.doe,id.vars = c("Longitude","Latitude","COUNTY_NAM"),variable.name="Year",value.name="Kills")
cty.deer.gunbuck1<-melt(cty.deer.gunbuck,id.vars = c("Longitude","Latitude","COUNTY_NAM"),variable.name="Year",value.name="Kills")
cty.deer.gundoe1<-melt(cty.deer.gundoe,id.vars = c("Longitude","Latitude","COUNTY_NAM"),variable.name="Year",value.name="Kills")
cty.deer.bowbuck1<-melt(cty.deer.bowbuck,id.vars = c("Longitude","Latitude","COUNTY_NAM"),variable.name="Year",value.name="Kills")
cty.deer.bowdoe1<-melt(cty.deer.bowdoe,id.vars = c("Longitude","Latitude","COUNTY_NAM"),variable.name="Year",value.name="Kills")
crash.gps1<-melt(crash.gps,id.vars = c("Longitude","Latitude","COUNTY_NAM"),variable.name="Year",value.name="Crashes")
adj.crash.gps1<-melt(adj.crash.gps,id.vars = c("Longitude","Latitude","COUNTY_NAM"),variable.name="Year",value.name="Crashes")
traffic.gps1<-melt(traffic.gps,id.vars = c("Longitude","Latitude","COUNTY_NAM"),variable.name="Year",value.name="Traffic")

#make year column
cty.deer1$Year<-as.numeric(sub('X', '', cty.deer1$Year))
cty.deer.buck1$Year<-as.numeric(sub('X', '', cty.deer.buck1$Year))
cty.deer.doe1$Year<-as.numeric(sub('X', '', cty.deer.doe1$Year))
cty.deer.gunbuck1$Year<-as.numeric(sub('X', '', cty.deer.gunbuck1$Year))
cty.deer.gundoe1$Year<-as.numeric(sub('X', '', cty.deer.gundoe1$Year))
cty.deer.bowbuck1$Year<-as.numeric(sub('X', '', cty.deer.bowbuck1$Year))
cty.deer.bowdoe1$Year<-as.numeric(sub('X', '', cty.deer.bowdoe1$Year))
crash.gps1$Year<-as.numeric(sub('X', '', crash.gps1$Year))
traffic.gps1$Year<-as.numeric(sub('X', '', traffic.gps1$Year))
adj.crash.gps1$Year<-as.numeric(sub('X', '', adj.crash.gps1$Year))
names(adj.crash.gps1)[names(adj.crash.gps1)=="Crashes"]<-"AdjDVC"

#merge files together
dat<-merge(cty.deer1,cty.deer.buck1,by=c("COUNTY_NAM","Year","Latitude","Longitude"))
names(dat)[names(dat)=="Kills.x"] <- "TotalKills"
names(dat)[names(dat)=="Kills.y"] <- "BuckKills"
dat<-merge(dat,cty.deer.doe1,by=c("COUNTY_NAM","Year","Latitude","Longitude"))
names(dat)[names(dat)=="Kills"] <- "DoeKills"
dat<-merge(dat,cty.deer.gunbuck1,by=c("COUNTY_NAM","Year","Latitude","Longitude"))
names(dat)[names(dat)=="Kills"] <- "GunBuckKills"
dat<-merge(dat,cty.deer.gundoe1,by=c("COUNTY_NAM","Year","Latitude","Longitude"))
names(dat)[names(dat)=="Kills"] <- "GunDoeKills"
dat$TotalGunKills<-dat$GunBuckKills+dat$GunDoeKills

#load DMU-Cty-Deer Range merged files
cty.flag<-"dez"
source('Code/DMU_Cty_DataStep.R')

cty.dmu.dr<-aggregate(AREA_SQM~UNIT_ID+COUNTY_NAM+Year,data=cty.dmu.dr,FUN=sum)
levels(cty.dmu.dr$COUNTY_NAM)[levels(cty.dmu.dr$COUNTY_NAM)=="Fond du Lac"] <- "Fond Du Lac"
levels(cty.dmu.dr$COUNTY_NAM)[levels(cty.dmu.dr$COUNTY_NAM)=="Saint Croix"] <- "St. Croix"

#calculate the amount of deer range by DMU, County and for DMUs used in the population estimates
deer.range<-(spread(aggregate(AREA_SQM~UNIT_ID+Year,data=cty.dmu.dr,FUN=sum),key = Year,value = AREA_SQM))#summarize habitat by DMU
deer.range1<-(spread(aggregate(AREA_SQM~COUNTY_NAM+Year,data=cty.dmu.dr,FUN=sum),key = Year,value = AREA_SQM))#summarize habitat by County
deer.range.wide<-deer.range[match(dens.est$Unit,deer.range$UNIT_ID),] #range to match wide format density and abundance
deer.range.wideH<-deer.range[match(hunter.dens$Unit,deer.range$UNIT_ID),] #range to match wide format density and abundance

#Run the next line to compute your own density estimates
if(dens.flag=="tom")
{
  dens.est2<-cbind(Unit=abun.est$Unit,abun.est[,2:dim(abun.est)[2]]/deer.range.wide[,2:dim(deer.range)[2]])
  dens.est2L<-melt(dens.est2,id.vars = c("Unit"),variable.name="Year",value.name="Dens")
  hunter.est<-cbind(Unit=hunter.abun$Unit,hunter.abun[,6:dim(hunter.abun)[2]]/deer.range.wideH[,2:dim(deer.range)[2]])
  hunter.est2L<-melt(hunter.est,id.vars = c("Unit"),variable.name="Year",value.name="HunterDens")
}

#Use the WIDNR density estimates
if(dens.flag=="dnr")
{
  dens.est2L<-melt(dens.est[,-2],id.vars = c("Unit"),variable.name="Year",value.name="Dens")
  hunter.est2L<-melt(hunter.dens[,-c(1,2,4,5)],id.vars = c("Unit"),variable.name="Year",value.name="HunterDens")
}
head(dens.est2L)
dens.est2L$Year<-as.numeric(sub('X', '', dens.est2L$Year))
head(hunter.est2L)
hunter.est2L$Year<-as.numeric(sub('X', '', hunter.est2L$Year))

#combine DMUs and Counties
dmu.cty.estL<-merge(dens.est2L,cty.dmu.dr,by.x=c("Unit","Year"),by.y=c("UNIT_ID","Year"),all=T)
dmu.cty.hunter.estL<-merge(hunter.est2L,cty.dmu.dr,by.x=c("Unit","Year"),by.y=c("UNIT_ID","Year"),all=T)

#generate abundance estimates by multiplying deer area by deer density
dmu.cty.estL$Abun<-dmu.cty.estL$Dens*dmu.cty.estL$AREA_SQM
dmu.cty.hunter.estL$Hunters<-dmu.cty.hunter.estL$HunterDens*dmu.cty.hunter.estL$AREA_SQM

#aggregate abundance by county, and merge with gps locations
cty.abun.estL<-aggregate(Abun~COUNTY_NAM+Year,FUN=sum,data=dmu.cty.estL)
cty.abun.estL<-merge(gps[,names(gps)%in%c("COUNTY_NAM","Longitude","Latitude")],cty.abun.estL,by.x=c("COUNTY_NAM"),by.y=c("COUNTY_NAM"),select.x=c("Latitude","Longitude"))

#aggregate hunters by county, and merge with gps locations
cty.hunter.estL<-aggregate(Hunters~COUNTY_NAM+Year,FUN=sum,data=dmu.cty.hunter.estL)
cty.hunter.estL<-merge(gps[,names(gps)%in%c("COUNTY_NAM","Longitude","Latitude")],cty.hunter.estL,by.x=c("COUNTY_NAM"),by.y=c("COUNTY_NAM"),select.x=c("Latitude","Longitude"),all=T)

#aggregate abundance for 2014-2016 data, combine with gps coordinates
cty.abun1.1416<-ddply(cty.abun.1416,c("County", "Year"), summarise, Abun=sum(PostHuntPop),.drop=FALSE)
cty.abun1.1416<-merge(gps[,names(gps)%in%c("COUNTY_NAM","Longitude","Latitude")],cty.abun1.1416,by.x=c("COUNTY_NAM"),by.y=c("County"),select.x=c("Latitude","Longitude"),all=T)

hunters1416L<-ddply(hunter.abun1416,c("County", "Year"), summarise, Hunters=sum(HunterAbun),.drop=FALSE)
hunters1416L<-merge(gps[,names(gps)%in%c("COUNTY_NAM","Longitude","Latitude")],hunters1416L,by.x=c("COUNTY_NAM"),by.y=c("County"),select.x=c("Latitude","Longitude"),all=T)

#combine 1981-2013 with 2014-2016 abundance
cty.abun.estL<-rbind(cty.abun.estL,cty.abun1.1416)
cty.abun.estL$Abun[is.na(cty.abun.estL$Abun)]<-0 #make two menominee NA values to 0

#combine 1981-2013 with 2014-2016 hunters
cty.hunter.estL<-rbind(cty.hunter.estL,hunters1416L)
#cty.hunter.estL$Hunters[is.na(cty.hunter.estL$Abun)]<-0 #make two menominee NA values to 0

#generate data in wide format
cty.abun.est<-spread(cty.abun.estL,key=Year,value=Abun)
hunter.abun.est<-spread(cty.hunter.estL,key=Year,value=Hunters)

#combine kill and abundance data into one dataframe
dat<-merge(dat,cty.abun.estL,by=c("COUNTY_NAM","Year","Latitude","Longitude"),all=T)

#combine with hunter data
dat<-merge(dat,cty.hunter.estL,by=c("COUNTY_NAM","Year","Latitude","Longitude"),all=T)

#combine with traffic and DVC data
dat<-merge(dat,crash.gps1,by=c("COUNTY_NAM","Year","Latitude","Longitude"),all=T)
dat<-merge(dat,traffic.gps1,by=c("COUNTY_NAM","Year","Latitude","Longitude"),all=T)
dat<-merge(dat,adj.crash.gps1,by=c("COUNTY_NAM","Year","Latitude","Longitude"),all=T)

#Combine with usda districts
dat<-merge(dat,usda,by.x=c("COUNTY_NAM"),by.y=c("County"),all=T)

#Drop Menominee county because it has very little data for anything (its an Indian reservation)
dat<-dat[!(dat$COUNTY_NAM=="Menominee"),]

#aggregate data by usda districts
usda.dat<-ddply(dat,.variables = c("Zone","Year"),summarise,TotalKills=sum(TotalKills),BuckKills=sum(BuckKills),DoeKills=sum(DoeKills),
                GunBuckKills=sum(GunBuckKills),GunDoeKills=sum(GunDoeKills),TotalGunKills=sum(TotalGunKills),Abun=sum(Abun),Hunters=sum(Hunters,na.rm=T),
                Crashes=sum(Crashes),Traffic=sum(Traffic),.drop=F)
#generate traffic-adjusted crashes
usda.dat$AdjDVC<-usda.dat$Crashes/usda.dat$Traffic
#combine with gps locations (centroid of USDA district)
usda.dat<-merge(usda.dat,usda.centroid,by="Zone",all=T)

usda.dat<-usda.dat[usda.dat$Year%in%c(minyear:maxyear),]

#LOAD NECESSARY FUNCTIONS
#source("c:/Users/reulab/Desktop/TomA/Rwork/wavelet_R_func/Fn_demean_standard_rows.R")
#source("c:/Users/reulab/Desktop/TomA/Rwork/Functions/Fn_flip_widetolong.R")

dat<-subset(dat,Year%in%c(minyear:maxyear))
cty.mat<-matrix(NA,length(unique(dat$COUNTY_NAM)),length(minyear:maxyear))
cty.list<-list()
for(j in 5:15){
  for(i in 1:length(unique(dat$COUNTY_NAM))){
    dat.sub<-dat[as.numeric(as.factor(as.character(dat$COUNTY_NAM)))==i,]
    dat.sub<-dat.sub[order(dat.sub$Year),]
    cty.mat[i,]<-dat.sub[,j]
  }
  cty.list[[j-4]]<-cty.mat
}
names(cty.list)<-colnames(dat)[5:15]
# dat.list<-list(total.kill=dat$TotalKills[dat$Year%in%c(minyear:maxyear)],
#                  buck.kill=dat$BuckKills[dat$Year%in%c(minyear:maxyear)],
#                  doe.kill=dat$DoeKills[dat$Year%in%c(minyear:maxyear)],
#                  gun.buck=dat$GunBuckKills[dat$Year%in%c(minyear:maxyear)],
#                  gun.doe=dat$GunDoeKills[dat$Year%in%c(minyear:maxyear)],
#                  total.gun=dat$TotalGunKills[dat$Year%in%c(minyear:maxyear)],
#                  AbunEst=dat$Abun[dat$Year%in%c(minyear:maxyear)],
#                  Hunters=dat$Hunters[dat$Year%in%c(minyear:maxyear)],
#                  DVC=dat$Crashes[dat$Year%in%c(minyear:maxyear)],
#                  Traffic=dat$Traffic[dat$Year%in%c(minyear:maxyear)])
# cty.list<-lapply(dat.list,function(x){x<-matrix(x,nrow=length(unique(dat$COUNTY_NAM)),ncol=length(minyear:maxyear),byrow=T);x})
openingday<-openingday[order(openingday$Year),]
openday.dat<-openingday$Jdate[openingday$Year%in%c(minyear:maxyear)]
openday.mat<-matrix(openday.dat,nrow=length(unique(dat$COUNTY_NAM)),ncol=length(openday.dat),byrow=T)
cty.list[[length(cty.list)+1]]<-openday.mat
names(cty.list)[12]<-"OpenDay"
cty.list<-lapply(cty.list,function(x){row.names(x)<-unique(dat$COUNTY_NAM);x})

#move data summarized by usda group into a list of matrices
# usda.list<-list(total.kill=usda.dat$TotalKills[usda.dat$Year%in%c(minyear:maxyear)],
#                buck.kill=usda.dat$BuckKills[usda.dat$Year%in%c(minyear:maxyear)],
#                doe.kill=usda.dat$DoeKills[usda.dat$Year%in%c(minyear:maxyear)],
#                gun.buck=usda.dat$GunBuckKills[usda.dat$Year%in%c(minyear:maxyear)],
#                gun.doe=usda.dat$GunDoeKills[usda.dat$Year%in%c(minyear:maxyear)],
#                total.gun=dat$TotalGunKills[dat$Year%in%c(minyear:maxyear)],
#                AbunEst=usda.dat$Abun[usda.dat$Year%in%c(minyear:maxyear)],
#                Hunters=usda.dat$Hunters[usda.dat$Year%in%c(minyear:maxyear)],
#                DVC=usda.dat$Crashes[usda.dat$Year%in%c(minyear:maxyear)],
#                Traffic=usda.dat$Traffic[usda.dat$Year%in%c(minyear:maxyear)])
# 
# usda.list<-lapply(usda.list,function(x){x<-matrix(x,nrow=length(unique(usda.dat$Zone)),ncol=length(minyear:maxyear),byrow=T);x})

usda.dat<-subset(usda.dat,Year%in%c(minyear:maxyear))
usda.mat<-matrix(NA,length(unique(usda.dat$Zone)),length(minyear:maxyear))
usda.list<-list()
for(j in 3:13){
  for(i in 1:length(unique(usda.dat$Zone))){
    dat.sub<-usda.dat[as.numeric(usda.dat$Zone)==i,]
    dat.sub<-dat.sub[order(dat.sub$Year),]
    usda.mat[i,]<-dat.sub[,j]
  }
  usda.list[[j-2]]<-usda.mat
}
names(usda.list)<-colnames(usda.dat)[3:13]
openingday<-openingday[order(openingday$Year),]
openday.dat<-openingday$Jdate[openingday$Year%in%c(minyear:maxyear)]
openday.mat<-matrix(openday.dat,nrow=length(unique(usda.dat$Zone)),ncol=length(openday.dat),byrow=T)
usda.list[[length(usda.list)+1]]<-openday.mat
crop<-crop[order(crop$Year,crop$Zone),]
crop.mat<-matrix(crop$PercentGrain/100,nrow=length(unique(usda.dat$Zone)),ncol=length(1997:maxyear),byrow=F)
usda.list[[length(usda.list)+1]]<-crop.mat
names(usda.list)[12:13]<-c("OpenDay","Crop")
usda.list<-lapply(usda.list,function(x){row.names(x)<-sort(unique(usda$Zone));x})

if(detr.flag=="norm"){
  cty.list<-lapply(cty.list,function(x){x<-CleanData(x[,colSums(is.na(x)) != nrow(x)],normalize=T)$cleandat;x})
  usda.list<-lapply(usda.list,function(x){x<-CleanData(x[,colSums(is.na(x)) != nrow(x)],normalize=T)$cleandat;x})
}

if(detr.flag=="detr"){
  cty.list<-lapply(cty.list,function(x){x<-CleanData(x[,colSums(is.na(x)) != nrow(x)],normalize=F,detrend=T)$cleandat;x})
  usda.list<-lapply(usda.list,function(x){x<-CleanData(x[,colSums(is.na(x)) != nrow(x)],normalize=F,detrend = T)$cleandat;x})
}
