#Deer and DVC Data Cleaning and Management Steps

#Load Data
gps <- read.csv("Data/WI_Cty_CentroidDD.csv")
usda.centroid<-read.csv("Data/USDA_centroids.csv")
state.totals <- read.csv("Data/State_Kill.csv")
crash <- read.csv("Data/DOT_Accidents1987_2016.csv")
adj.crash <- read.csv("Data/Adj_DOT_Accidents1987_2016.csv")
usda<-read.csv("Data/WI_USDA_Sections.csv")
traffic<-read.csv("Data/Cty_Traffic.csv")
abun.est<-read.csv("Data/DMU_DeerAbundance1981_2013.csv",na.strings = c(""))
dens.est<-read.csv("Data/DMU_DeerDensity1981_2013.csv",na.strings = c(""))
hunter.abun<-read.csv("Data/DMU_HunterAbundance1975_2013.csv",na.strings = c(""))
hunter.dens<-read.csv("Data/DMU_HunterDensity1975_2013.csv",na.strings = c(""))
hunter.dens1416<-read.csv("Data/hunterdens1416.csv",na.strings = c(""))
hunter.abun1416<-read.csv("Data/hunterabun1416.csv",na.strings = c(""))
cty.abun.1416 <- read.csv("Data/Cty_Abun_2014_2016.csv",na.strings = c(""))

#CLEAN DATA

#fix county names
levels(crash$County)[levels(crash$County)=="Sheboygon"] <- "Sheboygan"
levels(adj.crash$County)[levels(adj.crash$County)=="Sheboygon"] <- "Sheboygan"
levels(usda$County)[levels(usda$County)=="Saint Croix"] <- "St. Croix"
levels(usda$County)[levels(usda$County)=="Fond du Lac"] <- "Fond Du Lac"

#eliminate dashes in DMU values
abun.est$Unit<-noquote(gsub("-", "", abun.est$Unit))
dens.est$Unit<-noquote(gsub("-", "", dens.est$Unit))
names(crash)[names(crash)=="County"]<-"COUNTY_NAM"
names(adj.crash)[names(adj.crash)=="County"]<-"COUNTY_NAM"

names(traffic)[names(traffic)=="COUNTY"]<-"COUNTY_NAM"
hunter.abun<-hunter.abun[,!names(hunter.abun)%in%c("X60_64","X65_69","X70_74","X1975","X1976","X1977","X1978","X1979","X1980","X")] #eliminate years prior to 1981, and notes column
hunter.dens<-hunter.dens[,!names(hunter.dens)%in%c("X60_64","X65_69","X70_74","X1975","X1976","X1977","X1978","X1979","X1980")] #eliminate years prior to 1981
hunter.dens1416$County<-as.factor(toTitleCase(tolower(as.character(hunter.dens1416$County)))) #make names match in capitalization
hunter.abun1416$County<-as.factor(toTitleCase(tolower(as.character(hunter.abun1416$County)))) #make names match in capitalization

#generate lists of cwd counties and zones
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

#convert to long format
crashL<-melt(crash, id.vars = c("COUNTY_NAM"),variable.name="Year",value.name="Crashes")
adj.crashL<-melt(adj.crash,id.vars = c("COUNTY_NAM"),variable.name="Year",value.name="Crashes")
trafficL<-melt(traffic,id.vars = c("COUNTY_NAM"),variable.name="Year",value.name="Traffic")

#make year column
crashL$Year<-as.numeric(sub('X', '', crashL$Year))
trafficL$Year<-as.numeric(sub('X', '', trafficL$Year))
adj.crashL$Year<-as.numeric(sub('X', '', adj.crashL$Year))
names(adj.crashL)[names(adj.crashL)=="Crashes"]<-"AdjDVC"

#load DMU-Cty-Deer Range merged files
cty.flag<-"dez"
source('Code/DMU_Cty_DataStep.R')

cty.dmu.dr<-aggregate(AREA_SQM~UNIT_ID+COUNTY_NAM+Year,data=cty.dmu.dr,FUN=sum)
levels(cty.dmu.dr$COUNTY_NAM)[levels(cty.dmu.dr$COUNTY_NAM)=="Fond du Lac"] <- "Fond Du Lac"
levels(cty.dmu.dr$COUNTY_NAM)[levels(cty.dmu.dr$COUNTY_NAM)=="Saint Croix"] <- "St. Croix"

#calculate the amount of deer range by DMU, County and for DMUs used in the population estimates
deer.range <- as.data.frame(rbind(by(cty.dmu.dr$AREA_SQM,list(cty.dmu.dr$UNIT_ID,cty.dmu.dr$Year), sum,na.rm=T)))
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

#Add year column
dens.est2L$Year<-as.numeric(sub('X', '', dens.est2L$Year))
hunter.est2L$Year<-as.numeric(sub('X', '', hunter.est2L$Year))

#combine DMUs and Counties
dmu.cty.estL<-merge(dens.est2L,cty.dmu.dr,by.x=c("Unit","Year"),by.y=c("UNIT_ID","Year"),all=T)
dmu.cty.hunter.estL<-merge(hunter.est2L,cty.dmu.dr,by.x=c("Unit","Year"),by.y=c("UNIT_ID","Year"),all=T)

#generate abundance estimates by multiplying deer area by deer density
dmu.cty.estL$Abun<-dmu.cty.estL$Dens*dmu.cty.estL$AREA_SQM
dmu.cty.hunter.estL$Hunters<-dmu.cty.hunter.estL$HunterDens*dmu.cty.hunter.estL$AREA_SQM

#aggregate abundance by county
cty.abun.estL<-aggregate(Abun~COUNTY_NAM+Year,FUN=sum,data=dmu.cty.estL)

#aggregate hunters by county
cty.hunter.estL<-aggregate(Hunters~COUNTY_NAM+Year,FUN=sum,data=dmu.cty.hunter.estL)

#aggregate abundance for 2014-2016 data
abun1416 <- by(cty.abun.1416$PostHuntPop, INDICES=list(cty.abun.1416$County,cty.abun.1416$Year), FUN=sum)
cty.abun1.1416<-cbind(expand.grid(attributes(abun1416)$dimnames), as.vector(abun1416))
names(cty.abun1.1416)<-c("COUNTY_NAM","Year","Abun")

#aggregate hunters for 2014-2016 data
hunter1416 <- by(hunter.abun1416$HunterAbun, INDICES=list(hunter.abun1416$County, hunter.abun1416$Year), FUN=sum)
hunters1416L<-cbind(expand.grid(attributes(hunter1416)$dimnames), as.vector(hunter1416))
names(hunters1416L)<-c("COUNTY_NAM","Year","Hunters")

#combine 1981-2013 with 2014-2016 abundance
cty.abun.estL<-rbind(cty.abun.estL,cty.abun1.1416)
cty.abun.estL$Abun[is.na(cty.abun.estL$Abun)]<-0 #make two menominee NA values to 0

#combine 1981-2013 with 2014-2016 hunters
cty.hunter.estL<-rbind(cty.hunter.estL,hunters1416L)

#generate data in wide format
cty.abun.est<- as.data.frame(rbind(by(cty.abun.estL$Abun,list(cty.abun.estL$COUNTY_NAM,cty.abun.estL$Year), sum,na.rm=T)))
hunter.abun.est <- as.data.frame(rbind(by(cty.hunter.estL$Hunters,list(cty.hunter.estL$COUNTY_NAM,cty.hunter.estL$Year), sum,na.rm=T)))

#combine abudance data with hunter data
dat<-merge(cty.abun.estL,cty.hunter.estL,by=c("COUNTY_NAM","Year"),all=T)

#combine with traffic and DVC data
dat<-merge(dat,crashL,by=c("COUNTY_NAM","Year"),all=T)
dat<-merge(dat,trafficL,by=c("COUNTY_NAM","Year"),all=T)
dat<-merge(dat,adj.crashL,by=c("COUNTY_NAM","Year"),all=T)

#Combine with usda districts
dat<-merge(dat,usda,by.x=c("COUNTY_NAM"),by.y=c("County"),all=T)

#Drop Menominee county because it has very little data for anything (its an Indian reservation)
dat<-dat[!(dat$COUNTY_NAM=="Menominee"),]

rm(dmu.files,cty.abun.est,cty.abun.estL,cty.hunter.estL,hunter.abun.est,dmu.cty.estL,dmu.cty.hunter.estL)

#aggregate data by usda districts
usda.abun <- by(dat$Abun, INDICES=list(dat$Year, dat$Zone), FUN=sum)
usda.abun<-cbind(expand.grid(attributes(usda.abun)$dimnames), as.vector(usda.abun))
names(usda.abun)<-c("Year","Zone","Abun")

usda.hunters <- by(dat$Hunters, INDICES=list(dat$Year, dat$Zone), FUN=sum,na.rm=T)
usda.hunters<-cbind(expand.grid(attributes(usda.hunters)$dimnames), as.vector(usda.hunters))
names(usda.hunters)<-c("Year","Zone","Hunters")

usda.dvc <- by(dat$Crashes, INDICES=list(dat$Year, dat$Zone), FUN=sum)
usda.dvc<-cbind(expand.grid(attributes(usda.dvc)$dimnames), as.vector(usda.dvc))
names(usda.dvc)<-c("Year","Zone","Crashes")

usda.traf <- by(dat$Traffic, INDICES=list(dat$Year, dat$Zone), FUN=sum)
usda.traf<-cbind(expand.grid(attributes(usda.traf)$dimnames), as.vector(usda.traf))
names(usda.traf)<-c("Year","Zone","Traffic")

usda.dat<-merge(usda.abun,usda.hunters,by=c("Zone","Year"))
usda.dat<-merge(usda.dat,usda.dvc,by=c("Zone","Year"))
usda.dat<-merge(usda.dat,usda.traf,by=c("Zone","Year"))

#generate traffic-adjusted crashes
usda.dat$AdjDVC<-usda.dat$Crashes/usda.dat$Traffic

#Make sure only years need are selected
dat<-subset(dat,Year%in%c(minyear:maxyear))
usda.dat<-usda.dat[usda.dat$Year%in%c(minyear:maxyear),]

#convert to wide format
cty.mat<-matrix(NA,length(unique(dat$COUNTY_NAM)),length(minyear:maxyear))
cty.list<-list()
vars<-c("Abun","Hunters","Crashes","Traffic","AdjDVC")
for(j in vars){
  for(i in 1:length(unique(dat$COUNTY_NAM))){
    dat.sub<-dat[as.numeric(as.factor(as.character(dat$COUNTY_NAM)))==i,]
    dat.sub<-dat.sub[order(dat.sub$Year),]
    cty.mat[i,]<-dat.sub[,j]
  }
  cty.list[[j]]<-cty.mat
}
cty.list<-lapply(cty.list,function(x){row.names(x)<-unique(dat$COUNTY_NAM);x})
saveRDS(cty.list,"Results/cty.list.rds")

#Generate raw district data
usda.mat<-matrix(NA,length(unique(usda.dat$Zone)),length(minyear:maxyear))
usda.list<-list()
vars<-c("Abun","Hunters","Crashes","Traffic","AdjDVC")
for(j in vars){
  for(i in 1:length(unique(usda.dat$Zone))){
    dat.sub<-usda.dat[as.numeric(usda.dat$Zone)==i,]
    dat.sub<-dat.sub[order(dat.sub$Year),]
    usda.mat[i,]<-dat.sub[,j]
  }
  usda.list[[j]]<-usda.mat
}
usda.list<-lapply(usda.list,function(x){row.names(x)<-sort(unique(usda$Zone));x})
saveRDS(usda.list,"Results/usda.list.rds")
