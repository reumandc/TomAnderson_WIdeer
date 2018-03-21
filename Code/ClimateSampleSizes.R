samplesize<-matrix(NA,ncol=8,nrow=length(unique(wi.clim$STATION)))
month.list<-list()
all.sites.months<-list()
tmin.months<-vector()
tmax.months<-vector()
pcp.months<-vector()
snwd.months<-vector()

for(i in 1:length(unique(wi.clim$STATION)))
{
  #subset out each station
  stat<-subset(wi.clim,as.numeric(STATION)==i) 
  
  #calculate data series lengths, number of NAs, and total observations
  firstday<-min(stat$DATE) 
  lastday<-max(stat$DATE)
  datelength<-as.numeric(as.Date(as.character(lastday))-as.Date(as.character(firstday)))
  nobs<-dim(stat)[1]
  tmin.na<-sum(is.na(stat$TMIN))
  tmax.na<-sum(is.na(stat$TMAX))
  pcp.na<-sum(is.na(stat$PRCP))
  snwd.na<-sum(is.na(stat$SNWD))
  dat<-c(firstday,lastday,datelength,nobs,tmin.na,tmax.na,pcp.na,snwd.na)
  samplesize[i,]<-dat
  
  #extract monthly totals
  month.dat<-array(NA,c(length(min(stat$Year):max(stat$Year)),12,4))
  for(j in min(stat$Year):max(stat$Year))
  {
    for(k in 1:12)
    {
      tmin.months[k]<-sum(!is.na(stat$TMIN[stat$Year==j & stat$Month==k]))
      tmax.months[k]<-sum(!is.na(stat$TMAX[stat$Year==j & stat$Month==k]))
      pcp.months[k]<-sum(!is.na(stat$PRCP[stat$Year==j & stat$Month==k]))
      snwd.months[k]<-sum(!is.na(stat$SNWD[stat$Year==j & stat$Month==k]))
      
    }
    #month.dat[j-(min(stat$Year)-1),]<-c(tmin.months,tmax.months)
    month.dat[j-(min(stat$Year)-1),,1]<-tmin.months
    month.dat[j-(min(stat$Year)-1),,2]<-tmax.months
    month.dat[j-(min(stat$Year)-1),,3]<-pcp.months
    month.dat[j-(min(stat$Year)-1),,4]<-snwd.months
    #colnames(month.dat)<-rep(c("Ja","Fe","Ma","Ap","Ma","Jn","Jl","Au","Se","Oc","No","De"),1)
  }
  all.sites.months[[i]]<-month.dat
  
  samplesize<-data.frame(samplesize)
  samplesize$STATION<-as.factor(unique(sort(wi.clim$STATION)))
  names(samplesize)<-c("Begin","End","PotentialOBS","Nobs","TminNA","TmaxNA","PrcpNA","SnwdNA","STATION")
  samplesize$PotentialOBS<-as.numeric(as.character(samplesize$PotentialOBS))
  samplesize$Nobs<-as.numeric(as.character(samplesize$Nobs))
  samplesize$TminNA<-as.numeric(as.character(samplesize$TminNA))
  samplesize$TmaxNA<-as.numeric(as.character(samplesize$TmaxNA))
  samplesize$PrcpNA<-as.numeric(as.character(samplesize$PrcpNA))
  samplesize$SnwdNA<-as.numeric(as.character(samplesize$SnwdNA))
  samplesize$NmbYrs<-sapply(all.sites.months,function(x){dim(x)[1]})
  cty.samplesize<-merge(samplesize,noaa.coord1,by="STATION",all.x=T)
  names(all.sites.months)<-sort(unique(wi.temp$STATION))
  
  compl.stat<-samplesize$STATION[samplesize$NmbYrs>34 & samplesize$Nobs>12700] #41 stations/34 counties with complete data at 35 years
  
  aggregate(cbind(PotentialOBS,Nobs,TminNA,TmaxNA)~COUNTY_NAM,FUN=sum,data=cty.samplesize)#total observations by county
  aggregate(STATION~COUNTY_NAM,length,data=cty.samplesize)#number of stations per county
}
#trim climate data by stations with the most complete data
wi.clim.compl<-wi.clim[wi.clim$STATION %in% compl.stat,]
length(unique(wi.clim.compl$COUNTY_NAM)) # n = 46
length(unique(wi.clim.compl$STATION)) # n = 66
