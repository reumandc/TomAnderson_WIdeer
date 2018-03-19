#Code to generate abundance of deer for DMUs

#To use files where normal DMU boundaries are used
if(cty.flag=="norm"){
  drt<-"Data/DMU_Cty_DEZ/"
  temp = list.files(path=drt,pattern="*.csv")
  dmu.files = lapply(temp, function(x){read.csv(paste(drt,x,sep="/"),header=T)})
  years<-1:33
  dmu.list<-list()
  for(i in 1:length(dmu.files)){
    tmp<-dmu.files[[i]]
    tmp$Year<-rep(years[i],dim(tmp)[1])+1980
    dmu.list[[i]]<-tmp
  }
  for(i in 1:length(dmu.list)){
    if(i<22){
      names(dmu.list[[i]])<-names(dmu.list[[i]])
    }
    else{
      names(dmu.list[[i]])<-sub("^DEER_MGMT_$", "UNIT_ID", names(dmu.list[[i]]))
    }
  }
  lapply(dmu.list,names)
  dmu.names<-c('COUNTY_NAM', 'UNIT_ID', 'Year','AREA_SQM') 
  dmu.list<-lapply(dmu.list,function(x){x[,dmu.names]})
  lapply(dmu.list,dim)
  cty.dmu.dr<-Reduce(function(x, y) merge(x, y,all=TRUE), dmu.list)
  cty.dmu.dr$UNIT_ID<-noquote(gsub("-", "", cty.dmu.dr$UNIT_ID))
}

#To use files where deer eradication zone its own DMU
if(cty.flag=="dez"){
  drt<-"Data/DMU_Cty_DEZ/"
  temp = list.files(path=drt,pattern="*.csv")
  dmu.files = lapply(temp, function(x){read.csv(paste(drt,x,sep="/"),header=T)})
  years<-1:33
  dmu.list<-list()
  for(i in 1:length(dmu.files)){
    tmp<-dmu.files[[i]]
    tmp$Year<-rep(years[i],dim(tmp)[1])+1980
    dmu.list[[i]]<-tmp
  }
  for(i in 1:length(dmu.list)){
    if(i<22){
      names(dmu.list[[i]])<-names(dmu.list[[i]])
    }
    else{
      names(dmu.list[[i]])<-sub("^DEER_MGMT_$", "UNIT_ID", names(dmu.list[[i]]))
    }
  }
  #replace DMU Units IDs with SWDEZ or SEDEZ for eradication to units, to match with DMU population files
  for(i in 1:length(dmu.list)){
    if(dim(dmu.list[[i]])[2]<9){
      dmu.list[[i]]$UNIT_ID<-ifelse(dmu.list[[i]]$UNIT_DEZ=="SW-DEZ","SW-DEZ",as.character(dmu.list[[i]]$UNIT_ID))
      dmu.list[[i]]$UNIT_ID<-ifelse(dmu.list[[i]]$UNIT_DEZ=="SE-DEZ","SE-DEZ",as.character(dmu.list[[i]]$UNIT_ID))
    }
  }
  dmu.names<-c('COUNTY_NAM', 'UNIT_ID', 'Year','AREA_SQM') 
  dmu.list<-lapply(dmu.list,function(x){x[,dmu.names]})
  lapply(dmu.list,dim)
  cty.dmu.dr<-Reduce(function(x, y) merge(x, y,all=TRUE), dmu.list)
  cty.dmu.dr$UNIT_ID<-noquote(gsub("-", "", cty.dmu.dr$UNIT_ID))
}