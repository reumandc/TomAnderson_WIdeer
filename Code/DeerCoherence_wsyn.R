#Analysis with wsyn
library(wsyn)
bands<-rbind(c(3,7),c(3,4),c(4,7))

#Run coherence of climate indices and deer abundance
abun.dt<-cleandat(cty.list$Abun,clev=5,times=1981:2016)$cdat
climindex.dt<-lapply(climindex[c("WinterNAO","WinterPDO","WinterMEI","SummerNAO","SummerPDO","SummerMEI")],function(x){x<-cleandat(x,times=1981:2016,clev=5)$cdat;x})
clim.res<-list()
wlm_climabun<-list()
for(j in names(climindex.dt)){
  spatcoh.names<-paste(j,"Abun",sep=".")
  clim.res[[spatcoh.names]]<-coh(dat1=abun.dt,dat2=climindex.dt[[j]],times=minyear:maxyear,norm="powall",
                sigmethod="fast",nrand=10000,f0=1)
  mod<-wlm(list(abun.dt,climindex.dt[[j]]),times=minyear:maxyear,resp=1,pred=2,norm="powall")
  wlm_climabun[[spatcoh.names]]<-syncexpl(mod)
}
climabun_se<-list()
for(j in names(clim.res)){
  se_tmp<-matrix(NA,ncol=ncol(wlm_climabun[[1]])-2,nrow=nrow(bands))
  for(i in 1:nrow(bands)){
    clim.res[[j]]<-bandtest(clim.res[[j]],bands[i,])
    tmp<-wlm_climabun[[j]][wlm_climabun[[j]]$timescales>=bands[i,1] & wlm_climabun[[j]]$timescales<=bands[i,2],]
    se_tmp[i,]<-round(100*colMeans(tmp[,c(3:dim(tmp)[2])])/mean(tmp$sync),4)
  }
  climabun_se[[j]]<-se_tmp
}
clim.resP<-lapply(clim.res,get_bandp)
climabun_se<-lapply(climabun_se,function(x){colnames(x)<-colnames(tmp)[3:dim(tmp)[2]];x})
climabun_se<-lapply(climabun_se,function(x){row.names(x)<-paste(bands[,1],bands[,2],sep="-");x})


#Run coherence of winter weather and deer abundance
weath.res<-list()
wlm_weathabun<-list()
for(j in names(winter.clim)){
  winter.clim.dt.tmp<-winter.clim[[j]][!is.na(rowMeans(winter.clim[[j]])),]
  winter.clim.dt.tmp<-cleandat(winter.clim.dt.tmp,clev=5,times=minyear:maxyear)$cdat
  abun.dt.tmp<-abun.dt[!is.na(rowMeans(winter.clim[[j]])),]
  spatcoh.names<-paste(j,"Abun",sep=".")
  weath.res[[spatcoh.names]]<-coh(dat1=abun.dt.tmp,dat2=winter.clim.dt.tmp,times=minyear:maxyear,norm="powall",
                sigmethod="fast",nrand=10000,f0=1)
  mod<-wlm(list(abun.dt.tmp,winter.clim.dt.tmp),times=minyear:maxyear,resp=1,pred=2,norm="powall")
  wlm_weathabun[[spatcoh.names]]<-syncexpl(mod)
}
weathabun_se<-list()
for(j in names(weath.res)){
  se_tmp<-matrix(NA,ncol=ncol(wlm_weathabun[[1]])-2,nrow=nrow(bands))
  for(i in 1:nrow(bands)){
    weath.res[[j]]<-bandtest(weath.res[[j]],bands[i,])
    tmp<-wlm_weathabun[[j]][wlm_weathabun[[j]]$timescales>=bands[i,1] & wlm_weathabun[[j]]$timescales<=bands[i,2],]
    se_tmp[i,]<-round(100*colMeans(tmp[,c(3:dim(tmp)[2])])/mean(tmp$sync),4)
  }
  weathabun_se[[j]]<-se_tmp
}
weath.resP<-lapply(weath.res,get_bandp)
weathabun_se<-lapply(weathabun_se,function(x){colnames(x)<-colnames(tmp)[3:dim(tmp)[2]];x})
weathabun_se<-lapply(weathabun_se,function(x){row.names(x)<-paste(bands[,1],bands[,2],sep="-");x})

#Run coherence of climate indices and winter weather
TableS3<-matrix(NA,nrow=length(names(winter.clim)),ncol=length(names(climindex)))
for(j in 1:length(names(climindex))){
  for(i in 1:length(names(winter.clim))){
    name1<-names(climindex)[j]
    name2<-names(winter.clim)[i]
    climindex.tmp<-cleandat(climindex[[name1]],times=minyear:maxyear,clev=5)$cdat
    winter.clim.dt<-winter.clim[[name2]][!is.na(rowMeans(winter.clim[[name2]])),]
    climindex.tmp<-climindex.tmp[!is.na(rowMeans(winter.clim[[name2]])),]
    winter.clim.dt.tmp<-cleandat(winter.clim.dt,times=minyear:maxyear,clev=5)$cdat
    spatcoh.names<-paste(name1,name2,sep=".")
    weath.climind.res[[spatcoh.names]]<-coh(dat2=climindex.tmp,dat1=winter.clim.dt.tmp,times=minyear:maxyear,norm="powall",
                                            sigmethod="fast",nrand=10000,f0=1)
  }
}

for(j in names(weath.climind.res)){
  for(i in 1:nrow(bands)){
    weath.climind.res[[j]]<-bandtest(weath.climind.res[[j]],bands[i,])
  }
}
lapply(weath.climind.res,get_bandp)

#Run wavelet coherence among climate indices
TableS2<-matrix(NA,6,6)
for(i in 1:(length(climindex.dt)-1)){
  for(j in (i+1):length(climindex.dt)){
    name1<-names(climindex.dt)[i]
    name2<-names(climindex.dt)[j]
    spatcoh.names<-paste(name2,name1,sep=".")
    indices.res[[spatcoh.names]]<-coh(dat1=climindex.dt[[j]][1,],dat2=climindex.dt[[i]][1,],times=minyear:maxyear,norm="powall",
                                      sigmethod="fast",nrand=10000,f0=1)
  }
}
for(j in names(indices.res)){
  for(i in 1:nrow(bands)){
    indices.res[[j]]<-bandtest(indices.res[[j]],bands[i,])
  }
}
lapply(indices.res,get_bandp)

diag(TableS2)<-1
row.names(TableS2)<-names(climindex.dt)
colnames(TableS2)<-names(climindex.dt)
saveRDS(TableS2,file="Results/TabS2_results.rds")

## Run wavelet multiple regression model for deer abundance
# first consolidate data to dimensions of winter snow depth (60 counties)
snwd<-winter.clim$Snwd[!is.na(rowMeans(winter.clim$Snwd)),]
abun<-cty.list$Abun[!is.na(rowMeans(winter.clim$Snwd)),]
climindex.tmp<-lapply(climindex,function(x){x[!is.na(rowMeans(winter.clim$Snwd)),]})

#store data in list
all.dat<-list(abun,snwd,climindex.tmp$WinterMEI,climindex.tmp$WinterPDO)

#normalize data
norm.dat<-lapply(all.dat,function(x){x<-cleandat(x,clev=5,times=minyear:maxyear)$cdat;x})

#Run model and calculate expected synchrony
wlm_abun<-wlm(norm.dat,times=minyear:maxyear,resp=1,pred=2:4,norm="powall")
abun_se<-syncexpl(wlm_abun)
se_short<-abun_se[abun_se$timescales>=3 & abun_se$timescales<=7,]
round(100*colMeans(se_short[,c(3:12)])/mean(se_short$sync),4)

##Run coherence of hunters and abundance
#first filter data to match dimensions of hunter data
cty.list.dt<-lapply(cty.list[c('Abun','Hunters')],function(x){x<-x[,12:36];x})
cty.list.dt<-lapply(cty.list.dt,function(x){x[(!row.names(cty.list.dt$Abun)%in%cwd) & !is.na(rowMeans(cty.list.dt$Hunters)),]})
cty.list.dt<-lapply(cty.list.dt,function(x){x<-cleandat(x,clev=5,times=1992:2016)$cdat;x})

#run coherence
hunter.res<-coh(dat1=cty.list.dt$Abun,dat2=cty.list.dt$Hunters,times=1992:maxyear,norm="powall",
                sigmethod="fast",nrand=10000,f0=1)
hunter.bands<-rbind(c(3,7),c(2,2.5))

for(i in 1:nrow(hunter.bands)){
  hunter.res<-bandtest(hunter.res,hunter.bands[i,])
}
get_bandp(hunter.res)

##Run coherence of DVCs and abundance
#first filter data to match dimensions of DVC data
cty.list.dt<-lapply(cty.list[c('Abun','Crashes')],function(x){x<-cleandat(x[,!is.na(colSums(cty.list$Crashes))],clev = 5,times=1987:2016)$cdat;x})
dvc.res<-coh(dat1=cty.list.dt$Crashes,dat2=cty.list.dt$Abun,times=1987:maxyear,norm="powall",
             sigmethod="fast",nrand=10000,f0=1)
dvc.res<-bandtest(dvc.res,c(3,7))
get_bandp(dvc.res)
wlm_dvc<-wlm(list(cty.list.dt$Crashes,cty.list.dt$Abun),times=1987:maxyear,resp=1,pred=2,norm="powall",f0=1)
dvc_se<-syncexpl(wlm_dvc)
dvcabun_se<-dvc_se[dvc_se$timescales>=3 & dvc_se$timescales<=7,]
round(100*colMeans(dvcabun_se[,c(3:dim(dvcabun_se)[2])])/mean(dvcabun_se$sync),4)

#Run coherence of traffic-adjusted DVCs and abundance,again filtering to match dimensions of adjusted DVCs
cty.list.dt<-lapply(cty.list[c('Abun','AdjDVC')],function(x){x<-cleandat(x[,!is.na(colSums(cty.list$AdjDVC))],clev = 5,times=1988:2016)$cdat;x})
adjdvc.res<-coh(dat1=cty.list.dt$AdjDVC,dat2=cty.list.dt$Abun,times=1988:maxyear,norm="powall",
                sigmethod="fast",nrand=10000,f0=1)
adjdvc.res<-bandtest(adjdvc.res,c(3,7))
get_bandp(adjdvc.res)
wlm_adjdvc<-wlm(list(cty.list.dt$AdjDVC,cty.list.dt$Abun),times=1988:maxyear,resp=1,pred=2,norm="powall")
adjdvc_se<-syncexpl(wlm_adjdvc)
adjdvcabun_se<-adjdvc_se[adjdvc_se$timescales>=3 & adjdvc_se$timescales<=7,]
round(100*colMeans(adjdvcabun_se[,c(3:dim(adjdvcabun_se)[2])])/mean(adjdvcabun_se$sync),4)

Tab1<-cbind(rbind(do.call(rbind,weath.resP), 
            do.call(rbind,clim.resP),
            unlist(get_bandp(dvc.res)),
            unlist(get_bandp(adjdvc.res))),matrix(NA,nrow=94,ncol=4))
rbind(dvcabun_se,adjdvcabun_se)
for(j in 1:ncol(Tab1)){
  for(i in 1:nrow(Tab1)){
    if(Tab1[i,2]<0.05){
      Tab1[i,j]<-""
    }
  }
}


# USDA Analysis -----------------------------------------------------------

#Run coherence of traffic-adjusted DVCs and abundance,again filtering to match dimensions of adjusted DVCs

usda.list.dt<-lapply(usda.list[c('Abun','Crashes')],function(x){x<-cleandat(x[,!is.na(colSums(usda.list$Crashes))],clev = 5,times=1987:2016)$cdat;x})
usda.dvc.res<-coh(dat1=usda.list.dt$Crashes,dat2=usda.list.dt$Abun,times=1987:maxyear,norm="powall",
                  sigmethod="fast",nrand=10000,f0=1)
usda.dvc.res<-bandtest(usda.dvc.res,c(3,7))
get_bandp(usda.dvc.res)
usda.wlm_dvc<-wlm(list(usda.list.dt$Crashes,usda.list.dt$Abun),times=1987:maxyear,resp=1,pred=2,norm="powall")
usda.dvc_se<-syncexpl(usda.wlm_dvc)
usda.dvcabun_se<-usda.dvc_se[usda.dvc_se$timescales>=3 & usda.dvc_se$timescales<=7,]
round(100*colMeans(usda.dvcabun_se[,c(3:dim(usda.dvcabun_se)[2])])/mean(usda.dvcabun_se$sync),4)

#Run coherence of traffic-adjusted DVCs and abundance,again filtering to match dimensions of adjusted DVCs
usda.list.dt<-lapply(usda.list[c('Abun','AdjDVC')],function(x){x<-cleandat(x[,!is.na(colSums(cty.list$AdjDVC))],clev = 5,times=1988:2016)$cdat;x})
usda.adjdvc.res<-coh(dat1=usda.list.dt$AdjDVC,dat2=usda.list.dt$Abun,times=1988:maxyear,norm="powall",
                sigmethod="fast",nrand=10000,f0=1)
usda.adjdvc.res<-bandtest(adjdvc.res,c(3,7))
get_bandp(adjdvc.res)
usda.wlm_adjdvc<-wlm(list(usda.list.dt$AdjDVC,usda.list.dt$Abun),times=1988:maxyear,resp=1,pred=2,norm="powall")
usda.adjdvc_se<-syncexpl(usda.wlm_adjdvc)
usda.adjdvcabun_se<-usda.adjdvc_se[usda.adjdvc_se$timescales>=3 & usda.adjdvc_se$timescales<=7,]
usda.adjdvcabun_se<-round(100*colMeans(usda.adjdvcabun_se[,c(3:dim(usda.adjdvcabun_se)[2])])/mean(usda.adjdvcabun_se$sync),4)



