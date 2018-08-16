#Analysis with wsyn
bands<-rbind(c(3,7),c(3,4),c(4,7))

#Run coherence of climate indices and deer abundance
abun.dt<-cleandat(cty.list$Abun,clev=5,times=1981:2016)$cdat
climindex.dt<-lapply(climindex,function(x){x<-cleandat(x,times=1981:2016,clev=5)$cdat;x})
clim.res<-list()
wlm_climabun<-list()
for(j in names(climindex.dt)){
  spatcoh.names<-paste(j,"Abun",sep=".")
  clim.res[[spatcoh.names]]<-coh(dat1=abun.dt,dat2=climindex.dt[[j]],times=minyear:maxyear,norm="powall",
                sigmethod="fast",nrand=nsurrogs,f0=1)
  mod<-wlm(list(abun.dt,climindex.dt[[j]]),times=minyear:maxyear,resp=1,pred=2,norm="powall")
  wlm_climabun[[spatcoh.names]]<-syncexpl(mod)
}

#Calculate abundance p-values and synchrony explained by climate indices
climabun_se<-list()
for(j in names(clim.res)){
  se_mat<-matrix(NA,ncol=ncol(wlm_climabun[[1]])-2,nrow=nrow(bands))
  clim.res[[j]]<-bandtest(clim.res[[j]],bands[1,])
  tmp.se<-wlm_climabun[[j]][wlm_climabun[[j]]$timescales>=bands[1,1] & wlm_climabun[[j]]$timescales<=bands[1,2],]
  se_mat[1,]<-round(100*colMeans(tmp.se[,c(3:dim(tmp.se)[2])])/mean(tmp.se$sync),4)
  if(clim.res[[j]]$bandp['p_val']<0.06){
    for(i in 2:nrow(bands)){
      tmp.p<-bandtest(clim.res[[j]],bands[i,])
      clim.res[[j]]<-tmp.p
      tmp.se<-wlm_climabun[[j]][wlm_climabun[[j]]$timescales>=bands[i,1] & wlm_climabun[[j]]$timescales<=bands[i,2],]
      se_mat[i,]<-round(100*colMeans(tmp.se[,c(3:dim(tmp.se)[2])])/mean(tmp.se$sync),4)
    }
  }
  climabun_se[[j]]<-se_mat
}

clim.resP<-lapply(clim.res,get_bandp)
climabun_se<-lapply(climabun_se,function(x){colnames(x)<-colnames(tmp.se)[3:dim(tmp.se)[2]];x})
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
                sigmethod="fast",nrand=nsurrogs,f0=1)
  mod<-wlm(list(abun.dt.tmp,winter.clim.dt.tmp),times=minyear:maxyear,resp=1,pred=2,norm="powall")
  wlm_weathabun[[spatcoh.names]]<-syncexpl(mod)
}

#Calculate abundance synchrony explained by weather
weathabun_se<-list()
for(j in names(weath.res)){
  se_mat<-matrix(NA,ncol=ncol(wlm_weathabun[[1]])-2,nrow=nrow(bands))
  weath.res[[j]]<-bandtest(weath.res[[j]],bands[1,])
  tmp.se<-wlm_weathabun[[j]][wlm_weathabun[[j]]$timescales>=bands[1,1] & wlm_weathabun[[j]]$timescales<=bands[1,2],]
  se_mat[1,]<-round(100*colMeans(tmp.se[,c(3:dim(tmp.se)[2])])/mean(tmp.se$sync),4)
  if(weath.res[[j]]$bandp['p_val']<0.06){
    for(i in 2:nrow(bands)){
      tmp.p<-bandtest(weath.res[[j]],bands[i,])
      weath.res[[j]]<-tmp.p
      tmp.se<-wlm_weathabun[[j]][wlm_weathabun[[j]]$timescales>=bands[i,1] & wlm_weathabun[[j]]$timescales<=bands[i,2],]
      se_mat[i,]<-round(100*colMeans(tmp.se[,c(3:dim(tmp.se)[2])])/mean(tmp.se$sync),4)
    }
  }
  weathabun_se[[j]]<-se_mat
}

weath.resP<-lapply(weath.res,get_bandp)
weathabun_se<-lapply(weathabun_se,function(x){colnames(x)<-colnames(tmp.se)[3:dim(tmp.se)[2]];x})
weathabun_se<-lapply(weathabun_se,function(x){row.names(x)<-paste(bands[,1],bands[,2],sep="-");x})

#Run coherence of climate indices and winter weather
weath.climind.res<-list()
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
                                            sigmethod="fast",nrand=nsurrogs,f0=1)
  }
}

for(j in names(weath.climind.res)){
  weath.climind.res[[j]]<-bandtest(weath.climind.res[[j]],bands[1,])
}
TabS3<-do.call(rbind,lapply(weath.climind.res,get_bandp))
TabS3<-cbind(matrix(unlist(strsplit(row.names(TabS3),'.',fixed=TRUE)),ncol=2,byrow=T),TabS3)
TabS3$mn_phs<-as.numeric(as.character(ifelse(TabS3$p_val>0.06,"",TabS3$mn_phs)))
TableS3<-data.frame(TabS3[,1:2],paste(TabS3$ts_low_bd,TabS3$ts_hi_bd,sep="-"),TabS3[,c("p_val","mn_phs")])
colnames(TableS3)<-c("Predictor","Response","Timescale","Pvalue","MeanPhase")
saveRDS(TableS3,file="Results/TabS3_results.rds")

#Run wavelet coherence among climate indices
indices.res<-list()
for(i in 1:(length(climindex.dt)-1)){
  for(j in (i+1):length(climindex.dt)){
    name1<-names(climindex.dt)[i]
    name2<-names(climindex.dt)[j]
    spatcoh.names<-paste(name2,name1,sep=".")
    indices.res[[spatcoh.names]]<-coh(dat1=climindex.dt[[j]][1,],dat2=climindex.dt[[i]][1,],times=minyear:maxyear,norm="powall",
                                      sigmethod="fast",nrand=nsurrogs,f0=1)
  }
}
for(j in names(indices.res)){
  indices.res[[j]]<-bandtest(indices.res[[j]],bands[1,])
}

TabS2<-do.call(rbind,lapply(indices.res,get_bandp))
TabS2<-cbind(matrix(unlist(strsplit(row.names(TabS2),'.',fixed=TRUE)),ncol=2,byrow=T),TabS2)
TabS2$mn_phs<-as.numeric(as.character(ifelse(TabS2$p_val>0.06,"",TabS2$mn_phs)))
TableS2<-data.frame(TabS2[,1:2],paste(TabS2$ts_low_bd,TabS2$ts_hi_bd,sep="-"),TabS2[,c("p_val","mn_phs")])
colnames(TableS2)<-c("Predictor","Response","Timescale","Pvalue","MeanPhase")
TableS2
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
abunmod_se<-round(100*colMeans(se_short[,c(3:12)])/mean(se_short$sync),4)
saveRDS(abunmod_se[[1]],file="Results/abunmodelSyncExp.rds")

##Run coherence of hunters and abundance
#first filter data to match dimensions of hunter data
cty.list.dt<-lapply(cty.list[c('Abun','Hunters')],function(x){x<-x[,12:36];x})
cty.list.dt<-lapply(cty.list.dt,function(x){x[(!row.names(cty.list.dt$Abun)%in%cwd) & !is.na(rowMeans(cty.list.dt$Hunters)),]})
cty.list.dt<-lapply(cty.list.dt,function(x){x<-cleandat(x,clev=5,times=1992:2016)$cdat;x})

#run coherence
hunter.bands<-rbind(c(3,7),c(2,2.5))
hunter.res<-coh(dat1=cty.list.dt$Abun,dat2=cty.list.dt$Hunters,times=1992:maxyear,norm="powall",
                sigmethod="fast",nrand=nsurrogs,f0=1)
wlm_hunters<-wlm(list(cty.list.dt$Abun,cty.list.dt$Hunters),times=1992:maxyear,resp=1,pred=2,norm="powall",f0=1)
hunters_se<-syncexpl(wlm_hunters)

hunter_se<-matrix(NA,ncol=ncol(hunters_se)-2,nrow=nrow(hunter.bands))
for(i in 1:nrow(hunter.bands)){
  hunter.res<-bandtest(hunter.res,hunter.bands[i,])
  tmp<-hunters_se[hunters_se$timescales>=hunter.bands[i,1] & hunters_se$timescales<=hunter.bands[i,2],]
  hunter_se[i,]<-round(100*colMeans(tmp[,c(3:dim(tmp)[2])])/mean(tmp$sync),4)
  }
hunter.resP<-get_bandp(hunter.res)
colnames(hunter_se)<-colnames(tmp)[3:dim(tmp)[2]]
row.names(hunter_se)<-paste(hunter.bands[,1],hunter.bands[,2],sep="-")


##Run coherence of DVCs and abundance
#first filter data to match dimensions of DVC data
cty.list.dt<-lapply(cty.list[c('Abun','Crashes')],function(x){x<-cleandat(x[,!is.na(colSums(cty.list$Crashes))],clev = 5,times=1987:2016)$cdat;x})
dvc.res<-coh(dat1=cty.list.dt$Crashes,dat2=cty.list.dt$Abun,times=1987:maxyear,norm="powall",
             sigmethod="fast",nrand=nsurrogs,f0=1)
dvc.res<-bandtest(dvc.res,c(3,7))
dvc.resP<-get_bandp(dvc.res)
wlm_dvc<-wlm(list(cty.list.dt$Crashes,cty.list.dt$Abun),times=1987:maxyear,resp=1,pred=2,norm="powall",f0=1)
dvc_se<-syncexpl(wlm_dvc)
dvcabun_se<-dvc_se[dvc_se$timescales>=3 & dvc_se$timescales<=7,]
dvcabun_se37<-round(100*colMeans(dvcabun_se[,c(3:dim(dvcabun_se)[2])])/mean(dvcabun_se$sync),4)

#Run coherence of traffic-adjusted DVCs and abundance,again filtering to match dimensions of adjusted DVCs
cty.list.dt<-lapply(cty.list[c('Abun','AdjDVC')],function(x){x<-cleandat(x[,!is.na(colSums(cty.list$AdjDVC))],clev = 5,times=1988:2016)$cdat;x})
adjdvc.res<-coh(dat1=cty.list.dt$AdjDVC,dat2=cty.list.dt$Abun,times=1988:maxyear,norm="powall",
                sigmethod="fast",nrand=nsurrogs,f0=1)
adjdvc.res<-bandtest(adjdvc.res,c(3,7))
adjdvc.resP<-get_bandp(adjdvc.res)
wlm_adjdvc<-wlm(list(cty.list.dt$AdjDVC,cty.list.dt$Abun),times=1988:maxyear,resp=1,pred=2,norm="powall")
adjdvc_se<-syncexpl(wlm_adjdvc)
adjdvcabun_se<-adjdvc_se[adjdvc_se$timescales>=3 & adjdvc_se$timescales<=7,]
adjdvcabun_se37<-round(100*colMeans(adjdvcabun_se[,c(3:dim(adjdvcabun_se)[2])])/mean(adjdvcabun_se$sync),4)

#Run coherence of traffic and DVCs
cty.list.dt<-lapply(cty.list[c('Crashes','Traffic')],function(x){x<-cleandat(x[,!is.na(colSums(cty.list$Traffic))],clev = 5,times=1988:2016)$cdat;x})
traffic.dvc.res<-coh(dat1=cty.list.dt$Traffic,dat2=cty.list.dt$Crashes,times=1988:maxyear,norm="powall",
                sigmethod="fast",nrand=nsurrogs,f0=1)
traffic.dvc.res<-bandtest(traffic.dvc.res,c(3,7))
traffic.dvc.res<-bandtest(traffic.dvc.res,c(7,18))
traffic.dvc.resP<-get_bandp(traffic.dvc.res)
traffic.dvc.resP

#Build Table S1
tmp.list<-list()
for(i in names(clim.resP)){
  tmp.list[[i]]<-cbind(clim.resP[[i]],climabun_se[[i]])
}
TabS1<-do.call(rbind,tmp.list)
for(i in names(weath.resP)){
  tmp.list[[i]]<-cbind(weath.resP[[i]],weathabun_se[[i]])
}

##Need to add hunters to the table
TabS1<-rbind(do.call(rbind,tmp.list))
TabS1<-rbind(TabS1,cbind(do.call(cbind,hunter.resP),hunter_se))
row.names(TabS1)[(dim(TabS1)[1]-1):dim(TabS1)[1]]<-c("Hunter.Abun.1","Hunter.Abun.2")
TabS1<-rbind(TabS1,rbind(c(dvc.resP,dvcabun_se37),c(adjdvc.resP,adjdvcabun_se37)))
row.names(TabS1)[(dim(TabS1)[1]-1):dim(TabS1)[1]]<-c("Abun.DVC.1","Abun.AdjDVC.1")
TabS1.mat<-matrix(unlist(TabS1), ncol = length(TabS1), byrow = F)
TabS1.mat<-data.frame(matrix(unlist(strsplit(row.names(TabS1),'.',fixed=TRUE)),ncol=3,byrow=T)[,1:2],TabS1.mat)
colnames(TabS1.mat)<-c("Predictor","Response",colnames(TabS1))
TabS1.mat<-data.frame(TabS1.mat[,1:2],paste(TabS1.mat$ts_low_bd,TabS1.mat$ts_hi_bd,sep="-"),TabS1.mat[,c("p_val","mn_phs",'syncexpl',"crossterms","resids","pred1")])
colnames(TabS1.mat)<-c("Predictor","Response","Timescale","Pvalue","MeanPhase","SynchronyExplained","CrossTerms","Residuals","Pred1")
TabS1.mat<-TabS1.mat[,!colnames(TabS1.mat)%in%c("Residuals","Pred1")]
TableS1<-TabS1.mat[complete.cases(TabS1.mat),]
saveRDS(TableS1,file="Results/TableS1.rds")

# USDA Analysis -----------------------------------------------------------
if(scale.flag=="both"||scale.flag=="usda"){
library(wsyn)
bands<-rbind(c(3,7),c(3,4),c(4,7))
  
#Run coherence of climate indices and deer abundance
usda.clim.res<-list()
usda.wlm_climabun<-list()
usda.abun.dt<-cleandat(usda.list$Abun,clev=5,times=1981:2016)$cdat
climindex.usda<-climindex.usda[c("WinterNAO","WinterPDO","WinterMEI","SummerNAO","SummerPDO","SummerMEI")]
usda.climindex.dt<-lapply(climindex.usda,function(x){x<-cleandat(x,times=1981:2016,clev=5)$cdat;x})
usda.clim.res<-list()
wlm_climabun<-list()
for(j in names(usda.climindex.dt)){
  spatcoh.names<-paste(j,"Abun",sep=".")
  usda.clim.res[[spatcoh.names]]<-coh(dat1=usda.abun.dt,dat2=usda.climindex.dt[[j]],times=minyear:maxyear,norm="powall",
                                 sigmethod="fast",nrand=nsurrogs,f0=1)
  mod<-wlm(list(usda.abun.dt,usda.climindex.dt[[j]]),times=minyear:maxyear,resp=1,pred=2,norm="powall")
  usda.wlm_climabun[[spatcoh.names]]<-syncexpl(mod)
}

#Calculate abundance p-values and synchrony explained by climate indices
usda.climabun_se<-list()
for(j in names(usda.clim.res)){
  se_mat<-matrix(NA,ncol=ncol(usda.wlm_climabun[[1]])-2,nrow=nrow(bands))
  usda.clim.res[[j]]<-bandtest(usda.clim.res[[j]],bands[1,])
  tmp.se<-usda.wlm_climabun[[j]][usda.wlm_climabun[[j]]$timescales>=bands[1,1] & usda.wlm_climabun[[j]]$timescales<=bands[1,2],]
  se_mat[1,]<-round(100*colMeans(tmp.se[,c(3:dim(tmp.se)[2])])/mean(tmp.se$sync),4)
  if(usda.clim.res[[j]]$bandp['p_val']<0.06){
    for(i in 2:nrow(bands)){
      tmp.p<-bandtest(usda.clim.res[[j]],bands[i,])
      usda.clim.res[[j]]<-tmp.p
      tmp.se<-usda.wlm_climabun[[j]][usda.wlm_climabun[[j]]$timescales>=bands[i,1] & usda.wlm_climabun[[j]]$timescales<=bands[i,2],]
      se_mat[i,]<-round(100*colMeans(tmp.se[,c(3:dim(tmp)[2])])/mean(tmp.se$sync),4)
    }
  }
  usda.climabun_se[[j]]<-se_mat
}

usda.clim.resP<-lapply(usda.clim.res,get_bandp)
usda.climabun_se<-lapply(usda.climabun_se,function(x){colnames(x)<-colnames(tmp.se)[3:dim(tmp.se)[2]];x})
usda.climabun_se<-lapply(usda.climabun_se,function(x){row.names(x)<-paste(bands[,1],bands[,2],sep="-");x})

#Run coherence of winter weather and deer abundance
usda.weath.res<-list()
usda.wlm_weathabun<-list()
for(j in names(winter.clim.usda)){
  usda.winter.clim.dt<-cleandat(winter.clim.usda[[j]],clev=5,times=minyear:maxyear)$cdat
  spatcoh.names<-paste(j,"Abun",sep=".")
  usda.weath.res[[spatcoh.names]]<-coh(dat1=usda.abun.dt,dat2=usda.winter.clim.dt,times=minyear:maxyear,norm="powall",
                                  sigmethod="fast",nrand=nsurrogs,f0=1)
  mod<-wlm(list(usda.abun.dt,usda.winter.clim.dt),times=minyear:maxyear,resp=1,pred=2,norm="powall")
  usda.wlm_weathabun[[spatcoh.names]]<-syncexpl(mod)
}

#Calculate abundance synchrony explained by weather
usda.weathabun_se<-list()
for(j in names(usda.weath.res)){
  se_mat<-matrix(NA,ncol=ncol(usda.wlm_weathabun[[1]])-2,nrow=nrow(bands))
  usda.weath.res[[j]]<-bandtest(usda.weath.res[[j]],bands[1,])
  tmp.se<-usda.wlm_weathabun[[j]][usda.wlm_weathabun[[j]]$timescales>=bands[1,1] & usda.wlm_weathabun[[j]]$timescales<=bands[1,2],]
  se_mat[1,]<-round(100*colMeans(tmp.se[,c(3:dim(tmp.se)[2])])/mean(tmp.se$sync),4)
  if(usda.weath.res[[j]]$bandp['p_val']<0.06){
    for(i in 2:nrow(bands)){
      tmp.p<-bandtest(usda.weath.res[[j]],bands[i,])
      usda.weath.res[[j]]<-tmp.p
      tmp.se<-usda.wlm_weathabun[[j]][usda.wlm_weathabun[[j]]$timescales>=bands[i,1] & usda.wlm_weathabun[[j]]$timescales<=bands[i,2],]
      se_mat[i,]<-round(100*colMeans(tmp.se[,c(3:dim(tmp.se)[2])])/mean(tmp.se$sync),4)
    }
  }
  usda.weathabun_se[[j]]<-se_mat
}
usda.weath.res$Snwd.Abun<-bandtest(usda.weath.res$Snwd.Abun,c(3,4))
usda.weath.res$Snwd.Abun<-bandtest(usda.weath.res$Snwd.Abun,c(4,7))
tmp.se<-usda.wlm_weathabun$Snwd.Abun[usda.wlm_weathabun$Snwd.Abun$timescales>=bands[2,1] & usda.wlm_weathabun$Snwd.Abun$timescales<=bands[2,2],]
usda.weathabun_se$Snwd.Abun[2,]<-round(100*colMeans(tmp.se[,c(3:dim(tmp.se)[2])])/mean(tmp.se$sync),4)
tmp.se<-usda.wlm_weathabun$Snwd.Abun[usda.wlm_weathabun$Snwd.Abun$timescales>=bands[3,1] & usda.wlm_weathabun$Snwd.Abun$timescales<=bands[3,2],]
usda.weathabun_se$Snwd.Abun[3,]<-round(100*colMeans(tmp.se[,c(3:dim(tmp.se)[2])])/mean(tmp.se$sync),4)

usda.weath.resP<-lapply(usda.weath.res,get_bandp)
usda.weathabun_se<-lapply(usda.weathabun_se,function(x){colnames(x)<-colnames(tmp.se)[3:dim(tmp.se)[2]];x})
usda.weathabun_se<-lapply(usda.weathabun_se,function(x){row.names(x)<-paste(bands[,1],bands[,2],sep="-");x})

#Run coherence of climate indices and winter weather
usda.weath.climind.res<-list()
for(j in 1:length(names(climindex.usda))){
  for(i in 1:length(names(winter.clim.usda))){
    name1<-names(climindex.usda)[j]
    name2<-names(winter.clim.usda)[i]
    usda.climindex.tmp<-cleandat(climindex.usda[[name1]],times=minyear:maxyear,clev=5)$cdat
    usda.winter.clim.dt.tmp<-cleandat(winter.clim.usda[[name2]],times=minyear:maxyear,clev=5)$cdat
    spatcoh.names<-paste(name1,name2,sep=".")
    usda.weath.climind.res[[spatcoh.names]]<-coh(dat2=usda.climindex.tmp,dat1=usda.winter.clim.dt.tmp,times=minyear:maxyear,norm="powall",
                                            sigmethod="fast",nrand=nsurrogs,f0=1)
  }
}

for(j in names(usda.weath.climind.res)){
  usda.weath.climind.res[[j]]<-bandtest(usda.weath.climind.res[[j]],bands[1,])
}
TabS5<-do.call(rbind,lapply(weath.climind.res,get_bandp))
TabS5<-cbind(matrix(unlist(strsplit(row.names(TabS5),'.',fixed=TRUE)),ncol=2,byrow=T),TabS5)
TabS5$mn_phs<-as.numeric(as.character(ifelse(TabS5$p_val>0.06,"",TabS5$mn_phs)))
TableS5<-data.frame(TabS5[,1:2],paste(TabS5$ts_low_bd,TabS5$ts_hi_bd,sep="-"),TabS5[,c("p_val","mn_phs")])
colnames(TableS5)<-c("Predictor","Response","Timescale","Pvalue","MeanPhase")
saveRDS(TableS5,file="Results/TabS5_results.rds")

##Run coherence of hunters and abundance
#first filter data to match dimensions of hunter data
usda.list.dt<-lapply(usda.list[c('Abun','Hunters')],function(x){x<-x[,12:36];x})
usda.list.dt<-lapply(usda.list.dt,function(x){x[(!row.names(usda.list.dt$Abun)%in%cwd.usda) & !is.na(rowMeans(usda.list.dt$Hunters)),]})
usda.list.dt<-lapply(usda.list.dt,function(x){x<-cleandat(x,clev=5,times=1992:2016)$cdat;x})

#run coherence
hunter.bands<-rbind(c(3,7),c(2,2.5))
usda.hunter.res<-coh(dat1=usda.list.dt$Abun,dat2=usda.list.dt$Hunters,times=1992:maxyear,norm="powall",
                sigmethod="fast",nrand=nsurrogs,f0=1)
usda.wlm_hunters<-wlm(list(usda.list.dt$Abun,usda.list.dt$Hunters),times=1992:maxyear,resp=1,pred=2,norm="powall",f0=1)
usda.hunters_se<-syncexpl(usda.wlm_hunters)

usda.hunter_se<-matrix(NA,ncol=ncol(usda.hunters_se)-2,nrow=nrow(hunter.bands))
for(i in 1:nrow(hunter.bands)){
  usda.hunter.res<-bandtest(usda.hunter.res,hunter.bands[i,])
  tmp<-usda.hunters_se[usda.hunters_se$timescales>=hunter.bands[i,1] & usda.hunters_se$timescales<=hunter.bands[i,2],]
  usda.hunter_se[i,]<-round(100*colMeans(tmp[,c(3:dim(tmp)[2])])/mean(tmp$sync),4)
}
usda.hunter.resP<-get_bandp(usda.hunter.res)
colnames(usda.hunter_se)<-colnames(tmp)[3:dim(tmp)[2]]
row.names(usda.hunter_se)<-paste(hunter.bands[,1],hunter.bands[,2],sep="-")

#Run coherence of DVCs and traffic-adjusted DVCs with abundance,filtering to match temporal dimensions
usda.list.dt<-lapply(usda.list[c('Abun','Crashes')],function(x){x<-cleandat(x[,!is.na(colSums(usda.list$Crashes))],clev = 5,times=1987:2016)$cdat;x})
usda.dvc.res<-coh(dat1=usda.list.dt$Crashes,dat2=usda.list.dt$Abun,times=1987:maxyear,norm="powall",
                  sigmethod="fast",nrand=nsurrogs,f0=1)
usda.dvc.res<-bandtest(usda.dvc.res,c(3,7))
usda.dvc.resP<-get_bandp(usda.dvc.res)
usda.wlm_dvc<-wlm(list(usda.list.dt$Crashes,usda.list.dt$Abun),times=1987:maxyear,resp=1,pred=2,norm="powall")
usda.dvc_se<-syncexpl(usda.wlm_dvc)
usda.dvcabun_se<-usda.dvc_se[usda.dvc_se$timescales>=3 & usda.dvc_se$timescales<=7,]
usda.dvcabun_se37<-round(100*colMeans(usda.dvcabun_se[,c(3:dim(usda.dvcabun_se)[2])])/mean(usda.dvcabun_se$sync),4)

#Run coherence of traffic-adjusted DVCs and abundance,again filtering to match dimensions of adjusted DVCs
usda.list.dt<-lapply(usda.list[c('Abun','AdjDVC')],function(x){x<-cleandat(x[,!is.na(colSums(cty.list$AdjDVC))],clev = 5,times=1988:2016)$cdat;x})
usda.adjdvc.res<-coh(dat1=usda.list.dt$AdjDVC,dat2=usda.list.dt$Abun,times=1988:maxyear,norm="powall",
                sigmethod="fast",nrand=nsurrogs,f0=1)
usda.adjdvc.res<-bandtest(usda.adjdvc.res,c(3,7))
usda.adjdvc.resP<-get_bandp(usda.adjdvc.res)
usda.wlm_adjdvc<-wlm(list(usda.list.dt$AdjDVC,usda.list.dt$Abun),times=1988:maxyear,resp=1,pred=2,norm="powall")
usda.adjdvc_se<-syncexpl(usda.wlm_adjdvc)
usda.adjdvcabun_se<-usda.adjdvc_se[usda.adjdvc_se$timescales>=3 & usda.adjdvc_se$timescales<=7,]
usda.adjdvcabun_se37<-round(100*colMeans(usda.adjdvcabun_se[,c(3:dim(usda.adjdvcabun_se)[2])])/mean(usda.adjdvcabun_se$sync),4)

#Build Table S4
tmp.list<-list()
for(i in names(usda.clim.resP)){
  tmp.list[[i]]<-cbind(usda.clim.resP[[i]],usda.climabun_se[[i]])
}
TabS4<-do.call(rbind,tmp.list)
for(i in names(usda.weath.resP)){
  tmp.list[[i]]<-cbind(usda.weath.resP[[i]],usda.weathabun_se[[i]])
}

##Need to add hunters to the table
TabS4<-rbind(do.call(rbind,tmp.list))
TabS4<-rbind(TabS4,cbind(do.call(cbind,usda.hunter.resP),usda.hunter_se))
row.names(TabS4)[(dim(TabS4)[1]-1):dim(TabS4)[1]]<-c("Hunter.Abun.1","Hunter.Abun.2")
TabS4<-rbind(TabS4,rbind(c(usda.dvc.resP,usda.dvcabun_se37),c(usda.adjdvc.resP,usda.adjdvcabun_se37)))
row.names(TabS4)[(dim(TabS4)[1]-1):dim(TabS4)[1]]<-c("Abun.DVC.1","Abun.AdjDVC.1")
TabS4.mat<-matrix(unlist(TabS4), ncol = length(TabS4), byrow = F)
TabS4.mat<-data.frame(matrix(unlist(strsplit(row.names(TabS4),'.',fixed=TRUE)),ncol=3,byrow=T)[,1:2],TabS4.mat)
colnames(TabS4.mat)<-c("Predictor","Response",colnames(TabS4))
TabS4.mat<-data.frame(TabS4.mat[,1:2],paste(TabS4.mat$ts_low_bd,TabS4.mat$ts_hi_bd,sep="-"),TabS4.mat[,c("p_val","mn_phs",'syncexpl',"crossterms","resids","pred1")])
colnames(TabS4.mat)<-c("Predictor","Response","Timescale","Pvalue","MeanPhase","SynchronyExplained","CrossTerms","Residuals","Pred1")
TabS4.mat<-TabS4.mat[,!colnames(TabS4.mat)%in%c("Residuals","Pred1")]
TableS4<-TabS4.mat[complete.cases(TabS4.mat),]
saveRDS(TabS4.mat,file="Results/TableS4.rds")

}
