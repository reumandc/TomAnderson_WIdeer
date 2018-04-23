# Abundance Spatial Coherence Analysis ------------------------------------
climate.spcoh<-list()
climate.res<-list()
winter.spcoh<-list()
winter.res<-list()
hunter.spcoh<-list()
hunter.res<-list()
weath.climind.spcoh<-list()
weath.climind.res<-list()
indices.spcoh<-list()
indices.res<-list()
ranges=rbind(c(3,7),c(7,15))

#coherence of climate indices and abundance
abun.dt<-reumannplatz::CleanData(usda.list$Abun)$cleandat
climindex.dt<-lapply(climindex.usda[c(1:3,5:7)],function(x){x<-reumannplatz::CleanData(x)$cleandat;x})
for(j in names(climindex.dt)){
  spatcoh.names<-paste(j,"Abun",sep=".")
  climate.res[[spatcoh.names]]<-cohtestfast(dat2=climindex.dt[[j]],dat1=abun.dt,nsurrogs=nsurrogs,tsranges = ranges)
  climate.spcoh[[spatcoh.names]]<-swcoh(bio.dat=abun.dt,env.dat=climindex.dt[[j]],times=1981:2016)
}
lapply(climate.res,function(x){x$pvals})

#coherenc of winter weather and abundance
winter.clim.usda.dt<-lapply(winter.clim.usda,function(x){x<-reumannplatz::CleanData(x,normalize=T)$cleandat;x})
for(j in names(winter.clim.usda.dt)){
  spatcoh.names<-paste(j,"Abun",sep=".")
  winter.res[[spatcoh.names]]<-cohtestfast(dat2=winter.clim.usda.dt[[j]],dat1=abun.dt,nsurrogs=nsurrogs,tsranges=ranges)
  winter.spcoh[[spatcoh.names]]<-swcoh(bio.dat=abun.dt,env.dat=winter.clim.usda.dt[[j]],times=1981:2016)
}
lapply(winter.res,function(x){x$pvals})

#coherence of climate and weather
for(j in names(climindex.dt)){
  for(i in names(winter.clim.usda.dt)){
    spatcoh.names<-paste(j,i,sep=".")
    weath.climind.res[[spatcoh.names]]<-cohtestfast(dat2=climindex.dt[[j]],dat1=winter.clim.usda.dt[[i]],nsurrogs=nsurrogs,tsranges=ranges)
    weath.climind.spcoh[[spatcoh.names]]<-swcoh(bio.dat=winter.clim.usda.dt[[i]],env.dat=climindex.dt[[j]],times=1981:2016)
  }
}
lapply(weath.climind.res,function(x){x$pvals})

#coherence between climate indices
for(i in 1:(length(climindex.dt)-1)){
  for(j in (i+1):length(climindex.dt)){
    name1<-names(climindex.dt)[i]
    name2<-names(climindex.dt)[j]
    spatcoh.names<-paste(name2,name1,sep=".")
    indices.res[[spatcoh.names]]<-cohtestfast(dat1=climindex.dt[[j]],dat2=climindex.dt[[i]],nsurrogs=nsurrogs,tsranges=ranges)
    indices.spcoh[[spatcoh.names]]<-swcoh(bio.dat=climindex.dt[[j]],env.dat=climindex.dt[[i]],times=1981:2016)
  }
}
lapply(indices.res,function(x){x$pvals})

# Model Selection of Significant Coherence Pairs --------------------------

#store data in list
all.dat<-list(abun.dt,snow=winter.clim.usda.dt$Snwd,
              climindex.dt$WinterMEI,climindex.dt$WinterPDO)
abunmod<-wmrsig(indata=all.dat[c(1,2,3,4)],r=1,n=3,s=2,surr.test=F,n.surrog=1000)

# Sychrony explained --------------------------

#Synchrony explained in abundance by winter mei
win.mei.dt<-reumannplatz::CleanData(climindex.usda$WinterMEI,normalize=T)$cleandat
abun.dt<-reumannplatz::CleanData(usda.list$Abun,normalize=T)$cleandat
abunwinmei.es<-modelsyncexp(abun.dt,win.mei.dt,times=1981:2016,tsrange=c(3,7),plot=F) 
abunwinmei.es$avgsyncexp
abunwinmei.es$avgxterm

#synchrony explained in abundance by snow depth
snowd.dt<-reumannplatz::CleanData(winter.clim.usda$Snwd,normalize=T)$cleandat
abun.dt<-reumannplatz::CleanData(usda.list$Abun,normalize=T)$cleandat #truncate abundance by snow
abunsnwd.es<-modelsyncexp(abun.dt,snowd.dt,times=1981:2016,tsrange=c(3,5),plot=F)
abunsnwd.es$avgsyncexp
abunsnwd.es$avgxterm

#Synchrony explained in abundance by summer mei
sum.mei.dt<-reumannplatz::CleanData(climindex.usda$SummerMEI,normalize=T)$cleandat
abun.dt<-reumannplatz::CleanData(usda.list$Abun,normalize=T)$cleandat
abunsummei.es<-modelsyncexp(abun.dt,sum.mei.dt,times=1981:2016,tsrange=c(3,7),plot=F) 
abunsummei.es$avgsyncexp
abunsummei.es$avgxterm

#Synchrony explained in abundance by winter pdo
win.pdo.dt<-reumannplatz::CleanData(climindex.usda$WinterPDO,normalize=T)$cleandat
abun.dt<-reumannplatz::CleanData(usda.list$Abun,normalize=T)$cleandat
abunwinpdo.es<-modelsyncexp(abun.dt,win.pdo.dt,times=1981:2016,tsrange=c(3,7),plot=F) 
abunwinpdo.es$avgsyncexp
abunwinpdo.es$avgxterm

source("Functions/Fn_modelsyncexpwt.R")
source("Functions/Fn_syncexpplot.R")

abun.wt<-warray(all.dat[[1]],times=1981:2016)
dvc.wt<-warray(CleanData(usda.list$Crashes[,!is.na(colSums(usda.list$Crashes))])$cleandat,times=1987:2016)
model.es<-syncexpwt(abun.wt$wave.array,model = abunmod$pred.wt,times = abun.wt$times,timescales = abun.wt$timescales,tsrange = c(3,7),plot = F)
model.es$ave.syncexp

#coherence of hunters and abundance
abun<-usda.list$Abun[(!row.names(usda.list$Abun)%in%cwd.usda),12:36]
abun.dt<-reumannplatz::CleanData(abun,normalize=T)$cleandat
hunters<-usda.list$Hunters[(!row.names(usda.list$Hunters)%in%cwd.usda),12:36]
hunters.dt<-reumannplatz::CleanData(hunters,normalize=T)$cleandat
hunter.res<-cohtestfast(dat1=abun.dt,dat2=hunters.dt,nsurrogs=nsurrogs,tsranges=rbind(c(2,2.5)))
hunter.spcoh<-swcoh(bio.dat=abun.dt,env.dat=hunters.dt,times = 1981:2016)
hunter.res$pvals

#Synchrony explained in abundance by hunters
hunterabun.es<-modelsyncexp(abun.dt,hunters.dt,times=1992:2016,tsrange=c(2,2.5),plot=F) 
hunterabun.es$avgsyncexp
hunterabun.es$avgxterm

#coherence of DVCs and abundance
usda.list.dt<-lapply(usda.list[c(7,9)],function(x){x<-reumannplatz::CleanData(x[,!is.na(colSums(usda.list$Crashes))])$cleandat;x})
dvc.res<-cohtestfast(dat1=usda.list.dt$Crashes,dat2=usda.list.dt$Abun,nsurrogs=nsurrogs,tsranges=ranges)
dvc.spcoh<-swcoh(bio.dat=usda.list.dt$Crashes,env.dat=usda.list.dt$Abun,times = 1987:2016)
dvc.res$pvals

#Synchrony explained in DVCs by abundance
dvcabun.es<-modelsyncexp(usda.list.dt$Crashes,usda.list.dt$Abun,times=1987:2016,tsrange=c(3,7),plot=F) 
dvcabun.es$avgsyncexp
dvcabun.es$avgxterm

#coherence of adjusted DVCs and abundance
usda.list.dt<-lapply(usda.list[c(7,11)],function(x){x<-reumannplatz::CleanData(x[,!is.na(colSums(usda.list$AdjDVC))])$cleandat;x})
dvc.res1<-cohtestfast(dat1=usda.list.dt$AdjDVC,dat2=usda.list.dt$Abun,nsurrogs=nsurrogs,tsranges=ranges)
dvc.spcoh1<-swcoh(bio.dat=usda.list.dt$AdjDVC,env.dat=usda.list.dt$Abun,times = 1987:2016)
dvc.res1$pvals

#Synchrony explained in adjusted DVCs by abundance
dvcabun.es1<-modelsyncexp(usda.list.dt$AdjDVC,usda.list.dt$Abun,times=1988:2016,tsrange=c(3,7),plot=F) 
dvcabun.es1$avgsyncexp
dvcabun.es1$avgxterm

#calculate mean phases
source("Functions/Fn_phasemean.R")
phasemean(spatcoh = winter.spcoh$Snwd.Abun$empirical, timescales = winter.spcoh$Snwd.Abun$timescales,tsrange = c(3,5))/3.14
phasemean(spatcoh = climate.spcoh$WinterMEI.Abun$empirical, timescales = climate.spcoh$WinterMEI.Abun$timescales,tsrange=c(3,7))/3.14
phasemean(spatcoh = climate.spcoh$SummerMEI.Abun$empirical, timescales = climate.spcoh$SummerMEI.Abun$timescales,tsrange=c(3,7))/3.14
phasemean(spatcoh = climate.spcoh$WinterPDO.Abun$empirical, timescales = climate.spcoh$WinterPDO.Abun$timescales,tsrange=c(3,7))/3.14
phasemean(spatcoh = dvc.spcoh$empirical, timescales = dvc.spcoh$timescales,tsrange=c(3,7))/3.14
phasemean(spatcoh = dvc.spcoh1$empirical, timescales = dvc.spcoh1$timescales,tsrange=c(3,7))/3.14
phasemean(spatcoh = hunter.spcoh$empirical, timescales = hunter.spcoh$timescales,tsrange=c(2,2.5))/3.14
phasemean(spatcoh = indices.spcoh$WinterMEI.WinterPDO$empirical, timescales = indices.spcoh$WinterMEI.WinterPDO$timescales,tsrange=c(3,7))/3.14

