# Abundance Spatial Coherence Analysis ------------------------------------
climate.res<-list()
climate.spcoh<-list()
winter.spcoh<-list()
winter.res<-list()
hunter.spcoh<-list()
weath.climind.res<-list()
weath.climind.spcoh<-list()
indices.res<-list()
indices.spcoh<-list()
nsurrogs=nsurrogs
ranges=rbind(c(3,7),c(7,15))

#coherence of climate indices and abundance
abun.dt<-CleanData(cty.list$Abun,normalize=T)$cleandat
climindex.dt<-lapply(climindex[c(1:3,5:7)],function(x){x<-CleanData(x,normalize=T,each.ts = T)$cleandat;x})
for(j in names(climindex.dt)){
  spatcoh.names<-paste(j,"Abun",sep=".")
  climate.res[[spatcoh.names]]<-cohtestfast(dat2=climindex.dt[[j]],dat1=abun.dt,nsurrogs=nsurrogs,tsranges = ranges)
  climate.spcoh[[spatcoh.names]]<-swcoh(bio.dat=abun.dt,env.dat=climindex.dt[[j]],times=1981:2016)
}
lapply(climate.res,function(x){x$pvals})

#coherenc of winter weather and abundance
for(j in names(winter.clim)){
  winter.clim.dt.tmp<-winter.clim[[j]][!is.na(rowMeans(winter.clim[[j]])),]
  winter.clim.dt.tmp<-CleanData(winter.clim.dt.tmp,normalize=T)$cleandat
  abun.dt.tmp<-abun.dt[!is.na(rowMeans(winter.clim[[j]])),]
  spatcoh.names<-paste(j,"Abun",sep=".")
  winter.res[[spatcoh.names]]<-cohtestfast(dat2=winter.clim.dt.tmp,dat1=abun.dt.tmp,nsurrogs=nsurrogs,tsranges=ranges)
  winter.spcoh[[spatcoh.names]]<-swcoh(bio.dat=abun.dt.tmp,env.dat=winter.clim.dt.tmp,times=1981:2016)
}
lapply(winter.res,function(x){x$pvals})

#coherence of climate and weather
for(j in names(climindex)){
  for(i in names(winter.clim)){
    climindex.tmp<-CleanData(climindex[[j]],normalize=T)$cleandat
    winter.clim.dt<-winter.clim[[i]][!is.na(rowMeans(winter.clim[[i]])),]
    climindex.tmp<-climindex.tmp[!is.na(rowMeans(winter.clim[[i]])),]
    winter.clim.dt.tmp<-CleanData(winter.clim.dt,normalize=T)$cleandat
    spatcoh.names<-paste(j,i,sep=".")
    weath.climind.res[[spatcoh.names]]<-cohtestfast(dat2=climindex.tmp,dat1=winter.clim.dt.tmp,nsurrogs=nsurrogs,tsranges=ranges)
    weath.climind.spcoh[[spatcoh.names]]<-swcoh(bio.dat=winter.clim.dt.tmp,env.dat=climindex.tmp,times=1981:2016)
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
#consolidate data to winter temperature size (60 counties)
snwd<-winter.clim$Snwd[!is.na(rowMeans(winter.clim$Snwd)),]
abun<-cty.list$Abun[!is.na(rowMeans(winter.clim$Snwd)),]
climindex.tmp<-lapply(climindex,function(x){x[!is.na(rowMeans(winter.clim$Snwd)),]})

#store data in list
all.dat<-list(abun,snwd,climindex.tmp$WinterMEI,climindex.tmp$WinterPDO)

#normalize data
norm.dat<-lapply(all.dat,function(x){x<-CleanData(x,normalize=T)$cleandat;x})

#set up model
abunmod<-wmrsig(indata=norm.dat[c(1,2,3,4)],r=1,n=3,s=2,surr.test=F,n.surrog=1000)

# Sychrony explained --------------------------

#Synchrony explained in abundance by winter mei
win.mei.dt<-CleanData(climindex$WinterMEI,normalize=T)$cleandat
abun.dt<-CleanData(cty.list$Abun,normalize=T)$cleandat
abunwinmei.es<-modelsyncexp(abun.dt,win.mei.dt,times=1981:2016,tsrange=c(3,7),plot=F) 
abunwinmei.es$avgsyncexp
abunwinmei.es$avgxterm

#synchrony explained in abundance by snow depth
snowd.dt<-CleanData(winter.clim$Snwd[!is.na(rowSums(winter.clim$Snwd)),],normalize=T)$cleandat
abun.dt<-CleanData(cty.list$Abun[!is.na(rowSums(winter.clim$Snwd)),],normalize=T)$cleandat #truncate abundance by snow
abunsnwd.es<-modelsyncexp(abun.dt,snowd.dt,times=1981:2016,tsrange=c(3,5),plot=F)
abunsnwd.es$avgsyncexp
abunsnwd.es$avgxterm

#Synchrony explained in abundance by summer mei
sum.mei.dt<-CleanData(climindex$SummerMEI,normalize=T)$cleandat
abun.dt<-CleanData(cty.list$Abun,normalize=T)$cleandat
abunsummei.es<-modelsyncexp(abun.dt,sum.mei.dt,times=1981:2016,tsrange=c(3,7),plot=F) 
abunsummei.es$avgsyncexp
abunsummei.es$avgxterm

#Synchrony explained in abundance by winter pdo
win.pdo.dt<-CleanData(climindex$WinterPDO,normalize=T)$cleandat
abun.dt<-CleanData(cty.list$Abun,normalize=T)$cleandat
abunwinpdo.es<-modelsyncexp(abun.dt,win.pdo.dt,times=1981:2016,tsrange=c(3,7),plot=F) 
abunwinpdo.es$avgsyncexp
abunwinpdo.es$avgxterm

#calculate synchrony explained by a model of snow depth, winter mei and winter pdo
source("Functions/Fn_modelsyncexpwt.R")
source("Functions/Fn_syncexpplot.R")
abun.wt<-warray(norm.dat[[1]],times=1981:2016)
dvc.wt<-warray(CleanData(cty.list$Crashes[,!is.na(colSums(cty.list$Crashes))])$cleandat,times=1987:2016)
model.es<-syncexpwt(abun.wt$wave.array,model = abunmod$pred.wt,times = abun.wt$times,timescales = abun.wt$timescales,tsrange = c(3,7),plot = F)
model.es$ave.syncexp

#get phases from model
mod.coefs<-abunmod$coefs[abunmod$timescales<=7& abunmod$timescales>=3,1:3]

#Coherence of hunters and abundance
cty.list.dt<-lapply(cty.list[c(7,8)],function(x){x<-x[,12:36];x})
cty.list.dt<-lapply(cty.list.dt,function(x){x[(!row.names(cty.list.dt$Abun)%in%cwd) & !is.na(rowMeans(cty.list.dt$Hunters)),]})
cty.list.dt<-lapply(cty.list.dt,function(x){x<-CleanData(x)$cleandat;x})
hunter.res<-cohtestfast(dat1=cty.list.dt$Abun,dat2=cty.list.dt$Hunters,nsurrogs=nsurrogs,tsranges=rbind(c(3,7)))
hunter.spcoh<-swcoh(bio.dat=cty.list.dt$Abun,env.dat=cty.list.dt$Hunters,times = 1992:2016)
hunter.res$pvals

#Synchrony explained in abundance by hunters
hunterabun.es<-modelsyncexp(cty.list.dt$Abun,cty.list.dt$Hunters,times=1992:2016,tsrange=c(2,2.5),plot=F) 
hunterabun.es$avgsyncexp
hunterabun.es$avgxterm

#Coherence of DVCs and abundance
cty.list.dt<-lapply(cty.list[c(7,9)],function(x){x<-CleanData(x[,!is.na(colSums(cty.list$Crashes))])$cleandat;x})
dvc.res<-cohtestfast(dat1=cty.list.dt$Crashes,dat2=cty.list.dt$Abun,nsurrogs=nsurrogs,tsranges=ranges)
dvc.spcoh<-swcoh(bio.dat=cty.list.dt$Crashes,env.dat=cty.list.dt$Abun,times = 1981:2016)
dvc.res$pvals

#Synchrony explained in DVCs by abundance
dvcabun.es<-modelsyncexp(cty.list.dt$Crashes,cty.list.dt$Abun,times=1987:2016,tsrange=c(3,7),plot=F) 
dvcabun.es$avgsyncexp
dvcabun.es$avgxterm

#coherence of adjusted DVCs and abundance
cty.list.dt<-lapply(cty.list[c(7,11)],function(x){x<-CleanData(x[,!is.na(colSums(cty.list$AdjDVC))])$cleandat;x})
dvc.res1<-cohtestfast(dat1=cty.list.dt$AdjDVC,dat2=cty.list.dt$Abun,nsurrogs=nsurrogs,tsranges=ranges)
dvc.spcoh1<-swcoh(bio.dat=cty.list.dt$AdjDVC,env.dat=cty.list.dt$Abun,times = 1981:2016)
dvc.res1$pvals

#Synchrony explained in adjusted DVCs by abundance
dvcabun.es1<-modelsyncexp(cty.list.dt$AdjDVC,cty.list.dt$Abun,times=1988:2016,tsrange=c(3,7),plot=F) 
dvcabun.es1$avgsyncexp
dvcabun.es1$avgxterm

# Determine Mean Phase Relationships -------------------------------------------
source("Functions/Fn_phasemean.R")

phasemean(spatcoh = winter.spcoh$Snwd.Abun$empirical, timescales = winter.spcoh$Snwd.Abun$timescales,tsrange = c(3,7))/3.14
phasemean(spatcoh = climate.spcoh$WinterMEI.Abun$empirical, timescales = climate.spcoh$WinterMEI.Abun$timescales,tsrange=c(3,7))/3.14
phasemean(spatcoh = climate.spcoh$SummerMEI.Abun$empirical, timescales = climate.spcoh$SummerMEI.Abun$timescales,tsrange=c(3,7))/3.14
phasemean(spatcoh = climate.spcoh$WinterPDO.Abun$empirical, timescales = climate.spcoh$WinterPDO.Abun$timescales,tsrange=c(3,7))/3.14
phasemean(spatcoh = dvc.spcoh$empirical, timescales = dvc.spcoh$timescales,tsrange=c(3,7))/3.14
phasemean(spatcoh = dvc.spcoh1$empirical, timescales = dvc.spcoh1$timescales,tsrange=c(3,7))/3.14
phasemean(spatcoh = hunter.spcoh$empirical, timescales = hunter.spcoh$timescales,tsrange=c(2,2.5))/3.14
phasemean(spatcoh = indices.spcoh$WinterMEI.WinterPDO$empirical, timescales = indices.spcoh$WinterMEI.WinterPDO$timescales,tsrange=c(3,7))/3.14
phasemean(spatcoh = indices.spcoh$SummerMEI.WinterPDO$empirical, timescales = indices.spcoh$WinterMEI.WinterPDO$timescales,tsrange=c(3,7))/3.14
phasemean(spatcoh = indices.spcoh$SummerMEI.WinterMEI$empirical, timescales = indices.spcoh$WinterMEI.WinterPDO$timescales,tsrange=c(3,7))/3.14

# Do Statewide Analysis ---------------------------------------------------
ann.abun.dt<-CleanData(colSums(cty.list$Abun),normalize=T)$cleandat
ann.snow.dt<-CleanData(colSums(winter.clim$Snwd[!(is.na(rowSums(winter.clim$Snwd))),]),normalize=F)$cleandat
ann.hunter.dt<-CleanData(wi.deer$GunHunters[wi.deer$Year>1980])$cleandat

win.mei.dt<-CleanData(climindex$WinterMEI[1,])$cleandat
win.pdo.dt<-CleanData(climindex$WinterPDO[1,])$cleandat
sum.mei.dt<-CleanData(climindex$SummerMEI[1,])$cleandat

ann.abun.hunter<-cohtestfast(dat1=ann.abun.dt,dat2=ann.hunter.dt,tsranges=rbind(c(3,5)),nsurrogs = nsurrogs)
ann.abun.hunter$pvals

ann.abun.res<-cohtestfast(dat1=ann.abun.dt,dat2=ann.snow.dt,tsranges=ranges,nsurrogs = nsurrogs)
ann.abun.res$pvals

mod1<-cohtestfast(dat1=ann.abun.dt,dat2=win.mei.dt,tsranges=ranges,nsurrogs = nsurrogs)
mod2<-cohtestfast(dat1=ann.abun.dt,dat2=sum.mei.dt,tsranges=ranges,nsurrogs = nsurrogs)
mod3<-cohtestfast(dat1=ann.abun.dt,dat2=win.pdo.dt,tsranges=ranges,nsurrogs = nsurrogs)
mod1$pvals
mod2$pvals
mod3$pvals

mod4<-cohtestfast(dat1=ann.snow.dt,dat2=win.mei.dt,tsranges=ranges,nsurrogs = nsurrogs)
mod5<-cohtestfast(dat1=ann.snow.dt,dat2=sum.mei.dt,tsranges=ranges,nsurrogs = nsurrogs)
mod6<-cohtestfast(dat1=ann.snow.dt,dat2=win.pdo.dt,tsranges=ranges,nsurrogs = nsurrogs)
mod4$pvals
mod5$pvals
mod6$pvals

ann.dvc.dt<-CleanData(colSums(cty.list$Crashes[,!is.na(colSums(cty.list$Crashes))]))$cleandat
ann.abun.dt1<-CleanData(colSums(cty.list$Abun[,!is.na(colSums(cty.list$Crashes))]))$cleandat
ann.dvc.res<-cohtestfast(dat1=ann.dvc.dt,dat2=ann.abun.dt1,tsranges = rbind(c(3,5)),nsurrogs = nsurrogs)
ann.dvc.res$pvals

ann.adjdvc.dt<-CleanData(colSums(cty.list$AdjDVC[,!is.na(colSums(cty.list$AdjDVC))]))$cleandat
ann.abun.dt2<-CleanData(colSums(cty.list$Abun[,!is.na(colSums(cty.list$AdjDVC))]))$cleandat
ann.adjdvc.res<-cohtestfast(dat1=ann.adjdvc.dt,dat2=ann.abun.dt2,tsranges = rbind(c(3,7)),nsurrogs = nsurrogs)
ann.adjdvc.res$pvals

# Peak to trough distance ------------------------------------------------
ann.abun<-aggregate(Abun~Year,data=dat,FUN=sum)
ann.dvc<-aggregate(Crashes~Year,data=dat,FUN=sum)
par(mfrow=c(1,1))
plot(ann.dvc.dt~ann.dvc$Year,type="b",ylab="Box-Cox Normalized Abundance",xlab="Year")
lines(ann.abun$Year[ann.abun$Year>1986],ann.abun.dt1,col="red")
points(ann.abun$Year[ann.abun$Year>1986],ann.abun.dt1,col="red")
plot(Crashes~Year,data=ann.dvc,type="b")
dvc.resid<-residuals(lm(Crashes~poly(Year,3),data=ann.dvc))
abun.resid<-residuals(lm(Abun~Year,data=ann.abun))
plot(dvc.resid~ann.dvc$Year,type="b")

#values for deer abundance
peaks<-c(1989,1995,2000,2006,2012) #based on visual assessment
troughs<-c(1992,1997,2002,2009,2014) #based on visual assessment
diff<-matrix(NA,nrow=5,ncol=2)
for(i in 1:length(peaks)){
  for(i in 1:length(troughs)){
    diff[i,1]<-ann.abun$Abun[ann.abun$Year==peaks[i]]-ann.abun$Abun[ann.abun$Year==troughs[i]]
    diff[i,2]<-abun.resid[ann.abun$Year==peaks[i]]-abun.resid[ann.abun$Year==troughs[i]]
  }
}
colMeans(diff)

#values for dvcs
peaks1<-c(1990,1994,2003,2007,2012)
troughs1<-c(1991,1997,2005,2008,2014)
diff1<-matrix(NA,nrow=5,ncol=2)
for(i in 1:length(peaks1)){
  for(i in 1:length(troughs1)){
    diff1[i,1]<-ann.dvc$Crashes[ann.dvc$Year==peaks1[i]]-ann.dvc$Crashes[ann.dvc$Year==troughs1[i]]
    diff1[i,2]<-dvc.resid[ann.dvc$Year==peaks1[i]]-dvc.resid[ann.dvc$Year==troughs1[i]]
  }
}
colMeans(diff1)*2024
