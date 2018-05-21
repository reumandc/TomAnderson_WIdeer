# Abundance Spatial Coherence Analysis ------------------------------------
climate.res<-list()
climate.spcoh<-list()
winter.spcoh<-list()
winter.res<-list()
weath.climind.res<-list()
weath.climind.spcoh<-list()
indices.res<-list()
indices.spcoh<-list()
nsurrogs=nsurrogs
ranges=rbind(c(3,7))

climindex<-climindex[c("WinterNAO","WinterPDO","WinterMEI","SummerNAO","SummerPDO","SummerMEI")]#drop ENSO from all subsequent analysis, redundant to MEI

#coherence of climate indices and abundance
abun.dt<-Reumannplatz::CleanData(cty.list$Abun,normalize=T)$cleandat
climindex.dt<-lapply(climindex[c("WinterNAO","WinterPDO","WinterMEI","SummerNAO","SummerPDO","SummerMEI")],function(x){x<-Reumannplatz::CleanData(x,normalize=T,each.ts = T)$cleandat;x})
for(j in names(climindex.dt)){
  spatcoh.names<-paste(j,"Abun",sep=".")
  climate.res[[spatcoh.names]]<-cohtestfast(dat2=climindex.dt[[j]],dat1=abun.dt,nsurrogs=nsurrogs,tsranges = ranges)
  climate.spcoh[[spatcoh.names]]<-swcoh(bio.dat=abun.dt,env.dat=climindex.dt[[j]],times=1981:2016)
}
climpvals<-lapply(climate.res,function(x){x$pvals})

#coherence of winter weather and abundance
for(j in names(winter.clim)){
  winter.clim.dt.tmp<-winter.clim[[j]][!is.na(rowMeans(winter.clim[[j]])),]
  winter.clim.dt.tmp<-Reumannplatz::CleanData(winter.clim.dt.tmp,normalize=T)$cleandat
  abun.dt.tmp<-abun.dt[!is.na(rowMeans(winter.clim[[j]])),]
  spatcoh.names<-paste(j,"Abun",sep=".")
  winter.res[[spatcoh.names]]<-cohtestfast(dat2=winter.clim.dt.tmp,dat1=abun.dt.tmp,nsurrogs=nsurrogs,tsranges=ranges)
  winter.spcoh[[spatcoh.names]]<-swcoh(bio.dat=abun.dt.tmp,env.dat=winter.clim.dt.tmp,times=1981:2016)
}
weathpvals<-lapply(winter.res,function(x){x$pvals})

#rerun coherence for snow depth and abundance but for 3-5 year timescales
snow.tmp<-winter.clim$Snwd[!is.na(rowMeans(winter.clim$Snwd)),]
snow.dt<-Reumannplatz::CleanData(snow.tmp,normalize=T)$cleandat
abun.dt<-abun.dt[!is.na(rowMeans(winter.clim$Snwd)),]
snow.abun.spatcoh35<-cohtestfast(dat2=snow.dt,dat1=abun.dt,nsurrogs=nsurrogs,tsranges=rbind(c(3,5)))
snow_abun_pval3_5<-snow.abun.spatcoh35$pvals
saveRDS(snow_abun_pval3_5,file = "Results/snowabunpval3_5.rds")

#coherence of climate and weather
TableS3<-matrix(NA,nrow=length(names(winter.clim)),ncol=length(names(climindex)))
for(j in 1:length(names(climindex))){
  for(i in 1:length(names(winter.clim))){
    name1<-names(climindex)[j]
    name2<-names(winter.clim)[i]
    climindex.tmp<-Reumannplatz::CleanData(climindex[[name1]],normalize=T)$cleandat
    winter.clim.dt<-winter.clim[[name2]][!is.na(rowMeans(winter.clim[[name2]])),]
    climindex.tmp<-climindex.tmp[!is.na(rowMeans(winter.clim[[name2]])),]
    winter.clim.dt.tmp<-Reumannplatz::CleanData(winter.clim.dt,normalize=T)$cleandat
    spatcoh.names<-paste(name1,name2,sep=".")
    weath.climind.res[[spatcoh.names]]<-cohtestfast(dat2=climindex.tmp,dat1=winter.clim.dt.tmp,nsurrogs=500,tsranges=ranges)
    weath.climind.spcoh[[spatcoh.names]]<-swcoh(bio.dat=winter.clim.dt.tmp,env.dat=climindex.tmp,times=1981:2016)
    TableS3[i,j]<-weath.climind.res[[spatcoh.names]]$pvals
  }
}
saveRDS(TableS3,file="Results/TabS3_results.rds")
#coherence between climate indices
TableS2<-matrix(NA,6,6)
for(i in 1:(length(climindex.dt)-1)){
  for(j in (i+1):length(climindex.dt)){
    name1<-names(climindex.dt)[i]
    name2<-names(climindex.dt)[j]
    spatcoh.names<-paste(name2,name1,sep=".")
    indices.res[[spatcoh.names]]<-cohtestfast(dat1=climindex.dt[[j]],dat2=climindex.dt[[i]],nsurrogs=nsurrogs,tsranges=ranges)
    indices.spcoh[[spatcoh.names]]<-swcoh(bio.dat=climindex.dt[[j]],env.dat=climindex.dt[[i]],times=1981:2016)
    TableS2[j,i]<-indices.res[[spatcoh.names]]$pvals
  }
}
diag(TableS2)<-1
row.names(TableS2)<-names(climindex.dt)
colnames(TableS2)<-names(climindex.dt)
saveRDS(TableS2,file="Results/TabS2_results.rds")

# Model Selection of Significant Coherence Pairs --------------------------
#consolidate data to winter temperature size (60 counties)
snwd<-winter.clim$Snwd[!is.na(rowMeans(winter.clim$Snwd)),]
abun<-cty.list$Abun[!is.na(rowMeans(winter.clim$Snwd)),]
climindex.tmp<-lapply(climindex,function(x){x[!is.na(rowMeans(winter.clim$Snwd)),]})

#store data in list
all.dat<-list(abun,snwd,climindex.tmp$WinterMEI,climindex.tmp$WinterPDO)

#normalize data
norm.dat<-lapply(all.dat,function(x){x<-Reumannplatz::CleanData(x,normalize=T)$cleandat;x})

#set up model
abunmod<-wmrsig(indata=norm.dat[c(1,2,3,4)],r=1,n=3,s=2,surr.test=F,n.surrog=1000)

# Sychrony explained --------------------------

#Synchrony explained in abundance by winter mei
win.mei.dt<-Reumannplatz::CleanData(climindex$WinterMEI,normalize=T)$cleandat
abun.dt<-Reumannplatz::CleanData(cty.list$Abun,normalize=T)$cleandat
abunwinmei.es<-modelsyncexp(abun.dt,win.mei.dt,times=1981:2016,tsrange=c(3,7),plot=F) 
abun_wmei_syncexp<-abunwinmei.es$avgsyncexp
abun_wmei_xterms<-abunwinmei.es$avgxterm

#synchrony explained in abundance by snow depth
snowd.dt<-Reumannplatz::CleanData(winter.clim$Snwd[!is.na(rowSums(winter.clim$Snwd)),],normalize=T)$cleandat
abun.dt<-Reumannplatz::CleanData(cty.list$Abun[!is.na(rowSums(winter.clim$Snwd)),],normalize=T)$cleandat #truncate abundance by snow
abunsnwd.es37<-modelsyncexp(abun.dt,snowd.dt,times=1981:2016,tsrange=c(3,7),plot=F)
abunsnwd.es35<-modelsyncexp(abun.dt,snowd.dt,times=1981:2016,tsrange=c(3,5),plot=F)

abun_snwd_syncexp3_7<-abunsnwd.es37$avgsyncexp
abun_snwd_xterms3_7<-abunsnwd.es37$avgxterm
abun_snwd_syncexp3_5<-abunsnwd.es35$avgsyncexp
abun_snwd_xterms3_5<-abunsnwd.es35$avgxterm
saveRDS(abun_snwd_syncexp3_5,file="Results/abunsnowSyncExp3_5.rds")

#Synchrony explained in abundance by summer mei
sum.mei.dt<-Reumannplatz::CleanData(climindex$SummerMEI,normalize=T)$cleandat
abun.dt<-Reumannplatz::CleanData(cty.list$Abun,normalize=T)$cleandat
abunsummei.es<-modelsyncexp(abun.dt,sum.mei.dt,times=1981:2016,tsrange=c(3,7),plot=F) 
abun_smei_syncexp<-abunsummei.es$avgsyncexp
abun_smei_xterms<-abunsummei.es$avgxterm

#Synchrony explained in abundance by winter pdo
win.pdo.dt<-Reumannplatz::CleanData(climindex$WinterPDO,normalize=T)$cleandat
abun.dt<-Reumannplatz::CleanData(cty.list$Abun,normalize=T)$cleandat
abunwinpdo.es<-modelsyncexp(abun.dt,win.pdo.dt,times=1981:2016,tsrange=c(3,7),plot=F) 
abun_wpdo_syncexp<-abunwinpdo.es$avgsyncexp
abun_wpdo_xterms<-abunwinpdo.es$avgxterm

#calculate synchrony explained by a model of snow depth, winter mei and winter pdo
source("Functions/Fn_modelsyncexpwt.R")
source("Functions/Fn_syncexpplot.R")
abun.wt<-warray(norm.dat[[1]],times=1981:2016)
dvc.wt<-warray(Reumannplatz::CleanData(cty.list$Crashes[,!is.na(colSums(cty.list$Crashes))])$cleandat,times=1987:2016)
model.es<-syncexpwt(abun.wt$wave.array,model = abunmod$pred.wt,times = abun.wt$times,timescales = abun.wt$timescales,tsrange = c(3,7),plot = F)
abun_model_syncexp<-model.es$ave.syncexp
saveRDS(abun_model_syncexp,file="Results/abunmodelSyncExp.rds")

#get phases from model
mod.coefs<-abunmod$coefs[abunmod$timescales<=7& abunmod$timescales>=3,1:3]

#Coherence of hunters and abundance
cty.list.dt<-lapply(cty.list[c(7,8)],function(x){x<-x[,12:36];x})
cty.list.dt<-lapply(cty.list.dt,function(x){x[(!row.names(cty.list.dt$Abun)%in%cwd) & !is.na(rowMeans(cty.list.dt$Hunters)),]})
cty.list.dt<-lapply(cty.list.dt,function(x){x<-Reumannplatz::CleanData(x)$cleandat;x})
hunter.res3_7<-cohtestfast(dat1=cty.list.dt$Abun,dat2=cty.list.dt$Hunters,nsurrogs=nsurrogs,tsranges=ranges)
hunter.spcoh<-swcoh(bio.dat=cty.list.dt$Abun,env.dat=cty.list.dt$Hunters,times = 1992:2016)
hunterpval3_7<-hunter.res3_7$pvals

hunter.res2_2.5<-cohtestfast(dat1=cty.list.dt$Abun,dat2=cty.list.dt$Hunters,nsurrogs=nsurrogs,tsranges=rbind(c(2,2.5)))
hunterpval2_2.5<-hunter.res2_2.5$pvals
saveRDS(hunterpval2_2.5,file="Results/hunterabun_pval2_2.5.rds")

#Synchrony explained in abundance by hunters
hunterabun.es2_2.5<-modelsyncexp(cty.list.dt$Abun,cty.list.dt$Hunters,times=1992:2016,tsrange=c(2,2.5),plot=F) 
hunterabun_syncexp2_2.5<-hunterabun.es2_2.5$avgsyncexp
hunterabun.es2_2.5$avgxterm
saveRDS(hunterabun_syncexp2_2.5,file="Results/hunterabunSyncExp2_2.5.rds")

#Coherence of DVCs and abundance
cty.list.dt<-lapply(cty.list[c(7,9)],function(x){x<-Reumannplatz::CleanData(x[,!is.na(colSums(cty.list$Crashes))])$cleandat;x})
dvc.res<-cohtestfast(dat1=cty.list.dt$Crashes,dat2=cty.list.dt$Abun,nsurrogs=nsurrogs,tsranges=ranges)
dvc.spcoh<-swcoh(bio.dat=cty.list.dt$Crashes,env.dat=cty.list.dt$Abun,times = 1981:2016)
dvcpval<-dvc.res$pvals

#Synchrony explained in DVCs by abundance
dvcabun.es<-modelsyncexp(cty.list.dt$Crashes,cty.list.dt$Abun,times=1987:2016,tsrange=c(3,7),plot=F) 
dvc_abun_syncexp<-dvcabun.es$avgsyncexp
dvc_abun_xterms<-dvcabun.es$avgxterm

#coherence of adjusted DVCs and abundance
cty.list.dt<-lapply(cty.list[c(7,11)],function(x){x<-Reumannplatz::CleanData(x[,!is.na(colSums(cty.list$AdjDVC))])$cleandat;x})
adjdvc.res<-cohtestfast(dat1=cty.list.dt$AdjDVC,dat2=cty.list.dt$Abun,nsurrogs=nsurrogs,tsranges=ranges)
adjdvc.spcoh<-swcoh(bio.dat=cty.list.dt$AdjDVC,env.dat=cty.list.dt$Abun,times = 1981:2016)
adjdvcpval<-adjdvc.res$pvals

#Synchrony explained in adjusted DVCs by abundance
adjdvcabun.es<-modelsyncexp(cty.list.dt$AdjDVC,cty.list.dt$Abun,times=1988:2016,tsrange=c(3,7),plot=F) 
adjdvc_abun_syncexp<-adjdvcabun.es$avgsyncexp
adjdvc_abun_xterms<-adjdvcabun.es$avgxterm

#Make blank data frame to store results
tableS1.names<-c("Response","Predictor","P-value","Mean Phase","Synchrony Explained", "Average Cross Terms")
TableS1<-data.frame(matrix(NA, 14, 6,
                           dimnames=list(c(), tableS1.names)),
                    stringsAsFactors=F)

#store all p-values in vectors and add to Table S1
resp<-c(rep("Abundance",(nrow(TableS1)-2)),"DVCs","Adjusted DVCs")
preds<-c(unlist(names(winter.clim)),"Hunters",unlist(names(climindex.dt)),rep("Abundance",2))
pvals<-c(unlist(weathpvals),hunterpval3_7,unlist(climpvals),dvcpval,adjdvcpval)

TableS1$Response<-resp
TableS1$Predictor<-preds
TableS1$P.value<-pvals

#Store synchrony explained for significant variables
TableS1[TableS1$Response=="Abundance" & TableS1$Predictor=="Snwd",'Synchrony.Explained']<-abun_snwd_syncexp3_7
TableS1[TableS1$Response=="Abundance" & TableS1$Predictor=="WinterMEI",'Synchrony.Explained']<-abun_wmei_syncexp
TableS1[TableS1$Response=="Abundance" & TableS1$Predictor=="WinterPDO",'Synchrony.Explained']<-abun_wpdo_syncexp
TableS1[TableS1$Response=="Abundance" & TableS1$Predictor=="SummerMEI",'Synchrony.Explained']<-abun_smei_syncexp
TableS1[TableS1$Response=="DVCs" & TableS1$Predictor=="Abundance",'Synchrony.Explained']<-dvc_abun_syncexp
TableS1[TableS1$Response=="Adjusted DVCs" & TableS1$Predictor=="Abundance",'Synchrony.Explained']<-adjdvc_abun_syncexp

#Store average cross terms for 
TableS1[TableS1$Response=="Abundance" & TableS1$Predictor=="Snwd",'Average.Cross.Terms']<-abun_snwd_xterms3_7
TableS1[TableS1$Response=="Abundance" & TableS1$Predictor=="WinterMEI",'Average.Cross.Terms']<-abun_wmei_xterms
TableS1[TableS1$Response=="Abundance" & TableS1$Predictor=="WinterPDO",'Average.Cross.Terms']<-abun_wpdo_xterms
TableS1[TableS1$Response=="Abundance" & TableS1$Predictor=="SummerMEI",'Average.Cross.Terms']<-abun_smei_xterms
TableS1[TableS1$Response=="DVCs" & TableS1$Predictor=="Abundance",'Average.Cross.Terms']<-dvc_abun_xterms
TableS1[TableS1$Response=="Adjusted DVCs" & TableS1$Predictor=="Abundance",'Average.Cross.Terms']<-adjdvc_abun_xterms

# Determine Mean Phase Relationships -------------------------------------------
#compute mean phase and store in S1
source("Functions/Fn_phasemean.R")
TableS1[TableS1$Response=="Abundance" & TableS1$Predictor=="Snwd",'Mean.Phase']<-phasemean(spatcoh = winter.spcoh$Snwd.Abun$empirical, timescales = winter.spcoh$Snwd.Abun$timescales,tsrange = c(3,7))/3.14
TableS1[TableS1$Response=="Abundance" & TableS1$Predictor=="WinterMEI",'Mean.Phase']<-phasemean(spatcoh = climate.spcoh$WinterMEI.Abun$empirical, timescales = climate.spcoh$WinterMEI.Abun$timescales,tsrange=c(3,7))/3.14
TableS1[TableS1$Response=="Abundance" & TableS1$Predictor=="SummerMEI",'Mean.Phase']<-phasemean(spatcoh = climate.spcoh$SummerMEI.Abun$empirical, timescales = climate.spcoh$SummerMEI.Abun$timescales,tsrange=c(3,7))/3.14
TableS1[TableS1$Response=="Abundance" & TableS1$Predictor=="WinterPDO",'Mean.Phase']<-phasemean(spatcoh = climate.spcoh$WinterPDO.Abun$empirical, timescales = climate.spcoh$WinterPDO.Abun$timescales,tsrange=c(3,7))/3.14
TableS1[TableS1$Response=="DVCs" & TableS1$Predictor=="Abundance",'Mean.Phase']<-phasemean(spatcoh = dvc.spcoh$empirical, timescales = dvc.spcoh$timescales,tsrange=c(3,7))/3.14
TableS1[TableS1$Response=="Adjusted DVCs" & TableS1$Predictor=="Abundance",'Mean.Phase']<-phasemean(spatcoh = adjdvc.spcoh$empirical, timescales = adjdvc.spcoh$timescales,tsrange=c(3,7))/3.14
snowabun_phase3_5<-phasemean(spatcoh = winter.spcoh$Snwd.Abun$empirical, timescales = winter.spcoh$Snwd.Abun$timescales,tsrange = c(3,5))/3.14
hunterphase2_2.5<-phasemean(spatcoh = hunter.spcoh$empirical, timescales = hunter.spcoh$timescales,tsrange=c(2,2.5))/3.14
wpdo_wmei_phase<-phasemean(spatcoh = indices.spcoh$WinterMEI.WinterPDO$empirical, timescales = indices.spcoh$WinterMEI.WinterPDO$timescales,tsrange=c(3,7))/3.14
smei_wpdo_phase<-phasemean(spatcoh = indices.spcoh$SummerMEI.WinterPDO$empirical, timescales = indices.spcoh$WinterMEI.WinterPDO$timescales,tsrange=c(3,7))/3.14
smei_wmei_phase<-phasemean(spatcoh = indices.spcoh$SummerMEI.WinterMEI$empirical, timescales = indices.spcoh$WinterMEI.WinterPDO$timescales,tsrange=c(3,7))/3.14

saveRDS(TableS1,file = "Results/TableS1.rds")
saveRDS(hunterphase2_2.5,file="Results/hunterabunphase2_2.5")

# Do Statewide Analysis ---------------------------------------------------
ann.abun.dt<-Reumannplatz::CleanData(colSums(cty.list$Abun),normalize=T)$cleandat
ann.snow.dt<-Reumannplatz::CleanData(colSums(winter.clim$Snwd[!(is.na(rowSums(winter.clim$Snwd))),]),normalize=F)$cleandat
ann.hunter.dt<-Reumannplatz::CleanData(wi.deer$GunHunters[wi.deer$Year>1980])$cleandat

win.mei.dt<-Reumannplatz::CleanData(climindex$WinterMEI[1,])$cleandat
win.pdo.dt<-Reumannplatz::CleanData(climindex$WinterPDO[1,])$cleandat
sum.mei.dt<-Reumannplatz::CleanData(climindex$SummerMEI[1,])$cleandat

ann.abun.hunter<-cohtestfast(dat1=ann.abun.dt,dat2=ann.hunter.dt,tsranges=rbind(c(3,7)),nsurrogs = nsurrogs)
ann_abun_hunter_pval<-ann.abun.hunter$pvals

ann.abun.res<-cohtestfast(dat1=ann.abun.dt,dat2=ann.snow.dt,tsranges=ranges,nsurrogs = nsurrogs)
ann_abun_snow_pval<-ann.abun.res$pvals

ann.abun.wmei<-cohtestfast(dat1=ann.abun.dt,dat2=win.mei.dt,tsranges=ranges,nsurrogs = nsurrogs)
ann.abun.smei<-cohtestfast(dat1=ann.abun.dt,dat2=sum.mei.dt,tsranges=ranges,nsurrogs = nsurrogs)
ann.abun.wpdo<-cohtestfast(dat1=ann.abun.dt,dat2=win.pdo.dt,tsranges=ranges,nsurrogs = nsurrogs)
ann_abun_wmei_pval<-ann.abun.wmei$pvals
ann_abun_smei_pval<-ann.abun.smei$pvals
ann_abun_wpdo_pval<-ann.abun.wpdo$pvals

ann.snwd.wmei<-cohtestfast(dat1=ann.snow.dt,dat2=win.mei.dt,tsranges=ranges,nsurrogs = nsurrogs)
ann.snwd.smei<-cohtestfast(dat1=ann.snow.dt,dat2=sum.mei.dt,tsranges=ranges,nsurrogs = nsurrogs)
ann.snwd.wpdo<-cohtestfast(dat1=ann.snow.dt,dat2=win.pdo.dt,tsranges=ranges,nsurrogs = nsurrogs)
ann_snwd_wmei_pval<-ann.snwd.wmei$pvals
ann_snwd_smei_pval<-ann.snwd.smei$pvals
ann_snwd_wpdo_pval<-ann.snwd.wpdo$pvals

ann.dvc.dt<-Reumannplatz::CleanData(colSums(cty.list$Crashes[,!is.na(colSums(cty.list$Crashes))]))$cleandat
ann.abun.dt1<-Reumannplatz::CleanData(colSums(cty.list$Abun[,!is.na(colSums(cty.list$Crashes))]))$cleandat
ann.dvc.res<-cohtestfast(dat1=ann.dvc.dt,dat2=ann.abun.dt1,tsranges = rbind(c(3,5),c(3,7)),nsurrogs = nsurrogs)
ann_dvc_abun_pval3_5<-ann.dvc.res$pvals[1]
ann_dvc_abun_pval3_7<-ann.dvc.res$pvals[2]

ann.adjdvc.dt<-Reumannplatz::CleanData(colSums(cty.list$AdjDVC[,!is.na(colSums(cty.list$AdjDVC))]))$cleandat
ann.abun.dt2<-Reumannplatz::CleanData(colSums(cty.list$Abun[,!is.na(colSums(cty.list$AdjDVC))]))$cleandat
ann.adjdvc.res<-cohtestfast(dat1=ann.adjdvc.dt,dat2=ann.abun.dt2,tsranges = rbind(c(3,7)),nsurrogs = nsurrogs)
ann_adjdvc_abun_pval<-ann.adjdvc.res$pvals

# Peak to trough distance ------------------------------------------------
ann.abun<-aggregate(Abun~Year,data=dat,FUN=sum)
ann.dvc<-aggregate(Crashes~Year,data=dat,FUN=sum)
dvc.resid<-residuals(lm(Crashes~poly(Year,3),data=ann.dvc))
abun.resid<-residuals(lm(Abun~Year,data=ann.abun))

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
avg.deer.fluctuations<-colMeans(diff)[1]
saveRDS(avg.deer.fluctuations[1],file="Results/totaldeerfluctuations.rds")

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
avg.dvc.fluctuations<-colMeans(diff1)
avg.dvc.costLow<-avg.dvc.fluctuations*2024
avg.dvc.costHigh<-avg.dvc.fluctuations*8388
saveRDS(avg.dvc.fluctuations[1],file="Results/totaldvcfluctuations.rds")
saveRDS(avg.dvc.costLow[1],file="Results/totaldvccostsLow.rds")
saveRDS(avg.dvc.costHigh[1],file="Results/totaldvccostsHigh.rds")
