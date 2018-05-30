# Abundance Spatial Coherence Analysis ------------------------------------
#create lists to store results
climate.spcoh<-list()
climate.res<-list()
winter.spcoh<-list()
winter.res<-list()
weath.climind.spcoh<-list()
weath.climind.res<-list()
indices.spcoh<-list()
indices.res<-list()

#stores surrogates and timescales used in analyses
nsurrogs=nsurrogs
ranges=rbind(c(3,7))

#Run coherence of climate indices and abundance
abun.dt<-Reumannplatz::CleanData(usda.list$Abun)$cleandat
climindex.dt<-lapply(climindex.usda[c(1:3,5:7)],function(x){x<-Reumannplatz::CleanData(x)$cleandat;x})
for(j in names(climindex.dt)){
  spatcoh.names<-paste(j,"Abun",sep=".")
  climate.res[[spatcoh.names]]<-cohtestfast(dat2=climindex.dt[[j]],dat1=abun.dt,nsurrogs=nsurrogs,tsranges = ranges)
  climate.spcoh[[spatcoh.names]]<-swcoh(bio.dat=abun.dt,env.dat=climindex.dt[[j]],times=1981:2016)
}
usda_climpvals<-lapply(climate.res,function(x){x$pvals})

#Run coherence of winter weather and abundance
winter.clim.usda.dt<-lapply(winter.clim.usda,function(x){x<-Reumannplatz::CleanData(x,normalize=T)$cleandat;x})
for(j in names(winter.clim.usda.dt)){
  spatcoh.names<-paste(j,"Abun",sep=".")
  winter.res[[spatcoh.names]]<-cohtestfast(dat2=winter.clim.usda.dt[[j]],dat1=abun.dt,nsurrogs=nsurrogs,tsranges=ranges)
  winter.spcoh[[spatcoh.names]]<-swcoh(bio.dat=abun.dt,env.dat=winter.clim.usda.dt[[j]],times=1981:2016)
}
usda_weathpvals<-lapply(winter.res,function(x){x$pvals})
snow.res<-cohtestfast(dat2=winter.clim.usda.dt$Snwd,dat1=abun.dt,nsurrogs=nsurrogs,tsranges=rbind(c(3,5)))
usda_snowabun_pval3_5<-snow.res$pvals

#Run coherence of climate and winter weather
TableS5<-matrix(NA,nrow=length(names(winter.clim.usda)),ncol=length(names(climindex.dt)))
for(j in 1:length(names(climindex.dt))){
  for(i in 1:length(names(winter.clim.usda.dt))){
    name1<-names(climindex.dt)[j]
    name2<-names(winter.clim.usda.dt)[i]
    spatcoh.names<-paste(name1,name2,sep=".")
    weath.climind.res[[spatcoh.names]]<-cohtestfast(dat2=climindex.dt[[name1]],dat1=winter.clim.usda.dt[[name2]],nsurrogs=nsurrogs,tsranges=ranges)
    weath.climind.spcoh[[spatcoh.names]]<-swcoh(bio.dat=winter.clim.usda.dt[[name2]],env.dat=climindex.dt[[name1]],times=1981:2016)
    TableS5[i,j]<-weath.climind.res[[spatcoh.names]]$pvals
  }
}
colnames(TableS5)<-names(climindex.dt)
row.names(TableS5)<-names(winter.clim.usda.dt)

##Run wavelet multiple regression model of deer abundance with snow depth, winter MEI and winter PDO as predictors
#store data in list
all.dat<-list(abun.dt,snow=winter.clim.usda.dt$Snwd,
              climindex.dt$WinterMEI,climindex.dt$WinterPDO)
usda.abunmod<-wmrsig(indata=all.dat[c(1,2,3,4)],r=1,n=3,s=2,surr.test=F,n.surrog=1000)

#Determine synchrony explained in abundance by winter mei
win.mei.dt<-Reumannplatz::CleanData(climindex.usda$WinterMEI,normalize=T)$cleandat
abun.dt<-Reumannplatz::CleanData(usda.list$Abun,normalize=T)$cleandat
usda.abunwinmei.es<-modelsyncexp(abun.dt,win.mei.dt,times=1981:2016,tsrange=c(3,7),plot=F) 
usda_abunwmei_syncexp<-usda.abunwinmei.es$avgsyncexp
usda_abunwmei_xterms<-usda.abunwinmei.es$avgxterm

#Determine synchrony explained in abundance by snow depth over 3-5 year timescales
snowd.dt<-Reumannplatz::CleanData(winter.clim.usda$Snwd,normalize=T)$cleandat
abun.dt<-Reumannplatz::CleanData(usda.list$Abun,normalize=T)$cleandat #truncate abundance by snow
usda.abunsnwd.es<-modelsyncexp(abun.dt,snowd.dt,times=1981:2016,tsrange=c(3,5),plot=F)
usda_abunsnwd_syncexp<-usda.abunsnwd.es$avgsyncexp
usda_abunsnwd_xterms<-usda.abunsnwd.es$avgxterm

#Determine synchrony explained in abundance by summer mei
sum.mei.dt<-Reumannplatz::CleanData(climindex.usda$SummerMEI,normalize=T)$cleandat
abun.dt<-Reumannplatz::CleanData(usda.list$Abun,normalize=T)$cleandat
usda.abunsummei.es<-modelsyncexp(abun.dt,sum.mei.dt,times=1981:2016,tsrange=c(3,7),plot=F) 
usda_abunsmei_syncexp<-usda.abunsummei.es$avgsyncexp
usda_abunsmei_xterms<-usda.abunsummei.es$avgxterm

#Determine synchrony explained in abundance by winter pdo
win.pdo.dt<-Reumannplatz::CleanData(climindex.usda$WinterPDO,normalize=T)$cleandat
abun.dt<-Reumannplatz::CleanData(usda.list$Abun,normalize=T)$cleandat
usda.abunwinpdo.es<-modelsyncexp(abun.dt,win.pdo.dt,times=1981:2016,tsrange=c(3,7),plot=F) 
usda_abunwpdo_syncexp<-usda.abunwinpdo.es$avgsyncexp
usda_abunwpdo_xterms<-usda.abunwinpdo.es$avgxterm

#Determine synchrony explained by wavelet multiple regression model for abundance
source("Functions/Fn_modelsyncexpwt.R")
source("Functions/Fn_syncexpplot.R")
abun.wt<-warray(all.dat[[1]],times=1981:2016)
dvc.wt<-warray(CleanData(usda.list$Crashes[,!is.na(colSums(usda.list$Crashes))])$cleandat,times=1987:2016)
usda.model.es<-syncexpwt(abun.wt$wave.array,model = usda.abunmod$pred.wt,times = abun.wt$times,timescales = abun.wt$timescales,tsrange = c(3,7),plot = F)
usda_model_syncexp<-usda.model.es$ave.syncexp

#Filter data to match spatial and temporal dimensions, and run coherence of hunters and abundance
abun<-usda.list$Abun[(!row.names(usda.list$Abun)%in%cwd.usda),12:36]
abun.dt<-Reumannplatz::CleanData(abun,normalize=T)$cleandat
hunters<-usda.list$Hunters[(!row.names(usda.list$Hunters)%in%cwd.usda),12:36]
hunters.dt<-Reumannplatz::CleanData(hunters,normalize=T)$cleandat
usda.hunter.res37<-cohtestfast(dat1=abun.dt,dat2=hunters.dt,nsurrogs=nsurrogs,tsranges=ranges)
usda.hunter.spcoh37<-swcoh(bio.dat=abun.dt,env.dat=hunters.dt,times = 1981:2016)
usda_hunter_pval37<-usda.hunter.res37$pvals

usda.hunter.res2_2.5<-cohtestfast(dat1=abun.dt,dat2=hunters.dt,nsurrogs=nsurrogs,tsranges=rbind(c(2,2.5)))
usda.hunter.spcoh2_2.5<-swcoh(bio.dat=abun.dt,env.dat=hunters.dt,times = 1981:2016)
usda_hunter_pval2_2.5<-usda.hunter.res2_2.5$pvals

#Determine synchrony explained in abundance by hunters over 2-2.5 year timescales
usda.hunterabun.es<-modelsyncexp(abun.dt,hunters.dt,times=1992:2016,tsrange=c(2,2.5),plot=F) 
usda_hunterabun_syncexp2_2.5<-usda.hunterabun.es$avgsyncexp
usda_hunterabun_xterms<-usda.hunterabun.es$avgxterm

#Run coherence between DVCs and abundance
usda.list.dt<-lapply(usda.list[c(7,9)],function(x){x<-Reumannplatz::CleanData(x[,!is.na(colSums(usda.list$Crashes))])$cleandat;x})
usda.dvc.res<-cohtestfast(dat1=usda.list.dt$Crashes,dat2=usda.list.dt$Abun,nsurrogs=nsurrogs,tsranges=ranges)
usda.dvc.spcoh<-swcoh(bio.dat=usda.list.dt$Crashes,env.dat=usda.list.dt$Abun,times = 1987:2016)
usda_dvc_pval<-usda.dvc.res$pvals

#Determine synchrony explained in DVCs by abundance
usda.dvcabun.es<-modelsyncexp(usda.list.dt$Crashes,usda.list.dt$Abun,times=1987:2016,tsrange=c(3,7),plot=F) 
usda_dvc_syncexp<-usda.dvcabun.es$avgsyncexp
usda_dvc_xterms<-usda.dvcabun.es$avgxterm

#Run coherence of adjusted DVCs and abundance
usda.list.dt<-lapply(usda.list[c(7,11)],function(x){x<-Reumannplatz::CleanData(x[,!is.na(colSums(usda.list$AdjDVC))])$cleandat;x})
usda.adjdvc.res<-cohtestfast(dat1=usda.list.dt$AdjDVC,dat2=usda.list.dt$Abun,nsurrogs=nsurrogs,tsranges=ranges)
usda.adjdvc.spcoh<-swcoh(bio.dat=usda.list.dt$AdjDVC,env.dat=usda.list.dt$Abun,times = 1987:2016)
usda_adjdvc_pval<-usda.adjdvc.res$pvals

#Determine synchrony explained in adjusted DVCs by abundance
usda.adjdvcabun.es<-modelsyncexp(usda.list.dt$AdjDVC,usda.list.dt$Abun,times=1988:2016,tsrange=c(3,7),plot=F) 
usda_adjdvc_syncexp<-usda.adjdvcabun.es$avgsyncexp
usda_adjdvc_xterms<-usda.adjdvcabun.es$avgxterm

#Create empty data frame for Table S4
tableS4.names<-c("Response","Predictor","P-value","Mean Phase","Synchrony Explained", "Average Cross Terms")
TableS4<-data.frame(matrix(NA, 14, 6,
                           dimnames=list(c(), tableS4.names)),
                    stringsAsFactors=F)

#Store all p-values in vectors and add to Table S4
resp<-c(rep("Abundance",(nrow(TableS4)-2)),"DVCs","Adjusted DVCs")
preds<-c(unlist(names(winter.clim.usda)),"Hunters",unlist(names(climindex.dt)),rep("Abundance",2))
pvals<-c(unlist(usda_weathpvals),usda_hunter_pval37,unlist(usda_climpvals),usda_dvc_pval,usda_adjdvc_pval)
TableS4$Response<-resp
TableS4$Predictor<-preds
TableS4$P.value<-pvals

TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="Snwd",'Synchrony.Explained']<-usda_abunsnwd_syncexp
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="WinterMEI",'Synchrony.Explained']<-usda_abunwmei_syncexp
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="SummerMEI",'Synchrony.Explained']<-usda_abunsmei_syncexp
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="WinterPDO",'Synchrony.Explained']<-usda_abunwpdo_syncexp
TableS4[TableS4$Response=="DVCs" & TableS4$Predictor=="Abundance",'Synchrony.Explained']<-usda_dvc_syncexp
TableS4[TableS4$Response=="Adjusted DVCs" & TableS4$Predictor=="Abundance",'Synchrony.Explained']<-usda_adjdvc_syncexp

TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="Snwd",'Average.Cross.Terms']<-usda_abunsnwd_xterms
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="WinterMEI",'Average.Cross.Terms']<-usda_abunwmei_xterms
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="SummerMEI",'Average.Cross.Terms']<-usda_abunsmei_xterms
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="WinterPDO",'Average.Cross.Terms']<-usda_abunwpdo_xterms
TableS4[TableS4$Response=="DVCs" & TableS4$Predictor=="Abundance",'Average.Cross.Terms']<-usda_dvc_xterms
TableS4[TableS4$Response=="Adjusted DVCs" & TableS4$Predictor=="Abundance",'Average.Cross.Terms']<-usda_adjdvc_xterms

##Calculate mean phases between significant pairs of variables and store in Table S4
#load helper function
source("Functions/Fn_phasemean.R")
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="Snwd",'Mean.Phase']<-phasemean(spatcoh = winter.spcoh$Snwd.Abun$empirical, timescales = winter.spcoh$Snwd.Abun$timescales,tsrange = c(3,5))/3.14
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="WinterMEI",'Mean.Phase']<-phasemean(spatcoh = climate.spcoh$WinterMEI.Abun$empirical, timescales = climate.spcoh$WinterMEI.Abun$timescales,tsrange=c(3,7))/3.14
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="SummerMEI",'Mean.Phase']<-phasemean(spatcoh = climate.spcoh$SummerMEI.Abun$empirical, timescales = climate.spcoh$SummerMEI.Abun$timescales,tsrange=c(3,7))/3.14
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="WinterPDO",'Mean.Phase']<-phasemean(spatcoh = climate.spcoh$WinterPDO.Abun$empirical, timescales = climate.spcoh$WinterPDO.Abun$timescales,tsrange=c(3,7))/3.14
TableS4[TableS4$Response=="DVCs" & TableS4$Predictor=="Abundance",'Mean.Phase']<-phasemean(spatcoh = usda.dvc.spcoh$empirical, timescales = usda.dvc.spcoh$timescales,tsrange=c(3,7))/3.14
TableS4[TableS4$Response=="Adjusted DVCs" & TableS4$Predictor=="Abundance",'Mean.Phase']<-phasemean(spatcoh = usda.adjdvc.spcoh$empirical, timescales = usda.adjdvc.spcoh$timescales,tsrange=c(3,7))/3.14
usda_hunterphase2_2.5<-phasemean(spatcoh = usda.hunter.spcoh2_2.5$empirical, timescales = usda.hunter.spcoh2_2.5$timescales,tsrange=c(2,2.5))/3.14
usda_snwdabunphase_3_5<-phasemean(spatcoh = winter.spcoh$Snwd.Abun$empirical, timescales = winter.spcoh$Snwd.Abun$timescales,tsrange=c(3,5))/3.14

saveRDS(TableS4,file="Results/TabS4_results.rds")
