# Abundance Spatial Coherence Analysis ------------------------------------
#create lists to store results
usda.climate.spcoh<-list()
usda.climate.res<-list()
usda.winter.spcoh<-list()
usda.winter.res<-list()
usda.weath.climind.spcoh<-list()
usda.weath.climind.res<-list()
usda.indices.spcoh<-list()
usda.indices.res<-list()

#stores surrogates and timescales used in analyses
nsurrogs=nsurrogs
ranges=c(3,7)
snowranges=c(3,4)
climranges=c(4,7)

#Run coherence of climate indices and abundance
usda.abun.dt<-Reumannplatz::CleanData(usda.list$Abun)$cleandat
usda.climindex.dt<-lapply(climindex.usda[c(1:3,5:7)],function(x){x<-Reumannplatz::CleanData(x)$cleandat;x})
for(j in names(climindex.dt)){
  spatcoh.names<-paste(j,"Abun",sep=".")
  usda.climate.res[[spatcoh.names]]<-cohtestfast(dat2=usda.climindex.dt[[j]],dat1=usda.abun.dt,nsurrogs=nsurrogs,tsranges = rbind(ranges,snowranges,climranges))
  usda.climate.spcoh[[spatcoh.names]]<-swcoh(bio.dat=usda.abun.dt,env.dat=usda.climindex.dt[[j]],times=1981:2016)
}
usda_climpvals37<-lapply(usda.climate.res,function(x){x$pvals[1]})
usda_climpvals34<-lapply(usda.climate.res,function(x){x$pvals[2]})
usda_climpvals47<-lapply(usda.climate.res,function(x){x$pvals[3]})

#Run coherence of winter weather and abundance
winter.clim.usda.dt<-lapply(winter.clim.usda,function(x){x<-Reumannplatz::CleanData(x,normalize=T)$cleandat;x})
for(j in names(winter.clim.usda.dt)){
  spatcoh.names<-paste(j,"Abun",sep=".")
  usda.winter.res[[spatcoh.names]]<-cohtestfast(dat2=winter.clim.usda.dt[[j]],dat1=usda.abun.dt,nsurrogs=nsurrogs,tsranges=rbind(ranges,snowranges,climranges))
  usda.winter.spcoh[[spatcoh.names]]<-swcoh(bio.dat=usda.abun.dt,env.dat=winter.clim.usda.dt[[j]],times=1981:2016)
}
usda_weathpvals37<-lapply(usda.winter.res,function(x){x$pvals[1]})
usda_weathpvals34<-lapply(usda.winter.res,function(x){x$pvals[2]})
usda_weathpvals47<-lapply(usda.winter.res,function(x){x$pvals[3]})

#Run coherence of climate and winter weather
TableS5<-matrix(NA,nrow=length(names(winter.clim.usda)),ncol=length(names(usda.climindex.dt)))
for(j in 1:length(names(usda.climindex.dt))){
  for(i in 1:length(names(winter.clim.usda.dt))){
    name1<-names(usda.climindex.dt)[j]
    name2<-names(winter.clim.usda.dt)[i]
    spatcoh.names<-paste(name1,name2,sep=".")
    usda.weath.climind.res[[spatcoh.names]]<-cohtestfast(dat2=usda.climindex.dt[[name1]],dat1=winter.clim.usda.dt[[name2]],nsurrogs=nsurrogs,tsranges=ranges)
    usda.weath.climind.spcoh[[spatcoh.names]]<-swcoh(bio.dat=winter.clim.usda.dt[[name2]],env.dat=usda.climindex.dt[[name1]],times=1981:2016)
    TableS5[i,j]<-usda.weath.climind.res[[spatcoh.names]]$pvals
  }
}
colnames(TableS5)<-names(usda.climindex.dt)
row.names(TableS5)<-names(winter.clim.usda.dt)

##Run wavelet multiple regression model of deer abundance with snow depth, winter MEI and winter PDO as predictors
#store data in list
all.dat<-list(usda.abun.dt,snow=winter.clim.usda.dt$Snwd,
              usda.climindex.dt$WinterMEI,usda.climindex.dt$WinterPDO)
usda.abunmod<-wmrsig(indata=all.dat[c(1,2,3,4)],r=1,n=3,s=2,surr.test=F,n.surrog=1000)

#Determine synchrony explained in abundance by winter mei
win.mei.dt<-Reumannplatz::CleanData(climindex.usda$WinterMEI,normalize=T)$cleandat
usda.abun.dt<-Reumannplatz::CleanData(usda.list$Abun,normalize=T)$cleandat
usda.abunwinmei.es37<-modelsyncexp(usda.abun.dt,win.mei.dt,times=1981:2016,tsrange=c(3,7),plot=F) 
usda_abunwmei_syncexp37<-usda.abunwinmei.es37$avgsyncexp
usda_abunwmei_xterms37<-usda.abunwinmei.es37$avgxterm
usda.abunwinmei.es47<-modelsyncexp(usda.abun.dt,win.mei.dt,times=1981:2016,tsrange=c(4,7),plot=F) 
usda_abunwmei_syncexp47<-usda.abunwinmei.es47$avgsyncexp
usda_abunwmei_xterms47<-usda.abunwinmei.es47$avgxterm

#Determine synchrony explained in abundance by snow depth over 3-4 year timescales
usda.snowd.dt<-Reumannplatz::CleanData(winter.clim.usda$Snwd,normalize=T)$cleandat
usda.abun.dt<-Reumannplatz::CleanData(usda.list$Abun,normalize=T)$cleandat #truncate abundance by snow
usda.abunsnwd.es<-modelsyncexp(usda.abun.dt,usda.snowd.dt,times=1981:2016,tsrange=c(3,4),plot=F)
usda_abunsnwd_syncexp<-usda.abunsnwd.es$avgsyncexp
usda_abunsnwd_xterms<-usda.abunsnwd.es$avgxterm

#Determine synchrony explained in abundance by summer mei
usda.sum.mei.dt<-Reumannplatz::CleanData(climindex.usda$SummerMEI,normalize=T)$cleandat
usda.abun.dt<-Reumannplatz::CleanData(usda.list$Abun,normalize=T)$cleandat
usda.abunsummei.es<-modelsyncexp(usda.abun.dt,usda.sum.mei.dt,times=1981:2016,tsrange=c(3,7),plot=F) 
usda_abunsmei_syncexp<-usda.abunsummei.es$avgsyncexp
usda_abunsmei_xterms<-usda.abunsummei.es$avgxterm

#Determine synchrony explained in abundance by winter pdo
usda.win.pdo.dt<-Reumannplatz::CleanData(climindex.usda$WinterPDO,normalize=T)$cleandat
usda.abun.dt<-Reumannplatz::CleanData(usda.list$Abun,normalize=T)$cleandat
usda.abunwinpdo.es37<-modelsyncexp(usda.abun.dt,usda.win.pdo.dt,times=1981:2016,tsrange=c(3,7),plot=F) 
usda_abunwpdo_syncexp37<-usda.abunwinpdo.es37$avgsyncexp
usda_abunwpdo_xterms37<-usda.abunwinpdo.es37$avgxterm
usda.abunwinpdo.es47<-modelsyncexp(usda.abun.dt,usda.win.pdo.dt,times=1981:2016,tsrange=c(3,7),plot=F) 
usda_abunwpdo_syncexp47<-usda.abunwinpdo.es47$avgsyncexp
usda_abunwpdo_xterms47<-usda.abunwinpdo.es47$avgxterm

#Determine synchrony explained by wavelet multiple regression model for abundance
source("Functions/Fn_modelsyncexpwt.R")
source("Functions/Fn_syncexpplot.R")
usda.abun.wt<-warray(all.dat[[1]],times=1981:2016)
usda.dvc.wt<-warray(CleanData(usda.list$Crashes[,!is.na(colSums(usda.list$Crashes))])$cleandat,times=1987:2016)
usda.model.es<-syncexpwt(usda.abun.wt$wave.array,model = usda.abunmod$pred.wt,times = abun.wt$times,timescales = abun.wt$timescales,tsrange = c(3,7),plot = F)
usda_model_syncexp<-usda.model.es$ave.syncexp

#Filter data to match spatial and temporal dimensions, and run coherence of hunters and abundance
usda.abun<-usda.list$Abun[(!row.names(usda.list$Abun)%in%cwd.usda),12:36]
usda.abun.dt<-Reumannplatz::CleanData(usda.abun,normalize=T)$cleandat
usda.hunters<-usda.list$Hunters[(!row.names(usda.list$Hunters)%in%cwd.usda),12:36]
usda.hunters.dt<-Reumannplatz::CleanData(usda.hunters,normalize=T)$cleandat
usda.hunter.res37<-cohtestfast(dat1=usda.abun.dt,dat2=usda.hunters.dt,nsurrogs=nsurrogs,tsranges=rbind(ranges,c(2,2.5)))
usda.hunter.spcoh<-swcoh(bio.dat=usda.abun.dt,env.dat=usda.hunters.dt,times = 1981:2016)
usda_hunter_pval37<-usda.hunter.res37$pvals[1]
usda_hunter_pval2_2.5<-usda.hunter.res37$pvals[2]

#Determine synchrony explained in abundance by hunters over 2-2.5 year timescales
usda.hunterabun.es<-modelsyncexp(usda.abun.dt,usda.hunters.dt,times=1992:2016,tsrange=c(2,2.5),plot=F) 
usda_hunterabun_syncexp2_2.5<-usda.hunterabun.es$avgsyncexp
usda_hunterabun_xterms2_2.5<-usda.hunterabun.es$avgxterm

#Run coherence between DVCs and abundance
usda.list.dt<-lapply(usda.list[c("Abun","Crashes")],function(x){x<-Reumannplatz::CleanData(x[,!is.na(colSums(usda.list$Crashes))])$cleandat;x})
usda.dvc.res<-cohtestfast(dat1=usda.list.dt$Crashes,dat2=usda.list.dt$Abun,nsurrogs=nsurrogs,tsranges=ranges)
usda.dvc.spcoh<-swcoh(bio.dat=usda.list.dt$Crashes,env.dat=usda.list.dt$Abun,times = 1987:2016)
usda_dvc_pval<-usda.dvc.res$pvals

#Determine synchrony explained in DVCs by abundance
usda.dvcabun.es<-modelsyncexp(usda.list.dt$Crashes,usda.list.dt$Abun,times=1987:2016,tsrange=c(3,7),plot=F) 
usda_dvc_syncexp<-usda.dvcabun.es$avgsyncexp
usda_dvc_xterms<-usda.dvcabun.es$avgxterm

#Run coherence of adjusted DVCs and abundance
usda.list.dt<-lapply(usda.list[c("Abun","AdjDVC")],function(x){x<-Reumannplatz::CleanData(x[,!is.na(colSums(usda.list$AdjDVC))])$cleandat;x})
usda.adjdvc.res<-cohtestfast(dat1=usda.list.dt$AdjDVC,dat2=usda.list.dt$Abun,nsurrogs=nsurrogs,tsranges=ranges)
usda.adjdvc.spcoh<-swcoh(bio.dat=usda.list.dt$AdjDVC,env.dat=usda.list.dt$Abun,times = 1987:2016)
usda_adjdvc_pval<-usda.adjdvc.res$pvals

#Determine synchrony explained in adjusted DVCs by abundance
usda.adjdvcabun.es<-modelsyncexp(usda.list.dt$AdjDVC,usda.list.dt$Abun,times=1988:2016,tsrange=c(3,7),plot=F) 
usda_adjdvc_syncexp<-usda.adjdvcabun.es$avgsyncexp
usda_adjdvc_xterms<-usda.adjdvcabun.es$avgxterm

#Create empty data frame for Table S4
tableS4.names<-c("Response","Predictor","Timescale","P-value","Mean Phase","Synchrony Explained", "Average Cross Terms")
TableS4<-data.frame(matrix(NA, 21, 7,
                           dimnames=list(c(), tableS4.names)),
                    stringsAsFactors=F)

#Store all p-values in vectors and add to Table S4
resp<-c(rep("Abundance",(nrow(TableS4)-2)),"DVCs","Adjusted DVCs")
preds<-c(unlist(names(winter.clim.usda)),"Hunters",unlist(names(usda.climindex.dt)),rep(c("Snwd","WinterMEI","WinterPDO"),2),"Hunters",rep("Abundance",2))
pvals<-c(unlist(usda_weathpvals),usda_hunter_pval37,unlist(usda_climpvals),
         usda_weathpvals34$Snwd.Abun,usda_climpvals34$WinterMEI.Abun,usda_climpvals34$WinterPDO.Abun,
         usda_weathpvals47$Snwd.Abun,usda_climpvals47$WinterMEI.Abun,usda_climpvals47$WinterPDO.Abun,usda_hunter_pval2_2.5,
         usda_dvc_pval,usda_adjdvc_pval)
ts<-c(rep("3-7",12),rep("3-4",3),rep("4-7",3),"2-2.5",rep("3-7",2))
TableS4$Response<-resp
TableS4$Predictor<-preds
TableS4$P.value<-pvals
TableS4$Timescale<-ts

TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="WinterMEI"& TableS4$Timescale=="3-7",'Synchrony.Explained']<-usda_abunwmei_syncexp37
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="Snwd" & TableS4$Timescale=="3-4",'Synchrony.Explained']<-usda_abunsnwd_syncexp
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="WinterMEI"& TableS4$Timescale=="4-7",'Synchrony.Explained']<-usda_abunwmei_syncexp47
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="SummerMEI",'Synchrony.Explained']<-usda_abunsmei_syncexp
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="WinterPDO"& TableS4$Timescale=="3-7",'Synchrony.Explained']<-usda_abunwpdo_syncexp37
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="WinterPDO"& TableS4$Timescale=="4-7",'Synchrony.Explained']<-usda_abunwpdo_syncexp47
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="Hunters"& TableS4$Timescale=="2-2.5",'Synchrony.Explained']<-usda_hunterabun_syncexp2_2.5
TableS4[TableS4$Response=="DVCs" & TableS4$Predictor=="Abundance",'Synchrony.Explained']<-usda_dvc_syncexp
TableS4[TableS4$Response=="Adjusted DVCs" & TableS4$Predictor=="Abundance",'Synchrony.Explained']<-usda_adjdvc_syncexp

TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="Snwd" & TableS4$Timescale=="3-4",'Average.Cross.Terms']<-usda_abunsnwd_xterms
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="WinterMEI"& TableS4$Timescale=="3-7",'Average.Cross.Terms']<-usda_abunwmei_xterms37
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="WinterMEI"& TableS4$Timescale=="4-7",'Average.Cross.Terms']<-usda_abunwmei_xterms47
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="SummerMEI",'Average.Cross.Terms']<-usda_abunsmei_xterms
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="WinterPDO"& TableS4$Timescale=="3-7",'Average.Cross.Terms']<-usda_abunwpdo_xterms37
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="WinterPDO"& TableS4$Timescale=="4-7",'Average.Cross.Terms']<-usda_abunwpdo_xterms47
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="Hunters"& TableS4$Timescale=="2-2.5",'Average.Cross.Terms']<-usda_hunterabun_xterms2_2.5
TableS4[TableS4$Response=="DVCs" & TableS4$Predictor=="Abundance",'Average.Cross.Terms']<-usda_dvc_xterms
TableS4[TableS4$Response=="Adjusted DVCs" & TableS4$Predictor=="Abundance",'Average.Cross.Terms']<-usda_adjdvc_xterms

##Calculate mean phases between significant pairs of variables and store in Table S4
#load helper function
source("Functions/Fn_phasemean.R")
usda_snwdabun_phase34<-phasemean(spatcoh = usda.winter.spcoh$Snwd.Abun$empirical, timescales = usda.winter.spcoh$Snwd.Abun$timescales,tsrange=c(3,4))/3.14
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="Snwd"& TableS4$Timescale=="3-4",'Mean.Phase']<-usda_snwdabun_phase34
usda_wmeiabun_phase37<-phasemean(spatcoh = usda.climate.spcoh$WinterMEI.Abun$empirical, timescales = usda.climate.spcoh$WinterMEI.Abun$timescales,tsrange=c(3,7))/3.14
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="WinterMEI"& TableS4$Timescale=="3-7",'Mean.Phase']<-usda_wmeiabun_phase37
usda_wmeiabun_phase47<-phasemean(spatcoh = usda.climate.spcoh$WinterMEI.Abun$empirical, timescales = usda.climate.spcoh$WinterMEI.Abun$timescales,tsrange=c(4,7))/3.14
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="WinterMEI"& TableS4$Timescale=="4-7",'Mean.Phase']<-usda_wmeiabun_phase37
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="SummerMEI",'Mean.Phase']<-phasemean(spatcoh = usda.climate.spcoh$SummerMEI.Abun$empirical, timescales = usda.climate.spcoh$SummerMEI.Abun$timescales,tsrange=c(3,7))/3.14
usda_wpdoabun_phase37<-phasemean(spatcoh = usda.climate.spcoh$WinterPDO.Abun$empirical, timescales = usda.climate.spcoh$WinterPDO.Abun$timescales,tsrange=c(3,7))/3.14
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="WinterPDO"& TableS4$Timescale=="3-7",'Mean.Phase']<-usda_wpdoabun_phase37
usda_wpdoabun_phase47<-phasemean(spatcoh = usda.climate.spcoh$WinterPDO.Abun$empirical, timescales = usda.climate.spcoh$WinterPDO.Abun$timescales,tsrange=c(4,7))/3.14
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="WinterPDO"& TableS4$Timescale=="4-7",'Mean.Phase']<-usda_wpdoabun_phase47
TableS4[TableS4$Response=="DVCs" & TableS4$Predictor=="Abundance",'Mean.Phase']<-phasemean(spatcoh = usda.dvc.spcoh$empirical, timescales = usda.dvc.spcoh$timescales,tsrange=c(3,7))/3.14
TableS4[TableS4$Response=="Adjusted DVCs" & TableS4$Predictor=="Abundance",'Mean.Phase']<-phasemean(spatcoh = usda.adjdvc.spcoh$empirical, timescales = usda.adjdvc.spcoh$timescales,tsrange=c(3,7))/3.14
usda_hunterphase2_2.5<-phasemean(spatcoh = usda.hunter.spcoh2_2.5$empirical, timescales = usda.hunter.spcoh2_2.5$timescales,tsrange=c(2,2.5))/3.14
TableS4[TableS4$Response=="Abundance" & TableS4$Predictor=="Hunters"& TableS4$Timescale=="2-2.5",'Mean.Phase']<-usda_hunterphase2_2.5

saveRDS(TableS4,file="Results/TabS4_results.rds")
