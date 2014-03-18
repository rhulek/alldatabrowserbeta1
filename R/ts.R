ts<-function(records,centralValueType="median",whiskerValueType="5_95",transformationType="none") {
  
  casovani<-c()
  casovani<-c(casovani,"Pred nactenim knihovny",as.character(format(Sys.time(), "%H:%M:%OS3")))
  
  library(genasis)
  
  cenValue<-c()
  topValue<-c()
  botValue<-c()
  dateOfPoint<-c()
  nameOfSeries<-c()
  segment<-c()
  typeOfSeries<-c()
  globalUnit<-c()
  
  seriesDescription<-c()
  parameterDescription<-c()
  valueDescription<-c()
    
  ## Prvni opakovani cyklu - primarni casove rady
  # i cyklus bezi pres sites
  for (i in 1:length(records)) {
    loca<-as.character(records[[i]]$rowLabel)
    value         <-c()
    loqValue      <-c()
    loqMethodCode <-c()
    unit          <-c()
    dateTime      <-c()
    dateTimeString<-c()
    timeLength    <-c()
    
    casovani<-c(casovani,paste0("1. cyklus, ",i,". iterace"),as.character(format(Sys.time(), "%H:%M:%OS3")))
    
    for (j in 1:length(records[[i]]$values)) {
      value         <-c(value,         records[[i]]$values[[j]]$value)
      loqValue      <-c(loqValue,      records[[i]]$values[[j]]$loqValue)
      loqMethodCode <-c(loqMethodCode, records[[i]]$values[[j]]$loqMethodCode)
      unit          <-c(unit,          records[[i]]$values[[j]]$unit)
      dateTime      <-c(dateTime,      records[[i]]$values[[j]]$dateTime)
      dateTimeString<-c(dateTimeString,records[[i]]$values[[j]]$dateTimeString)
      timeLength    <-c(timeLength,    records[[i]]$values[[j]]$timeLength)
    }
    
    dateTime<-as.Date(dateTimeString)
    
    value         <-value[order(dateTimeString)]
    loqValue      <-loqValue[order(dateTimeString)]
    loqMethodCode <-loqMethodCode[order(dateTimeString)]
    unit          <-unit[order(dateTimeString)]
    dateTime      <-dateTime[order(dateTimeString)]
    timeLength    <-timeLength[order(dateTimeString)]
    dateTimeString<-dateTimeString[order(dateTimeString)]
    
    # Nahrada LoQ (v promenne valu budou hodnoty vstupujici do vypoctu)
    valu<-value
    valu[which(is.na(valu)&loqMethodCode=="INS")]<-loqValue[which(is.na(valu)&loqMethodCode=="INS")]*1/2
    
    # Logaritmicka transformace
    if (transformationType=="log") {
      valu<-log(valu)
    }
   
    if (length(dateTime)>1) {
      hole<-3*mean(dateTime[-1]-dateTime[-length(records[[i]]$values)],trim=0.05)
    } else {
      hole<-0
    }
    
    # k udava poradi segmentu jedne casove rady
    k<-1
    
    # j cyklus bezi pres jednotliva mereni
    for (j in 1:length(records[[i]]$values)) {
      if (j!=1) {
        if ((dateTime[j]-dateTime[j-1])>hole) {
          k<-k+1
        }
      }
      
      # Vystupni promenne
      cenValue    <-c(cenValue,valu[j])
      botValue    <-c(botValue,NA)
      topValue    <-c(topValue,NA)
      dateOfPoint <-c(dateOfPoint,dateTimeString[j])
      nameOfSeries<-c(nameOfSeries,loca)
      segment     <-c(segment,k)
      typeOfSeries<-c(typeOfSeries,"prim")
      globalUnit  <-c(globalUnit,as.character(unique(unit)))
    }
    
    # Popis primarnich casovych rad v 1. cyklu    
    if (k==1) {
      res<-genstatistic(valu,dateTime)$res
      
      parameterNames<-c("delta",
                        "mannKendall",
                        "mannKendallP",
                        "daniels",
                        "danielsP",
                        "mean",
                        "sd",
                        "geomean",
                        "gsd",
                        "median",
                        "min",
                        "max",
                        "perc5",
                        "perc25",
                        "perc75",
                        "perc95",
                        "geoMean95CIUpperBound",
                        "geoMean95CILowerBound",
                        "mean95CIUpperBound",
                        "mean95CILowerBound")
      
      parameterValues<-c(res$delta,
                         res$"Mann-Kendall",
                         res$MKp,
                         res$Daniels,
                         res$Dp,
                         res$mean,
                         res$sd,
                         res$"geom. mean",
                         res$"geom. sd",
                         res$median,
                         res$min,
                         res$max,
                         quantile05(valu),
                         quantile25(valu),
                         quantile75(valu),
                         quantile95(valu),
                         res$"geom. mean"*res$"geom. sd"^qnorm(0.975),
                         res$"geom. mean"*res$"geom. sd"^qnorm(0.025),
                         res$mean+qnorm(0.975)*res$sd,
                         res$mean+qnorm(0.025)*res$sd)
      
      seriesDescription<-c(seriesDescription,rep(loca,length(parameterNames)))
      parameterDescription<-c(parameterDescription,parameterNames)
      valueDescription<-c(valueDescription,parameterValues)
    }
  }
  
    
  ## Druhe opakovani cyklu - vypocet rad agregaci
  
  # Vyber agregacnich funkci
  whisk<-c("5_95","25_75","min_max","2iq","ci")
  whisl<-c("quantile05","quantile25","min","iql","cil")
  whisu<-c("quantile95","quantile75","max","iqu","ciu")
  
  centv<-c("mean","median","geomean")
  centf<-c("arimean","quantile50","geomean")
  
  f1<-centf[which(centv==centralValueType)]
  f2<-whisl[which(whisk==whiskerValueType)]
  f3<-whisu[which(whisk==whiskerValueType)]
  
  # i cyklus bezi pres sites
  for (i in 1:length(records)) {
    loca<-as.character(records[[i]]$rowLabel)
    value         <-c()
    loqValue      <-c()
    loqMethodCode <-c()
    unit          <-c()
    dateTime      <-c()
    dateTimeString<-c()
    timeLength    <-c()
    
    casovani<-c(casovani,paste0("2. cyklus, ",i,". iterace"),as.character(format(Sys.time(), "%H:%M:%OS3")))
    
    
    for (j in 1:length(records[[i]]$values)) {
      value         <-c(value,         records[[i]]$values[[j]]$value)
      loqValue      <-c(loqValue,      records[[i]]$values[[j]]$loqValue)
      loqMethodCode <-c(loqMethodCode, records[[i]]$values[[j]]$loqMethodCode)
      unit          <-c(unit,          records[[i]]$values[[j]]$unit)
      dateTime      <-c(dateTime,      records[[i]]$values[[j]]$dateTime)
      dateTimeString<-c(dateTimeString,records[[i]]$values[[j]]$dateTimeString)
      timeLength    <-c(timeLength,    records[[i]]$values[[j]]$timeLength)
    }
    
    dateTime<-as.Date(dateTimeString)
    
    value         <-value[order(dateTimeString)]
    loqValue      <-loqValue[order(dateTimeString)]
    loqMethodCode <-loqMethodCode[order(dateTimeString)]
    unit          <-unit[order(dateTimeString)]
    dateTime      <-dateTime[order(dateTimeString)]
    timeLength    <-timeLength[order(dateTimeString)]
    dateTimeString<-dateTimeString[order(dateTimeString)]
    
    # Nahrada LoQ (v promenne valu budou hodnoty vstupujici do vypoctu)
    valu<-value
    valu[which(is.na(valu)&loqMethodCode=="INS")]<-loqValue[which(is.na(valu)&loqMethodCode=="INS")]*1/2
    
    # Logaritmicka transformace
    if (transformationType=="log") {
      valu<-log(valu)
    }
    
    # Jednotky musi byt stejne
    if (length(unique(unit))>1) {
      stop("Units of some records differ!")
    } else {
      unit<-unique(unit)
    }
    
    # Rocni agregace do noveho data.frame aggr
    year<-gendate(substr(dateTimeString,1,4))
    aggr<-data.frame(aggregate(valu,by=list(year),FUN=f1)[,2],
                     as.character(unit),
                     centralValueType,
                     aggregate(valu,by=list(year),FUN=f3)[,2],
                     aggregate(valu,by=list(year),FUN=f2)[,2],
                     whiskerValueType,
                     gendate(aggregate(dateTime,by=list(year),FUN=mean)[,1]),
                     as.character(aggregate(dateTime,by=list(year),FUN=mean)[,1]),
                     as.character(aggregate(valu,by=list(year),FUN=length)[,2]),
                     as.character(aggregate(value,by=list(year),FUN=loqlength)[,2]))
    colnames(aggr)<-c("centralValue","unit","centralValueType","whiskerTopValue","whiskerBottomValue","whiskerType","dateTime","dateTimeString","n","nUnderLOQ")
    
    return(aggr$dateTimeString[j])
    
    if (nrow(aggr)>1) {
      hole<-3*mean(aggr$dateTime[-1]-aggr$dateTime[-nrow(aggr)],trim=0.05)
    } else {
      hole<-0
    }
    
    # k udava poradi segmentu jedne casove rady
    k<-1
    
    # j cyklus bezi pres jednotliva mereni
    for (j in 1:nrow(aggr)) {
      if (j!=1) {
        if ((aggr$dateTime[j]-aggr$dateTime[j-1])>hole) {
          k<-k+1
        }
      }
      # Vystupni promenne
      cenValue    <-c(cenValue,aggr$centralValue[j])
      botValue    <-c(botValue,aggr$whiskerBottomValue[j])
      topValue    <-c(topValue,aggr$whiskerTopValue[j])
      dateOfPoint <-c(dateOfPoint,aggr$dateTimeString[j])
      nameOfSeries<-c(nameOfSeries,loca)
      segment     <-c(segment,k)
      typeOfSeries<-c(typeOfSeries,"aggr")
      globalUnit  <-c(globalUnit,as.character(unique(aggr$unit)))
    }
    
    # Popis agregovanych casovych rad ve 2. cyklu
    if (k==1) {
      res<-genstatistic(aggr$centralValue,aggr$dateTime)$res
      
      parameterNames<-c("delta",
                        "mannKendall",
                        "mannKendallP",
                        "daniels",
                        "danielsP",
                        "mean",
                        "sd",
                        "geomean",
                        "gsd",
                        "median",
                        "min",
                        "max",
                        "perc5",
                        "perc25",
                        "perc75",
                        "perc95",
                        "geoMean95CIUpperBound",
                        "geoMean95CILowerBound",
                        "mean95CIUpperBound",
                        "mean95CILowerBound")
      
      parameterValues<-c(res$delta,
                         res$"Mann-Kendall",
                         res$MKp,
                         res$Daniels,
                         res$Dp,
                         res$mean,
                         res$sd,
                         res$"geom. mean",
                         res$"geom. sd",
                         res$median,
                         res$min,
                         res$max,
                         quantile05(aggr$centralValue),
                         quantile25(aggr$centralValue),
                         quantile75(aggr$centralValue),
                         quantile95(aggr$centralValue),
                         res$"geom. mean"*res$"geom. sd"^qnorm(0.975),
                         res$"geom. mean"*res$"geom. sd"^qnorm(0.025),
                         res$mean+qnorm(0.975)*res$sd,
                         res$mean+qnorm(0.025)*res$sd)
      
      seriesDescription<-c(seriesDescription,rep(loca,length(parameterNames)))
      parameterDescription<-c(parameterDescription,parameterNames)
      valueDescription<-c(valueDescription,parameterValues)
    }
  }


  ## Treti opakovani cyklu - vypocet trendu primarnich rad
  
  # i cyklus bezi pres sites
  for (i in 1:length(records)) {
    loca<-as.character(records[[i]]$rowLabel)
    value         <-c()
    loqValue      <-c()
    loqMethodCode <-c()
    unit          <-c()
    dateTime      <-c()
    dateTimeString<-c()
    timeLength    <-c()
    
    casovani<-c(casovani,paste0("3. cyklus, ",i,". iterace"),as.character(format(Sys.time(), "%H:%M:%OS3")))
    
    
    for (j in 1:length(records[[i]]$values)) {
      value         <-c(value,         records[[i]]$values[[j]]$value)
      loqValue      <-c(loqValue,      records[[i]]$values[[j]]$loqValue)
      loqMethodCode <-c(loqMethodCode, records[[i]]$values[[j]]$loqMethodCode)
      unit          <-c(unit,          records[[i]]$values[[j]]$unit)
      dateTime      <-c(dateTime,      records[[i]]$values[[j]]$dateTime)
      dateTimeString<-c(dateTimeString,records[[i]]$values[[j]]$dateTimeString)
      timeLength    <-c(timeLength,    records[[i]]$values[[j]]$timeLength)
    }
    
    dateTime<-as.Date(dateTimeString)
    
    value         <-value[order(dateTimeString)]
    loqValue      <-loqValue[order(dateTimeString)]
    loqMethodCode <-loqMethodCode[order(dateTimeString)]
    unit          <-unit[order(dateTimeString)]
    dateTime      <-dateTime[order(dateTimeString)]
    timeLength    <-timeLength[order(dateTimeString)]
    dateTimeString<-dateTimeString[order(dateTimeString)]
    
    # Nahrada LoQ (v promenne valu budou hodnoty vstupujici do vypoctu)
    valu<-value
    valu[which(is.na(valu)&loqMethodCode=="INS")]<-loqValue[which(is.na(valu)&loqMethodCode=="INS")]*1/2
    
    if (length(dateTime)>1) {
      hole<-3*mean(dateTime[-1]-dateTime[-length(records[[i]]$values)],trim=0.05)
    } else {
      hole<-0
    }
    
    if (max(c(hole,dateTime[-1]-dateTime[-length(records[[i]]$values)]))>hole|(length(dateTime)<3)) { # Hole added to the vector to avoid problems with vector of length 1.
      series<-NA
    } else {
      curve<-data.frame(as.Date(as.numeric(genplot(valu,dateTime,n=20,distr="lnorm",plot=FALSE)$belt[1,]),origin="1970-01-01"),
                        as.numeric(genplot(valu,dateTime,n=20,distr="lnorm",plot=FALSE)$line[1,]),
                        as.numeric(genplot(valu,dateTime,n=20,distr="lnorm",plot=FALSE)$lower[1,]),
                        as.numeric(genplot(valu,dateTime,n=20,distr="lnorm",plot=FALSE)$upper[1,]))
      colnames(curve)<-c("belt","line","lower","upper")
      
      # Logaritmace v pripade log transformace (trend bude linearni)
      if (transformationType=="log") {
        curve$line<-log(curve$line)
        curve$lower<-log(curve$lower)
        curve$upper<-log(curve$upper)
      }
      
      # j cyklus bezi pres jednotlive body krivky
      for (j in 1:nrow(curve)) {
        
        # Vystupni promenne
        cenValue    <-c(cenValue,curve$line[j])
        botValue    <-c(botValue,curve$lower[j])
        topValue    <-c(topValue,curve$upper[j])
        dateOfPoint <-c(dateOfPoint,as.character(curve$belt[j]))
        nameOfSeries<-c(nameOfSeries,loca)
        segment     <-c(segment,1)
        typeOfSeries<-c(typeOfSeries,"prim_trend") 
        globalUnit  <-c(globalUnit,as.character(unique(unit)))
      }
    }
    
    # Popis trendovych krivek v 3. cyklu.    
    if (max(c(hole,dateTime[-1]-dateTime[-length(records[[i]]$values)]))<=hole) {
      parameterNames<-c("slope",
                        "intercept")
      
      parameterValues<-c(genplot(valu,dateTime,n=20,distr="lnorm",plot=FALSE)$slope,
                         genplot(valu,dateTime,n=20,distr="lnorm",plot=FALSE)$intercept)
      
      seriesDescription<-c(seriesDescription,rep(loca,length(parameterNames)))
      parameterDescription<-c(parameterDescription,parameterNames)
      valueDescription<-c(valueDescription,parameterValues)
    }
  }
  
  ## Ctvrte opakovani cyklu - vypocet trendu agregovanych rad
  
  # Vyber agregacnich funkci
  whisk<-c("5_95","25_75","min_max","2iq","ci")
  whisl<-c("quantile05","quantile25","min","iql","cil")
  whisu<-c("quantile95","quantile75","max","iqu","ciu")
  
  centv<-c("mean","median","geomean")
  centf<-c("arimean","quantile50","geomean")
  
  f1<-centf[which(centv==centralValueType)]
  f2<-whisl[which(whisk==whiskerValueType)]
  f3<-whisu[which(whisk==whiskerValueType)]
  
  # i cyklus bezi pres sites
  for (i in 1:length(records)) {
    loca<-as.character(records[[i]]$rowLabel)
    value        <-c()
    loqValue     <-c()
    loqMethodCode<-c()
    unit         <-c()
    dateTime     <-c()
    dateTimeString<-c()
    timeLength   <-c()
    
    casovani<-c(casovani,paste0("4. cyklus, ",i,". iterace"),as.character(format(Sys.time(), "%H:%M:%OS3")))
    
    
    for (j in 1:length(records[[i]]$values)) {
      value         <-c(value,         records[[i]]$values[[j]]$value)
      loqValue      <-c(loqValue,      records[[i]]$values[[j]]$loqValue)
      loqMethodCode <-c(loqMethodCode, records[[i]]$values[[j]]$loqMethodCode)
      unit          <-c(unit,          records[[i]]$values[[j]]$unit)
      dateTime      <-c(dateTime,      records[[i]]$values[[j]]$dateTime)
      dateTimeString<-c(dateTimeString,records[[i]]$values[[j]]$dateTimeString)
      timeLength    <-c(timeLength,    records[[i]]$values[[j]]$timeLength)
    }
    
    dateTime<-as.Date(dateTimeString)
    
    value         <-value[order(dateTimeString)]
    loqValue      <-loqValue[order(dateTimeString)]
    loqMethodCode <-loqMethodCode[order(dateTimeString)]
    unit          <-unit[order(dateTimeString)]
    dateTime      <-dateTime[order(dateTimeString)]
    timeLength    <-timeLength[order(dateTimeString)]
    dateTimeString<-dateTimeString[order(dateTimeString)]
    
    # Nahrada LoQ (v promenne valu budou hodnoty vstupujici do vypoctu)
    valu<-value
    valu[which(is.na(valu)&loqMethodCode=="INS")]<-loqValue[which(is.na(valu)&loqMethodCode=="INS")]*1/2
    
    # Jednotky musi byt stejne
    if (length(unique(unit))>1) {
      stop("Units of some records differ!")
    } else {
      unit<-unique(unit)
    }
    
    # Rocni agregace do noveho data.frame aggr
    year<-gendate(substr(dateTimeString,1,4))
    
    aggr<-data.frame(aggregate(valu,by=list(year),FUN=f1)[,2],
                     as.character(unit),
                     centralValueType,
                     aggregate(valu,by=list(year),FUN=f3)[,2],
                     aggregate(valu,by=list(year),FUN=f2)[,2],
                     whiskerValueType,
                     gendate(aggregate(dateTime,by=list(year),FUN=mean)[,1]),
                     as.character(aggregate(dateTime,by=list(year),FUN=mean)[,1]),
                     as.character(aggregate(valu,by=list(year),FUN=length)[,2]),
                     as.character(aggregate(value,by=list(year),FUN=loqlength)[,2]))
    colnames(aggr)<-c("centralValue","unit","centralValueType","whiskerTopValue","whiskerBottomValue","whiskerType","dateTime","dateTimeString","n","nUnderLOQ")
    
    if (nrow(aggr)>1) {
      hole<-3*mean(aggr$dateTime[-1]-aggr$dateTime[-nrow(aggr)],trim=0.05)
    } else {
      hole<-0
    }
    
    if ((max(c(hole,aggr$dateTime[-1]-aggr$dateTime[-nrow(aggr)]))>hole)|(nrow(aggr)<3)) {
      series<-NA
    } else {
      curve<-data.frame(as.Date(as.numeric(genplot(aggr$centralValue,aggr$dateTime,n=20,distr="lnorm",plot=FALSE)$belt[1,]),origin="1970-01-01"),
                        as.numeric(genplot(aggr$centralValue,aggr$dateTime,n=20,distr="lnorm",plot=FALSE)$line[1,]),
                        as.numeric(genplot(aggr$centralValue,aggr$dateTime,n=20,distr="lnorm",plot=FALSE)$lower[1,]),
                        as.numeric(genplot(aggr$centralValue,aggr$dateTime,n=20,distr="lnorm",plot=FALSE)$upper[1,]))
      colnames(curve)<-c("belt","line","lower","upper")
      
      # Logaritmace v pripade log transformace (trend bude linearni)
      if (transformationType=="log") {
        curve$line<-log(curve$line)
        curve$lower<-log(curve$lower)
        curve$upper<-log(curve$upper)
      }
      
      # j cyklus bezi pres jednotlive body krivky
      for (j in 1:nrow(curve)) {
        # Vystupni promenne
        cenValue    <-c(cenValue,curve$line[j])
        botValue    <-c(botValue,curve$lower[j])
        topValue    <-c(topValue,curve$upper[j])
        dateOfPoint <-c(dateOfPoint,as.character(curve$belt[j]))
        nameOfSeries<-c(nameOfSeries,loca)
        segment     <-c(segment,1)
        typeOfSeries<-c(typeOfSeries,"aggr_trend")
        globalUnit  <-c(globalUnit,as.character(unique(aggr$unit)))
      }
    }
    
    # Popis trendovych agregovanych krivek ve 4. cyklu.    
    if (max(c(hole,aggr$dateTime[-1]-aggr$dateTime[-nrow(aggr)]))<=hole) {
      parameterNames<-c("slope",
                        "intercept")
      
      parameterValues<-c(genplot(aggr$centralValue,aggr$dateTime,n=20,distr="lnorm",plot=FALSE)$slope,
                         genplot(aggr$centralValue,aggr$dateTime,n=20,distr="lnorm",plot=FALSE)$intercept)
      
      seriesDescription<-c(seriesDescription,rep(loca,length(parameterNames)))
      parameterDescription<-c(parameterDescription,parameterNames)
      valueDescription<-c(valueDescription,parameterValues)
    }
  }
  

  ## Vypocet prostorove agregovane rady z jednotlivych rocnich agregaci (jen jednou pro cely datovy soubor)
  valu<-as.numeric(cenValue)
  data<-as.Date(dateOfPoint)
  unit<-globalUnit
  
  unit<-unit[order(data)]
  valu<-valu[order(data)]
  data<-data[order(data)]
  
  # Jednotky musi byt stejne
  if (length(unique(unit))>1) {
    stop("Units of some records differ!")
  } else {
    unit<-unique(unit)
  }
  
  whisk<-c("5_95","25_75","min_max","2iq","ci")
  whisl<-c("quantile05","quantile25","min","iql","cil")
  whisu<-c("quantile95","quantile75","max","iqu","ciu")
  
  centv<-c("mean","median","geomean")
  centf<-c("arimean","quantile50","geomean")
  
  f1<-centf[which(centv==centralValueType)]
  f2<-whisl[which(whisk==whiskerValueType)]
  f3<-whisu[which(whisk==whiskerValueType)]
  
  # Rocni agregace do noveho data.frame aggr
  year<-gendate(substr(data,1,4))
  aggr<-data.frame(aggregate(valu,by=list(year),FUN=f1)[,2],
                   as.character(unit),
                   centralValueType,
                   aggregate(valu,by=list(year),FUN=f3)[,2],
                   aggregate(valu,by=list(year),FUN=f2)[,2],
                   whiskerValueType,
                   gendate(aggregate(data,by=list(year),FUN=mean)[,1]),
                   as.character(aggregate(data,by=list(year),FUN=mean)[,1]),
                   as.character(aggregate(valu,by=list(year),FUN=length)[,2]),
                   as.character(aggregate(valu,by=list(year),FUN=loqlength)[,2]))
  colnames(aggr)<-c("centralValue","unit","centralValueType","whiskerTopValue","whiskerBottomValue","whiskerType","dateTime","dateTimeString","n","nUnderLOQ")
  
  if (nrow(aggr)>1) {
    hole<-3*mean(aggr$dateTime[-1]-aggr$dateTime[-nrow(aggr)],trim=0.05)
  } else {
    hole<-0
  }
  
  # k udava poradi segmentu jedne casove rady
  k<-1
  
  # j cyklus bezi pres jednotliva mereni
  for (j in 1:nrow(aggr)) {
    # Vystupni promenne
    cenValue    <-c(cenValue,aggr$centralValue[j])
    botValue    <-c(botValue,aggr$whiskerBottomValue[j])
    topValue    <-c(topValue,aggr$whiskerTopValue[j])
    dateOfPoint <-c(dateOfPoint,aggr$dateTimeString[j])
    nameOfSeries<-c(nameOfSeries,"Total")
    segment     <-c(segment,k)
    typeOfSeries<-c(typeOfSeries,"whol")
    globalUnit  <-c(globalUnit,as.character(unique(aggr$unit)))
    
    if (j!=1) {
      if ((aggr$dateTime[j]-aggr$dateTime[j-1])>hole) {
        k<-k+1
      }
    }
  }

  # Popis celkove agregovane rady
  if (k==1) {
    res<-genstatistic(aggr$centralValue,aggr$dateTime)$res
    
    parameterNames<-c("delta",
                      "mannKendall",
                      "mannKendallP",
                      "daniels",
                      "danielsP",
                      "mean",
                      "sd",
                      "geomean",
                      "gsd",
                      "median",
                      "min",
                      "max",
                      "perc5",
                      "perc25",
                      "perc75",
                      "perc95",
                      "geoMean95CIUpperBound",
                      "geoMean95CILowerBound",
                      "mean95CIUpperBound",
                      "mean95CILowerBound")
    
    parameterValues<-c(res$delta,
                          res$"Mann-Kendall",
                          res$MKp,
                          res$Daniels,
                          res$Dp,
                          res$mean,
                          res$sd,
                          res$"geom. mean",
                          res$"geom. sd",
                          res$median,
                          res$min,
                          res$max,
                          quantile05(aggr$centralValue),
                          quantile25(aggr$centralValue),
                          quantile75(aggr$centralValue),
                          quantile95(aggr$centralValue),
                          res$"geom. mean"*res$"geom. sd"^qnorm(0.975),
                          res$"geom. mean"*res$"geom. sd"^qnorm(0.025),
                          res$mean+qnorm(0.975)*res$sd,
                          res$mean+qnorm(0.025)*res$sd)
    
    seriesDescription<-c(seriesDescription,rep("Total",length(parameterNames)))
    parameterDescription<-c(parameterDescription,parameterNames)
    valueDescription<-c(valueDescription,parameterValues)
  }
  
  casovani<-c(casovani,"Konec",as.character(format(Sys.time(), "%H:%M:%OS3")))
  
  return(list(cenValue,botValue,topValue,dateOfPoint,nameOfSeries,segment,typeOfSeries,seriesDescription,parameterDescription,valueDescription))
}