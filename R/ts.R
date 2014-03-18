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
      dateOfPoint <-c(dateOfPoint,dateTime[j])
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
  
    
  
  
  dateOfPoint<-as.character(as.Date(dateOfPoint,origin="1970-01-01"))
  
  casovani<-c(casovani,"Konec",as.character(format(Sys.time(), "%H:%M:%OS3")))
  
  return(list(cenValue,botValue,topValue,dateOfPoint,nameOfSeries,segment,typeOfSeries,seriesDescription,parameterDescription,valueDescription))
}