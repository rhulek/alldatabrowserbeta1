cil<-function(x) {
  il<-mean(x,na.rm=TRUE)+qnorm(0.025)*sd(x,na.rm=TRUE)/sqrt(length(which(!is.na(x))))
  return(il)
}



