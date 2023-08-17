
BINSEG = function(sumstat, pen = 0, cost_func = "norm.mean", shape = 1, minseglen = 2,  Q=5){
  
  n=length(sumstat[,1])-1
  if(n<2){stop('Data must have atleast 2 observations to fit a changepoint model.')}
  if(Q>((n/2)+1)){stop(paste('Q is larger than the maximum number of segments',(n/2)+1))}
  if(Q<=0){stop(paste('Q is the maximum number of changepoints so should be greater than 0'))}
  
  storage.mode(sumstat) = 'double'
  
  cptsout=rep(0,Q) # sets up null vector for changepoint answer
  likeout=rep(0,Q) # sets up null vector for likelihood of changepoints in cptsout
  storage.mode(cptsout)='integer'
  storage.mode(likeout)='double'
  op_cps=0
  
  #on.exit(.C("FreeBinSeg",answer[[6]],PACKAGE='changepoint'))
  # answer=.C('PELT',cost_func, y3, y2,y,as.integer(n),as.double(pen),cptsout,as.integer(error),as.double(shape))
  answer=.C('binseg',cost_func, sumstat,as.integer(n),as.double(pen),as.integer(Q),cptsout, as.integer(minseglen), likeout, as.integer(op_cps), as.double(shape))
  if(answer[[9]]==Q){warning('The number of changepoints identified is Q, it is advised to increase Q to make sure changepoints have not been missed.')}
  if(answer[[9]]==0){cpts=n}
  else{cpts=c(sort(answer[[6]][1:answer[[9]]]),n)}
  return(list(cps=rbind(answer[[6]],2*answer[[8]]),cpts=cpts,op.cpts=answer[[9]],pen=pen))
  ##answer[6] is cptsout, answer[8] is likeout ("beta value")
}
CROPS <- function(data, penalty="CROPS", pen.value, method="PELT", test.stat="Normal", class=TRUE, param.est=TRUE, minseglen, shape, func){
  if(method != "PELT"){stop('CROPS is a valid penalty choice only if method="PELT", please change your method or your penalty.')}
  mu <- mean(data)
  sumstat=cbind(c(0,cumsum(coredata(data))),c(0,cumsum(coredata(data)^2)),cumsum(c(0,(coredata(data)-mu)^2)))
  
  switch(test.stat,
         "Normal" = {stat = "norm"},
         "Exponential" = {stat = "exp"},
         "Gamma" = {stat = "gamma"},
         "Poisson" = {stat = "poisson"},
         {stop("Only Normal, Exponential, Gamma and Poisson are valid test statistics")}
  )
  costfunc = paste0(func, ".", stat)
  
  out = range_of_penalties(sumstat, cost=costfunc, min_pen=pen.value[1], max_pen=pen.value[2], minseglen=minseglen)
  
  if(func=="var"){
    cpttype="variance"
  }else if(func=="meanvar"){
    cpttype="mean and variance"
  }else{
    cpttype="mean"
  }
  
  if(class==TRUE){
    ans = class_input(data=data,cpttype=cpttype, method="PELT", test.stat=test.stat, penalty=penalty, pen.value=pen.value, minseglen=minseglen, param.estimates=param.est, out=out,shape=shape)
    if(func=="var"){
      param.est(ans)=c(param.est(ans),mean=mu)
    }
    return(ans)
  }else{return(out)}
}
PELT = function(sumstat,pen=0, cost_func = "norm.mean", shape = 1, minseglen = 1){
  # function that uses the PELT method to calculate changes in mean where the segments in the data are assumed to be Normal
  n = length(sumstat[,1]) - 1
  if(n<2){stop('Data must have atleast 2 observations to fit a changepoint model.')}
  
  storage.mode(sumstat) = 'double'
  error=0
  
  lastchangelike = array(0,dim = n+1)
  lastchangecpts = array(0,dim = n+1)
  numchangecpts = array(0,dim = n+1)
  
  cptsout=rep(0,n) # sets up null vector for changepoint answer
  storage.mode(cptsout)='integer'
  
  answer=list()
  answer[[6]]=1
  on.exit(.C("FreePELT",answer[[6]]))
  
  storage.mode(lastchangelike) = 'double'
  storage.mode(lastchangecpts) = 'integer'
  storage.mode(numchangecpts) = 'integer'
  
  # answer=.C('PELT',cost_func, y3, y2,y,as.integer(n),as.double(pen),cptsout,as.integer(error),as.double(shape))
  answer=.C('PELT',cost_func, sumstat,as.integer(n),as.double(pen),cptsout,as.integer(error),as.double(shape), as.integer(minseglen), lastchangelike, lastchangecpts,numchangecpts)
  
  if(answer[[6]]>0){
    stop("C code error:",answer[[6]],call.=F)
  }
  return(list(lastchangecpts=answer[[10]],cpts=sort(answer[[5]][answer[[5]]>0]), lastchangelike=answer[[9]], ncpts=answer[[11]]))
  
}
class_input <- function(data, cpttype, method, test.stat, penalty, pen.value, minseglen, param.estimates, out=list(), Q=NA, shape=NA){
  if(method=="BinSeg" || penalty=="CROPS"){
    ans=new("cpt.range")
  }else{
    ans=new("cpt")
  }
  
  data.set(ans)=data;cpttype(ans)=cpttype;method(ans)=method; test.stat(ans)=test.stat;pen.type(ans)=penalty;pen.value(ans)=pen.value;minseglen(ans)=minseglen;ans@date=date();
  if(penalty!="CROPS"){ # crops is only one that doesn't give a single set of cpts
    cpts(ans)=out[[2]]
    
    if(param.estimates==TRUE){
      if(test.stat == "Gamma"){
        ans=param(ans, shape)
      }else{
        ans=param(ans)
      }
    }
  }
  
  if(method=="PELT"){
    ncpts.max(ans)=Inf
  }
  else if(method=="AMOC"){
    ncpts.max(ans)=1
  }
  else{
    ncpts.max(ans)=Q
  }
  
  if(method=="BinSeg"){
    l=list()
    for(i in 1:(length(out$cps)/2)){
      l[[i]] = out$cps[1,1:i] 
    }
    m = t(sapply(l, '[', 1:max(sapply(l, length))))
    
    cpts.full(ans)=m
    pen.value.full(ans)=out$cps[2,]
  }else if(penalty=="CROPS"){
    m = t(sapply(out[[2]], '[', 1:max(sapply(out[[2]], length))))
    
    cpts.full(ans) = m
    pen.value.full(ans) = out[[1]][1,]
    if(test.stat=="Gamma"){param.est(ans)$shape=shape}
  }
  
  return(ans)
}
cpt.mean=function(data,penalty="MBIC",pen.value=0,method="AMOC",Q=5,test.stat="Normal",class=TRUE,param.estimates=TRUE,minseglen=1){
  checkData(data)
  if(minseglen<1){minseglen=1;warning('Minimum segment length for a change in mean is 1, automatically changed to be 1.')}
  if(test.stat=="CUSUM"){stop("CUSUM test has moved to the changepoint.np package.")}
  else if(!(test.stat=="Normal")){ stop("Invalid test statistic, must be Normal.") }
  
  if(penalty == "CROPS"){
    # browser()
    if(is.numeric(pen.value)){
      if(length(pen.value) == 2){
        if(pen.value[2] < pen.value[1]){
          pen.value = rev(pen.value)
        }
        #run range of penalties
        return(CROPS(data=data, method=method, pen.value=pen.value, test.stat=test.stat, class=class, param.est=param.estimates, minseglen=minseglen, func="mean"))
      }else{
        stop('The length of pen.value must be 2')
      }
    }else{
      stop('For CROPS, pen.value must be supplied as a numeric vector and must be of length 2')
    }
  }
  
  if(test.stat=="Normal"){
    if(method=="AMOC"){
      return(single.mean.norm(data,penalty,pen.value,class,param.estimates,minseglen))
    }
    else if(method=="PELT" || method=="BinSeg"){
      
      return(multiple.mean.norm(data,mul.method=method,penalty,pen.value,Q,class,param.estimates,minseglen)) 
    }
    else{
      stop("Invalid Method, must be AMOC, PELT or BinSeg.")
    }
  }
}
#cpt.reg=function(data,penalty="MBIC",pen.value=0,method="AMOC",Q=5,test.stat="Normal",class=TRUE,param.estimates=TRUE){
#	if(test.stat !="Normal"){ stop("Invalid test statistic, must be Normal") }
#	if(method=="AMOC"){
#		return(single.reg.norm(data,penalty,pen.value,class,param.estimates))
#	}
#	else if(method=="PELT" || method=="BinSeg"){
#		return(multiple.reg.norm(data,mul.method=method,penalty,pen.value,Q,class,param.estimates))
#	}
#	else if(method=="SegNeigh"){
#		warning("SegNeigh is computationally slow, use PELT instead")
#		return(multiple.reg.norm(data,mul.method=method,penalty,pen.value,Q,class,param.estimates))
#	}
#	else{
#		stop("Invalid Method, must be AMOC, PELT, SegNeigh or BinSeg")
#	}
#}
cpt.var=function(data,penalty="MBIC",pen.value=0,know.mean=FALSE, mu=NA,method="AMOC",Q=5,test.stat="Normal",class=TRUE,param.estimates=TRUE,minseglen=2){
  checkData(data)
  if(minseglen<2){minseglen=2;warning('Minimum segment length for a change in variance is 2, automatically changed to be 2.')}
  
  if(penalty == "CROPS"){
    # browser()
    if(is.numeric(pen.value)){
      if(length(pen.value) == 2){
        if(pen.value[2] < pen.value[1]){
          pen.value = rev(pen.value)
        }
        #run range of penalties
        return(CROPS(data=data, method=method,pen.value=pen.value, test.stat=test.stat, class=class, param.est=param.estimates, minseglen=minseglen, func="var"))
      }else{
        stop('The length of pen.value must be 2')
      }
    }else{
      stop('For CROPS, pen.value must be supplied as a numeric vector and must be of length 2')
    }
  }
  
  if(test.stat =="Normal"){
    
    if(method=="AMOC"){
      return(single.var.norm(data,penalty,pen.value,know.mean,mu,class,param.estimates,minseglen))
    }
    else if(method=="PELT" || method=="BinSeg"){
      
      return(multiple.var.norm(data,mul.method=method,penalty,pen.value,Q,know.mean,mu,class,param.estimates,minseglen))
    }
    else{
      stop("Invalid Method, must be AMOC, PELT, BinSeg.")
    }
  }
  else if(test.stat=="CSS"){
    stop("CSS test has moved to the changepoint.np package.")
  }
  else{
    stop("Invalid test statistic, must be Normal.")
  }
}
cpt.meanvar=function(data,penalty="MBIC",pen.value=0,method="AMOC",Q=5,test.stat="Normal",class=TRUE,param.estimates=TRUE,shape=1,minseglen=2){
  checkData(data)
  if(minseglen<2){
    if(!(minseglen==1 & (test.stat=="Poisson"|test.stat=="Exponential"))){
      minseglen=2;warning('Minimum segment length for a change in mean and variance is 2, automatically changed to be 2.')}
  }
  if(penalty == "CROPS"){
    if(is.numeric(pen.value)){
      if(length(pen.value) == 2){
        if(pen.value[2] < pen.value[1]){
          pen.value = rev(pen.value)
        }
        #run range of penalties
        return(CROPS(data=data, method=method,pen.value=pen.value, test.stat=test.stat, class=class, param.est=param.estimates, minseglen=minseglen, shape=shape, func="meanvar"))
      }else{
        stop('The length of pen.value must be 2')
      }
    }else{
      stop('For CROPS, pen.value must be supplied as a numeric vector and must be of length 2')
    }
  }
  if(test.stat=="Normal"){
    
    if(method=="AMOC"){
      return(single.meanvar.norm(data,penalty,pen.value,class,param.estimates,minseglen))
    }
    else if(method=="PELT" || method=="BinSeg"){
      
      return(multiple.meanvar.norm(data,mul.method=method,penalty,pen.value,Q,class,param.estimates,minseglen))
    }
    else{
      stop("Invalid Method, must be AMOC, PELT or BinSeg.")
    }
  }
  else if(test.stat=="Gamma"){
    if(method=="AMOC"){
      return(single.meanvar.gamma(data,shape,penalty,pen.value,class,param.estimates,minseglen))
    }
    else if(method=="PELT" || method=="BinSeg"){
      return(multiple.meanvar.gamma(data,shape,mul.method=method,penalty,pen.value,Q,class,param.estimates,minseglen))
    }
    else{
      stop("Invalid Method, must be AMOC, PELT or BinSeg.")
    }
  }
  else if(test.stat=="Exponential"){
    if(method=="AMOC"){
      return(single.meanvar.exp(data,penalty,pen.value,class,param.estimates,minseglen))
    }
    else if(method=="PELT" || method=="BinSeg"){
      
      return(multiple.meanvar.exp(data,mul.method=method,penalty,pen.value,Q,class,param.estimates,minseglen))
    }
    else{
      stop("Invalid Method, must be AMOC, PELT or BinSeg.")
    }
  }
  else if(test.stat=="Poisson"){
    if(method=="AMOC"){
      return(single.meanvar.poisson(data,penalty,pen.value,class,param.estimates,minseglen))
    }
    else if(method=="PELT" || method=="BinSeg"){
      
      return(multiple.meanvar.poisson(data,mul.method=method,penalty,pen.value,Q,class,param.estimates,minseglen))
    }
    else{
      stop("Invalid Method, must be AMOC, PELT or BinSeg.")
    }
  }
  else{
    stop("Invalid test statistic, must be Normal, Gamma, Exponential or Poisson.")
  }
}
checkData = function(data){
  if(!is.numeric(data)){
    stop("Only numeric data allowed")
  }  
  if(anyNA(data)){stop("Missing value: NA is not allowed in the data as changepoint methods assume regularly spaced data.")}
  
}
setClass("cpt",slots=list(data.set="ts", cpttype="character", method="character", 	test.stat="character",pen.type="character",pen.value="numeric",minseglen="numeric",cpts="numeric",ncpts.max="numeric",param.est="list",date="character",version="character"),prototype=prototype(cpttype="Not Set",date=date(),version=as(packageVersion("changepoint"),'character')))

setClass("cpt.reg",slots=list(data.set="matrix", cpttype="character", method="character", test.stat="character",pen.type="character",pen.value="numeric",minseglen="numeric",cpts="numeric",ncpts.max="numeric",param.est="list",date="character",version="character"),prototype=prototype(cpttype="regression",date=date(),version=as(packageVersion("changepoint"),"character")))

#   setClass("cpt", representation(), prototype())
# # cpts is the optimal segementation
#   
setClass("cpt.range",slots=list(cpts.full="matrix", pen.value.full="numeric"), prototype=prototype(), contains="cpt")
# cpts.full is the entire matrix
# pen.value.full (beta) values as an extra slot (vector)

# retrival functions for slots
if(!isGeneric("data.set")) {
  if (is.function("data.set")){
    fun <- data.set
  }
  else {fun <- function(object){
    standardGeneric("data.set")
  }
  }
  setGeneric("data.set", fun)
}
setMethod("data.set","cpt",function(object) coredata(object@data.set))
setMethod("data.set","cpt.reg",function(object) coredata(object@data.set))

if(!isGeneric("data.set.ts")) {
  if (is.function("data.set.ts")){
    fun <- data.set.ts
  }
  else {fun <- function(object){
    standardGeneric("data.set.ts")
  }
  }
  setGeneric("data.set.ts", fun)
}
setMethod("data.set.ts","cpt",function(object) object@data.set)

if(!isGeneric("cpttype")) {
  if (is.function("cpttype")){
    fun <- cpttype
  }
  else {fun <- function(object){
    standardGeneric("cpttype")
  }
  }
  setGeneric("cpttype", fun)
}
setMethod("cpttype","cpt",function(object) object@cpttype)
setMethod("cpttype","cpt.reg",function(object) object@cpttype)

if(!isGeneric("method")) {
  if (is.function("method")){
    fun <- method
  }
  else {fun <- function(object){
    standardGeneric("method")
  }
  }
  setGeneric("method", fun)
}
setMethod("method","cpt",function(object) object@method)
setMethod("method","cpt.reg",function(object) object@method)

# distribution remains for backwards compatability, changed to test.stat version 1.0
if(!isGeneric("distribution")) {
  if (is.function("distribution")){
    fun <- distribution
  }
  else {fun <- function(object){
    standardGeneric("distribution")
  }
  }
  setGeneric("distribution", fun)
}
setMethod("distribution","cpt",function(object) object@test.stat)
setMethod("distribution","cpt.reg",function(object) object@test.stat)

if(!isGeneric("test.stat")) {
  if (is.function("test.stat")){
    fun <- test.stat
  }
  else {fun <- function(object){
    standardGeneric("test.stat")
  }
  }
  setGeneric("test.stat", fun)
}
setMethod("test.stat","cpt",function(object) object@test.stat)
setMethod("test.stat","cpt.reg",function(object) object@test.stat)

if(!isGeneric("pen.type")) {
  if (is.function("pen.type")){
    fun <- pen.type
  }
  else {fun <- function(object){
    standardGeneric("pen.type")
  }
  }
  setGeneric("pen.type", fun)
}
setMethod("pen.type","cpt",function(object) object@pen.type)
setMethod("pen.type","cpt.reg",function(object) object@pen.type)

if(!isGeneric("pen.value")) {
  if (is.function("pen.value")){
    fun <- pen.value
  }
  else {fun <- function(object){
    standardGeneric("pen.value")
  }
  }
  setGeneric("pen.value", fun)
}
setMethod("pen.value","cpt",function(object) object@pen.value)
setMethod("pen.value","cpt.reg",function(object) object@pen.value)

if(!isGeneric("pen.value.full")) {
  if (is.function("pen.value.full")){
    fun <- pen.value.full
  }
  else {fun <- function(object){
    standardGeneric("pen.value.full")
  }
  }
  setGeneric("pen.value.full", fun)
}
setMethod("pen.value.full","cpt.range",function(object) object@pen.value.full)

if(!isGeneric("minseglen")) {
  if (is.function("minseglen")){
    fun <- minseglen
  }
  else {fun <- function(object){
    standardGeneric("minseglen")
  }
  }
  setGeneric("minseglen", fun)
}
setMethod("minseglen","cpt",function(object) object@minseglen)

if(!isGeneric("cpts")) {
  if (is.function("cpts")){
    fun <- cpts
  }
  else {fun <- function(object){
    standardGeneric("cpts")
  }
  }
  setGeneric("cpts", fun)
}
setMethod("cpts","cpt",function(object) object@cpts[-length(object@cpts)])
setMethod("cpts","cpt.reg",function(object) object@cpts[-length(object@cpts)])

if(!isGeneric("cpts.full")) {
  if (is.function("cpts.full")){
    fun <- cpts.full
  }
  else {fun <- function(object){
    standardGeneric("cpts.full")
  }
  }
  setGeneric("cpts.full", fun)
}
setMethod("cpts.full","cpt.range",function(object) object@cpts.full)

if(!isGeneric("cpts.ts")) {
  if (is.function("cpts.ts")){
    fun <- cpts.ts
  }
  else {fun <- function(object){
    standardGeneric("cpts.ts")
  }
  }
  setGeneric("cpts.ts", fun)
}
setMethod("cpts.ts","cpt",function(object) index(data.set.ts(object))[cpts(object)] )

if(!isGeneric("ncpts.max")) {
  if (is.function("ncpts.max")){
    fun <- ncpts.max
  }
  else {fun <- function(object){
    standardGeneric("ncpts.max")
  }
  }
  setGeneric("ncpts.max", fun)
}
setMethod("ncpts.max","cpt",function(object) object@ncpts.max)
setMethod("ncpts.max","cpt.reg",function(object) object@ncpts.max)

if(!isGeneric("param.est")) {
  if (is.function("param.est")){
    fun <- param.est
  }
  else {fun <- function(object){
    standardGeneric("param.est")
  }
  }
  setGeneric("param.est", fun)
}
setMethod("param.est","cpt",function(object) object@param.est)
setMethod("param.est","cpt.reg",function(object) object@param.est)


setMethod("coef","cpt",function(object) object@param.est)
setMethod("coef","cpt.reg",function(object) object@param.est)

# ncpts function
if(!isGeneric("ncpts")) {
  if (is.function("ncpts")){
    fun <- ncpts
  }
  else {fun <- function(object){
    standardGeneric("ncpts")
  }
  }
  setGeneric("ncpts", fun)
}
setMethod("ncpts","cpt",function(object) length(cpts(object)))
setMethod("ncpts","cpt.reg",function(object) length(cpts(object)))

# seg.len function
if(!isGeneric("seg.len")) {
  if (is.function("seg.len")){
    fun <- seg.len
  }
  else {fun <- function(object){
    standardGeneric("seg.len")
  }
  }
  setGeneric("seg.len", fun)
}
setMethod("seg.len","cpt",function(object){object@cpts-c(0,object@cpts[-length(object@cpts)])})
setMethod("seg.len","cpt.reg",function(object){object@cpts-c(0,object@cpts[-length(object@cpts)])})
#i.e. if there is a changepoint in the data, return segment length. If not, return length of the data

# nseg function
if(!isGeneric("nseg")) {
  if (is.function("nseg")){
    fun <- nseg
  }
  else {fun <- function(object){
    standardGeneric("nseg")
  }
  }
  setGeneric("nseg", fun)
}
setMethod("nseg","cpt",function(object){ncpts(object)+1})
setMethod("nseg","cpt.reg",function(object){ncpts(object)+1})


# replacement functions for slots
setGeneric("data.set<-", function(object, value) standardGeneric("data.set<-"))
setReplaceMethod("data.set", "cpt", function(object, value) {
  if(is.ts(value)){object@data.set <- value}else{object@data.set <- ts(value)}
  return(object)
})
setReplaceMethod("data.set", "cpt.reg", function(object, value) {
  object@data.set <- value
  return(object)
})

setGeneric("cpttype<-", function(object, value) standardGeneric("cpttype<-"))
setReplaceMethod("cpttype", "cpt", function(object, value) {
  object@cpttype <- value
  return(object)
})
setReplaceMethod("cpttype", "cpt.reg", function(object, value) {
  object@cpttype <- value
  return(object)
})

setGeneric("method<-", function(object, value) standardGeneric("method<-"))
setReplaceMethod("method", "cpt", function(object, value) {
  object@method <- value
  return(object)
})
setReplaceMethod("method", "cpt.reg", function(object, value) {
  object@method <- value
  return(object)
})

# distribution remains for backwards compatability, changed to test.stat version 1.0
setGeneric("distribution<-", function(object, value) standardGeneric("distribution<-"))
setReplaceMethod("distribution", "cpt", function(object, value) {
  object@test.stat <- value
  return(object)
})
setReplaceMethod("distribution", "cpt.reg", function(object, value) {
  object@test.stat <- value
  return(object)
})

setGeneric("test.stat<-", function(object, value) standardGeneric("test.stat<-"))
setReplaceMethod("test.stat", "cpt", function(object, value) {
  object@test.stat <- value
  return(object)
})
setReplaceMethod("test.stat", "cpt.reg", function(object, value) {
  object@test.stat <- value
  return(object)
})

setGeneric("pen.type<-", function(object, value) standardGeneric("pen.type<-"))
setReplaceMethod("pen.type", "cpt", function(object, value) {
  object@pen.type <- value
  return(object)
})
setReplaceMethod("pen.type", "cpt.reg", function(object, value) {
  object@pen.type <- value
  return(object)
})

setGeneric("pen.value<-", function(object, value) standardGeneric("pen.value<-"))
setReplaceMethod("pen.value", "cpt", function(object, value) {
  object@pen.value <- value
  return(object)
})
setReplaceMethod("pen.value", "cpt.reg", function(object, value) {
  object@pen.value <- value
  return(object)
})

setGeneric("minseglen<-", function(object, value) standardGeneric("minseglen<-"))
setReplaceMethod("minseglen", "cpt", function(object, value) {
  object@minseglen <- value
  return(object)
})
setReplaceMethod("minseglen", "cpt.range", function(object, value) {
  object@minseglen <- value
  return(object)
})
setReplaceMethod("minseglen", "cpt.reg", function(object, value) {
  object@minseglen <- value
  return(object)
})

setGeneric("cpts<-", function(object, value) standardGeneric("cpts<-"))
setReplaceMethod("cpts", "cpt", function(object, value) {
  if((cpttype(object)=="meanar")|(cpttype(object)=="trendar")){
    n=length(object@data.set)-1
  }
  else{n=length(object@data.set)}
  
  if(value[length(value)]==n){object@cpts <- value}
  else{		object@cpts <- c(value,n)  }
  return(object)
})
setReplaceMethod("cpts", "cpt.reg", function(object, value) {
  if(value[length(value)]==nrow(object@data.set)){object@cpts <- value}
  else{  	object@cpts <- c(value,nrow(object@data.set))  }
  return(object)
})

setGeneric("ncpts.max<-", function(object, value) standardGeneric("ncpts.max<-"))
setReplaceMethod("ncpts.max", "cpt", function(object, value) {
  object@ncpts.max <- value
  return(object)
})
setReplaceMethod("ncpts.max", "cpt.reg", function(object, value) {
  object@ncpts.max <- value
  return(object)
})

setGeneric("param.est<-", function(object, value) standardGeneric("param.est<-"))
setReplaceMethod("param.est", "cpt", function(object, value) {
  object@param.est <- value
  return(object)
})
setReplaceMethod("param.est", "cpt.reg", function(object, value) {
  object@param.est <- value
  return(object)
})

setGeneric("cpts.full<-", function(object, value) standardGeneric("cpts.full<-"))
setReplaceMethod("cpts.full", "cpt.range", function(object, value) {
  object@cpts.full <- value
  return(object)
})
setGeneric("pen.value.full<-", function(object, value) standardGeneric("pen.value.full<-"))
setReplaceMethod("pen.value.full", "cpt.range", function(object, value) {
  object@pen.value.full <- value
  return(object)
})
setGeneric("pen.value.input<-", function(object, value) standardGeneric("pen.value.input<-"))
setReplaceMethod("pen.value.input", "cpt", function(object, value) {
  object@pen.value.input <- value
  return(object)
})


# parameter functions
setGeneric("param", function(object,...) standardGeneric("param"))
setMethod("param", "cpt", function(object,shape,...) {			
  param.mean=function(object){
    cpts=c(0,object@cpts)
    #nseg=length(cpts)-1
    data=data.set(object)
    tmpmean=NULL
    for(j in 1:nseg(object)){
      tmpmean[j]=mean(data[(cpts[j]+1):(cpts[j+1])])
    }
    return(tmpmean)
  }
  param.var=function(object){
    cpts=c(0,object@cpts)
    #nseg=length(cpts)-1
    data=data.set(object)
    seglen=seg.len(object)
    tmpvar=NULL
    for(j in 1:nseg(object)){
      tmpvar[j]=var(data[(cpts[j]+1):(cpts[j+1])])
    }
    tmpvar=tmpvar*(seglen-1)/seglen # correctly for the fact that the MLE estimate is /n but the var function is /n-1
    return(tmpvar)
  }
  param.scale=function(object,shape){
    cpts=c(0,object@cpts)
    #nseg=length(cpts)-1
    data=data.set(object)
    y=c(0,cumsum(data))
    tmpscale=NULL
    for(j in 1:nseg(object)){
      tmpscale[j]=(y[(cpts[j+1]+1)]-y[(cpts[j]+1)])/((cpts[j+1]-cpts[j])*shape)
    }
    return(tmpscale)			
  }
  param.trend=function(object){
    cpts=c(0,object@cpts)
    seglen=seg.len(object)
    data=data.set(object)
    n=length(data)
    sumstat=cbind(cumsum(c(0,data)),cumsum(c(0,data*c(1:n))))
    cptsumstat=matrix(sumstat[object@cpts+1,]-sumstat[c(0,cpts(object))+1,],ncol=2)
    cptsumstat[,2]=cptsumstat[,2]-cptsumstat[,1]*c(0,cpts(object)) # i.e. creating newx3
    
    thetaS=(2*cptsumstat[,1]*(2*seglen + 1) - 6*cptsumstat[,2]) / (2*seglen*(2*seglen + 1) - 3*seglen*(seglen+1))
    thetaT=(6*cptsumstat[,2])/((seglen+1)*(2*seglen+1)) + (thetaS * (1-((3*seglen)/((2*seglen)+1))))
    return(cbind(thetaS,thetaT))
  }
  param.meanar=function(object){
    seglen=seg.len(object)
    data=data.set(object)
    n=length(data)-1
    sumstat=cbind(cumsum(c(0,data[-1])),cumsum(c(0,data[-(n+1)])),cumsum(c(0,data[-1]*data[-(n+1)])),cumsum(c(0,data[-1]^2)),cumsum(c(0,data[-(n+1)]^2)))
    cptsumstat=matrix(sumstat[object@cpts+1,]-sumstat[c(0,cpts(object))+1,],ncol=5)
    beta2=(2*seglen*cptsumstat[,3]-cptsumstat[,1]*cptsumstat[,2])/(2*seglen*cptsumstat[,5]*(1-cptsumstat[,2]^2));
    beta1=(2*cptsumstat[,1]-beta2*cptsumstat[,2])/(2*seglen);
    
    return(cbind(beta1,beta2))
  }
  param.trendar=function(object){
    seglen=seg.len(object)
    data=data.set(object)
    n=length(data)-1
    sumstat=cbind(cumsum(c(0,data[-1])),cumsum(c(0,data[-(n+1)])),cumsum(c(0,data[-1]*data[-(n+1)])),cumsum(c(0,data[-1]*c(1:n))),cumsum(c(0,data[-(n+1)]*c(0:(n-1)))),cumsum(c(0,data[-1]^2)),cumsum(c(0,data[-(n+1)]^2)))
    cptsumstat=matrix(sumstat[object@cpts+1,]-sumstat[c(0,cpts(object))+1,],ncol=7)
    cptsumstat[,4]=cptsumstat[,4]-cptsumstat[,1]*c(0,cpts(object)) # i.e. creating newx4
    cptsumstat[,5]=cptsumstat[,5]-cptsumstat[,2]*c(0,cpts(object)) # i.e. creating newx5
    betatop=seglen*(seglen-1)*(seglen*(seglen-1)*cptsumstat[,3] + 2*(2*seglen+1)*cptsumstat[,1]*(cptsumstat[,5]-seglen*cptsumstat[,2]) + 6*cptsumstat[,4]*(cptsumstat[,2]-cptsumstat[,5]))
    betabottom=seglen*(seglen-1)*cptsumstat[,7] + 2*(2*seglen+1)*cptsumstat[,2]*(seglen*cptsumstat[,2]-cptsumstat[,5]) + 6*cptsumstat[,5]*(cptsumstat[,5]-cptsumstat[,2]);
    beta=betatop/betabottom;  
    thetajpo=(6*(seglen+2)*(cptsumstat[,4]-beta*cptsumstat[,5]))/((seglen+1)*(2*seglen+1)) - 2*(cptsumstat[,1]-beta*cptsumstat[,2])
    thetaj=(2*(2*seglen+1)*(cptsumstat[,1]-beta*cptsumstat[,2])-6*(cptsumstat[,4]-beta*cptsumstat[,5]))/(seglen-1)
    
    return(cbind(beta,thetajpo,thetaj))
  }
  if(cpttype(object)=="mean"){
    param.est(object)<-list(mean=param.mean(object))
  }
  else if(cpttype(object)=="variance"){
    param.est(object)<-list(variance=param.var(object))
  }
  else if(cpttype(object)=="mean and variance"){
    if(test.stat(object)=="Normal"){
      param.est(object)<-list(mean=param.mean(object),variance=param.var(object))
    }
    else if(test.stat(object)=="Gamma"){
      param.est(object)<-list(scale=param.scale(object,shape=shape),shape=shape)
    }
    else if(test.stat(object)=="Exponential"){
      param.est(object)<-list(rate=1/param.mean(object))
    }
    else if(test.stat(object)=="Poisson"){
      param.est(object)<-list(lambda=param.mean(object))
    }
    else{
      stop("Unknown test statistic for a change in mean and variance")
    }
  }
  else if(cpttype(object)=="trend"){
    if(test.stat(object)=="Normal"){
      tmp=param.trend(object)
      param.est(object)<-list(thetaS=tmp[,1],thetaT=tmp[,2])
    }
    else{
      stop("Unknown test statistic for a change in trend")
    }
  }
  else if(cpttype(object)=="trendar"){
    if(test.stat(object)=="Normal"){
      tmp=param.trendar(object)
      param.est(object)<-list(beta=tmp[,1],thetajpo=tmp[,2],thetaj=tmp[,3])
    }
    else{
      stop("Unknown test statistic for a change in trend+ar")
    }
  }
  else if(cpttype(object)=="meanar"){
    if(test.stat(object)=="Normal"){
      tmp=param.meanar(object)
      param.est(object)<-list(beta1=tmp[,1],beta2=tmp[,2])
    }
    else{
      stop("Unknown test statistic for a change in mean+ar")
    }
  }
  else{
    stop("Unknown changepoint type, must be 'mean', 'variance', 'mean and variance', 'trend', 'meanar' or 'trendar'.")
  }
  return(object)
})

setMethod("param", "cpt.range", function(object,ncpts=NA,shape,...) {
  if(is.na(ncpts)){
    cpts=c(0,object@cpts)
  }
  else{
    ncpts.full=apply(cpts.full(object),1,function(x){sum(x>0,na.rm=TRUE)})
    row=try(which(ncpts.full==ncpts),silent=TRUE)
    if(class(row)=='try-error'){
      stop("Your input object doesn't have a segmentation with the requested number of changepoints.")
    }
    cpts=c(0,cpts.full(object)[row,1:ncpts],length(data.set(object)))
  }
  
  param.mean=function(object,cpts){
    nseg=length(cpts)-1
    data=data.set(object)
    tmpmean=NULL
    for(j in 1:nseg){
      tmpmean[j]=mean(data[(cpts[j]+1):(cpts[j+1])])
    }
    return(tmpmean)
  }
  param.var=function(object,cpts){
    nseg=length(cpts)-1
    data=data.set(object)
    seglen=cpts[-1]-cpts[-length(cpts)]
    tmpvar=NULL
    for(j in 1:nseg){
      tmpvar[j]=var(data[(cpts[j]+1):(cpts[j+1])])
    }
    tmpvar=tmpvar*(seglen-1)/seglen
    return(tmpvar)
  }
  param.scale=function(object,cpts,shape){
    nseg=length(cpts)-1
    data=data.set(object)
    y=c(0,cumsum(data))
    tmpscale=NULL
    for(j in 1:nseg){
      tmpscale[j]=(y[(cpts[j+1]+1)]-y[(cpts[j]+1)])/((cpts[j+1]-cpts[j])*shape)
    }
    return(tmpscale)			
  }
  param.trend=function(object,cpts){
    seglen=cpts[-1]-cpts[-length(cpts)]
    data=data.set(object)
    n=length(data)
    sumstat=cbind(cumsum(c(0,data)),cumsum(c(0,data*c(1:n))))
    cptsumstat=matrix(sumstat[object@cpts+1,]-sumstat[c(0,cpts(object))+1,],ncol=2)
    cptsumstat[,2]=cptsumstat[,2]-cptsumstat[,1]*c(0,cpts(object)) # i.e. creating newx3
    
    thetaS=(2*cptsumstat[,1]*(2*seglen + 1) - 6*cptsumstat[,2]) / (2*seglen*(2*seglen + 1) - 3*seglen*(seglen+1))
    thetaT=(6*cptsumstat[,2])/((seglen+1)*(2*seglen+1)) + (thetaS * (1-((3*seglen)/((2*seglen)+1))))
    return(cbind(thetaS,thetaT))
  }
  param.meanar=function(object,cpts){
    seglen=cpts[-1]-cpts[-length(cpts)]
    data=data.set(object)
    n=length(data)-1
    sumstat=cbind(cumsum(c(0,data[-1])),cumsum(c(0,data[-(n+1)])),cumsum(c(0,data[-1]*data[-(n+1)])),cumsum(c(0,data[-1]^2)),cumsum(c(0,data[-(n+1)]^2)))
    cptsumstat=matrix(sumstat[object@cpts+1,]-sumstat[c(0,cpts(object))+1,],ncol=5)
    beta2=(2*seglen*cptsumstat[,3]-cptsumstat[,1]*cptsumstat[,2])/(2*seglen*cptsumstat[,5]*(1-cptsumstat[,2]^2));
    beta1=(2*cptsumstat[,1]-beta2*cptsumstat[,2])/(2*seglen);
    
    return(cbind(beta1,beta2))
  }
  param.trendar=function(object,cpts){
    seglen=cpts[-1]-cpts[-length(cpts)]
    data=data.set(object)
    n=length(data)-1
    sumstat=cbind(cumsum(c(0,data[-1])),cumsum(c(0,data[-(n+1)])),cumsum(c(0,data[-1]*data[-(n+1)])),cumsum(c(0,data[-1]*c(1:n))),cumsum(c(0,data[-(n+1)]*c(0:(n-1)))))
    cptsumstat=matrix(sumstat[object@cpts+1,]-sumstat[c(0,cpts(object))+1,],ncol=7)
    cptsumstat[,4]=cptsumstat[,4]-cptsumstat[,1]*c(0,cpts(object)) # i.e. creating newx4
    cptsumstat[,5]=cptsumstat[,5]-cptsumstat[,2]*c(0,cpts(object)) # i.e. creating newx5
    betatop=seglen*(seglen-1)*(seglen*(seglen-1)*cptsumstat[,3] + 2*(2*seglen+1)*cptsumstat[,1]*(cptsumstat[,5]-seglen*cptsumstat[,2]) + 6*cptsumstat[,4]*(cptsumstat[,2]-cptsumstat[,5]))
    betabottom=seglen*(seglen-1)*cptsumstat[,7] + 2*(2*seglen+1)*cptsumstat[,2]*(seglen*cptsumstat[,2]-cptsumstat[,5]) + 6*cptsumstat[,5]*(cptsumstat[,5]-cptsumstat[,2]);
    beta=betatop/betabottom;  
    thetajpo=(6*(seglen+2)*(cptsumstat[,4]-beta*cptsumstat[,5]))/((seglen+1)*(2*seglen+1)) - 2*(cptsumstat[,1]-beta*cptsumstat[,2])
    thetaj=(2*(2*seglen+1)*(cptsumstat[,1]-beta*cptsumstat[,2])-6*(cptsumstat[,4]-beta*cptsumstat[,5]))/(seglen-1)
    
    return(cbind(beta,thetajpo,thetaj))
  }
  
  if(cpttype(object)=="mean"){
    param.est<-list(mean=param.mean(object,cpts))
  }
  else if(cpttype(object)=="variance"){
    param.est<-list(variance=param.var(object,cpts))
  }
  else if(cpttype(object)=="mean and variance"){
    if(test.stat(object)=="Normal"){
      param.est<-list(mean=param.mean(object,cpts),variance=param.var(object,cpts))
    }
    else if(test.stat(object)=="Gamma"){
      param.est<-list(scale=param.scale(object,cpts,shape=shape),shape=shape)
    }
    else if(test.stat(object)=="Exponential"){
      param.est<-list(rate=1/param.mean(object,cpts))
    }
    else if(test.stat(object)=="Poisson"){
      param.est<-list(lambda=param.mean(object,cpts))
    }
    else{
      stop("Unknown test statistic for a change in mean and variance")
    }
  }
  else if(cpttype(object)=="trend"){
    if(test.stat(object)=="Normal"){
      tmp=param.trend(object)
      param.est(object)<-list(thetaS=tmp[,1],thetaT=tmp[,2])
    }
    else{
      stop("Unknown test statistic for a change in trend")
    }
  }
  else if(cpttype(object)=="trendar"){
    if(test.stat(object)=="Normal"){
      tmp=param.trendar(object)
      param.est(object)<-list(beta=tmp[,1],thetajpo=tmp[,2],thetaj=tmp[,3])
    }
    else{
      stop("Unknown test statistic for a change in trend+ar")
    }
  }
  else if(cpttype(object)=="meanar"){
    if(test.stat(object)=="Normal"){
      tmp=param.meanar(object)
      param.est(object)<-list(beta1=tmp[,1],beta2=tmp[,2])
    }
    else{
      stop("Unknown test statistic for a change in mean+ar")
    }
  }
  else{
    stop("Unknown changepoint type, must be 'mean', 'variance', 'mean and variance', 'trend', 'meanar' or 'trendar'")
  }
  if(is.na(ncpts)){
    param.est(object)=param.est
    return(object)
  }
  out=new('cpt.range')
  param.est(out)=param.est
  return(out)
})

setMethod("param", "cpt.reg", function(object,shape,...) {
  param.norm=function(object){
    cpts=c(0,object@cpts)
    #	nseg=length(cpts)-1 #nseg(object)
    data=data.set(object)
    p=ncol(data)-1
    tmpbeta=matrix(NA,ncol=p,nrow=nseg(object))
    tmpsigma=rep(NA,nseg(object))
    for(j in 1:nseg(object)){
      formula=paste('-1+data[',cpts[j]+1,':',cpts[j+1],',2]',sep='')
      if(p>1){
        for(i in 2:p){
          formula=paste(formula,'+data[',(cpts[j]+1),':',cpts[j+1],',',i+1,']',sep='')
        }
      }
      tmpfit=eval(parse(text=paste('lm(data[',(cpts[j]+1),':',cpts[j+1],',1]~',formula,')',sep='')))
      tmpbeta[j,]=tmpfit$coefficients
      tmpsigma[j]=sum(tmpfit$residuals^2)/(length(tmpfit$residuals)-length(tmpfit$coefficients)) ##var(tmpfit$residuals)
    }
    return(list(beta=tmpbeta,sig2=tmpsigma))
  }
  if(test.stat(object)=="Normal"){
    param.est(object)<-param.norm(object)
  }
  else{
    stop("Unknown test statistic, must be 'Normal'")
  }
  return(object)
})

# summary functions
setMethod("summary","cpt",function(object){
  cat("Created Using changepoint version",object@version,'\n')
  cat("Changepoint type      : Change in",cpttype(object),'\n')
  cat("Method of analysis    :",method(object),"\n")
  cat("Test Statistic  :", test.stat(object),"\n")
  cat("Type of penalty       :", pen.type(object), "with value,",pen.value(object),"\n")
  cat("Minimum Segment Length :", minseglen(object),"\n")
  cat("Maximum no. of cpts   :", ncpts.max(object),"\n")
  if(length(cpts(object))<=20){cat("Changepoint Locations :",cpts(object),"\n")}
  else{cat("Number of changepoints:", ncpts(object),"\n")}
})

setMethod("summary","cpt.range",function(object){
  cat("Created Using changepoint version",object@version,'\n')
  cat("Changepoint type      : Change in",cpttype(object),'\n')
  cat("Method of analysis    :",method(object),"\n")
  cat("Test Statistic  :", test.stat(object),"\n")
  cat("Type of penalty       :", pen.type(object), "with value,",pen.value(object),"\n")
  cat("Minimum Segment Length :", minseglen(object),"\n")
  cat("Maximum no. of cpts   :", ncpts.max(object),"\n")
  if(length(cpts(object))<=20){cat("Changepoint Locations :",cpts(object),"\n")}
  else{cat("Number of changepoints:", ncpts(object),"\n")}
  if((nrow(cpts.full(object))<=5)&(ncol(cpts.full(object)<=20))){cat("Range of segmentations:\n");print(cpts.full(object));cat("\n For penalty values:", pen.value.full(object),"\n")}
  else{cat("Number of segmentations recorded:", nrow(cpts.full(object)), " with between ", sum(cpts.full(object)[nrow(cpts.full(object)),]>0,na.rm=T), " and ", sum(cpts.full(object)[1,]>0,na.rm=T), "changepoints.\n Penalty value ranges from:",min(pen.value.full(object))," to ",max(pen.value.full(object)))}
})

setMethod("summary","cpt.reg",function(object){
  cat("Created Using changepoint version",object@version,'\n')
  cat("Changepoint type     : Change in",cpttype(object),'\n')
  cat("Method of analysis   :",method(object),"\n")
  cat("Test Statistic :", test.stat(object),"\n")
  cat("Type of penalty      :", pen.type(object), "with value,",pen.value(object),"\n")
  cat("Maximum no. of cpts   :", ncpts.max(object),"\n")
  if(length(cpts(object))<=20){cat("Changepoint Locations :",cpts(object),"\n")}
  else{cat("Number of changepoints:", ncpts(object),"\n")}
})

# show functions
setMethod("show","cpt",function(object){
  cat("Class 'cpt' : Changepoint Object\n")
  cat("       ~~   : S4 class containing", length(attributes(object))-1, "slots with names\n")
  cat("             ", names(attributes(object))[1:(length(attributes(object))-1)], "\n\n")
  cat("Created on  :", object@date, "\n\n")
  cat("summary(.)  :\n----------\n")
  summary(object)
})
setMethod("show","cpt.reg",function(object){
  cat("Class 'cpt.reg' : Changepoint Regression Object\n")
  cat("       ~~   : S4 class containing", length(attributes(object))-1, "slots with names\n")
  cat("             ", names(attributes(object))[1:(length(attributes(object))-1)], "\n\n")
  cat("Created on  :", object@date, "\n\n")
  cat("summary(.)  :\n----------\n")
  summary(object)
})

# plot functions
setMethod("plot","cpt",function(x,cpt.col='red',cpt.width=1,cpt.style=1,...){
  if(length(param.est(x))==0){# i.e. parameter.estimates=FALSE in call
    cat('Calculating parameter estimates...')
    object=param(x)
    cat('done.\n')
  }
  plot(data.set.ts(x),...)
  if(cpttype(x)=="variance" || cpttype(x)=="nonparametric (empirical_distribution)"){
    abline(v=index(data.set.ts(x))[cpts(x)],col=cpt.col,lwd=cpt.width,lty=cpt.style)
  }
  else if(cpttype(x)=="mean"  ||  cpttype(x)=="mean and variance"){
    #nseg=length(cpts(x))+1
    cpts=c(0,x@cpts)
    if((test.stat(x)=="Normal")||(test.stat(x)=="CUSUM")){
      means=param.est(x)$mean
    }
    else if(test.stat(x)=="Gamma"){
      means=param.est(x)$scale*param.est(x)$shape
    }
    else if(test.stat(x)=="Exponential"){
      means=1/param.est(x)$rate
    }
    else if(test.stat(x)=="Poisson"){
      means=param.est(x)$lambda
    }
    else{
      stop('Invalid Changepoint test statistic')
    }
    for(i in 1:nseg(x)){
      segments(index(data.set.ts(x))[cpts[i]+1],means[i],index(data.set.ts(x))[cpts[i+1]],means[i],col=cpt.col,lwd=cpt.width,lty=cpt.style)
    }
  }
  else if(cpttype(x)=="trend"){
    cpts=c(0,x@cpts)
    intercept=rep(param.est(x)$thetaS,x@cpts-c(0,cpts(x)))
    slope=rep(param.est(x)$thetaT-param.est(x)$thetaS,x@cpts-c(0,cpts(x)))/rep(x@cpts-c(0,cpts(x)),x@cpts-c(0,cpts(x)))
    cptn=rep(c(0,cpts(x)),x@cpts-c(0,cpts(x)))
    n=length(data.set(x))
    means=intercept+slope*((1:n)-cptn)
    for(i in 1:nseg(x)){
      segments(index(data.set.ts(x))[cpts[i]+1],means[cpts[i]+1],index(data.set.ts(x))[cpts[i+1]],means[cpts[i+1]],col=cpt.col,lwd=cpt.width,lty=cpt.style)
    }
  }
  else{
    stop('Invalid Changepoint Type for plotting.\n Can only plot mean, variance, mean and variance')
  }
})

setMethod("plot","cpt.range",function(x,ncpts=NA,diagnostic=FALSE,cpt.col='red',cpt.width=1,cpt.style=1,...){
  if(diagnostic==TRUE){
    n.changepoints = apply(cpts.full(x), 1, function(x) sum(x > 0, na.rm = TRUE))
    penalty.values = pen.value.full(x)
    if (is.null(list(...)$type)) {
      # By default, the type of the diagnostic plots is "lines".
      plot(x = n.changepoints, y = penalty.values, xlab = 'Number of Changepoints', ylab = 'Penalty Value', type = "l", ...)
    } else {
      plot(x = n.changepoints, y = penalty.values, xlab = 'Number of Changepoints', ylab = 'Penalty Value', ...)
    }
    return(invisible(NULL))
  }
  plot(data.set.ts(x),...)
  if(is.na(ncpts)){
    if(pen.type(x)=="CROPS"){
      stop('CROPS does not supply an optimal set of changepoints, set ncpts to the desired segmentation to plot or use diagnostic=TRUE to identify an appropriate number of changepoints')
    }
    cpts.to.plot=cpts(x)
    param.est=x
  }
  else{
    ncpts.full=apply(cpts.full(x),1,function(x){sum(x>0,na.rm=TRUE)})
    row=which(ncpts.full==ncpts)
    if(length(row)==0){
      stop(paste("Your input object doesn't have a segmentation with the requested number of changepoints.\n Possible ncpts are: "),paste(ncpts.full,collapse=','))
    }
    cpts.to.plot=cpts.full(x)[row,1:ncpts]
    if(test.stat(x)=="Gamma"){
      param.est=param(x,ncpts,shape=param.est(x)$shape)
    }
    else{
      param.est=param(x,ncpts)
    }
  }
  if(cpttype(x)=="variance"){
    abline(v=index(data.set.ts(x))[cpts.to.plot],col=cpt.col,lwd=cpt.width,lty=cpt.style)
  }
  else if(cpttype(x)=="mean"  ||  cpttype(x)=="mean and variance"){
    if((test.stat(x)=="Normal")||(test.stat(x)=="CUSUM")){
      means=param.est(param.est)$mean
    }
    else if(test.stat(x)=="Gamma"){
      means=param.est(param.est)$scale*param.est(param.est)$shape
    }
    else if(test.stat(x)=="Exponential"){
      means=1/param.est(param.est)$rate
    }
    else if(test.stat(x)=="Poisson"){
      means=param.est(param.est)$lambda
    }
    else{
      stop('Invalid Changepoint test statistic')
    }
    nseg=length(means)
    cpts.to.plot=c(0,cpts.to.plot,length(data.set(x)))
    for(i in 1:nseg){
      segments(index(data.set.ts(x))[cpts.to.plot[i]+1],means[i],index(data.set.ts(x))[cpts.to.plot[i+1]],means[i],col=cpt.col,lwd=cpt.width,lty=cpt.style)
    }
  }
  else{
    stop('Invalid Changepoint Type for plotting.\n Can only plot mean, variance, mean and variance')
  }
})

setMethod("plot","cpt.reg",function(x,cpt.col='red',cpt.width=1,cpt.style=1,...){
  if(length(param.est(x))==0){# i.e. parameter.estimates=FALSE in call
    cat('Calculating parameter estimates...')
    object=param(x)
    cat('done.\n')
  }
  plot(data.set(x)[,1],type='l',...)
  if(test.stat(x)=="Normal"){
    cpts=c(0,x@cpts)
    betas=param.est(x)$beta
    for(i in 1:nseg(x)){
      lines((cpts[i]+1):cpts[i+1],betas[i,]%*%t(data.set(x)[(cpts[i]+1):cpts[i+1],-1]),col=cpt.col,lwd=cpt.width,lty=cpt.style)
    }
  }
  else{
    stop('Invalid Changepoint test statistic')
  }
})

# likelihood functions
setMethod("logLik", "cpt", function(object) {
  if(length(param.est(object))==0){# i.e. parameter.estimates=FALSE in call
    cat('Calculating parameter estimates...')
    object=param(object)
    cat('done.\n')
  }
  if(test.stat(object)=="Normal"){
    if(cpttype(object)=="mean"){
      means=rep(param.est(object)$mean,object@cpts-c(0,cpts(object)))
      rss=sum((data.set(object)-means)^2)
      n=length(data.set(object))
      like=n*(log(2*pi)+log(rss/n)+1) # -2*loglik
      cpts=c(0,object@cpts)
      if(pen.type(object)=="MBIC"){
        like=c(like, like+(nseg(object)-2)*pen.value(object)+sum(log(cpts[-1]-cpts[-(nseg(object)+1)])))
      }else{
        like=c(like,like+(nseg(object)-1)*pen.value(object))
      }
    }
    else if(cpttype(object)=="variance"){
      rss=c(0,cumsum((data.set(object)-param.est(object)$mean)^2))
      cpts=c(0,object@cpts)
      n=length(data.set(object))
      seglen=seg.len(object)
      sigmas=(rss[cpts[-1]+1]-rss[cpts[-length(cpts)]+1])/seglen
      like=n*log(2*pi)+sum(seglen*log(sigmas))+n
      if(pen.type(object)=="MBIC"){
        like=c(like, like+(nseg(object)-2)*pen.value(object)+sum(log(seglen)))
      }else{
        like=c(like,like+(nseg(object)-1)*pen.value(object))
      }
    }
    else if(cpttype(object)=="mean and variance"){
      means=rep(param.est(object)$mean,object@cpts-c(0,cpts(object)))
      rss=sum((data.set(object)-means)^2)
      n=length(data.set(object))
      cpts=c(0,object@cpts)
      seglen=seg.len(object)
      sigmas=param.est(object)$variance
      like=n*log(2*pi)+sum(seglen*log(sigmas))+n
      if(pen.type(object)=="MBIC"){
        like=c(like,like+(nseg(object)-2)*pen.value(object)+sum(log(seglen)))
      }else{
        like=c(like,like+(nseg(object)-1)*pen.value(object))
      }
    }
    else if(cpttype(object)=="trend"){
      intercept=rep(param.est(object)$thetaS,object@cpts-c(0,cpts(object)))
      slope=rep(param.est(object)$thetaT-param.est(object)$thetaS,object@cpts-c(0,cpts(object)))/rep(object@cpts-c(0,cpts(object)),object@cpts-c(0,cpts(object)))
      cptn=rep(c(0,cpts(object)),object@cpts-c(0,cpts(object)))
      n=length(data.set(object))
      means=intercept+slope*((1:n)-cptn)
      rss=sum((data.set(object)-means)^2)
      like=n*(log(2*pi)+log(rss/n)+1) # -2*loglik
      if(pen.type(object)=="MBIC"){
        like=c(like, like+(nseg(object)-2)*pen.value(object)+sum(log(seg.len(object))))
      }else{
        like=c(like,like+(nseg(object)-1)*pen.value(object))
      }
    }
    else if(cpttype(object)=="trendar"){
      seglen=seg.len(object)
      intercept=rep(param.est(object)$thetaj,seglen)
      slope=rep(param.est(object)$thetajpo-param.est(object)$thetaj,seglen)/rep(seglen,seglen)
      ar=rep(param.est(object)$beta,seglen)
      cptn=rep(c(0,cpts(object)),seglen)
      n=length(data.set(object))
      means=NULL;means[1]=0
      for(i in 2:n){means[i]=intercept+slope*((1:n)-cptn)+ar*means[i-1]}
      means=means[-1]
      rss=sum((data.set(object)[-1]-means)^2)
      like=n*(log(2*pi)+log(rss/n)+1) # -2*loglik
      if(pen.type(object)=="MBIC"){
        like=c(like, like+(nseg(object)-2)*pen.value(object)+sum(log(seg.len(object))))
      }else{
        like=c(like,like+(nseg(object)-1)*pen.value(object))
      }
    }
    else if(cpttype(object)=="meanar"){
      seglen=seg.len(object)
      intercept=rep(param.est(object)$beta1,seglen)
      ar=rep(param.est(object)$beta2,seglen)
      cptn=rep(c(0,cpts(object)),seglen)
      n=length(data.set(object))
      means[1]=0;for(i in 2:n){means[i]=intercept+ar*means[i-1]}
      means=means[-1]
      rss=sum((data.set(object)[-1]-means)^2)
      like=n*(log(2*pi)+log(rss/n)+1) # -2*loglik
      if(pen.type(object)=="MBIC"){
        like=c(like, like+(nseg(object)-2)*pen.value(object)+sum(log(seg.len(object))))
      }else{
        like=c(like,like+(nseg(object)-1)*pen.value(object))
      }
    }
    else{
      stop("Unknown changepoint type, must be 'mean', 'variance', 'mean and variance', 'trend', 'meanar' or 'trendar'")
    }
  }
  else if(test.stat(object)=="Gamma"){
    if(cpttype(object)!="mean and variance"){
      stop("Unknown changepoint type for test.stat='Gamma', must be 'mean and variance'")
    }
    else{
      warning("Not changed to be -2*logLik")
      mll.meanvarg=function(x,n,shape){
        return(n*shape*log(n*shape)-n*shape*log(x))
      }
      y=c(0,cumsum(data.set(object)))
      shape=param.est(object)$shape
      cpts=c(0,object@cpts)
      #nseg=length(cpts)-1
      tmplike=0
      for(j in 1:nseg(object)){
        tmplike=tmplike+mll.meanvarg(y[cpts[j+1]+1]-y[cpts[j]+1],cpts[j+1]-cpts[j],shape)
      }
      if(pen.type(object)=="MBIC"){
        like=c(tmplike, tmplike+(nseg(object)-2)*pen.value(object)+sum(log(seg.len(object))))
      }else{
        like=c(tmplike,tmplike+(nseg(object)-1)*pen.value(object))
      }
    }
  }
  else if(test.stat(object)=="Exponential"){
    if(cpttype(object)!="mean and variance"){
      stop("Unknown changepoint type for test.stat='Exponential', must be 'mean and variance'")
    }
    else{
      warning("Not changed to be -2*logLik")
      mll.meanvare=function(x,n){
        return(n*log(n)-n*log(x))
      }
      y=c(0,cumsum(data.set(object)))
      cpts=c(0,object@cpts)
      #nseg=length(cpts)-1
      tmplike=0
      for(j in 1:nseg(object)){
        tmplike=tmplike+mll.meanvare(y[cpts[j+1]+1]-y[cpts[j]+1],cpts[j+1]-cpts[j])
      }
      if(pen.type(object)=="MBIC"){
        like=c(tmplike, tmplike+(nseg(object)-2)*pen.value(object)+sum(log(seg.len(object))))
      }else{
        like=c(tmplike,tmplike+(nseg(object)-1)*pen.value(object))
      }
    }
  }
  else if(test.stat(object)=="Poisson"){
    if(cpttype(object)!="mean and variance"){
      stop("Unknown changepoint type for test.stat='Poisson', must be 'mean and variance'")
    }
    else{
      warning("Not changed to be -2*logLik")
      mll.meanvarp=function(x,n){
        return(x*log(x)-x*log(n))
      }
      y=c(0,cumsum(data.set(object)))
      cpts=c(0,object@cpts)
      #nseg=length(cpts)-1
      tmplike=0
      for(j in 1:nseg(object)){
        tmplike=tmplike+mll.meanvarp(y[cpts[j+1]+1]-y[cpts[j]+1],cpts[j+1]-cpts[j])
      }
      if(pen.type(object)=="MBIC"){
        like=c(tmplike, tmplike+(nseg(object)-2)*pen.value(object)+sum(log(seg.len(object))))
      }else{
        like=c(tmplike,tmplike+(nseg(object)-1)*pen.value(object))
      }
    }
  }
  else{stop("logLik is only valid for distributional assumptions, not CUSUM or CSS")}
  names(like)=c("-2*logLik","-2*Loglike+pen")
  return(like)
})

setMethod("logLik", "cpt.range", function(object,ncpts=NA) {
  # warning("Not changed to be -2*logLik")
  if(is.na(ncpts)){
    if(pen.type(object)=="CROPS"){
      stop('CROPS does not supply an optimal set of changepoints, set ncpts argument to the desired segmentation to plot or use diagnostic=TRUE to identify an appropriate number of changepoints')
    }
    cpts=c(0,object@cpts)
    pen.value=pen.value(object)
  }
  else{
    ncpts.full=apply(cpts.full(object),1,function(x){sum(x>0,na.rm=TRUE)})
    row=which(ncpts.full==ncpts)
    if(length(row)==0){
      stop(paste("Your input object doesn't have a segmentation with the requested number of changepoints.\n Possible ncpts are: "),paste(ncpts.full,collapse=','))
    }
    cpts=c(0,cpts.full(object)[row,1:ncpts],length(data.set(object)))
    pen.value=pen.value.full(object)[row]
  }
  nseg=length(cpts)-1
  
  if(test.stat(object)=="Normal"){
    if(cpttype(object)=="mean"){
      mll.mean=function(x2,x,n){
        return( x2-(x^2)/n)
      }
      y2=c(0,cumsum(data.set(object)^2))
      y=c(0,cumsum(data.set(object)))
      tmplike=0
      for(j in 1:nseg){
        tmplike=tmplike+mll.mean(y2[cpts[j+1]+1]-y2[cpts[j]+1],y[cpts[j+1]+1]-y[cpts[j]+1],cpts[j+1]-cpts[j])
      }
      ##c(tmplike, tmplike+(nseg-2)*pen.value(object)+sum(log(cpts[-1]-cpts[-(nseg+1)])))
      if(pen.type(object)=="MBIC"){
        like=c(tmplike, tmplike+(nseg-2)*pen.value+sum(log(seg.len(object))))
      }else{
        like=c(tmplike,tmplike+(nseg-1)*pen.value)
      }
      names(like)=c("-like","-likepen")
    }
    else if(cpttype(object)=="variance"){
      mll.var=function(x,n){
        neg=x<=0
        x[neg==TRUE]=0.00000000001    
        return( n*(log(2*pi)+log(x/n)+1))
      }
      y2=c(0,cumsum((data.set(object)-param.est(object)$mean)^2))
      tmplike=0
      for(j in 1:nseg){
        tmplike=tmplike+mll.var(y2[cpts[j+1]+1]-y2[cpts[j]+1],cpts[j+1]-cpts[j])
      }
      if(pen.type(object)=="MBIC"){
        like=c(tmplike, tmplike+(nseg-2)*pen.value+sum(log(seg.len(object))))
      }else{
        like=c(tmplike,tmplike+(nseg-1)*pen.value)
      }
      names(like)=c("-like","-likepen")
    }
    else if(cpttype(object)=="mean and variance"){
      mll.meanvar=function(x2,x,n){
        sigmasq=(1/n)*(x2-(x^2)/n)
        neg=sigmasq<=0
        sigmasq[neg==TRUE]=0.00000000001
        return( n*(log(2*pi)+log(sigmasq)+1))
      }
      y2=c(0,cumsum(data.set(object)^2))
      y=c(0,cumsum(data.set(object)))
      tmplike=0
      for(j in 1:nseg){
        tmplike=tmplike+mll.meanvar(y2[cpts[j+1]+1]-y2[cpts[j]+1],y[cpts[j+1]+1]-y[cpts[j]+1],cpts[j+1]-cpts[j])
      }
      if(pen.type(object)=="MBIC"){
        like=c(tmplike, tmplike+(nseg-2)*pen.value+sum(log(seg.len(object))))
      }else{
        like=c(tmplike,tmplike+(nseg-1)*pen.value)
      }
      names(like)=c("-like","-likepen")
    }
    else{
      stop("Unknown changepoint type, must be 'mean', 'variance' or 'mean and variance'")
    }
  }
  else if(test.stat(object)=="Gamma"){
    if(cpttype(object)!="mean and variance"){
      stop("Unknown changepoint type for test.stat='Gamma', must be 'mean and variance'")
    }
    else{
      mll.meanvarg=function(x,n,shape){
        return(n*shape*log(n*shape)-n*shape*log(x))
      }
      y=c(0,cumsum(data.set(object)))
      shape=param.est(object)$shape
      cpts=c(0,object@cpts)
      #nseg=length(cpts)-1
      tmplike=0
      for(j in 1:nseg){
        tmplike=tmplike+mll.meanvarg(y[cpts[j+1]+1]-y[cpts[j]+1],cpts[j+1]-cpts[j],shape)
      }
      if(pen.type(object)=="MBIC"){
        like=c(tmplike, tmplike+(nseg-2)*pen.value+sum(log(seg.len(object))))
      }else{
        like=c(tmplike,tmplike+(nseg-1)*pen.value)
      }
      names(like)=c("-like","-likepen")
    }
  }
  else if(test.stat(object)=="Exponential"){
    if(cpttype(object)!="mean and variance"){
      stop("Unknown changepoint type for test.stat='Exponential', must be 'mean and variance'")
    }
    else{
      mll.meanvare=function(x,n){
        return(n*log(n)-n*log(x))
      }
      y=c(0,cumsum(data.set(object)))
      cpts=c(0,object@cpts)
      #nseg=length(cpts)-1
      tmplike=0
      for(j in 1:nseg){
        tmplike=tmplike+mll.meanvare(y[cpts[j+1]+1]-y[cpts[j]+1],cpts[j+1]-cpts[j])
      }
      if(pen.type(object)=="MBIC"){
        like=c(tmplike, tmplike+(nseg-2)*pen.value+sum(log(seg.len(object))))
      }else{
        like=c(tmplike,tmplike+(nseg-1)*pen.value)
      }
      names(like)=c("-like","-likepen")
    }
  }
  else if(test.stat(object)=="Poisson"){
    if(cpttype(object)!="mean and variance"){
      stop("Unknown changepoint type for test.stat='Poisson', must be 'mean and variance'")
    }
    else{
      mll.meanvarp=function(x,n){
        return(x*log(x)-x*log(n))
      }
      y=c(0,cumsum(data.set(object)))
      cpts=c(0,object@cpts)
      #nseg=length(cpts)-1
      tmplike=0
      for(j in 1:nseg){
        tmplike=tmplike+mll.meanvarp(y[cpts[j+1]+1]-y[cpts[j]+1],cpts[j+1]-cpts[j])
      }
      if(pen.type(object)=="MBIC"){
        like=c(tmplike, tmplike+(nseg-2)*pen.value+sum(log(seg.len(object))))
      }else{
        like=c(tmplike,tmplike+(nseg-1)*pen.value)
      }
      names(like)=c("-like","-likepen")
    }
  }
  else{stop("logLik is only valid for distributional assumptions, not CUSUM or CSS")}
  return(like)
})

setMethod("logLik", "cpt.reg", function(object) {
  if(length(param.est(object))==0){# i.e. parameter.estimates=FALSE in call
    cat('Calculating parameter estimates...')
    object=param(object)
    cat('done.\n')
  }
  if(test.stat(object)=="Normal"){
    cpts=c(0,object@cpts)
    seglen=seg.len(object)
    data=data.set(object)
    beta=param.est(object)$beta
    sigmas=param.est(object)$sig2
    rss=NULL
    for(i in 1:length(seglen)){
      rss[i]=sum((data[(cpts[i]+1):cpts[i+1],1]-data[(cpts[i]+1):cpts[i+1],-1]%*%beta[i,])^2)
    }
    like=sum(seglen*log(2*pi*sigmas))+sum(rss/sigmas)
    if(pen.type(object)=="MBIC"){
      like=c(like, like+(nseg(object)-2)*pen.value(object)+sum(log(seg.len(object))))
    }else{
      like=c(like,like+(nseg(object)-1)*pen.value(object))
    }
  }
  else{stop("logLik is only valid for Normal distributional assumption.")}
  return(like)
})

setGeneric("likelihood", function(object) standardGeneric("likelihood"))
setMethod("likelihood", "cpt", function(object) {
  return(logLik(object))
})

# acf functions
# 	setGeneric("acf", function(object,...) standardGeneric("acf"))
# 	setMethod("acf", "cpt", function(object,lag.max=NULL,...) {
#     cpts=c(0,object@cpts)
#     nseg=nseg(object)
#     data=data.set(object)
#     for(i in 1:nseg){
#       stats::acf(data[(cpts[i]+1):cpts[i+1]],main=paste("Series part:",(cpts[i]+1),":",cpts[i+1]),...)
#     }
# 	})

# setMethod("acf", "cpt.reg", function(object,lag.max=NULL,...) {
#   cpts=c(0,object@cpts)
#   nseg=nseg(object)
#   data=data.set(object)[,1]
#   for(i in 1:nseg){
#     stats::acf(data[(cpts[i]+1):cpts[i+1]],main=paste("Series part:",(cpts[i]+1),"-",cpts[i+1]),...)
#   }
# })
# 

data_input <- function(data, method, pen.value, costfunc, minseglen, Q, var=0, shape=1){
  if(var !=0){
    mu<-var
  }else{
    mu <- mean(data)
  }
  sumstat=cbind(c(0,cumsum(coredata(data))),c(0,cumsum(coredata(data)^2)),cumsum(c(0,(coredata(data)-mu)^2)))
  if(method=="PELT"){
    #out=PELT.meanvar.norm(coredata(data),pen.value)
    out=PELT(sumstat,pen=pen.value,cost_func = costfunc,minseglen=minseglen, shape=shape)  ## K NEW ##
    #cpts=out[[2]]
  }
  else if(method=="BinSeg"){
    out=BINSEG(sumstat,pen=pen.value,cost_func = costfunc,minseglen=minseglen,Q=Q, shape=shape)  ## K NEW ##
    #cpts=out[[2]]
    #   		out=binseg.meanvar.norm(coredata(data),Q,pen.value)
    # 			if(out$op.cpts==0){cpts=n}
    # 			else{cpts=c(sort(out$cps[1,1:out$op.cpts]),n)}
    # the above is now inside the BINSEG function
  }
  else {
    stop('Unknown method, should be either PELT or BinSeg.')
  }
  return(out)
  
}
decision<-function(tau,null,alt=NA,penalty="MBIC",n=0,diffparam=1,pen.value=0){
  if(sum(is.na(alt))){
    if(length(tau)!=length(null)){
      stop("Lengths of tau and null do not match")
    }
  }
  else{
    if((length(tau)!=length(null))||(length(tau)!=length(alt))){
      stop("Lengths of tau, null and alt do not match")
    }
  }
  single.decision=function(tau,null,alt,n=0,diffparam=1,pen.value=0){
    if(is.na(alt)){teststat=null}
    else{teststat=null-alt}
    if(teststat>=pen.value){return(tau)}
    else{return(n)}
  }
  if(length(tau)==1){
    out=single.decision(tau,null,alt,n,diffparam,pen.value)
    names(out)="cpt"
    return(list(cpt=out,pen=pen.value))
  }
  else{
    rep=length(tau)
    out=NULL
    for(i in 1:rep){
      out[i]=single.decision(tau[i],null[i],alt[i],n,diffparam,pen.value)
    }
    names(out)=rep("cpt",rep)
    return(list(cpt=out,pen=pen.value))
  }
}
single.meanvar.exp.calc = function(data,extrainf=TRUE,minseglen){
    singledim=function(data,extrainf=TRUE,minseglen){
      n=length(data)
      y=c(0,cumsum(data))
      null=2*n*log(y[n+1])-2*n*log(n)
      taustar=minseglen:(n-minseglen)
      tmp=2*taustar*log(y[taustar+1]) -2*taustar*log(taustar) + 2*(n-taustar)*log((y[n+1]-y[taustar+1]))-2*(n-taustar)*log(n-taustar)
      
      tau=which(tmp==min(tmp,na.rm=T))[1]
      taulike=tmp[tau]
      tau=tau+minseglen-1 # correcting for the fact that we are starting at minseglen
      if(extrainf==TRUE){
        out=c(tau,null,taulike)
        names(out)=c('cpt','null','alt')
        return(out)
      }
      else{
        return(tau)
      }
    }
    
    if(is.null(dim(data))==TRUE){
      # single data set
      cpt=singledim(data,extrainf,minseglen)
      return(cpt)
    }
    else{
      rep=nrow(data)
      n=ncol(data)
      cpt=NULL
      if(extrainf==FALSE){
        for(i in 1:rep){
          cpt[i]=singledim(data[i,],extrainf,minseglen)
        }
      }
      else{
        cpt=matrix(0,ncol=3,nrow=rep)
        for(i in 1:rep){
          cpt[i,]=singledim(data[i,],extrainf,minseglen)
        }
        colnames(cpt)=c('cpt','null','alt')
      }
      return(cpt)
    }
  }
single.meanvar.exp<-function(data,penalty="MBIC",pen.value=0,class=TRUE,param.estimates=TRUE,minseglen){
  if(sum(data<0)>0){stop('Exponential test statistic requires positive data')}
  if(is.null(dim(data))==TRUE){
    # single dataset
    n=length(data)
  }
  else{
    n=ncol(data)
  }
  if(n<4){stop('Data must have atleast 4 observations to fit a changepoint model.')}
  if(n<(2*minseglen)){stop('Minimum segment legnth is too large to include a change in this data')}
  
  penalty_decision(penalty, pen.value, n, diffparam=1, asymcheck="meanvar.exp", method="AMOC")
  if(is.null(dim(data))==TRUE){
    tmp=single.meanvar.exp.calc(coredata(data),extrainf=TRUE,minseglen)
    if(penalty=="MBIC"){
      tmp[3]=tmp[3]+log(tmp[1])+log(n-tmp[1]+1)
    }
    ans=decision(tmp[1],tmp[2],tmp[3],penalty,n,diffparam=1,pen.value)
    if(class==TRUE){
      return(class_input(data, cpttype="mean and variance", method="AMOC", test.stat="Exponential", penalty=penalty, pen.value=ans$pen, minseglen=minseglen, param.estimates=param.estimates, out=c(0,ans$cpt)))
    }
    else{ 
      an=(2*log(log(n)))^(1/2)
      bn=2*log(log(n))+(1/2)*log(log(log(n)))-(1/2)*log(pi)
      out=c(ans$cpt,exp(-2*exp(-an*sqrt(abs(tmp[2]-tmp[3]))+an*bn))-exp(-2*exp(an*bn)))  # Chen & Gupta (2000) pg149
      names(out)=c('cpt','p value')
      return(out)
    }
  }
  else{ 
    tmp=single.meanvar.exp.calc(data,extrainf=TRUE,minseglen)
    if(penalty=="MBIC"){
      tmp[,3]=tmp[,3]+log(tmp[,1])+log(n-tmp[,1]+1)
    }
    ans=decision(tmp[,1],tmp[,2],tmp[,3],penalty,n,diffparam=1,pen.value)
    if(class==TRUE){
      rep=nrow(data)
      out=list()
      for(i in 1:rep){
        out[[i]]=class_input(data[i,], cpttype="mean and variance", method="AMOC", test.stat="Exponential", penalty=penalty, pen.value=ans$pen, minseglen=minseglen, param.estimates=param.estimates, out=c(0,ans$cpt[i]))
      }
      return(out)
    }
    else{ 
      an=(2*log(log(n)))^(1/2)
      bn=2*log(log(n))+(1/2)*log(log(log(n)))-(1/2)*log(pi)
      out=cbind(ans$cpt,exp(-2*exp(-an*sqrt(abs(tmp[,2]-tmp[,3]))+bn))-exp(-2*exp(bn)))  # Chen & Gupta (2000) pg149
      colnames(out)=c('cpt','p value')
      rownames(out)=NULL
      return(out)
    }
  }
}


# PELT.meanvar.exp=function(data,pen=0,nprune=FALSE){
#   mll.meanvar.EFK=function(x,n){
#     return( 2*n*log(x)-2*n*log(n))
#   }
#   n=length(data)
#   y=c(0,cumsum(data))
# 
#   lastchangecpts=matrix(NA,nrow=n,ncol=2)
#   lastchangelike=matrix(NA,nrow=n,ncol=2)
#   checklist=NULL
#   lastchangelike[1,]=c(mll.meanvar.EFK(y[2],1),mll.meanvar.EFK(y[n+1]-y[2],n-1)+pen)
#   lastchangecpts[1,]=c(0,1)
#   lastchangelike[2,]=c(mll.meanvar.EFK(y[3],2),mll.meanvar.EFK(y[n+1]-y[3],n-2)+pen)
#   lastchangecpts[2,]=c(0,2)
#   lastchangelike[3,]=c(mll.meanvar.EFK(y[4],3),mll.meanvar.EFK(y[n+1]-y[4],n-3)+pen)
#   lastchangecpts[3,]=c(0,3)
#   noprune=NULL
#   for(tstar in 4:n){
#     tmplike=NULL
#     tmpt=c(checklist, tstar-2)
#     tmplike=lastchangelike[tmpt,1]+mll.meanvar.EFK(y[tstar+1]-y[tmpt+1],tstar-tmpt)+pen
#     if(tstar==n){
#       lastchangelike[tstar,]=c(min(c(tmplike,mll.meanvar.EFK(y[tstar+1]-y[1],tstar)),na.rm=TRUE),0)
#     }
#     else{
#       lastchangelike[tstar,]=c(min(c(tmplike,mll.meanvar.EFK(y[tstar+1]-y[1],tstar)),na.rm=TRUE),mll.meanvar.EFK(y[n+1]-y[tstar+1],n-tstar)+pen)
#     }
#     if(lastchangelike[tstar,1]==mll.meanvar.EFK(y[tstar+1]-y[1],tstar)){
#       lastchangecpts[tstar,]=c(0,tstar)
#     }
#     else{
#       cpt=tmpt[tmplike==lastchangelike[tstar,1]][1]
#       lastchangecpts[tstar,]=c(cpt,tstar)
#     }
#     checklist=tmpt[tmplike<=lastchangelike[tstar,1]+pen]
#     if(nprune==TRUE){
#       noprune=c(noprune,length(checklist))
#     }
#   }
#   if(nprune==TRUE){
#     return(nprune=noprune)
#   }
#   else{
#     fcpt=NULL
#     last=n
#     while(last!=0){
# 	fcpt=c(fcpt,lastchangecpts[last,2])
# 	last=lastchangecpts[last,1]
#     }
#     return(cpt=sort(fcpt))
#   }
# }

# binseg.meanvar.exp=function(data,Q=5,pen=0){
#   mll.meanvar=function(x,n){
#     return(n*log(n)-n*log(x))
#   }
#   n=length(data)
#   y=c(0,cumsum(data))
#   tau=c(0,n)
#   cpt=matrix(0,nrow=2,ncol=Q)
#   oldmax=1000
# 
#   for(q in 1:Q){
#     lambda=rep(0,n-1)
#     i=1
#     st=tau[1]+1;end=tau[2]
#     null=mll.meanvar(y[end+1]-y[st],end-st+1)
#     for(j in 1:(n-1)){
#       if(j==end){
#         st=end+1;i=i+1;end=tau[i+1]
#         null=mll.meanvar(y[end+1]-y[st],end-st+1)
#       }else{
# 	if((j-st)<2){lambda[j]=-1*10^(100)}
# 	else if((end-j)<2){lambda[j]=-1*10^(100)}
# 	else{lambda[j]=mll.meanvar(y[j+1]-y[st],j-st+1)+mll.meanvar(y[end+1]-y[j+1],end-j)-null}
#       }
#     }
#     k=which.max(lambda)[1]
#     cpt[1,q]=k;cpt[2,q]=min(oldmax,max(lambda))
#     oldmax=min(oldmax,max(lambda))
#     tau=sort(c(tau,k))
#   }
#   op.cps=NULL
#   p=1:(Q-1)
#   for(i in 1:length(pen)){
#     criterion=(2*cpt[2,])>=pen[i]
#     if(sum(criterion)==0){
#       op.cps=0
#     }
#     else{
#       op.cps=c(op.cps,max(which((criterion)==TRUE)))
#     }
#   }
#   return(list(cps=cpt,op.cpts=op.cps,pen=pen))
# }

multiple.meanvar.exp=function(data,mul.method="PELT",penalty="MBIC",pen.value=0,Q=5,class=TRUE,param.estimates=TRUE,minseglen){
  if(sum(data<0)>0){stop('Exponential test statistic requires positive data')}
  if(!((mul.method=="PELT")||(mul.method=="BinSeg"))){
    stop("Multiple Method is not recognised, must be PELT or BinSeg.")
  }
  costfunc = "meanvar.exp"
  if(penalty=="MBIC"){
    costfunc = "meanvar.exp.mbic"
  }
  
  diffparam=1
  if(is.null(dim(data))==TRUE){
    # single dataset
    n=length(data)
  }
  else{
    n=ncol(data)
  }
  if(n<(2*minseglen)){stop('Minimum segment legnth is too large to include a change in this data')}
  
  pen.value = penalty_decision(penalty, pen.value, n, diffparam=1, asymcheck=costfunc, method=mul.method)
  
  if(is.null(dim(data))==TRUE){
    # single dataset
    out = data_input(data=data,method=mul.method,pen.value=pen.value,costfunc=costfunc,minseglen=minseglen,Q=Q)
    
    if(class==TRUE){
      return(class_input(data, cpttype="mean and variance", method=mul.method, test.stat="Exponential", penalty=penalty, pen.value=pen.value, minseglen=minseglen, param.estimates=param.estimates, out=out, Q=Q))
    }
    else{ return(out[[2]])}
  }
  else{
    rep=nrow(data)
    out=list()
    for(i in 1:rep){
      out[[i]]=data_input(data[i,],method=mul.method,pen.value=pen.value,costfunc=costfunc,minseglen=minseglen,Q=Q)
    }
    
    cpts=lapply(out, '[[', 2)
    
    if(class==TRUE){
      ans=list()
      for(i in 1:rep){
        ans[[i]]=class_input(data[i,], cpttype="mean and variance", method=mul.method, test.stat="Exponential", penalty=penalty, pen.value=pen.value, minseglen=minseglen, param.estimates=param.estimates, out=out[[i]], Q=Q)
      }
      return(ans)
    }
    else{return(cpts)}
  }
}
single.meanvar.gamma.calc = function(data,shape=1,extrainf=TRUE,minseglen){
    singledim=function(data,shape,extrainf=TRUE,minseglen){
      n=length(data)
      y=c(0,cumsum(data))
      null=2*n*shape*log(y[n+1])-2*n*shape*log(n*shape)
      taustar=minseglen:(n-minseglen)
      tmp=2*taustar*shape*log(y[taustar+1]) -2*taustar*shape*log(taustar*shape) + 2*(n-taustar)*shape*log((y[n+1]-y[taustar+1]))-2*(n-taustar)*shape*log((n-taustar)*shape)
      
      tau=which(tmp==min(tmp,na.rm=T))[1]
      taulike=tmp[tau]
      tau=tau+minseglen-1 # correcting for the fact that we are starting at minseglen
      if(extrainf==TRUE){
        out=c(tau,null,taulike)
        names(out)=c('cpt','null','alt')
        return(out)
      }
      else{
        return(tau)
      }
    }
    
    
    if(is.null(dim(data))==TRUE){
      # single data set
      cpt=singledim(data,shape,extrainf,minseglen)
      return(cpt)
    }
    else{
      rep=nrow(data)
      n=ncol(data)
      cpt=NULL
      if(length(shape)==1){
        shape=rep(shape,rep)
      }
      if(extrainf==FALSE){
        for(i in 1:rep){
          cpt[i]=singledim(data[i,],shape[i],extrainf,minseglen)
        }
      }
      else{
        cpt=matrix(0,ncol=3,nrow=rep)
        for(i in 1:rep){
          cpt[i,]=singledim(data[i,],shape[i],extrainf,minseglen)
        }
        colnames(cpt)=c('cpt','null','alt')
      }
      return(cpt)
    }
  }
single.meanvar.gamma<-function(data,shape=1,penalty="MBIC",pen.value=0,class=TRUE,param.estimates=TRUE,minseglen){
  if(sum(data<=0)>0){stop('Gamma test statistic requires positive data')}
  if(is.null(dim(data))==TRUE){
    # single dataset
    n=length(data)
  }
  else{
    n=ncol(data)
  }
  if(n<4){stop('Data must have atleast 4 observations to fit a changepoint model.')}
  if(n<(2*minseglen)){stop('Minimum segment legnth is too large to include a change in this data')}
  
  pen.value = penalty_decision(penalty, pen.value, n, diffparam=1, asymcheck="meanvar.gamma", method="AMOC")
  if(is.null(dim(data))==TRUE){
    tmp=single.meanvar.gamma.calc(coredata(data),shape,extrainf=TRUE,minseglen)
    if(penalty=="MBIC"){
      tmp[3]=tmp[3]+log(tmp[1])+log(n-tmp[1]+1)
    }
    ans=decision(tmp[1],tmp[2],tmp[3],penalty,n,diffparam=1,pen.value)
    if(class==TRUE){
      return(class_input(data, cpttype="mean and variance", method="AMOC", test.stat="Gamma", penalty=penalty, pen.value=ans$pen, minseglen=minseglen, param.estimates=param.estimates, out=c(0,ans$cpt), shape=shape))
    }
    else{ 
      return(ans$cpt)
    }
  }
  else{ 
    tmp=single.meanvar.gamma.calc(data,shape,extrainf=TRUE,minseglen)
    if(penalty=="MBIC"){
      tmp[,3]=tmp[,3]+log(tmp[,1])+log(n-tmp[,1]+1)
    }
    ans=decision(tmp[,1],tmp[,2],tmp[,3],penalty,n,diffparam=1,pen.value)
    if(class==TRUE){
      rep=nrow(data)
      out=list()
      for(i in 1:rep){
        out[[i]]=class_input(data[i,], cpttype="mean and variance", method="AMOC", test.stat="Gamma", penalty=penalty, pen.value=ans$pen, minseglen=minseglen, param.estimates=param.estimates, out=c(0,ans$cpt[i]), shape=shape)
      }
      return(out)
    }
    else{ return(ans$cpt)}
  }
}

# PELT.meanvar.gamma=function(data,shape=1,pen=0,nprune=FALSE){
#   mll.meanvar.EFK=function(x,n,shape){
#     return( 2*n*shape*log(x)-2*n*shape*log(n*shape))
#   }
#   n=length(data)
#   y=c(0,cumsum(data))
# 
#   lastchangecpts=matrix(NA,nrow=n,ncol=2)
#   lastchangelike=matrix(NA,nrow=n,ncol=2)
#   checklist=NULL
#   lastchangelike[1,]=c(mll.meanvar.EFK(y[2],1,shape),mll.meanvar.EFK(y[n+1]-y[2],n-1,shape)+pen)
#   lastchangecpts[1,]=c(0,1)
#   lastchangelike[2,]=c(mll.meanvar.EFK(y[3],2,shape),mll.meanvar.EFK(y[n+1]-y[3],n-2,shape)+pen)
#   lastchangecpts[2,]=c(0,2)
#   lastchangelike[3,]=c(mll.meanvar.EFK(y[4],3,shape),mll.meanvar.EFK(y[n+1]-y[4],n-3,shape)+pen)
#   lastchangecpts[3,]=c(0,3)
#   noprune=NULL
#   for(tstar in 4:n){
#     tmplike=NULL
#     tmpt=c(checklist, tstar-2)
#     tmplike=lastchangelike[tmpt,1]+mll.meanvar.EFK(y[tstar+1]-y[tmpt+1],tstar-tmpt,shape)+pen
#     if(tstar==n){
#       lastchangelike[tstar,]=c(min(c(tmplike,mll.meanvar.EFK(y[tstar+1]-y[1],tstar,shape)),na.rm=TRUE),0)
#     }
#     else{
#       lastchangelike[tstar,]=c(min(c(tmplike,mll.meanvar.EFK(y[tstar+1]-y[1],tstar,shape)),na.rm=TRUE),mll.meanvar.EFK(y[n+1]-y[tstar+1],n-tstar,shape)+pen)
#     }
#     if(lastchangelike[tstar,1]==mll.meanvar.EFK(y[tstar+1]-y[1],tstar,shape)){
#       lastchangecpts[tstar,]=c(0,tstar)
#     }
#     else{
#       cpt=tmpt[tmplike==lastchangelike[tstar,1]][1]
#       lastchangecpts[tstar,]=c(cpt,tstar)
#     }
#     checklist=tmpt[tmplike<=lastchangelike[tstar,1]+pen]
#     if(nprune==TRUE){
#       noprune=c(noprune,length(checklist))
#     }
#   }
#   if(nprune==TRUE){
#     return(nprune=noprune)
#   }
#   else{
#     fcpt=NULL
#     last=n
#     while(last!=0){
# 	fcpt=c(fcpt,lastchangecpts[last,2])
# 	last=lastchangecpts[last,1]
#     }
#     return(cpt=sort(fcpt))
#   }
# }


# binseg.meanvar.gamma=function(data,shape=1,Q=5,pen=0){
#   mll.meanvar=function(x,n,shape){
#     return(n*shape*log(n*shape)-n*shape*log(x))
#   }
#   n=length(data)
#   y=c(0,cumsum(data))
#   tau=c(0,n)
#   cpt=matrix(0,nrow=2,ncol=Q)
#   oldmax=1000
# 
#   for(q in 1:Q){
#     lambda=rep(0,n-1)
#     i=1
#     st=tau[1]+1;end=tau[2]
#     null=mll.meanvar(y[end+1]-y[st],end-st+1,shape)
#     for(j in 1:(n-1)){
#       if(j==end){
#         st=end+1;i=i+1;end=tau[i+1]
#         null=mll.meanvar(y[end+1]-y[st],end-st+1,shape)
#       }else{
# 	if((j-st)<2){lambda[j]=-1*10^(100)}
# 	else if((end-j)<2){lambda[j]=-1*10^(100)}
# 	else{lambda[j]=mll.meanvar(y[j+1]-y[st],j-st+1,shape)+mll.meanvar(y[end+1]-y[j+1],end-j,shape)-null}
#       }
#     }
#     k=which.max(lambda)[1]
#     cpt[1,q]=k;cpt[2,q]=min(oldmax,max(lambda))
#     oldmax=min(oldmax,max(lambda))
#     tau=sort(c(tau,k))
#   }
#   op.cps=NULL
#   p=1:(Q-1)
#   for(i in 1:length(pen)){
#     criterion=(2*cpt[2,])>=pen[i]
#     if(sum(criterion)==0){
#       op.cps=0
#     }
#     else{
#       op.cps=c(op.cps,max(which((criterion)==TRUE)))
#     }
#   }
#   return(list(cps=cpt,op.cpts=op.cps,pen=pen))
# }


multiple.meanvar.gamma=function(data,shape=1,mul.method="PELT",penalty="MBIC",pen.value=0,Q=5,class=TRUE,param.estimates=TRUE,minseglen){
  if(sum(data<=0)>0){stop('Gamma test statistic requires positive data')}
  if(!((mul.method=="PELT")||(mul.method=="BinSeg"))){
    stop("Multiple Method is not recognised, must be PELT or BinSeg.")
  }
  costfunc = "meanvar.gamma"
  if(penalty=="MBIC"){
    costfunc = "meanvar.gamma.mbic"
  }
  diffparam=1
  if(is.null(dim(data))==TRUE){
    # single dataset
    n=length(data)
    shape=shape[1]
  }
  else{
    n=ncol(data)
  }
  if(n<(2*minseglen)){stop('Minimum segment legnth is too large to include a change in this data')}
  
  pen.value = penalty_decision(penalty, pen.value, n, diffparam, asymcheck = costfunc, method=mul.method)
  
  if(is.null(dim(data))==TRUE){
    # single dataset
    out = data_input(data=data,method=mul.method,pen.value=pen.value,costfunc=costfunc,minseglen=minseglen,Q=Q,shape=shape)
    
    if(class==TRUE){
      return(class_input(data, cpttype="mean and variance", method=mul.method, test.stat="Gamma", penalty=penalty, pen.value=pen.value, minseglen=minseglen, param.estimates=param.estimates, out=out, Q=Q, shape=shape))
    }
    else{ return(out[[2]])}
  }
  else{
    rep=nrow(data)
    out=list()
    if(length(shape)!=rep){
      shape=rep(shape,rep)
    }
    for(i in 1:rep){
      out[[i]]=data_input(data[i,],method=mul.method,pen.value=pen.value,costfunc=costfunc,minseglen=minseglen,Q=Q, shape=shape)
    }
    
    cpts=lapply(out, '[[', 2)
    
    if(class==TRUE){
      ans=list()
      for(i in 1:rep){
        ans[[i]]=class_input(data[i,], cpttype="mean and variance", method=mul.method, test.stat="Gamma", penalty=penalty, pen.value=pen.value, minseglen=minseglen, param.estimates=param.estimates, out=out[[2]], Q=Q, shape=shape[[i]])
      }
      return(ans)
    }
    else{return(cpts)}
  }
}
# PELT.var.norm=function(data,pen=0,know.mean=FALSE,mu=NA,nprune=FALSE){
#   mll.var.EFK=function(x,n){
#     neg=x<=0
#     x[neg==TRUE]=0.00000000001    
#     return( n*(log(2*pi)+log(x/n)+1))
#   }
#   if((know.mean==FALSE)&(is.na(mu))){
# 	mu=mean(data)
#   }
#   n=length(data)
#   y2=c(0,cumsum((data-mu)^2))
# 
#   lastchangecpts=matrix(NA,nrow=n,ncol=2)
#   lastchangelike=matrix(NA,nrow=n,ncol=2)
#   checklist=NULL
#   lastchangelike[1,]=c(mll.var.EFK(y2[2],1),mll.var.EFK(y2[n+1]-y2[2],n-1)+pen)
#   lastchangecpts[1,]=c(0,1)
#   lastchangelike[2,]=c(mll.var.EFK(y2[3],2),mll.var.EFK(y2[n+1]-y2[3],n-2)+pen)
#   lastchangecpts[2,]=c(0,2)
#   lastchangelike[3,]=c(mll.var.EFK(y2[4],3),mll.var.EFK(y2[n+1]-y2[4],n-3)+pen)
#   lastchangecpts[3,]=c(0,3)
#   noprune=NULL
#   for(tstar in 4:n){
#     tmplike=NULL
#     tmpt=c(checklist, tstar-2)
#     tmplike=lastchangelike[tmpt,1]+mll.var.EFK(y2[tstar+1]-y2[tmpt+1],tstar-tmpt)+pen
#     if(tstar==n){
#       lastchangelike[tstar,]=c(min(c(tmplike,mll.var.EFK(y2[tstar+1]-y2[1],tstar)),na.rm=TRUE),0)
#     }
#     else{
#       lastchangelike[tstar,]=c(min(c(tmplike,mll.var.EFK(y2[tstar+1]-y2[1],tstar)),na.rm=TRUE),mll.var.EFK(y2[n+1]-y2[tstar+1],n-tstar)+pen)
#     }
#     if(lastchangelike[tstar,1]==mll.var.EFK(y2[tstar+1]-y2[1],tstar)){
#       lastchangecpts[tstar,]=c(0,tstar)
#     }
#     else{
#       cpt=tmpt[tmplike==lastchangelike[tstar,1]][1]
#       lastchangecpts[tstar,]=c(cpt,tstar)
#     }
#     checklist=tmpt[tmplike<=lastchangelike[tstar,1]+pen]
#     if(nprune==TRUE){
#       noprune=c(noprune,length(checklist))
#     }
#   }
#   if(nprune==TRUE){
#     return(nprune=noprune)
#   }
#   else{
#     fcpt=NULL
#     last=n
#     while(last!=0){
# 	fcpt=c(fcpt,lastchangecpts[last,2])
# 	last=lastchangecpts[last,1]
#     }
#     return(cpt=sort(fcpt))
#   }
# }


# PELT.mean.norm=function(data,pen=0,nprune=FALSE){
#   mll.mean.EFK=function(x2,x,n){
#     return( x2-(x^2)/n)
#   }
#   n=length(data)
#   y2=c(0,cumsum(data^2))
#   y=c(0,cumsum(data))
# 
#   lastchangecpts=matrix(NA,nrow=n,ncol=2)
#   lastchangelike=matrix(NA,nrow=n,ncol=2)
#   checklist=NULL
#   lastchangelike[1,]=c(mll.mean.EFK(y2[2],y[2],1),mll.mean.EFK(y2[n+1]-y2[2],y[n+1]-y[2],n-1)+pen)
#   lastchangecpts[1,]=c(0,1)
#   lastchangelike[2,]=c(mll.mean.EFK(y2[3],y[3],2),mll.mean.EFK(y2[n+1]-y2[3],y[n+1]-y[3],n-2)+pen)
#   lastchangecpts[2,]=c(0,2)
#   lastchangelike[3,]=c(mll.mean.EFK(y2[4],y[4],3),mll.mean.EFK(y2[n+1]-y2[4],y[n+1]-y[4],n-3)+pen)
#   lastchangecpts[3,]=c(0,3)
#   noprune=NULL
#   for(tstar in 4:n){
#     tmplike=NULL
#     tmpt=c(checklist, tstar-2)
#     tmplike=lastchangelike[tmpt,1]+mll.mean.EFK(y2[tstar+1]-y2[tmpt+1],y[tstar+1]-y[tmpt+1],tstar-tmpt)+pen
#     if(tstar==n){
#       lastchangelike[tstar,]=c(min(c(tmplike,mll.mean.EFK(y2[tstar+1],y[tstar+1],tstar)),na.rm=TRUE),0)
#     }
#     else{
#       lastchangelike[tstar,]=c(min(c(tmplike,mll.mean.EFK(y2[tstar+1],y[tstar+1],tstar)),na.rm=TRUE),mll.mean.EFK(y2[n+1]-y2[tstar+1],y[n+1]-y[tstar+1],n-tstar)+pen)
#     }
#     if(lastchangelike[tstar,1]==mll.mean.EFK(y2[tstar+1],y[tstar+1],tstar)){
#       lastchangecpts[tstar,]=c(0,tstar)
#     }
#     else{
#       cpt=tmpt[tmplike==lastchangelike[tstar,1]][1]
#       lastchangecpts[tstar,]=c(cpt,tstar)
#     }
#     checklist=tmpt[tmplike<=lastchangelike[tstar,1]+pen]
#     if(nprune==TRUE){
#       noprune=c(noprune,length(checklist))
#     }
#   }
#   if(nprune==TRUE){
#     return(noprune)
#   }
#   else{
#     fcpt=NULL
#     last=n
#     while(last!=0){
# 	fcpt=c(fcpt,lastchangecpts[last,2])
# 	last=lastchangecpts[last,1]
#     }
#     return(cpt=sort(fcpt))
#   }
# }

# PELT.meanvar.norm=function(data,pen=0,nprune=FALSE){
#   mll.meanvar.EFK=function(x2,x,n){
#     sigmasq=(1/n)*(x2-(x^2)/n)
#     neg=sigmasq<=0
#     sigmasq[neg==TRUE]=0.00000000001
#     return(n*(log(2*pi)+log(sigmasq)+1))
#   }
#   n=length(data)
#   y2=c(0,cumsum(data^2))
#   y=c(0,cumsum(data))
# 
#   lastchangecpts=matrix(NA,nrow=n,ncol=2)
#   lastchangelike=matrix(NA,nrow=n,ncol=2)
#   checklist=NULL
#   lastchangelike[1,]=c(mll.meanvar.EFK(y2[2],y[2],1),mll.meanvar.EFK(y2[n+1]-y2[2],y[n+1]-y[2],n-1)+pen)
#   lastchangecpts[1,]=c(0,1)
#   lastchangelike[2,]=c(mll.meanvar.EFK(y2[3],y[3],2),mll.meanvar.EFK(y2[n+1]-y2[3],y[n+1]-y[3],n-2)+pen)
#   lastchangecpts[2,]=c(0,2)
#   lastchangelike[3,]=c(mll.meanvar.EFK(y2[4],y[4],3),mll.meanvar.EFK(y2[n+1]-y2[4],y[n+1]-y[4],n-3)+pen)
#   lastchangecpts[3,]=c(0,3)
#   noprune=NULL
#   for(tstar in 4:n){
#     tmplike=NULL
#     tmpt=c(checklist, tstar-2)
#     tmplike=lastchangelike[tmpt,1]+mll.meanvar.EFK(y2[tstar+1]-y2[tmpt+1],y[tstar+1]-y[tmpt+1],tstar-tmpt)+pen
#     if(tstar==n){
#       lastchangelike[tstar,]=c(min(c(tmplike,mll.meanvar.EFK(y2[tstar+1],y[tstar+1],tstar)),na.rm=TRUE),0)
#     }
#     else{
#       lastchangelike[tstar,]=c(min(c(tmplike,mll.meanvar.EFK(y2[tstar+1],y[tstar+1],tstar)),na.rm=TRUE),mll.meanvar.EFK(y2[n+1]-y2[tstar+1],y[n+1]-y[tstar+1],n-tstar)+pen)
#     }
#     if(lastchangelike[tstar,1]==mll.meanvar.EFK(y2[tstar+1],y[tstar+1],tstar)){
#       lastchangecpts[tstar,]=c(0,tstar)
#     }
#     else{
#       cpt=tmpt[tmplike==lastchangelike[tstar,1]][1]
#       lastchangecpts[tstar,]=c(cpt,tstar)
#     }
#     checklist=tmpt[tmplike<=lastchangelike[tstar,1]+pen]
#     if(nprune==TRUE){
#       noprune=c(noprune,length(checklist))
#     }
#   }
#   if(nprune==TRUE){
#     return(noprune)
#   }
#   else{
#     fcpt=NULL
#     last=n
#     while(last!=0){
# 	fcpt=c(fcpt,lastchangecpts[last,2])
# 	last=lastchangecpts[last,1]
#     }
#     return(cpt=sort(fcpt))
#   }
# }



#binseg.var.norm=function(data,Q=5,pen=0,know.mean=FALSE,mu=NA){
#  mll.var=function(x,n){
#    neg=x<=0
#    x[neg==TRUE]=0.00000000001    
#    return( -0.5*n*(log(2*pi)+log(x/n)+1))
#  }
#  n=length(data)
#  if((know.mean==FALSE)&(is.na(mu))){
#	mu=mean(data)
#  }
#  y2=c(0,cumsum((data-mu)^2))
#  tau=c(0,n)
#  cpt=matrix(0,nrow=2,ncol=Q)
#  oldmax=1000
#
#  for(q in 1:Q){
#    lambda=rep(0,n-1)
#    i=1
#    st=tau[1]+1;end=tau[2]
#    null=mll.var(y2[end+1]-y2[st],end-st+1)
#    for(j in 1:(n-1)){
#      if(j==end){
#        st=end+1;i=i+1;end=tau[i+1]
#        null=mll.var(y2[end+1]-y2[st],end-st+1)
#      }else{
#        lambda[j]=mll.var(y2[j+1]-y2[st],j-st+1)+mll.var(y2[end+1]-y2[j+1],end-j)-null
#      }
#    }
#    k=which.max(lambda)[1]
#    cpt[1,q]=k;cpt[2,q]=min(oldmax,max(lambda))
#    oldmax=min(oldmax,max(lambda))
#    tau=sort(c(tau,k))
#  }
#  op.cps=NULL
#  p=1:(Q-1)
#  for(i in 1:length(pen)){
#    criterion=(2*cpt[2,])>=pen[i]
#    if(sum(criterion)==0){
#      op.cps=0
#    }
#    else{
#      op.cps=c(op.cps,max(which((criterion)==TRUE)))
#    }
#  }
#  return(list(cps=cpt,op.cpts=op.cps,pen=pen))
#}


# binseg.mean.norm=function(data,Q=5,pen=0){
#   mll.mean=function(x2,x,n){
#     return( -0.5*(x2-(x^2)/n))
#   }
#   n=length(data)
#   y2=c(0,cumsum(data^2))
#   y=c(0,cumsum(data))
#   tau=c(0,n)
#   cpt=matrix(0,nrow=2,ncol=Q)
#   oldmax=1000
# 
#   for(q in 1:Q){
#     lambda=rep(0,n-1)
#     i=1
#     st=tau[1]+1;end=tau[2]
#     null=mll.mean(y2[end+1]-y2[st],y[end+1]-y[st],end-st+1)
#     for(j in 1:(n-1)){
#       if(j==end){
#         st=end+1;i=i+1;end=tau[i+1]
#         null=mll.mean(y2[end+1]-y2[st],y[end+1]-y[st],end-st+1)
#       }else{
#         lambda[j]=mll.mean(y2[j+1]-y2[st],y[j+1]-y[st],j-st+1)+mll.mean(y2[end+1]-y2[j+1],y[end+1]-y[j+1],end-j)-null
#       }
#     }
#     k=which.max(lambda)[1]
#     cpt[1,q]=k;cpt[2,q]=min(oldmax,max(lambda)) # done so that when we do the decision later we can take the max(which(criterion==T)), rather than min(which(criterion==F))-1
#     oldmax=min(oldmax,max(lambda))
#     tau=sort(c(tau,k))
#   }
#   op.cps=NULL
#   p=1:(Q-1)
#   for(i in 1:length(pen)){
#     criterion=(2*cpt[2,])>=pen[i]
#     if(sum(criterion)==0){
#       op.cps=0
#     }
#     else{
#       op.cps=c(op.cps,max(which((criterion)==TRUE)))
#     }
#   }
#   return(list(cps=cpt,op.cpts=op.cps,pen=pen))
# }

# binseg.meanvar.norm=function(data,Q=5,pen=0){
#   mll.meanvar=function(x2,x,n){
#     sigmasq=(1/n)*(x2-(x^2)/n)
#     neg=sigmasq<=0
#     sigmasq[neg==TRUE]=0.00000000001
#     return(-(n/2)*(log(2*pi)+log(sigmasq)+1))
#   }
#   n=length(data)
#   y2=c(0,cumsum(data^2))
#   y=c(0,cumsum(data))
#   tau=c(0,n)
#   cpt=matrix(0,nrow=2,ncol=Q)
#   oldmax=1000
# 
#   for(q in 1:Q){
#     lambda=rep(0,n-1)
#     i=1
#     st=tau[1]+1;end=tau[2]
#     null=mll.meanvar(y2[end+1]-y2[st],y[end+1]-y[st],end-st+1)
#     for(j in 1:(n-1)){
#       if(j==end){
#         st=end+1;i=i+1;end=tau[i+1]
#         null=mll.meanvar(y2[end+1]-y2[st],y[end+1]-y[st],end-st+1)
#       }else{
# 	if((j-st)<2){lambda[j]=-1*10^(100)}
# 	else if((end-j)<2){lambda[j]=-1*10^(100)}
# 	else{lambda[j]=mll.meanvar(y2[j+1]-y2[st],y[j+1]-y[st],j-st+1)+mll.meanvar(y2[end+1]-y2[j+1],y[end+1]-y[j+1],end-j)-null}
#       }
#     }
#     k=which.max(lambda)[1]
#     cpt[1,q]=k;cpt[2,q]=min(oldmax,max(lambda))
#     oldmax=min(oldmax,max(lambda))
#     tau=sort(c(tau,k))
#   }
#   op.cps=NULL
#   p=1:(Q-1)
#   for(i in 1:length(pen)){
#     criterion=(2*cpt[2,])>=pen[i]
#     if(sum(criterion)==0){
#       op.cps=0
#     }
#     else{
#       op.cps=c(op.cps,max(which((criterion)==TRUE)))
#     }
#   }
#   return(list(cps=cpt,op.cpts=op.cps,pen=pen))
# }
multiple.var.norm=function(data,mul.method="PELT",penalty="MBIC",pen.value=0,Q=5,know.mean=FALSE,mu=NA,class=TRUE,param.estimates=TRUE, minseglen=2){
  if(!((mul.method=="PELT")||(mul.method=="BinSeg"))){
    stop("Multiple Method is not recognised, must be PELT or BinSeg.")
  }
  costfunc = "var.norm"
  if(penalty =="MBIC"){
    costfunc = "var.norm.mbic"
  }
  diffparam=1
  if(is.null(dim(data))==TRUE){
    # single dataset
    n=length(data)
    mu=mu[1]
  }
  else{
    n=ncol(data)
  }
  if(n<4){stop('Data must have atleast 4 observations to fit a changepoint model.')}
  if(n<(2*minseglen)){stop('Minimum segment legnth is too large to include a change in this data')}
  
  pen.value = penalty_decision(penalty, pen.value, n, diffparam, asymcheck=costfunc, method=mul.method)
  
  if(is.null(dim(data))==TRUE){
    # single dataset
    if((know.mean==FALSE)&(is.na(mu))){
      mu=mean(coredata(data))
    }
    out = data_input(data=data,method=mul.method,pen.value=pen.value,costfunc=costfunc,minseglen=minseglen,Q=Q,var=mu)
    
    if(class==TRUE){
      out=class_input(data, cpttype="variance", method=mul.method, test.stat="Normal", penalty=penalty, pen.value=pen.value, minseglen=minseglen, param.estimates=param.estimates, out=out, Q=Q)
      param.est(out)=c(param.est(out),mean=mu)
      return(out)
    }
    else{ return(out[[2]])}
  }
  else{
    rep=nrow(data)
    out=list()
    if(length(mu)!=rep){
      mu=rep(mu,rep)
    }
    
    for(i in 1:rep){
      if((know.mean==FALSE)&(is.na(mu[i]))){
        mu=mean(coredata(data[i,]))
      }
      out[[i]]=data_input(data[i,],method=mul.method,pen.value=pen.value,costfunc=costfunc,minseglen=minseglen,Q=Q,var=mu)
    }
    
    cpts=lapply(out, '[[', 2)
    
    if(class==TRUE){
      ans=list()
      for(i in 1:rep){
        ans[[i]]=class_input(data[i,], cpttype="variance", method=mul.method, test.stat="Normal", penalty=penalty, pen.value=pen.value, minseglen=minseglen, param.estimates=param.estimates, out=out[[i]], Q=Q)
        param.est(ans[[i]])=c(param.est(ans[[i]]),mean=mu[i])
      }
      return(ans)
    }
    else{return(cpts)}
  }
}
multiple.mean.norm=function(data,mul.method="PELT",penalty="MBIC",pen.value=0,Q=5,class=TRUE,param.estimates=TRUE,minseglen){
  if(!((mul.method=="PELT")||(mul.method=="BinSeg"))){
    stop("Multiple Method is not recognised, must be PELT or BinSeg")
  }
  costfunc = "mean.norm"
  if(penalty=="MBIC"){
    costfunc = "mean.norm.mbic"
  }
  diffparam=1
  if(is.null(dim(data))==TRUE){
    # single dataset
    n=length(data) # still works if data is of class ts
  }
  else{
    n=ncol(data)
  }
  if(n<(2*minseglen)){stop('Minimum segment legnth is too large to include a change in this data')}
  
  pen.value = penalty_decision(penalty, pen.value, n, diffparam, asymcheck = costfunc, method=mul.method)
  
  if(is.null(dim(data))==TRUE){
    # single dataset
    out = data_input(data=data,method=mul.method,pen.value=pen.value,costfunc=costfunc,minseglen=minseglen,Q=Q)
    
    if(class==TRUE){
      return(class_input(data, cpttype="mean", method=mul.method, test.stat="Normal", penalty=penalty, pen.value=pen.value, minseglen=minseglen, param.estimates=param.estimates, out=out, Q=Q))
    }
    else{ return(out[[2]])}
  }
  else{
    rep=nrow(data)
    out=list()
    if(class==TRUE){cpts=list()}
    for(i in 1:rep){
      out[[i]]=data_input(data[i,],method=mul.method,pen.value=pen.value,costfunc=costfunc,minseglen=minseglen,Q=Q)
    }
    
    cps=lapply(out, '[[', 2)
    
    if(class==TRUE){
      ans=list()
      for(i in 1:rep){
        ans[[i]]=class_input(data[i,], cpttype="mean", method=mul.method, test.stat="Normal", penalty=penalty, pen.value=pen.value, minseglen=minseglen, param.estimates=param.estimates, out=out[[i]], Q=Q)
      }
      return(ans)
    }
    else{return(cps)}
  }
}
multiple.meanvar.norm=function(data,mul.method="PELT",penalty="MBIC",pen.value=0,Q=5,class=TRUE,param.estimates=TRUE,minseglen){
  if(!((mul.method=="PELT")||(mul.method=="BinSeg"))){
    stop("Multiple Method is not recognised, must be PELT or BinSeg.")
  }
  costfunc = "meanvar.norm"
  if(penalty=="MBIC"){
    costfunc = "meanvar.norm.mbic"
  }
  diffparam=2
  if(is.null(dim(data))==TRUE){
    # single dataset
    n=length(data)
  }
  else{
    n=ncol(data)
  }
  if(n<(2*minseglen)){stop('Minimum segment legnth is too large to include a change in this data')}
  
  pen.value = penalty_decision(penalty, pen.value, n, diffparam, asymcheck = costfunc, method=mul.method)
  
  if(is.null(dim(data))==TRUE){
    # single dataset
    out = data_input(data=data,method=mul.method,pen.value=pen.value,costfunc=costfunc,minseglen=minseglen,Q=Q)
    
    if(class==TRUE){
      return(class_input(data, cpttype="mean and variance", method=mul.method, test.stat="Normal", penalty=penalty, pen.value=pen.value, minseglen=minseglen, param.estimates=param.estimates, out=out, Q=Q))
    }
    else{ return(out[[2]])}
  }
  else{
    rep=nrow(data)
    out=list()
    for(i in 1:rep){
      out[[i]]=data_input(data[i,],method=mul.method,pen.value=pen.value,costfunc=costfunc,minseglen=minseglen,Q=Q)
    }
    
    cps=lapply(out, '[[', 2)
    
    if(class==TRUE){
      ans=list()
      for(i in 1:rep){
        ans[[i]] = class_input(data[i,], cpttype="mean and variance", method=mul.method, test.stat="Normal", penalty=penalty, pen.value=pen.value, minseglen=minseglen, param.estimates=param.estimates, out=out[[i]], Q=Q)
      }
      return(ans)
    }
    else{return(cps)}
  }
}
penalty_decision = function(penalty, pen.value, n, diffparam, asymcheck, method){
  
  
  if((penalty=="SIC0") || (penalty=="BIC0")){
    pen.return=diffparam*log(n)
  }
  else if((penalty=="SIC") || (penalty=="BIC")){
    pen.return=(diffparam+1)*log(n)
  }
  else if(penalty=="MBIC"){
    pen.return=(diffparam+2)*log(n)
  }
  else if((penalty=="SIC1") || (penalty=="BIC1")){
    stop("SIC1 and BIC1 have been depreciated, use SIC or BIC for the same result.")
  }
  else if(penalty=="AIC0"){
    pen.return=2*diffparam
  }
  else if(penalty=="AIC"){
    pen.return=2*(diffparam+1)
  }
  else if(penalty=="AIC1"){
    stop("AIC1 has been depreciated, use AIC for the same result.")
  }
  else if(penalty=="Hannan-Quinn0"){
    pen.return=2*diffparam*log(log(n))
  }
  else if(penalty=="Hannan-Quinn"){
    pen.return=2*(diffparam+1)*log(log(n))
  }
  else if(penalty=="Hannan-Quinn1"){
    stop("Hannan-Quinn1 has been depreciated, use Hannan-Quinn for the same result.")
  }
  else if(penalty=="None"){
    pen.return=0
  }
  else if((penalty!="Manual")&&(penalty!="Asymptotic")){
    stop('Unknown Penalty')
  }
  if((penalty=="Manual")&&(is.numeric(pen.value)==FALSE)){
    pen.value=try(eval(parse(text=paste(pen.value))),silent=TRUE)
    if(class(pen.value)=='try-error'){
      stop('Your manual penalty cannot be evaluated')
    }else{
      pen.return=pen.value
    }
  }
  
  if((penalty=="Manual")&&(is.numeric(pen.value)==TRUE)){
    pen.return=pen.value
  }
  
  
  if(penalty=="Asymptotic"){
    if((pen.value <= 0) || (pen.value > 1)){
      stop('Asymptotic penalty values must be > 0 and <= 1')
    }
    if(method != "AMOC"){
      warning('Asymptotic penalty value is not accurate for multiple changes, it should be treated the same as a manual penalty choice.')
    }
    if(asymcheck == "mean.norm"){
      alpha=pen.value
      alogn=(2*log(log(n)))^(-(1/2))
      blogn=(alogn^(-1))+(1/2)*alogn*log(log(log(n)))
      pen.return=(-alogn*log(log((1-alpha+exp(-2*(pi^(1/2))*exp(blogn/alogn)))^(-1/(2*(pi^(1/2))))))+blogn)^2
    }else if(asymcheck == "var.norm"){ 
      alpha=pen.value
      alogn=sqrt(2*log(log(n)))
      blogn=2*log(log(n))+ (log(log(log(n))))/2 - log(gamma(1/2))
      pen.return=(-(log(log((1-alpha+exp(-2*exp(blogn)))^(-1/2))))/alogn + blogn/alogn)^2
    }else if(asymcheck == "meanvar.norm"){
      alpha=pen.value
      alogn=sqrt(2*log(log(n)))
      blogn=2*log(log(n))+ log(log(log(n)))
      pen.return=(-(log(log((1-alpha+exp(-2*exp(blogn)))^(-1/2))))/alogn + blogn/alogn)^2
    }else if(asymcheck == "reg.norm"){
      alpha=pen.value
      top=-(log(log((1 - alpha + exp(-2*exp(2*(log(log(n)))+(diffparam/2)*(log(log(log(n))))- log(gamma(diffparam/2)))))^(-1/2))))  +  2*(log(log(n)))+(diffparam/2)*(log(log(log(n))))- log(gamma(diffparam/2))
      bottom=(2*log(log(n)))^(1/2)
      pen.return=(top/bottom)^2
    }else if(asymcheck == "var.css"){
      if(pen.value==0.01){pen.return=1.628}
      else if(pen.value==0.05){pen.return=1.358}
      else if(pen.value==0.1){pen.return=1.224}
      else if(pen.value==0.25){pen.return=1.019}
      else if(pen.value==0.5){pen.return=0.828}
      else if(pen.value==0.75){pen.return=0.677}
      else if(pen.value==0.9){pen.return=0.571}
      else if(pen.value==0.95){pen.return=0.520}
      else{stop('Only alpha values of 0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95 are valid for CSS')}
    }else if(asymcheck == "mean.cusum"){
      stop('Asymptotic penalties have not been implemented yet for CUSUM')
    }else if(asymcheck == "meanvar.gamma"){
      stop('Asymptotic penalties for the Gamma test statistic are not defined, please choose an alternative penalty type')
    }else if(asymcheck == "meanvar.exp"){
      alpha=pen.value
      an=(2*log(log(n)))^(1/2)
      bn=2*log(log(n))+(1/2)*log(log(log(n)))-(1/2)*log(pi)
      pen.return=(-1/an)*log(-0.5*log(1-alpha))+bn
      if(alpha==1){pen.return=1.42417} # value of 1 gives log(0), this is alpha=0.99999999999999993
    }else if(asymcheck == "meanvar.poisson"){
      stop('Asymptotic penalties for the Poisson test statistic are not available yet, please choose an alternative penalty type')
    }
  }
  #if(method=="AMOC"){
  #  pen.return=pen.value
  #}
  if(pen.return < 0){
    stop('pen.value cannot be negative, please change your penalty value')
  }else{
    return(pen.return)
  }
}
single.meanvar.poisson.calc = function(data,extrainf=TRUE,minseglen){
    singledim=function(data,extrainf=TRUE,minseglen){
      n=length(data)
      y=c(0,cumsum(data))
      if(y[n+1]==0){
        null=Inf
      }
      else{
        null=2*y[n+1]*log(n) - 2*y[n+1]*log(y[n+1])
      }
      taustar=minseglen:(n-minseglen)
      tmp=2*log(taustar)*y[taustar+1] -2*y[taustar+1]*log(y[taustar+1]) + 2*log(n-taustar)*(y[n+1]-y[taustar+1])-2*(y[n+1]-y[taustar+1])*log((y[n+1]-y[taustar+1]))
      if(sum(is.na(tmp))!=0){
        tmp[which(is.na(tmp))]=Inf
      }
      tau=which(tmp==min(tmp,na.rm=T))[1]
      taulike=tmp[tau]
      tau=tau+minseglen-1 # correcting for the fact that we are starting at minseglen
      if(extrainf==TRUE){
        out=c(tau,null,taulike)
        names(out)=c('cpt','null','alt')
        return(out)
      }
      else{
        return(tau)
      }
    }
    
    if(is.null(dim(data))==TRUE){
      # single data set
      cpt=singledim(data,extrainf,minseglen)
      return(cpt)
    }
    else{
      rep=nrow(data)
      n=ncol(data)
      cpt=NULL
      if(extrainf==FALSE){
        for(i in 1:rep){
          cpt[i]=singledim(data[i,],extrainf,minseglen)
        }
      }
      else{
        cpt=matrix(0,ncol=3,nrow=rep)
        for(i in 1:rep){
          cpt[i,]=singledim(data[i,],extrainf,minseglen)
        }
        colnames(cpt)=c('cpt','null','alt')
      }
      return(cpt)
    }
  }
single.meanvar.poisson<-function(data,penalty="MBIC",pen.value=0,class=TRUE,param.estimates=TRUE,minseglen){
  if((sum(data<0)>0)){stop('Poisson test statistic requires positive data')}
  if(sum(as.integer(data)==data)!=length(data)){stop('Poisson test statistic requires integer data')}
  if(is.null(dim(data))==TRUE){
    # single dataset
    n=length(data)
  }
  else{
    n=ncol(data)
  }
  if(n<4){stop('Data must have atleast 4 observations to fit a changepoint model.')}
  if(n<(2*minseglen)){stop('Minimum segment legnth is too large to include a change in this data')}
  
  pen.value = penalty_decision(penalty, pen.value, n, diffparam=1, asymcheck="meanvar.poisson", method="AMOC")   
  if(is.null(dim(data))==TRUE){
    tmp=single.meanvar.poisson.calc(coredata(data),extrainf=TRUE,minseglen)
    if(penalty=="MBIC"){
      tmp[3]=tmp[3]+log(tmp[1])+log(n-tmp[1]+1)
    }
    ans=decision(tmp[1],tmp[2],tmp[3],penalty,n,diffparam=1,pen.value)
    if(class==TRUE){
      return(class_input(data, cpttype="mean and variance", method="AMOC", test.stat="Poisson", penalty=penalty, pen.value=ans$pen, minseglen=minseglen, param.estimates=param.estimates, out=c(0,ans$cpt)))
    }
    else{ return(ans$cpt)}
  }
  else{ 
    tmp=single.meanvar.poisson.calc(data,extrainf=TRUE,minseglen)
    if(penalty=="MBIC"){
      tmp[,3]=tmp[,3]+log(tmp[,1])+log(n-tmp[,1]+1)
    }
    ans=decision(tmp[,1],tmp[,2],tmp[,3],penalty,n,diffparam=1,pen.value)
    if(class==TRUE){
      rep=nrow(data)
      out=list()
      for(i in 1:rep){
        out[[i]]=class_input(data[i,], cpttype="mean and variance", method="AMOC", test.stat="Poisson", penalty=penalty, pen.value=ans$pen, minseglen=minseglen, param.estimates=param.estimates, out=c(0,ans$cpt[i]))
      }
      return(out)
    }
    else{ return(ans$cpt)}
  }
}
multiple.meanvar.poisson=function(data,mul.method="PELT",penalty="MBIC",pen.value=0,Q=5,class=TRUE,param.estimates=TRUE,minseglen){
  if((sum(data<0)>0)){stop('Poisson test statistic requires positive data')}
  if(sum(as.integer(data)==data)!=length(data)){stop('Poisson test statistic requires integer data')}
  if(!((mul.method=="PELT")||(mul.method=="BinSeg"))){
    stop("Multiple Method is not recognised, must be PELT or BinSeg.")
  }
  costfunc = "meanvar.poisson"
  if(penalty=="MBIC"){
    costfunc = "meanvar.poisson.mbic"
  }
  
  diffparam=1
  if(is.null(dim(data))==TRUE){
    # single dataset
    n=length(data)
  }
  else{
    n=ncol(data)
  }
  if(n<(2*minseglen)){stop('Minimum segment legnth is too large to include a change in this data')}
  
  pen.value = penalty_decision(penalty, pen.value, n, diffparam=1, asymcheck=costfunc, method=mul.method)
  if(is.null(dim(data))==TRUE){
    # single dataset
    out = data_input(data=data,method=mul.method,pen.value=pen.value,costfunc=costfunc,minseglen=minseglen,Q=Q)
    
    if(class==TRUE){
      return(class_input(data, cpttype="mean and variance", method=mul.method, test.stat="Poisson", penalty=penalty, pen.value=pen.value, minseglen=minseglen, param.estimates=param.estimates, out=out, Q=Q))
    }
    else{ return(out[[2]])}
  }
  else{
    rep=nrow(data)
    out=list()
    for(i in 1:rep){
      out[[i]]=data_input(data[i,],method=mul.method,pen.value=pen.value,costfunc=costfunc,minseglen=minseglen,Q=Q)  
    }
    
    cpts=lapply(out, '[[', 2)
    
    if(class==TRUE){
      ans=list()
      for(i in 1:rep){
        ans[[i]]=class_input(data[i,], cpttype="mean and variance", method=mul.method, test.stat="Poisson", penalty=penalty, pen.value=pen.value, minseglen=minseglen, param.estimates=param.estimates, out=out[[i]], Q=Q)
      }
      return(ans)
    }
    else{return(cpts)}
  }
}
######## Function to run PELT for a Range Of Penalty values.  ############
range_of_penalties <- function(sumstat,cost = "mean.norm",PELT = T,min_pen=log(length(sumstat)/3-1),max_pen=10*log(length(sumstat)/3-1), shape = 1, minseglen) {
  
  NCALC=0
  pen_interval <- c(min_pen,max_pen)
  n = length(sumstat)/3 - 1
  
  test_penalties <- NULL
  numberofchangepoints <- NULL
  penal <- NULL
  overall_cost <- array() 
  segmentations <- NULL
  b_between <- array() 
  
  ##### Want to store and use Func, M and CP in PELT 
  
  count <- 0
  
  while (length(pen_interval) > 0){
    
    new_numcpts <- array()
    new_penalty <- array()
    new_cpts <- array() 
    
    for (b in 1:length(pen_interval)) {
      
      ans<- PELT(sumstat,pen=pen_interval[b], cost_func = cost , shape = shape, minseglen = minseglen)
      resultingcpts <- ans[[2]]
      new_numcpts[b] <- length(resultingcpts)
      new_cpts[b] <- list(resultingcpts[-length(resultingcpts)])
      new_penalty[b] <- ans[[3]][n+1]-(ans[[4]][n+1]-1)*pen_interval[b]
    }
    
    if (count == 0){
      print(paste("Maximum number of runs of algorithm = ", new_numcpts[1] - new_numcpts[2] + 2, sep = ""))
      count <- count + length(new_numcpts)
      print(paste("Completed runs = ", count, sep = ""))
    }
    
    else{
      count <- count + length(new_numcpts)
      print(paste("Completed runs = ", count, sep = ""))
    }
    
    
    ## Add the values calculated to the already stored values 
    test_penalties <- unique((sort(c(test_penalties,pen_interval))))
    new_numcpts <- c(numberofchangepoints,new_numcpts)
    new_penalty <- c(penal,new_penalty)
    
    new_cpts <- c(segmentations,new_cpts)
    numberofchangepoints <- -sort(-new_numcpts) ##can use sort to re-order 
    penal <- sort(new_penalty)
    
    ls <- array()
    
    for (l in 1:length(new_cpts)){
      ls[l] <- length(new_cpts[[l]])    
    } 
    
    
    ls1 <- sort(ls,index.return = T, decreasing = T)
    ls1 <- ls1$ix
    
    
    segmentations <- new_cpts[c(ls1)]
    
    pen_interval <- NULL
    tmppen_interval <- NULL
    
    for (i in 1:(length(test_penalties)-1)){
      if(abs(numberofchangepoints[i]-numberofchangepoints[i+1])>1){ ##only need to add a beta if difference in cpts>1
        j <- i+1
        tmppen_interval <- (penal[j] - penal[i]) * (((numberofchangepoints[i]) - (numberofchangepoints[j]))^-1)
        pen_interval <- c(pen_interval, tmppen_interval )
      }
    }
    
    if(length(pen_interval)>0){
      for(k in length(pen_interval):1){ 
        index <- which.min(abs(pen_interval[k]-test_penalties)) 
        if (isTRUE(all.equal(pen_interval[k], test_penalties[index]))){
          pen_interval=pen_interval[-k]
        }
      }
    }
  }
  
  
  ##PRUNE VALUES WITH SAME num_cp
  for(j in length(test_penalties):2){
    if(numberofchangepoints[j]==numberofchangepoints[j-1]){
      numberofchangepoints=numberofchangepoints[-j]
      test_penalties=test_penalties[-j]
      penal=penal[-j]
      segmentations = segmentations[-j]
    }
  }
  
  
  
  ###calculate beta intervals
  nb=length(test_penalties)
  beta.int=rep(0,nb)
  beta.e=rep(0,nb)
  for(k in 1:nb){
    if(k==1){
      beta.int[1]=test_penalties[1]
    }else{
      beta.int[k]=beta.e[k-1]   
    }
    if(k==nb){
      beta.e[k]=test_penalties[k]
    }else{
      beta.e[k]=(penal[k]-penal[k+1])/(numberofchangepoints[k+1]-numberofchangepoints[k])
    }
    
  }
  
  return(list(cpt.out = rbind(beta_interval = beta.int,numberofchangepoints,penalised_cost = penal),changepoints = segmentations))
  #segmentations is output matrix
  #beta.int
}
single.mean.norm.calc = function(data,extrainf=TRUE,minseglen){
    singledim=function(data,extrainf=TRUE,minseglen){
      n=length(data)
      y=c(0,cumsum(data))
      y2=c(0,cumsum(data^2))
      null=y2[n+1]-y[n+1]^2/n
      taustar=minseglen:(n-minseglen+1)
      tmp=y2[taustar+1]-y[taustar+1]^2/taustar + (y2[n+1]-y2[taustar+1]) - ((y[n+1]-y[taustar+1])^2)/(n-taustar)
      
      tau=which(tmp==min(tmp,na.rm=T))[1]
      taulike=tmp[tau]
      tau=tau+minseglen-1 # correcting for the fact that we are starting at minseglen
      if(extrainf==TRUE){
        out=c(tau,null,taulike)
        names(out)=c('cpt','null','alt')
        return(out)
      }
      else{
        return(tau)
      }
    }
    
    
    if(is.null(dim(data))==TRUE){
      # single data set
      cpt=singledim(data,extrainf,minseglen)
      return(cpt)
    }
    else{
      rep=nrow(data)
      n=ncol(data)
      cpt=NULL
      if(extrainf==FALSE){
        for(i in 1:rep){
          cpt[i]=singledim(data[i,],extrainf,minseglen)
        }
      }
      else{
        cpt=matrix(0,ncol=3,nrow=rep)
        for(i in 1:rep){
          cpt[i,]=singledim(data[i,],extrainf,minseglen)
        }
        colnames(cpt)=c('cpt','null','alt')
      }
      return(cpt)
    }
  }
single.mean.norm<-function(data,penalty="MBIC",pen.value=0,class=TRUE,param.estimates=TRUE,minseglen){
  if(is.null(dim(data))==TRUE){
    # single dataset
    n=length(data)
  }
  else{
    n=ncol(data)
  }
  if(n<2){stop('Data must have atleast 2 observations to fit a changepoint model.')}
  if(n<(2*minseglen)){stop('Minimum segment legnth is too large to include a change in this data')}
  
  pen.value = penalty_decision(penalty, pen.value, n, diffparam=1, asymcheck="mean.norm", method="AMOC")
  if(is.null(dim(data))==TRUE){ # single dataset
    tmp=single.mean.norm.calc(coredata(data),extrainf=TRUE,minseglen)
    if(penalty=="MBIC"){
      tmp[3]=tmp[3]+log(tmp[1])+log(n-tmp[1]+1)
    }
    ans=decision(tmp[1],tmp[2],tmp[3],penalty,n,diffparam=1,pen.value)
    
    if(class==TRUE){
      return(class_input(data, cpttype="mean", method="AMOC", test.stat="Normal", penalty=penalty, pen.value=ans$pen, minseglen=minseglen, param.estimates=param.estimates, out=c(0, ans$cpt)))
    }
    else{ 
      alogn=(2*log(log(n)))^(-(1/2))
      blogn=(alogn^(-1))+(1/2)*alogn*log(log(log(n)))  # Chen & Gupta (2000) pg10
      out=c(ans$cpt,exp(-2*(pi^(1/2))*exp(-alogn*sqrt(abs(tmp[2]-tmp[3]))+(alogn^{-1})*blogn))-exp(-2*(pi^(1/2))*exp((alogn^{-1})*blogn)))
      names(out)=c('cpt','conf.value')
      return(out)
    }
  }
  else{ 
    tmp=single.mean.norm.calc(data,extrainf=TRUE,minseglen)
    if(penalty=="MBIC"){
      tmp[,3]=tmp[,3]+log(tmp[,1])+log(n-tmp[,1]+1)
    }
    ans=decision(tmp[,1],tmp[,2],tmp[,3],penalty,n,diffparam=1,pen.value)
    if(class==TRUE){
      rep=nrow(data)
      out=list()
      for(i in 1:rep){
        out[[i]] = class_input(data, cpttype="mean", method="AMOC", test.stat="Normal", penalty=penalty, pen.value=ans$pen, minseglen=minseglen, param.estimates=param.estimates, out=c(0, ans$cpt[i]))
      }
      return(out)
    }
    else{ 
      alogn=(2*log(log(n)))^(-(1/2))
      blogn=(alogn^(-1))+(1/2)*alogn*log(log(log(n)))  # Chen & Gupta (2000) pg10
      out=cbind(ans$cpt,exp(-2*(pi^(1/2))*exp(-alogn*sqrt(abs(tmp[,2]-tmp[,3]))+(alogn^{-1})*blogn))-exp(-2*(pi^(1/2))*exp((alogn^{-1})*blogn)))
      colnames(out)=c('cpt','conf.value')
      rownames(out)=NULL
      return(out)
    }
  }
}
single.var.norm.calc <- function(data,mu,extrainf=TRUE,minseglen){
  n=length(data)
  y=c(0,cumsum((data-mu)^2))
  null=n*log(y[n+1]/n)
  taustar=minseglen:(n-minseglen+1)
  sigma1=y[taustar+1]/taustar
  neg=sigma1<=0
  sigma1[neg==TRUE]=1*10^(-10)
  sigman=(y[n+1]-y[taustar+1])/(n-taustar)
  neg=sigman<=0
  sigman[neg==TRUE]=1*10^(-10)
  tmp=taustar*log(sigma1) + (n-taustar)*log(sigman)
  
  tau=which(tmp==min(tmp,na.rm=T))[1]
  taulike=tmp[tau]
  tau=tau+minseglen-1 # correcting for the fact that we are starting at minseglen
  if(extrainf==TRUE){
    out=c(tau,null,taulike)
    names(out)=c('cpt','null','alt')
    return(out)
  }
  else{
    return(tau)
  }
}    
single.var.norm<-function(data,penalty="MBIC",pen.value=0,know.mean=FALSE,mu=NA,class=TRUE,param.estimates=TRUE,minseglen){
  if(is.null(dim(data))==TRUE){
    # single dataset
    n=length(data)
    mu=mu[1]
  }
  else{
    n=ncol(data)
  }
  if(n<4){stop('Data must have atleast 4 observations to fit a changepoint model.')}
  if(n<(2*minseglen)){stop('Minimum segment legnth is too large to include a change in this data')}
  
  pen.value = penalty_decision(penalty, pen.value, n, diffparam=1, asymcheck="var.norm", method="AMOC")
  
  if(is.null(dim(data))==TRUE){
    if((know.mean==FALSE)&(is.na(mu))){
      mu=mean(coredata(data))
    }
    tmp=single.var.norm.calc(coredata(data),mu,extrainf=TRUE,minseglen)
    if(penalty=="MBIC"){
      tmp[3]=tmp[3]+log(tmp[1])+log(n-tmp[1]+1)
    }
    ans=decision(tmp[1],tmp[2],tmp[3],penalty,n,diffparam=1,pen.value)
    if(class==TRUE){
      out=class_input(data, cpttype="variance", method="AMOC", test.stat="Normal", penalty=penalty, pen.value=pen.value, minseglen=minseglen, param.estimates=param.estimates, out=c(0, ans$cpt))
      param.est(out)=c(param.est(out),mean=mu)
      return(out)
    }
    else{ 
      alogn=sqrt(2*log(log(n)))
      blogn=2*log(log(n))+ (log(log(log(n))))/2 - log(gamma(1/2))
      out=c(ans$cpt,exp(-2*exp(-alogn*sqrt(abs(tmp[2]-tmp[3]))+blogn))-exp(-2*exp(blogn)))  # Chen & Gupta (2000) pg27
      names(out)=c('cpt','conf.value')
      return(out)
    }
  }
  else{ 
    rep=nrow(data)
    tmp=matrix(0,ncol=3,nrow=rep)
    if(length(mu)!=rep){
      mu=rep(mu,rep)
    }
    for(i in 1:rep){
      if((know.mean==FALSE)&(is.na(mu[i]))){
        mu=mean(coredata(data[i,]))
      }
      tmp[i,]=single.var.norm.calc(data[i,],mu[i],extrainf=TRUE,minseglen)
    }
    
    if(penalty=="MBIC"){
      tmp[,3]=tmp[,3]+log(tmp[,1])+log(n-tmp[,1]+1)
    }
    ans=decision(tmp[,1],tmp[,2],tmp[,3],penalty,n,diffparam=1,pen.value)
    if(class==TRUE){
      out=list()
      for(i in 1:rep){
        out[[i]] = class_input(data, cpttype="variance", method="AMOC", test.stat="Normal", penalty=penalty, pen.value=ans$pen, minseglen=minseglen, param.estimates=param.estimates, out=c(0, ans$cpt[i]))
        param.est(out[[i]])=c(param.est(out[[i]]),mean=mu[i])
      }
      return(out)
    }
    else{ 
      alogn=sqrt(2*log(log(n)))
      blogn=2*log(log(n))+ (log(log(log(n))))/2 - log(gamma(1/2))
      out=cbind(ans$cpt,exp(-2*exp(-alogn*sqrt(abs(tmp[,2]-tmp[,3]))+blogn))-exp(-2*exp(blogn)))  # Chen & Gupta (2000) pg27
      colnames(out)=c('cpt','conf.value')
      rownames(out)=NULL
      return(out)
    }
  }
}
single.meanvar.norm.calc = function(data,extrainf=TRUE,minseglen){
    singledim=function(data,extrainf=TRUE,minseglen){
      n=length(data)
      y=c(0,cumsum(data))
      y2=c(0,cumsum((data)^2))
      null=n*log((y2[n+1]-(y[n+1]^2/n))/n)
      taustar=minseglen:(n-minseglen+1)
      sigma1=((y2[taustar+1]-(y[taustar+1]^2/taustar))/taustar)
      neg=sigma1<=0
      sigma1[neg==TRUE]=1*10^(-10)
      sigman=((y2[n+1]-y2[taustar+1])-((y[n+1]-y[taustar+1])^2/(n-taustar)))/(n-taustar)
      neg=sigman<=0
      sigman[neg==TRUE]=1*10^(-10)
      tmp=taustar*log(sigma1) + (n-taustar)*log(sigman)
      
      tau=which(tmp==min(tmp,na.rm=T))[1]
      taulike=tmp[tau]
      tau=tau+minseglen-1 # correcting for the fact that we are starting at minseglen
      if(extrainf==TRUE){
        out=c(tau,null,taulike)
        names(out)=c('cpt','null','alt')
        return(out)
      }
      else{
        return(tau)
      }
    }
    
    
    if(is.null(dim(data))==TRUE){
      # single data set
      cpt=singledim(data,extrainf,minseglen)
      return(cpt)
    }
    else{
      rep=nrow(data)
      n=ncol(data)
      cpt=NULL
      if(extrainf==FALSE){
        for(i in 1:rep){
          cpt[i]=singledim(data[i,],extrainf,minseglen)
        }
      }
      else{
        cpt=matrix(0,ncol=3,nrow=rep)
        for(i in 1:rep){
          cpt[i,]=singledim(data[i,],extrainf,minseglen)
        }
        colnames(cpt)=c('cpt','null','alt')
      }
      return(cpt)
    }
  }
single.meanvar.norm<-function(data,penalty="MBIC",pen.value=0,class=TRUE,param.estimates=TRUE,minseglen){
  if(is.null(dim(data))==TRUE){
    # single dataset
    n=length(data)
  }
  else{
    n=ncol(data)
  } 
  if(n<4){stop('Data must have atleast 4 observations to fit a changepoint model.')}
  if(n<(2*minseglen)){stop('Minimum segment legnth is too large to include a change in this data')}
  
  pen.value = penalty_decision(penalty, pen.value, n, diffparam=1, asymcheck="meanvar.norm", method="AMOC")
  
  if(is.null(dim(data))==TRUE){
    tmp=single.meanvar.norm.calc(coredata(data),extrainf=TRUE,minseglen)
    if(penalty=="MBIC"){
      tmp[3]=tmp[3]+log(tmp[1])+log(n-tmp[1]+1)
    }
    ans=decision(tmp[1],tmp[2],tmp[3],penalty,n,diffparam=2,pen.value)
    if(class==TRUE){
      return(class_input(data, cpttype="mean and variance", method="AMOC", test.stat="Normal", penalty=penalty, pen.value=ans$pen, minseglen=minseglen, param.estimates=param.estimates, out=c(0,ans$cpt)))
    }
    else{ 
      alogn=sqrt(2*log(log(n)))
      blogn=2*log(log(n))+ (log(log(log(n))))
      out=c(ans$cpt,exp(-2*exp(-alogn*sqrt(abs(tmp[2]-tmp[3]))+blogn))-exp(-2*exp(blogn)))   # Chen & Gupta (2000) pg54
      names(out)=c('cpt','conf.value')
      return(out)
    }
  }
  else{ 
    tmp=single.meanvar.norm.calc(data,extrainf=TRUE,minseglen)
    if(penalty=="MBIC"){
      tmp[,3]=tmp[,3]+log(tmp[,1])+log(n-tmp[,1]+1)
    }
    ans=decision(tmp[,1],tmp[,2],tmp[,3],penalty,n,diffparam=2,pen.value)
    if(class==TRUE){
      rep=nrow(data)
      out=list()
      for(i in 1:rep){
        out[[i]]= class_input(data[i,], cpttype="mean and variance", method="AMOC", test.stat="Normal", penalty=penalty, pen.value=ans$pen, minseglen=minseglen, param.estimates=param.estimates, out=c(0,ans$cpt[i]))
      }
      return(out)
    }
    else{ 
      alogn=sqrt(2*log(log(n)))
      blogn=2*log(log(n))+ (log(log(log(n))))
      out=cbind(ans$cpt,exp(-2*exp(-alogn*sqrt(abs(tmp[,2]-tmp[,3]))+blogn))-exp(-2*exp(blogn)))   # Chen & Gupta (2000) pg54
      colnames(out)=c('cpt','conf.value')
      rownames(out)=NULL
      return(out)
    }
  }
}
.onAttach = function(libname, pkgname) {
  f <- read.dcf(file.path(libname, pkgname, "DESCRIPTION"),
                c("Version", "Date"))
  packageStartupMessage('Successfully loaded changepoint package version ',
                        f[1,1],'\n NOTE: Predefined penalty values changed in version 2.2.  Previous penalty values with a postfix 1 i.e. SIC1 are now without i.e. SIC ', 
                        'and previous penalties without a postfix i.e. SIC are now with a postfix 0 i.e. SIC0. See NEWS and help files for further details.')
}


