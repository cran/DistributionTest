zk.test <-
function(x,distr,para,N=1000){
  DNAME <- deparse(substitute(x))
  x <- sort(x[complete.cases(x)])
  n <- length(x)
  if (n < 8) 
    stop("sample size must be greater than 7")
  Distr<-c("unif","exp","norm","lognorm","gamma","weibull","t")
  if (!is.element(distr,Distr))
    stop("unsupported distribution")
  if (missing(para)) 
    para <- NULL
  if(distr=="unif"){
    if(is.null(para)){
      # Moment estimation
      a<-mean(x)-sqrt(3)*sd(x)
      b<-mean(x)+sqrt(3)*sd(x)
    }else{
      a<-para$min
      b<-para$max
    }
    p<-punif(x,min = a, max =b,lower.tail = TRUE, log.p = FALSE )
    for(i in 1:n){
      if(p[i]==0){p[i]=p[i]+0.0000000001}
      if(p[i]==1){p[i]=p[i]-0.0000000001}
    }
  }
  if(distr=="exp"){
    if (any(x<=0)) {
      stop("sample must be greater than 0")
    }
    if(is.null(para)){
      a<-1/mean(x)# Maximum likelihood estimation
    }else{
      a<-para$rate
    }
    p<-pexp(x,rate = a,lower.tail = TRUE, log.p = FALSE )
    for(i in 1:n){
      if(p[i]==0){p[i]=p[i]+0.0000000001}
      if(p[i]==1){p[i]=p[i]-0.0000000001}
    }
  }
  if(distr=="norm"){
    if(is.null(para)){
      p <- pnorm((x - mean(x))/sd(x),mean = 0,sd = 1,lower.tail = TRUE,log.p = FALSE)
    }else{
      a<-para$mean
      b<-para$sd
      p <- pnorm(x,mean = a,sd = b,lower.tail = TRUE,log.p = FALSE)
    }
  }
  if(distr=="lognorm"){
    if (any(x<=0)) {
      stop("sample must be greater than 0")
    }
    if(is.null(para)){
      p<-pnorm((log(x)-mean(log(x)))/sd(log(x)), mean = 0,sd = 1, lower.tail = TRUE, log.p = FALSE)
    }else{
      a<-para$mean
      b<-para$sd
      p <- pnorm(log(x),mean = a,sd = b,lower.tail = TRUE,log.p = FALSE)
    }
  }
  if(distr=="gamma"){
    if (any(x<=0)) {
      stop("sample must be greater than 0")
    }
    if(is.null(para)){
      # Maximum likelihood estimation
      H<-log(mean(x)/((cumprod(x)[n])^(1/n)))
      a<-0
      if(H>0&H<=0.5772){
        a<-(0.500876+0.1648852*H-0.054427*H^2)/H
      }
      if(H>0.6772&H<=17){
        a<-(8.898919+9.05995*H+0.9775373*H^2)/(H*(17.79728+11.968477*H+H^2))
      }
      b<-mean(x)/a
    }else{
      a<-para$shape
      b<-para$scale
    }
    p<-pgamma(x,shape = a,scale = b,lower.tail = TRUE, log.p = FALSE)
  }
  if(distr=="weibull"){
    if (any(x<=0)) {
      stop("sample must be greater than 0")
    }
    if(is.null(para)){
      # Estimation based on extreme value distribution
      a<-1.283/sd(log(x))
      b<-exp(mean(log(x))+0.5772*0.7797*sd(log(x)))
    }else{
      a<-para$shape
      b<-para$scale
    }
    p<-pweibull(x,shape = a,scale = b,lower.tail = TRUE, log.p = FALSE)
  }
  if(distr=="t"){
    if(is.null(para)){
      s<-fitdistr(x,densfun = "t")
      a<-unname(s$estimate)[3]
    }else{
      a<-para$df
    }
    p<-pt(x,df = a,lower.tail = TRUE, log.p = FALSE)
  }
  
  logp1 <- log((seq(1:n)-0.5)/(n*p))
  logp2 <- log((n-seq(1:n)+0.5)/(n*(1-p)))
  h<-(seq(1:n)-0.5)*logp1+(n-seq(1:n)+0.5)*logp2
  ZK <- max(h) 
  
  # Monte Carlo simulation
  S<-0
  if(distr=="unif"){
    for (j in 1:N) {
      R<-runif(n,a,b)
      R<-sort(R)
      p<-punif(R,min = a, max =b,lower.tail = TRUE, log.p = FALSE )
      for(i in 1:n){
        if(p[i]==0){p[i]=p[i]+0.0000000001}
        if(p[i]==1){p[i]=p[i]-0.0000000001}
      }
      logp1 <- log((seq(1:n)-0.5)/(n*p))
      logp2 <- log((n-seq(1:n)+0.5)/(n*(1-p)))
      h<-(seq(1:n)-0.5)*logp1+(n-seq(1:n)+0.5)*logp2
      zk <- max(h)
      S <- S + (zk > ZK)
    }
  }
  if(distr=="exp"){
    for (j in 1:N) {
      R<-rexp(n,a)
      R<-sort(R)
      p<-pexp(R,rate=a,lower.tail = TRUE, log.p = FALSE )
      logp1 <- log((seq(1:n)-0.5)/(n*p))
      logp2 <- log((n-seq(1:n)+0.5)/(n*(1-p)))
      h<-(seq(1:n)-0.5)*logp1+(n-seq(1:n)+0.5)*logp2
      zk <- max(h)
      S <- S + (zk > ZK)
    }
  }
  if(distr=="norm"||distr=="lognorm"){
    if(is.null(para)){
      for (j in 1:N) {
        R<-rnorm(n)
        R<-sort(R)
        p <- pnorm((R-mean(R))/sd(R),mean = 0,sd=1,lower.tail = TRUE,log.p = FALSE)
        logp1 <- log((seq(1:n)-0.5)/(n*p))
        logp2 <- log((n-seq(1:n)+0.5)/(n*(1-p)))
        h<-(seq(1:n)-0.5)*logp1+(n-seq(1:n)+0.5)*logp2
        zk <- max(h)
        S <- S + (zk > ZK)
      }
    }else{
      for (j in 1:N) {
        R<-rnorm(n,a,b)
        R<-sort(R)
        p <- pnorm((R-mean(R))/sd(R),mean = 0,sd=1,lower.tail = TRUE,log.p = FALSE)
        logp1 <- log((seq(1:n)-0.5)/(n*p))
        logp2 <- log((n-seq(1:n)+0.5)/(n*(1-p)))
        h<-(seq(1:n)-0.5)*logp1+(n-seq(1:n)+0.5)*logp2
        zk <- max(h)
        S <- S + (zk > ZK)
      }
    }
  }
  if(distr=="gamma"){
    for (j in 1:N) {
      R<-rgamma(n,shape = a,scale = b)
      R<-sort(R)
      p<-pgamma(R,shape = a,scale = b,lower.tail = TRUE, log.p = FALSE)
      logp1 <- log((seq(1:n)-0.5)/(n*p))
      logp2 <- log((n-seq(1:n)+0.5)/(n*(1-p)))
      h<-(seq(1:n)-0.5)*logp1+(n-seq(1:n)+0.5)*logp2
      zk <- max(h)
      S <- S + (zk > ZK)
    }
  }
  if(distr=="weibull"){
    for (j in 1:N) {
      R<-rweibull(n,shape = a,scale = b)
      R<-sort(R)
      p<-pweibull(R,shape = a,scale = b,lower.tail = TRUE, log.p = FALSE)
      logp1 <- log((seq(1:n)-0.5)/(n*p))
      logp2 <- log((n-seq(1:n)+0.5)/(n*(1-p)))
      h<-(seq(1:n)-0.5)*logp1+(n-seq(1:n)+0.5)*logp2
      zk <- max(h)
      S <- S + (zk > ZK)
    }
  }
  if(distr=="t"){
    for (j in 1:N) {
      R<-rt(n,df = a)
      R<-sort(R)
      p<-pt(R,df = a,lower.tail = TRUE, log.p = FALSE)
      logp1 <- log((seq(1:n)-0.5)/(n*p))
      logp2 <- log((n-seq(1:n)+0.5)/(n*(1-p)))
      h<-(seq(1:n)-0.5)*logp1+(n-seq(1:n)+0.5)*logp2
      zk <- max(h)
      S <- S + (zk > ZK)
    }
  }
  pval<-S/N
  RVAL<-list(statistic = c(ZK = ZK),p.value = pval,method = "ZK test for given distribution",data.name = DNAME)
  class(RVAL) <- "htest"    
  return(RVAL)
}
