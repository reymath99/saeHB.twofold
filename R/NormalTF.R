#' @title Small Area Estimation Using Hierarchical Bayesian Method under Twofold Subarea Level Model with Normal distribution
#'
#' @description
#' \itemize{
#'  \item {This function is implemented to variable of interest \eqn{y} that assumed to be a Normal Distribution. The range of data is \eqn{-\infty <y<\infty}}
#'  \item {This function gives estimation of subarea and area means simultaneously under Twofold Subarea Level Small Area Estimation Model Using Hierarchical Bayesian Method with Normal distribution}
#' }
#'
#' @param formula Formula that describe the fitted model
#' @param vardir Sampling variances of direct estimations on each subarea
#' @param area Index that describes the code relating to area in each subarea.This should be defined for aggregation to get area estimator. Index start from 1 until m
#' @param weight Vector contain proportion units or proportion of population on each subarea. \eqn{w_{ij}}
#' @param iter.update Number of updates perform in Gibbs Sampling with default \code{3}
#' @param iter.mcmc Number of total iteration per chain perform in Gibbs Sampling with default \code{2000}
#' @param thin Thinning rate perform in Gibbs Sampling and it must be a positive integer with default \code{1}
#' @param burn.in Number of burn in period in Gibbs Sampling with default \code{1000}
#' @param data The data frame
#' @param coef Initial value for mean on coefficient's prior distribution or \eqn{\beta}'s prior distribution
#' @param var.coef Initial value for varians on coefficient's prior distribution or \eqn{\beta}'s prior distribution
#'
#' @return This function returns a list with following objects:
#' \describe{
#'   \item{Est_sub}{A dataframe that contains the values, standar deviation, and quantile of Subarea mean Estimates using Twofold Subarea level model under Hierarchical Bayes method}
#'   \item{Est_area}{A dataframe that contains the values, standar deviation, and quantile of Area mean Estimates using Twofold Subarea level model under Hierarchical Bayes method}
#'   \item{refVar}{A dataframe that contains estimated subarea and area random effect variance \eqn{(\sigma_{u}^{2}} and \eqn{\sigma_{v}^{2})}}
#'   \item{coefficient}{A dataframe that contains the estimated model coefficient \eqn{\beta}}
#'   \item{plot}{Trace, Density, Autocorrelation Function Plot of coefficient}
#' }
#'
#' @import stringr
#' @import coda
#' @import rjags
#' @import stats
#' @import grDevices
#' @import graphics
#' @import data.table
#' @import utils
#'
#' @export NormalTF
#'
#' @examples
#' ##load dataset for data without any nonsampled subarea
#' data(dataTwofold)
#'
#' #formula of fitted model
#' formula=y~x1+x2
#'
#' #model fitting
#' mod=NormalTF(formula,vardir="vardir",area="codearea",weight="w",data=dataTwofold)
#'
#' #estimate
#' mod$Est_sub #Subarea mean estimate
#' mod$Est_area #area mean estimate
#' mod$coefficient #coefficient estimate
#' mod$refVar #random effect subarea and area estimates
#'
#' ##for dataset with nonsampled subarea use dataTwofoldNS
NormalTF<-function (formula, vardir, area, weight,iter.update = 3, iter.mcmc = 2000, thin = 1, burn.in = 1000, data,coef,var.coef){
  opar <- par("mar")
  on.exit(par(opar))
  result <- list()
  formuladata <- model.frame(formula, data, na.action = NULL)
  if (any(is.na(formuladata[, -1]))) {
    stop("Auxiliary Variables contains NA values")
  }
  auxVar <- as.matrix(formuladata[, -1])
  p<-ncol(auxVar)
  nvar <- p + 1
  formuladata <- data.frame(formuladata, vardir = data[,vardir],codearea=data[,area], weight=data[,weight])
  if (iter.update < 3){
    stop("the number of iteration updates at least 3 times")
  }
  for (i in 1:nrow(formuladata)) {
    if (!is.na(formuladata[i, 1])) {
      if (is.na(formuladata[i, (nvar + 1)])) {
        stop(formula[2], "[", i, "] is not NA but vardir is NA")
      }else if(is.na(formuladata[i, (nvar + 2)])){
        stop(formula[2], "[", i, "] is not NA but codearea is NA")
      }else if(is.na(formuladata[i, (nvar + 3)])){
        stop(formula[2], "[", i, "] is not NA but weight is NA")
      }
    }
  }

  if (!missing(var.coef)){
    if( length(var.coef) != nvar ){
      stop("length of vector var.coef does not match the number of regression coefficients, the length must be ",nvar)
    }
    tau.b.value = 1/var.coef
  }else{
    tau.b.value = 1/rep(1,nvar)
  }

  if (!missing(coef)){
    if( length(coef) != nvar ){
      stop("length of vector coef does not match the number of regression coefficients, the length must be ",nvar)
    }
    mu.b.value = coef
  } else {
    mu.b.value = rep(0,nvar)
  }

  if (!any(is.na(formuladata[, 1]))) {
    formuladata <- as.matrix(na.omit(formuladata))
    x <- model.matrix(formula, data = as.data.frame(formuladata))
    n <- nrow(formuladata)
    m <- length(unique(formuladata[,(nvar+2)]))
    for (i in 1:n) {
      if (formuladata[i, (nvar + 1)] == 0) {
        stop("Vardir for ", formula[2], "[", i, "] must not be 0")
      }
    }
    mu.b=mu.b.value
    tau.b=tau.b.value
    tau.ua=tau.ub=tau.va=tau.vb=var_area=var_sub=1
    Iter=iter.update
    for(i in 1:Iter){
      mydata <- list(n = n,
                     m = m,
                     nvar = nvar,
                     y = formuladata[, 1],
                     x = as.matrix(x[, -1]),
                     mu.b = mu.b,
                     tau.b = tau.b,
                     tau.ua = tau.ua, tau.ub = tau.ub,
                     tau.va = tau.va, tau.vb = tau.vb,
                     vardir = formuladata[,(nvar + 1)],
                     state2 = formuladata[,(nvar + 2)]
      )
      inits<-list(u = rep(0,n), v=rep(0,m), beta =mu.b, sigma_u2=1, sigma_v2=1)
      cat(
        "model {
                  #N observations
                  for (i in 1:n) {
                    y[i] ~ dnorm(theta[i], tau[i])
                    tau[i]<- 1/vardir[i]
                    theta[i]<- beta[1] + inprod(beta[2:nvar],x[i,]) + v[state2[i]]+u[i]
                    u[i]~dnorm(0, sigma_u2)
                  }

                  # M states
                  for (j in 1:m) {
                    v[j] ~ dnorm(0,sigma_v2)
                  }

                  # Priors for beta
                  for (k in 1:nvar){
                    beta[k]~dnorm(mu.b[k],tau.b[k])
                  }

                  # Priors for variance over area and subarea

                  sigma_u2 ~ dgamma(tau.ua, tau.ub)
                  sigma_v2 ~ dgamma(tau.va, tau.vb)
                  var_sub <- 1/sigma_u2
                  var_area <- 1/sigma_v2
                  }",file="TwoFoldNormalRev8.txt"
      )
      jags.m <- jags.model(file = "TwoFoldNormalRev8.txt",
                           data = mydata, inits = inits,
                           n.chains = 1, n.adapt = 500)
      file.remove("TwoFoldNormalRev8.txt")
      params <- c("theta","beta","sigma_u2","sigma_v2","var_sub","var_area")
      samps <- coda.samples(jags.m, params, n.iter = iter.mcmc, thin = thin)
      samps1 <- window(samps, start = burn.in + 1, end = iter.mcmc)
      result_samps = summary(samps1)
      var_area= tail(result_samps$statistics,2)[1,1]
      var_sub= tail(result_samps$statistics,2)[2,1]

      beta = result_samps$statistics[1:(nvar), 1:2]
      for (i in 1:nvar) {
        mu.b[i] = beta[i, 1]
        tau.b[i] = 1/(beta[i, 2]^2)
      }

      taures = result_samps$statistics[(nvar+1):(nvar+2), 1:2]
      tau.ua = taures[1, 1]^2/taures[1, 2]^2
      tau.ub = taures[1, 1]/taures[1, 2]^2
      tau.va = taures[2, 1]^2/taures[2, 2]^2
      tau.vb = taures[2, 1]/taures[2, 2]^2

    }

    result_samps <- summary(samps1)
    b.varnames <- list()
    for (i in 1:(nvar)) {
      idx.b.varnames <- as.character(i - 1)
      b.varnames[i] <- str_replace_all(paste("beta[",
                                             idx.b.varnames, "]"), pattern = " ",
                                       replacement = "")
    }
    result_mcmc <- samps1[, c(1:(nvar))]
    colnames(result_mcmc[[1]]) <- b.varnames
    var_area<- tail(result_samps$statistics,2)[1,1]
    var_sub<- tail(result_samps$statistics,2)[2,1]
    refVari<-data.frame(var_area,var_sub)
    beta <- result_samps$statistics[1:(nvar), 1:2]
    rownames(beta) <- b.varnames
    theta <- result_samps$statistics[(nvar+3):(nvar+3+n-1),1:2]
    Estimation1 <- data.frame(theta)
    colnames(Estimation1) <- c("mean", "sd")
    Quantiles <- as.data.frame(result_samps$quantiles[c( (1:(nvar)),((nvar+3):(nvar+3+n-1)) ),])
    q_beta <- Quantiles[1:(nvar), ]
    q_mu<-Quantiles[-c(1:nvar),]
    rownames(q_beta) <- b.varnames
    beta <- cbind(beta, q_beta)
    Estimation <- data.frame(Estimation1,q_mu)
    colnames(Estimation) <- c("Mean","SD","2.5%","25%","50%","75%","97.5%")
    w<-gr<-0
    Estimation_area<-data.frame(cbind(Estimation1,w=data[,weight],code=data[,area]))
    Estimation_area[,"wm"]=Estimation_area$w*Estimation_area$mean
    Estimation_area[,"ws"]=(Estimation_area$w)^2*(Estimation_area$sd)^2
    Est_area<-aggregate(wm~code,data=Estimation_area,FUN = sum)
    sdarea<-aggregate(ws~code,data=Estimation_area,FUN=sum)
    result_mcmc_area<-data.frame(t(samps1[[1]][, c((nvar+3):(nvar+3+n-1))]))
    dtarea<-data.table(result_mcmc_area,w=data[,weight],gr=data[,area])
    Quantilesdt<-dtarea[, lapply(.SD, function(x, w) sum(x*w), w=w), by=gr][, w := NULL]
    Quantiles_area<-apply(Quantilesdt[,-1],MARGIN = 1,FUN=function(x){
      quantile(x,probs = c(0.025,0.25,0.50,0.75,0.975))
    })
    Est_area2<-data.frame(cbind(Est_area$wm,sqrt(sdarea$ws),t(Quantiles_area)))
    colnames(Est_area2)<-c("Mean","SD","2.5%","25%","50%","75%","97.5%")
  }else{

    formuladata <- as.data.frame(formuladata)
    x <- as.matrix(formuladata[, 2:nvar])
    n <- nrow(formuladata)
    m <- length(unique(formuladata[,(nvar+2)]))
    formuladata$idx <- rep(1:n)
    data_sampled <- na.omit(formuladata)
    data_nonsampled <- formuladata[-data_sampled$idx, ]
    r = data_nonsampled$idx
    n1 = nrow(data_sampled)
    n2 = nrow(data_nonsampled)

    for (i in 1:n1) {
      if (data_sampled[i, (nvar + 1)] == 0) {
        stop("Vardir for ", formula[2], "[",
             data_sampled$idx[i], "] must not be 0")
      }
    }

    mu.b=mu.b.value
    tau.b=tau.b.value
    tau.ua=tau.ub=tau.va=tau.vb=var_area=var_sub=1
    Iter=iter.update

    for(i in 1:Iter){
      mydata <- list(n1 = n1,
                     n2 = n2,
                     m = m,
                     nvar = nvar,
                     y_sampled = data_sampled[, 1],
                     x_sampled = data_sampled[,2:nvar],
                     x_nonsampled=data_nonsampled[,2:nvar],
                     mu.b = mu.b,
                     tau.b = tau.b,
                     tau.ua = tau.ua, tau.ub = tau.ub,
                     tau.va = tau.va, tau.vb = tau.vb,
                     vardir = data_sampled[,(nvar + 1)],
                     state2_sampled = data_sampled[,(nvar + 2)],
                     state2_nonsampled=data_nonsampled[,(nvar+2)]
      )
      inits<-list(u1 = rep(0,n1),u2=rep(0,n2), v=rep(0,m), beta =mu.b, sigma_u2=1, sigma_v2=1)
      cat(
        "model {
                  #N observations sampled
                  for (i in 1:n1) {
                    y_sampled[i] ~ dnorm(theta[i], tau[i])
                    tau[i]<- 1/vardir[i]
                    theta[i]<- beta[1] + inprod(beta[2:nvar],x_sampled[i,]) + v[state2_sampled[i]]+u1[i]
                    u1[i]~dnorm(0, sigma_u2)
                  }

                  #N observations non_sampled y=xb+vi (Torabi,2014)
                  for (j in 1:n2) {
                    theta.nonsampled[j] <- mu.b[1] + inprod(mu.b[2:nvar],x_nonsampled[j,]) + v[state2_nonsampled[j]]+u2[j]
                    u2[j]~dnorm(0,sigma_u2)
                  }


                  # M states
                  for (l in 1:m) {
                    v[l] ~ dnorm(0,sigma_v2)
                  }

                  # Priors for beta
                  for (k in 1:nvar){
                    beta[k]~dnorm(mu.b[k],tau.b[k])
                  }

                  # Priors for variance over area and subarea

                  sigma_u2 ~ dgamma(tau.ua, tau.ub)
                  sigma_v2 ~ dgamma(tau.va, tau.vb)
                  var_sub <- 1/sigma_u2
                  var_area <- 1/sigma_v2
                  }",file="TwoFoldNormalRev8.txt"
      )
      jags.m <- jags.model(file = "TwoFoldNormalRev8.txt",
                           data = mydata, inits = inits,
                           n.chains = 1, n.adapt = 500)
      file.remove("TwoFoldNormalRev8.txt")
      params <- c("theta","theta.nonsampled","beta","sigma_u2","sigma_v2","var_sub","var_area")
      samps <- coda.samples(jags.m, params, n.iter = iter.mcmc, thin = thin)
      samps1 <- window(samps, start = burn.in + 1, end = iter.mcmc)
      result_samps = summary(samps1)
      var_area= tail(result_samps$statistics,2)[1,1]
      var_sub= tail(result_samps$statistics,2)[2,1]
      beta = result_samps$statistics[1:(nvar), 1:2]
      for (i in 1:nvar) {
        mu.b[i] = beta[i, 1]
        tau.b[i] = 1/(beta[i, 2]^2)
      }
      taures = result_samps$statistics[(nvar+1):(nvar+2), 1:2]
      tau.ua = taures[1, 1]^2/taures[1, 2]^2
      tau.ub = taures[1, 1]/taures[1, 2]^2
      tau.va = taures[2, 1]^2/taures[2, 2]^2
      tau.vb = taures[2, 1]/taures[2, 2]^2
    }
    result_samps <- summary(samps1)
    b.varnames <- list()
    for (i in 1:(nvar)) {
      idx.b.varnames <- as.character(i - 1)
      b.varnames[i] <- str_replace_all(paste("beta[",
                                             idx.b.varnames, "]"), pattern = " ",
                                       replacement = "")
    }
    result_mcmc <- samps1[, c(1:(nvar))]
    colnames(result_mcmc[[1]]) <- b.varnames
    var_area<- tail(result_samps$statistics,2)[1,1]
    var_sub<- tail(result_samps$statistics,2)[2,1]
    refVari<-data.frame(var_area,var_sub)
    beta <- result_samps$statistics[1:(nvar), 1:2]
    rownames(beta) <- b.varnames
    theta <- result_samps$statistics[(nvar+3):(nvar+3+n1-1),1:2]
    theta_ns<-result_samps$statistics[(nvar+3+n1):(nvar+3+n-1),1:2]
    Estimation1 <- matrix(rep(0,n),n,2)
    Estimation1[r,]<-theta_ns
    Estimation1[-r,]<-theta
    Estimation1 <- as.data.frame(Estimation1)
    colnames(Estimation1) <- c("mean", "sd")
    Quantiles <- as.data.frame(result_samps$quantiles)
    q_beta <- Quantiles[1:(nvar), ]
    q_mu<-Quantiles[(nvar+3):(nvar+3+n1-1),]
    q_mu.nonsampled<-Quantiles[(nvar+3+n1):(nvar+3+n-1),]
    q_Estimation <- matrix(0,n,5)
    for (i in 1:5){
      q_Estimation[r,i] <- q_mu.nonsampled[,i]
      q_Estimation[-r,i] <- q_mu[,i]
    }
    rownames(q_beta) <- b.varnames
    beta <- cbind(beta, q_beta)
    Estimation <- data.frame(Estimation1,q_Estimation)
    colnames(Estimation) <- c("Mean","SD","2.5%","25%","50%","75%","97.5%")
    w<-gr<-0
    Estimation_area<-data.frame(cbind(Estimation1,w=data[,weight],code=data[,area]))
    Estimation_area[,"wm"]=Estimation_area$w*Estimation_area$mean
    Estimation_area[,"ws"]=(Estimation_area$w)^2*(Estimation_area$sd)^2
    Est_area<-aggregate(wm~code,data=Estimation_area,FUN = sum)
    sdarea<-aggregate(ws~code,data=Estimation_area,FUN=sum)
    result_mcmc_area_s<-data.frame(t(samps1[[1]][, c((nvar+3):(nvar+3+n1-1))]))
    result_mcmc_area_ns<-data.frame(t(samps1[[1]][, c((nvar+3+n1):(nvar+3+n-1))]))
    result_mcmc_area <- matrix(0,n,ncol(result_mcmc_area_s))
    for(i in 1:ncol(result_mcmc_area_s)){
      result_mcmc_area[r,i]<-result_mcmc_area_ns[,i]
      result_mcmc_area[-r,i]<-result_mcmc_area_s[,i]
    }
    result_mcmc_area<-data.frame(result_mcmc_area)
    dtarea<-data.table(result_mcmc_area,w=data[,weight],gr=data[,area])
    Quantilesdt<-dtarea[, lapply(.SD, function(x, w) sum(x*w), w=w), by=gr][, w := NULL]
    Quantiles_area<-apply(Quantilesdt[,-1],MARGIN = 1,FUN=function(x){
      quantile(x,probs = c(0.025,0.25,0.50,0.75,0.975))})
    Est_area2<-data.frame(cbind(Est_area$wm,sqrt(sdarea$ws),t(Quantiles_area)))
    colnames(Est_area2)<-c("Mean","SD","2.5%","25%","50%","75%","97.5%")
  }

  result$Est_sub = Estimation
  result$Est_area = Est_area2
  result$refvar = refVari
  result$coefficient = beta
  result$plot = list(graphics.off(), par(mar = c(2, 2, 2, 2)),
                     autocorr.plot(result_mcmc, col = "brown2", lwd = 2),
                     plot(result_mcmc, col = "brown2", lwd = 2))
  return(result)
}
