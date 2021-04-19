#' @export
velocitiesState <- function(statesLong, vaccinesLong, stateAbbrev, stateInterventions = NULL, minCases = 100, endDate = "2099-01-01") {

  if(!is.null(stateInterventions)) {
    intervDate <- stateInterventions$interventionDate[stateInterventions$stateAbbreviation == stateAbbrev]
    if(intervDate == "") intervDate <- "2099-01-01"
  }
  stateLong <- statesLong[which(statesLong$state == stateAbbrev),]
  # print(head(stateLong$date))
  stateLong <- stateLong[which(stateLong$date %in% vaccinesLong$date), ]
  # print(stateLong$date)
  # print(vaccinesLong)
  vaccineLong <- vaccinesLong[which(vaccinesLong$location == stateAbbrev),]
  vaccineLong <- vaccineLong[which(vaccineLong$date %in% stateLong$date), ]
  # n.day <- nrow(stateLong)
  n.day <- nrow(vaccineLong)
  # print(nrow(vaccineLong))
  # print(nrow(stateLong))
  daysAll <- NULL
  # first_y <- as.numeric(substr(stateLong$date[i],1,4))
  # first_m <- as.numeric(substr(stateLong$date[i],5,6))
  # first_d <- as.numeric(substr(stateLong$date[i],7,8))
  for(i in 1:n.day) {
    daysAll[i] <-(paste(substr(stateLong$date[i],1,4),
                        substr(stateLong$date[i],5,6),
                        substr(stateLong$date[i],7,8), sep = "-") )
  }
  daysAll <- as.Date(daysAll)

  stateLong <- stateLong[order(daysAll, decreasing = F),]
  # print(colnames(stateLong))
  vaccineLong <- vaccineLong[order(daysAll, decreasing = F),]
  daysAll   <- daysAll[order(daysAll, decreasing = F)]
        
  stateAll <- cbind(stateLong$positive, stateLong$death, stateLong$hospitalizedCurrently, vaccineLong$people_vaccinated)
  # print(stateAll)
  # ts_pos <- ts(stateLong[,1], frequency=365, start=1)
  # ts_death <- ts(stateLong[,2], frequency=365, start=1)
  # ts_hosp <- ts(stateLong[,3], frequency=365, start=1)
  # ts_vacc <- ts(stateLong[,4], frequency=365, start=1)
  state <- stateAll[which(stateAll[,1] >= minCases),]
  # print(state)
  days  <- as.Date(daysAll[which(stateAll[,1] >= minCases)])
  
  postIntervention <- 1 * (days > intervDate)
  y <- log(state[!is.na(state[1,]) & (days <= endDate), 1])
  x <- which(!is.na(state[1,]) & (days <= endDate))
  splineCases <- smooth.spline(x = x[y > - Inf], y = y[y > -Inf])
  derivCases <- predict(splineCases, deriv = 1)
  # print(derivCases)
  y <- log(state[!is.na(state[,1]) & !is.na(state[,4]) & (days <= endDate), 1])
  # print(y)
  x <- seq(1, length(y), by=1)
  # print(x)
  v <- log(state[!is.na(state[,1]) & !is.na(state[,4]) & (days <= endDate), 4])
  # print(v)
  ys <- diff(y)/diff(x)
  vs <- diff(v)/diff(x)
  # print(lapply(as.numeric(vs), as.integer))
  # x <- seq(1, length(ycs.prime), by=1)
  ts_v <- ts(vs, frequency=365, start=1)
  ts_y <- ts(ys, frequency=365, start=1)
  xreg <- as.matrix(ts_v)
  # print(xreg)
  # print(ycs.prime)
  fit_basic1 <- auto.arima(ts_y, xreg=xreg)
  # print(fit_basic1)
  # checkresiduals(fit_basic1)
  # print(as.data.frame(fitted(fit_basic1)))
  # forecast_1 <- forecast(fit_basic1,xreg = xreg)
  # print(forecast_1)
  # splineCases <- smooth.spline(x = x[y > - Inf], y = y[y > -Inf])
  # derivCases <- predict(splineCases, deriv = 1)
  # print(derivCases)

  # stateAll <- rbind(stateLong$positive, stateLong$death, stateLong$hospitalizedCurrently)
  # state <- stateAll[,which(stateAll[1,] >= minCases)]
  # days  <- as.Date(daysAll[which(stateAll[1,] >= minCases)])
  
  # postIntervention <- 1 * (days > intervDate)
  
  # y <- log(state[1, !is.na(state[1,]) & (days <= endDate)])
  # x <- which(!is.na(state[1,]) & (days <= endDate))
  # splineCases <- smooth.spline(x = x[y > - Inf], y = y[y > -Inf])
  # derivCases <- predict(splineCases, deriv = 1)
  # print("derivCases")
  # print(derivCases)

  ts_y <- ts(derivCases$y, frequency=365, start=1)
  x <- seq(1, nrow(as.data.frame(fitted.values(fit_basic1))), by=1)
  plot(fit_basic1$x,col="red", ylab = "d/dt Log Cumulative Cases",   xlab = "Days Since 100+ Cases")
  lines(fitted(fit_basic1),col="blue")
  lines(ts_y, col="green")

  
  yd <- log(state[!is.na(state[,2]) & !is.na(state[,4]) & (days <= endDate), 1])
  # print(y)
  xd <- seq(1, length(yd), by=1)
  # print(x)
  vd <- log(state[!is.na(state[,2]) & !is.na(state[,4]) & (days <= endDate), 4])
  # print(v)
  yds <- diff(yd)/diff(xd)
  vds <- diff(vd)/diff(xd)
  # x <- seq(1, length(ycs.prime), by=1)
  ts_vd <- ts(vds, frequency=365, start=1)
  ts_yd <- ts(yds, frequency=365, start=1)
  xregd <- as.matrix(ts_vd)
  # print(xreg)
  # print(ycs.prime)
  fit_basic2 <- auto.arima(ts_yd, xreg=xregd)
  # print(fit_basic2)
  # checkresiduals(fit_basic1)
  xd <- seq(1, nrow(as.data.frame(fitted.values(fit_basic2))), by=1)
  # print(as.data.frame(fitted.values(fit_basic2)))

  # y <- log(state[2, !is.na(state[2,]) & (days <= endDate)])
  # x <- which(!is.na(state[2,]) & (days <= endDate))
  # splineDeaths <- smooth.spline(x = x[y > - Inf], y = y[y > -Inf])
  # derivDeaths <- predict(splineDeaths, deriv = 1)

  return(list(cases = data.frame(y = as.numeric(fitted(fit_basic1)),
                                 t = x,
                                 u = state[x, 1],
                                 postIntervention = postIntervention[x],
                                 deaths = state[x,2]),
              deaths = data.frame(y = as.numeric(fitted(fit_basic2)),
                                  t = xd,
                                  u = state[xd, 2],
                                  hosps = state[xd, 3],
                                 postIntervention = postIntervention[xd]),
              days = days,
              intervDate = intervDate))
}

#' @export
riasJagsModel <- function(){
  # Likelihood:
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], tau[i])

    mu[i] <-  a[loc[i]] + b[loc[i]] * t[i] + (g[loc[i]] + d[loc[i]] * t[i]) * postIntervention[i] 
    # + v[(vacc[i])*(vacc[i])] + w[(vacc[i])*(vacc[i])] * t[i]
    tau[i] <- exp(alpha[loc[i]] + beta[loc[i]] * t[i])
  }
  
  # Priors:
  for(j in 1:nLoc) {
    a[j] ~ dnorm(mu_a, 1) # random intercept for location
    b[j] ~ dnorm(mu_b, 10) # random slope for location
    g[j] ~ dnorm(mu_g, 1) # random effect for post-intervention
    d[j] ~ dnorm(mu_d, 10) # random slope for post-intervention
    v[j] ~ dnorm(mu_v, 1) # random effect for vacc
    w[j] ~ dnorm(mu_w, 10) # random slope for vacc
    alpha[j] ~ dnorm(mu_alpha, 1)
       beta[j]  ~ dnorm(mu_beta, 10)
  }
  mu_a  ~ dnorm(-1, 1) # random intercept mean
  mu_b  ~ dnorm(0, 10) # random slope mean
  mu_g  ~ dnorm(0, 1) # random effect mean
  mu_d  ~ dnorm(-.05, 10) # random slope mean
  mu_v  ~ dnorm(0, 1) # random vacc mean
  mu_w  ~ dnorm(-.05, 10) # random vacc slope mean
  mu_alpha ~ dnorm(0,1)
  mu_beta  ~ dnorm(0, 10)
}

#' @export
constErrLogNorm <- function(x, t, u, beta, alpha) {
  err <- (u - exp((1/c(beta)) * exp(c(beta) * t + c(alpha)) + c(x))) ^2
  return(sum(err))
}

#' @export
#' @importFrom foreach %dopar% %:% 
caseModelConstant <- function(velocityPosterior, intervention = 1) {
  if(intervention == 0) {
    beta = velocityPosterior$BUGSoutput$sims.list$b
    alpha = velocityPosterior$BUGSoutput$sims.list$a
  } else if(intervention == 1) {
    beta = velocityPosterior$BUGSoutput$sims.list$b +
           velocityPosterior$BUGSoutput$sims.list$d
    alpha = velocityPosterior$BUGSoutput$sims.list$a +
            velocityPosterior$BUGSoutput$sims.list$g
  }
  
  constants <- foreach::foreach(k = 1:nrow(velocityPosterior$BUGSoutput$sims.matrix), .combine = 'rbind') %:%
     foreach::foreach(i = 1:velocLogCasesList$nLoc, .combine = 'c') %dopar% {
       optim(6, constErrLogNorm, u = velocLogCases$u[(velocLogCases$loc == i) & (velocLogCases$postIntervention == intervention)],
            t = velocLogCases$t[(velocLogCases$loc == i) & (velocLogCases$postIntervention == intervention)],
         beta = beta[k,i],
         alpha = alpha[k,i],
       method = "SANN")$par
    }
  return(constants)
}
