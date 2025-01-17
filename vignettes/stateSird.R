## ----load_libraries-----------------------------------------------------------
# install.packages("covidStateSird", repos = "http://cran.us.r-project.org")
# install.packages("rjags", repos = "http://cran.us.r-project.org")
# library(rjags)
# install.packages("randomForest", repo="http://cran.r-project.org", dep=T)
library(randomForest)
library(covidStateSird)
library(foreach)
library(dplyr)
# install.packages("ggplot2", repo="http://cran.r-project.org", dep=T)
# install.packages("forecast", repo="http://cran.r-project.org", dep=T)
# install.packages("astsa", repo="http://cran.r-project.org", dep=T)
# install.packages("nlme", repo="http://cran.r-project.org", dep=T)
library(ggplot2)
library(forecast)
library(astsa)
library(nlme)


# the last day of data to use
endDate <- "2021-03-07"
minCase <- 100

set.seed(525600)

covidDir <- ".."

doParallel::registerDoParallel(cores=5)

outputPath <- file.path(covidDir, "Output", Sys.Date())
dir.create(outputPath)
dir.create(file.path(outputPath, "Tables/"))
dir.create(file.path(outputPath, "Plots/"))
dir.create(file.path(outputPath, "Data/"))

stateInterventions <- read.csv(paste0(covidDir, "/Data/StateInterventionDates.csv"),
                               stringsAsFactors = FALSE,
                               header = TRUE)

stateDataCSV <- "https://raw.githubusercontent.com/COVID19Tracking/covid-public-api/master/v1/states/daily.csv"

vaccineData <- read.csv(paste0(covidDir, "/Data/us_state_vaccinations.csv"),
                               stringsAsFactors = FALSE,
                               header = TRUE)
stateCovidData <- read.csv(stateDataCSV)

covariates <- read.csv(paste0(covidDir, "/Data/COVID-19 Location Covariates - analyticStates.csv"))
covariates[, c("heartDisease", "lungCancer", "diabetes", "copd")] <-
  covariates[, c("heartDisease", "lungCancer", "diabetes", "copd")] / 100000

states <- c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA",
            "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD",
            "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH",
            "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC",
            "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")

velocLogCases <- velocLogDeaths <- data.frame()
loc <- 0
for(i in 1:length(states)) {
  loc <- loc + 1
  # vacc <- vaccineData[which(vaccineData$location == stateAbbrev),]
  # vacc <- vacc$
  velocLoc <- velocitiesState(stateCovidData, vaccineData, states[i], stateInterventions, minCases = minCase, endDate = endDate)
  population <- stateInterventions$statePopulation[stateInterventions$stateAbbreviation == states[i]]
  
  velocLogCases <- rbind(velocLogCases, cbind(velocLoc$cases,  loc, row.names = NULL))
  # print(colnames(velocLogCases))
}
print(head(velocLogCases))
velocLogCasesList <- as.list(velocLogCases)
velocLogCasesList$N <- nrow(velocLogCases)
velocLogCasesList$nLoc <- length(unique(velocLogCasesList$loc))
velocLogCasesList$y[velocLogCasesList$y <= 0] <- NA
print(names(velocLogCasesList))
params <- c("mu_a", "mu_b", "tau", "a", "b", "g", "d", "mu_g", "mu_d", "mu", "alpha", "beta", "mu_alpha", "mu_beta")
 
start <- Sys.time()

velocModel <- R2jags::jags(data = velocLogCasesList, inits = NULL, parameters.to.save = params,
  model.file = riasJagsModel, n.chains = 3, n.iter = 100, n.burnin = 10, n.thin = 10, DIC = F)

timeElapsed <- (Sys.time() - start)

plotVelocityFit(velocModel$BUGSoutput$mean,
                stateCovidData, vaccineData, states, stateInterventions,
                fileName = paste0(outputPath, "/Plots/velocityModelFit.pdf"))

c_S_pre <- caseModelConstant(velocModel, intervention = 0)
c_S_post <- caseModelConstant(velocModel, intervention = 1)

rownames(c_S_pre) <- NULL
rownames(c_S_post) <- NULL

posteriorSamples <- velocModel$BUGSoutput$sims.list

posteriorSamples[["c_pre"]] <- c_S_pre
posteriorSamples[["c_post"]] <- c_S_post

save(posteriorSamples,
  file = paste0(outputPath, "/CasePosteriorSamples", endDate, ".Rdata"))

# run state SIRD models
randomForestDeathModel <- deathForest(stateCovidData, vaccineData, states, covariates, 21, fileOut = paste0(outputPath, "/randomForestDeathModel.Rdata"))

load(paste0(outputPath, "/CasePosteriorSamples", endDate, ".Rdata"))

# stateSird("CA", covariates, stateInterventions, stateCovidData, randomForestDeathModel,
# posteriorSamples, rfError = T, plots = T)

# foreach(i = 1:length(states)) %dopar% {
#   stateSird(states[i], covariates, stateInterventions, stateCovidData, randomForestDeathModel,
#   posteriorSamples, rfError = T)
# }
 for(i in 1:length(states)) {
   stateSird(states[i], covariates, stateInterventions, stateCovidData, vaccineData, randomForestDeathModel,
   posteriorSamples, rfError = T)
 }

