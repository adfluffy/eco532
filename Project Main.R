################################################################################
# Final Project
#     Author: Devan Arnold
#       Date: 4/23/24
#     Course:
#       ECO 532 - Applied Time Series Econometrics
#       Dr. Thomas Wiesen
#       University of Maine, Spring 2024
################################################################################

# Import the relevant libraries for time series data analysis
library(readxl) 
library(tseries) 
library(astsa)
library(forecast)
library(modelsummary)
library(urca)

################################################################################
# Sets the working directory to the parent of the current file location
## Choose the current working file
this_script <- file.choose()

## Determines if the chosen file path uses the Windows or Mac style file path
mac_pc_test <- grepl("/",this_script,fixed = F)

## Break the file path into the parent path and the child path based on the
##  location of the 'scripts' project folder. The logic checks to see if the
##  file path is from a Mac or Windows system and adjusts accordingly
if (mac_pc_test) {
  print("Mac True Selected")
  this_script_broken <- unlist(strsplit(this_script,"scripts"))
} else {
  print("Mac False Selected")
  this_script_broken <- unlist(strsplit(this_script,"scripts"))
  
}

## Use the above split to determine the parent file path of the script
script_parent <- this_script_broken[1]

## Sets the working directory based on the above logic
setwd(script_parent)

################################################################################
# Load in GDP data and store in master data frame object
nom_gdp_data = read_excel(file.choose())

# Attach the file for more ease of data calling for the next few lines
attach(nom_gdp_data)
gdp_nobs <- nrow(nom_gdp_data)
# Create time series objects for GDP (Y), Consumption (C), Investment (I), Government
#   Expenditures (G), and Net Exports (NX)
Y_nom_ts <- ts(data = `Gross Domestic Product (Y)`,start = c(1947,1),end = c(2023,4),frequency = 4)
C_nom_ts <- ts(data = `Personal Consumption Expenditures (C)`,start = c(1947,1),end = c(2023,4),frequency = 4)
I_nom_ts <- ts(data = `Gross Private Domestic Investment (I)`,start = c(1947,1),end = c(2023,4),frequency = 4)
G_nom_ts <- ts(data = `Government Consumption Expenditures and Gross Investment (G)`,start = c(1947,1),end = c(2023,4),frequency = 4)
NX_nom_ts <- ts(data = `Net Exports (NX)`,start = c(1947,1),end = c(2023,4),frequency = 4)


# Load in Real GDP data and store in master data frame object
real_gdp_data = read_excel(file.choose())

# Attach the file for more ease of data calling for the next few lines
attach(real_gdp_data)
real_gdp_nobs <- nrow(real_gdp_data)
# Create time series objects for GDP (Y), Consumption (C), Investment (I), Government
#   Expenditures (G), Net Exports (NX), and the residuals from converting the
#   GDP components to 2017 dollars
Y_real_ts <- ts(data = `Gross Domestic Product (Y)`,start = c(1947,1),end = c(2023,4),frequency = 4)
C_real_ts <- ts(data = `Personal Consumption Expenditures (C)`,start = c(1947,1),end = c(2023,4),frequency = 4)
I_real_ts <- ts(data = `Gross Private Domestic Investment (I)`,start = c(1947,1),end = c(2023,4),frequency = 4)
G_real_ts <- ts(data = `Government Consumption Expenditures and Gross Investment (G)`,start = c(1947,1),end = c(2023,4),frequency = 4)
NX_real_ts <- ts(data = `Net Exports (NX)`,start = c(1947,1),end = c(2023,4),frequency = 4)
chained_residual_ts <- ts(data = Residual,start = c(1947,1),end = c(2023,4),frequency = 4)


################################################################################
# Begin Box-Jenkins Model Selection Process
#   Identification Stage
################################################################################
# Create plots to determine if there is a time trend in the data
# Plots for nominal GDP data
composite_nom_ts <- ts.union(C_nom_ts, I_nom_ts, G_nom_ts, NX_nom_ts)
plot.ts(composite_nom_ts,main = "Nominal GDP Components (Quarterly, 1947-2023)")
ts.plot(composite_nom_ts,main = "Nominal GDP Components (Quarterly, 1947-2023)")

plot.ts(Y_nom_ts,main = "Nominal GDP (Quarterly, 1947-2023)")

# Plots for real GDP data
composite_real_ts <- ts.union(C_real_ts, I_real_ts, G_real_ts, NX_real_ts)
plot.ts(composite_real_ts,main = "Real GDP Components (Quarterly, 1947-2023)")
ts.plot(composite_real_ts,main = "Real GDP Components (Quarterly, 1947-2023)")

plot.ts(Y_real_ts,main="Real GDP (Quarterly, 1947-2023)")

plot.ts(chained_residual_ts,main="Real GDP Residual from Chain 2017 Conversion")

# Upon initial graphical inspection, it appears that both real and nominal GDP
#   and their components have some type of time trend in the data. Given this,
#   I suspect that differencing would be appropriate. I will continue by performing
#   a Dickey-Fuller test on this data to confirm those suspicion
################################################################################
# Dickey-Fuller Tests

# Creates data frames that store the relevant tau and tau critical values based on 
# the number of observations in the data. Labels the row and column names accordingly
tau_critical <- data.frame(col = c("tau1","tau2","tau3"))
tau_critical <- cbind(tau_critical,
                      c(-2.58,-3.44,-3.98),
                      c(-1.95,-2.87,-3.42),
                      c(-1.62,-2.57,-3.13))
colnames(tau_critical) <- c("Test","OnePCT","FivePCT","TenPCT")


phi_critical <- data.frame(col = c("phi1","phi2","phi3"))
phi_critical <- cbind(phi_critical,
                      c(6.47,6.15,8.34),
                      c(4.61,4.71,6.30),
                      c(3/79,4.05,5.36))
colnames(phi_critical) <- c("Test","OnePCT","FivePCT","TenPCT")


# Performs Dickey-Fuller Tests on nominal GDP for each category of possible model
Y_nom_test_trend <- ur.df(Y_nom_ts,selectlags ='BIC', type='trend')
Y_nom_test_drift <- ur.df(Y_nom_ts,selectlags ='BIC', type='drift')
Y_nom_test_none <- ur.df(Y_nom_ts,selectlags ='BIC', type='none')

# Summarizes the first Dickey-Fuller Test to console
summary(Y_nom_test_trend)

# Summarizes the second Dickey-Fuller Test to console
summary(Y_nom_test_drift)

# Summarizes the final Dickey-Fuller Test to console
summary(Y_nom_test_none)
# Records the test statistic values for later calculation
Y_nom_tstats <- list("tau1" = 11.5295,
                     "tau2" = 8.721, 
                     "tau3" = 3.4189,
                     "phi1" = 66.6145, 
                     "phi2" = 44.2641, 
                     "phi3" = 37.9034)

# Initializes Data Frame for storing the results of test statistic analysis
Y_nom_results <- data.frame(col=c("phi","tau"))

# Iterates through all tests and compares them to the t-stats that were taken from
#   the various Dickey-Fuller tests that were performed
for(i in 1:length(tau_critical$Test)){
  # Gets a phi and tau test to begin with
  current_phi <- phi_critical$Test[i]
  current_tau <- tau_critical$Test[i]
  
  # Logic for phi comparisons. Compares the associated phiX test statistic from
  #   the Y_tstat list with the critical values for our given sample size. Stores
  #   the result of the comparison as a string, either confirming a rejection of
  #   the null or a failure to reject, and if the null was able to be rejected
  #   the level of confidence is also recorded. 
  if (Y_nom_tstats[current_phi] >= phi_critical$OnePCT[i]){
    df_phi_test_result <- "Reject (1%)"
  } else if (Y_nom_tstats[current_phi] >= phi_critical$FivePCT[i]){
    df_phi_test_result <- "Reject (5%)"
  } else if (Y_nom_tstats[current_phi] >= phi_critical$TenPCT[i]){
    df_phi_test_result <- "Reject (10%)"
  } else {
    df_phi_test_result <- "Fail"
  }
  
  # Performs above logic again, but for the tauX test
  if (Y_nom_tstats[current_tau] <= tau_critical$OnePCT[i]){
    df_tau_test_result <- "Reject (1%)"
  } else if (Y_nom_tstats[current_tau] <= tau_critical$FivePCT[i]){
    df_tau_test_result <- "Reject (5%)"
  } else if (Y_nom_tstats[current_tau] <= tau_critical$TenPCT[i]){
    df_tau_test_result <- "Reject (10%)"
  } else {
    df_tau_test_result <- "Fail"
  }
  
  # Records the results of the above logic in a new column of the Y_results data frame
  Y_nom_results <- cbind(Y_nom_results,c(df_phi_test_result,df_tau_test_result))
}
# Appropriately names the columns based on the phiX/tauX test performed
colnames(Y_nom_results) <- c("Test",1:length(phi_critical$Test))
# Names the rows to match the phi or tau tests, respectively
 
# Outputs the results to the console
Y_nom_results

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Performs the above analysis again, but this time using the real GDP data
# Performs Dickey-Fuller Tests on nominal GDP for each category of possible model
Y_real_test_trend <- ur.df(Y_real_ts,selectlags ='BIC', type='trend')
Y_real_test_drift <- ur.df(Y_real_ts,selectlags ='BIC', type='drift')
Y_real_test_none <- ur.df(Y_real_ts,selectlags ='BIC', type='none')

# Summarizes the first Dickey-Fuller Test to console
summary(Y_real_test_trend)

# Summarizes the second Dickey-Fuller Test to console
summary(Y_real_test_drift)

# Summarizes the final Dickey-Fuller Test to console
summary(Y_real_test_none)
# Records the test statistic values for later calculation
Y_real_tstats <- list("tau1" = 8.6073,
                      "tau2" = 3.4598, 
                      "tau3" = -0.9695,
                      "phi1" = 38.7108, 
                      "phi2" = 26.8724, 
                      "phi3" = 7.3942)

# Initializes Data Frame for storing the results of test statistic analysis
Y_real_results <- data.frame(col=c("phi","tau"))

# Iterates through all tests and compares them to the t-stats that were taken from
#   the various Dickey-Fuller tests that were performed
for(i in 1:length(tau_critical$Test)){
  # Gets a phi and tau test to begin with
  current_phi <- phi_critical$Test[i]
  current_tau <- tau_critical$Test[i]
  
  # Logic for phi comparisons. Compares the associated phiX test statistic from
  #   the Y_real_tstat list with the critical values for our given sample size. Stores
  #   the result of the comparison as a string, either confirming a rejection of
  #   the null or a failure to reject, and if the null was able to be rejected
  #   the level of confidence is also recorded. 
  if (Y_real_tstats[current_phi] >= phi_critical$OnePCT[i]){
    df_phi_test_result <- "Reject (1%)"
  } else if (Y_real_tstats[current_phi] >= phi_critical$FivePCT[i]){
    df_phi_test_result <- "Reject (5%)"
  } else if (Y_real_tstats[current_phi] >= phi_critical$TenPCT[i]){
    df_phi_test_result <- "Reject (10%)"
  } else {
    df_phi_test_result <- "Fail"
  }
  
  # Performs above logic again, but for the tauX test
  if (Y_real_tstats[current_tau] <= tau_critical$OnePCT[i]){
    df_tau_test_result <- "Reject (1%)"
  } else if (Y_real_tstats[current_tau] <= tau_critical$FivePCT[i]){
    df_tau_test_result <- "Reject (5%)"
  } else if (Y_real_tstats[current_tau] <= tau_critical$TenPCT[i]){
    df_tau_test_result <- "Reject (10%)"
  } else {
    df_tau_test_result <- "Fail"
  }
  
  # Records the results of the above logic in a new column of the Y_real_results data frame
  Y_real_results <- cbind(Y_real_results,c(df_phi_test_result,df_tau_test_result))
}
# Appropriately names the columns based on the phiX/tauX test performed
colnames(Y_real_results) <- c("Test",1:length(phi_critical$Test))
# Names the rows to match the phi or tau tests, respectively

# Outputs the results to the console
Y_real_results



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Given the above analysis I will also perform DF tests on differenced GDP data
#   for both nominal and real GDP data

## Creates a time-series object of the 1st differenced Real and Nominal GDP 
Y_d1_real_ts <- diff(Y_real_ts)
Y_d1_nom_ts <- diff(Y_nom_ts)

## Plots the Real and Nominal 1st differenced GDP
plot.ts(Y_d1_real_ts,main = "Real GDP, 1st Diff (Quarterly, 1947-2023)")
plot.ts(Y_d1_nom_ts,main = "Nominal GDP, 1st Diff (Quarterly, 1947-2023)")

# Looking at these two charts, an issue that I had suspected becomes clear - real
#   GDP removes the influence of inflation on GDP, whereas nominal GDP also includes
#   inflationary trends. Given this, and dispite the residual issue when converting
#   to 2017 dollar, real GDP will be the superior measure

# Performs Dickey-Fuller Tests on real 1st differenced GDP for each category of 
#   possible model
Y_d1_real_test_trend <- ur.df(Y_d1_real_ts,selectlags ='BIC', type='trend')
Y_d1_real_test_drift <- ur.df(Y_d1_real_ts,selectlags ='BIC', type='drift')
Y_d1_real_test_none <- ur.df(Y_d1_real_ts,selectlags ='BIC', type='none')

# Summarizes the first Dickey-Fuller Test to console
summary(Y_d1_real_test_trend)

# Summarizes the second Dickey-Fuller Test to console
summary(Y_d1_real_test_drift)

# Summarizes the final Dickey-Fuller Test to console
summary(Y_d1_real_test_none)
# Records the test statistic values for later calculation
Y_d1_real_tstats <- list("tau1" = -10.0256,
                         "tau2" = -13.118, 
                         "tau3" = -13.94,
                         "phi1" = 86.0498, 
                         "phi2" = 64.7827, 
                         "phi3" = 97.1651)

# Initializes Data Frame for storing the results of test statistic analysis
Y_d1_real_results <- data.frame(col=c("phi","tau"))

# Iterates through all tests and compares them to the t-stats that were taken from
#   the various Dickey-Fuller tests that were performed
for(i in 1:length(tau_critical$Test)){
  # Gets a phi and tau test to begin with
  current_phi <- phi_critical$Test[i]
  current_tau <- tau_critical$Test[i]
  
  # Logic for phi comparisons. Compares the associated phiX test statistic from
  #   the Y_tstat list with the critical values for our given sample size. Stores
  #   the result of the comparison as a string, either confirming a rejection of
  #   the null or a failure to reject, and if the null was able to be rejected
  #   the level of confidence is also recorded. 
  if (Y_d1_real_tstats[current_phi] >= phi_critical$OnePCT[i]){
    df_phi_test_result <- "Reject (1%)"
  } else if (Y_d1_real_tstats[current_phi] >= phi_critical$FivePCT[i]){
    df_phi_test_result <- "Reject (5%)"
  } else if (Y_d1_real_tstats[current_phi] >= phi_critical$TenPCT[i]){
    df_phi_test_result <- "Reject (10%)"
  } else {
    df_phi_test_result <- "Fail"
  }
  
  # Performs above logic again, but for the tauX test
  if (Y_d1_real_tstats[current_tau] <= tau_critical$OnePCT[i]){
    df_tau_test_result <- "Reject (1%)"
  } else if (Y_d1_real_tstats[current_tau] <= tau_critical$FivePCT[i]){
    df_tau_test_result <- "Reject (5%)"
  } else if (Y_d1_real_tstats[current_tau] <= tau_critical$TenPCT[i]){
    df_tau_test_result <- "Reject (10%)"
  } else {
    df_tau_test_result <- "Fail"
  }
  
  # Records the results of the above logic in a new column of the Y_results data frame
  Y_d1_real_results <- cbind(Y_d1_real_results,c(df_phi_test_result,df_tau_test_result))
}
# Appropriately names the columns based on the phiX/tauX test performed
colnames(Y_d1_real_results) <- c("Test",1:length(phi_critical$Test))
# Names the rows to match the phi or tau tests, respectively

# Outputs the results to the console
Y_d1_real_results


################################################################################
# ACf and PACF Tests on Y

# Now that I have run the Dicker-Fuller tests and determined that a 1st differenced
#   real GDP measure is the appropriate one for my analysis, now I will take a look
#   at the ACf and PACF of the data to try to determine what model will be best
#   going forward.
acf2(Y_d1_real_ts, max.lag=30, plot=TRUE) 
# The ACF and PACF are really not very convincing in telling me what model would best
#   describe this data. I am going to use the Auto ARIMA function to see where that
#   lands me
auto.arima(Y_d1_real_ts, trace=TRUE, stepwise = FALSE, ic="bic", seasonal=FALSE, stationary = FALSE)
auto.arima(Y_d1_real_ts, trace=TRUE, stepwise = FALSE, ic="aic", seasonal=FALSE, stationary = FALSE)

auto.arima(Y_real_ts, trace=TRUE, stepwise = FALSE, ic="bic", seasonal=FALSE, stationary = FALSE)
auto.arima(Y_real_ts, trace=TRUE, stepwise = FALSE, ic="aic", seasonal=FALSE, stationary = FALSE)


# Looks like a 2nd difference of the data is warranted. Let me see what the DF tests look 
#   like on that type of transformation...
Y_d2_real_ts <- diff(Y_d1_real_ts)
# Performs Dickey-Fuller Tests on real 1st differenced GDP for each category of 
#   possible model
Y_d2_real_test_trend <- ur.df(Y_d2_real_ts,selectlags ='BIC', type='trend')
Y_d2_real_test_drift <- ur.df(Y_d2_real_ts,selectlags ='BIC', type='drift')
Y_d2_real_test_none <- ur.df(Y_d2_real_ts,selectlags ='BIC', type='none')

# Summarizes the first Dickey-Fuller Test to console
summary(Y_d2_real_test_trend)

# Summarizes the second Dickey-Fuller Test to console
summary(Y_d2_real_test_drift)

# Summarizes the final Dickey-Fuller Test to console
summary(Y_d2_real_test_none)
# Records the test statistic values for later calculation
Y_d2_real_tstats <- list("tau1" = -22.7034,
                         "tau2" = -22.6668, 
                         "tau3" = -22.6297,
                         "phi1" = 256.8924, 
                         "phi2" = 170.7017, 
                         "phi3" = 256.0525)

# Initializes Data Frame for storing the results of test statistic analysis
Y_d2_real_results <- data.frame(col=c("phi","tau"))

# Iterates through all tests and compares them to the t-stats that were taken from
#   the various Dickey-Fuller tests that were performed
for(i in 1:length(tau_critical$Test)){
  # Gets a phi and tau test to begin with
  current_phi <- phi_critical$Test[i]
  current_tau <- tau_critical$Test[i]
  
  # Logic for phi comparisons. Compares the associated phiX test statistic from
  #   the Y_tstat list with the critical values for our given sample size. Stores
  #   the result of the comparison as a string, either confirming a rejection of
  #   the null or a failure to reject, and if the null was able to be rejected
  #   the level of confidence is also recorded. 
  if (Y_d2_real_tstats[current_phi] >= phi_critical$OnePCT[i]){
    df_phi_test_result <- "Reject (1%)"
  } else if (Y_d2_real_tstats[current_phi] >= phi_critical$FivePCT[i]){
    df_phi_test_result <- "Reject (5%)"
  } else if (Y_d2_real_tstats[current_phi] >= phi_critical$TenPCT[i]){
    df_phi_test_result <- "Reject (10%)"
  } else {
    df_phi_test_result <- "Fail"
  }
  
  # Performs above logic again, but for the tauX test
  if (Y_d2_real_tstats[current_tau] <= tau_critical$OnePCT[i]){
    df_tau_test_result <- "Reject (1%)"
  } else if (Y_d2_real_tstats[current_tau] <= tau_critical$FivePCT[i]){
    df_tau_test_result <- "Reject (5%)"
  } else if (Y_d2_real_tstats[current_tau] <= tau_critical$TenPCT[i]){
    df_tau_test_result <- "Reject (10%)"
  } else {
    df_tau_test_result <- "Fail"
  }
  
  # Records the results of the above logic in a new column of the Y_results data frame
  Y_d2_real_results <- cbind(Y_d2_real_results,c(df_phi_test_result,df_tau_test_result))
}
# Appropriately names the columns based on the phiX/tauX test performed
colnames(Y_d2_real_results) <- c("Test",1:length(phi_critical$Test))
# Names the rows to match the phi or tau tests, respectively

# Outputs the results to the console
Y_d2_real_results

# Given all of that above analysis, I am going to move forward with an AARIMA(0,2,2)
#   model performed on the real GDP measure. Given my goal of forecast coherence,
#   I should also use the same model for the real GDP components. Let's see how that works

acf2(Y_d2_real_ts, max.lag=30, plot=TRUE)


################################################################################
# Dickey-Fuller Tests on Real GDP Components
# Consumption
C_real_test_trend <- ur.df(C_real_ts,selectlags ='BIC', type='trend')
C_real_test_drift <- ur.df(C_real_ts,selectlags ='BIC', type='drift')
C_real_test_none <- ur.df(C_real_ts,selectlags ='BIC', type='none')

# Summarizes the first Dickey-Fuller Test to console
summary(C_real_test_trend)

# Summarizes the second Dickey-Fuller Test to console
summary(C_real_test_drift)

# Summarizes the final Dickey-Fuller Test to console
summary(C_real_test_none)
# Records the test statistic values for later calculation
C_real_tstats <- list("tau1" = 8.7703,
                         "tau2" = 3.9753, 
                         "tau3" = -0.7418,
                         "phi1" = 39.6422, 
                         "phi2" = 27.484, 
                         "phi3" = 9.3052)

# Investment
I_real_test_trend <- ur.df(I_real_ts,selectlags ='BIC', type='trend')
I_real_test_drift <- ur.df(I_real_ts,selectlags ='BIC', type='drift')
I_real_test_none <- ur.df(I_real_ts,selectlags ='BIC', type='none')

# Summarizes the first Dickey-Fuller Test to console
summary(I_real_test_trend)

# Summarizes the second Dickey-Fuller Test to console
summary(I_real_test_drift)

# Summarizes the final Dickey-Fuller Test to console
summary(I_real_test_none)
# Records the test statistic values for later calculation
I_real_tstats <- list("tau1" = 3.0886,
                      "tau2" = 1.1973, 
                      "tau3" = -1.5754,
                      "phi1" = 5.151, 
                      "phi2" = 4.8703, 
                      "phi3" = 2.8243)

# Government Expenditures
G_real_test_trend <- ur.df(G_real_ts,selectlags ='BIC', type='trend')
G_real_test_drift <- ur.df(G_real_ts,selectlags ='BIC', type='drift')
G_real_test_none <- ur.df(G_real_ts,selectlags ='BIC', type='none')

# Summarizes the first Dickey-Fuller Test to console
summary(G_real_test_trend)

# Summarizes the second Dickey-Fuller Test to console
summary(G_real_test_drift)

# Summarizes the final Dickey-Fuller Test to console
summary(G_real_test_none)
# Records the test statistic values for later calculation
G_real_tstats <- list("tau1" = 4.6197,
                      "tau2" = -0.3854, 
                      "tau3" = -2.2882,
                      "phi1" = 14.1087, 
                      "phi2" = 11.2299, 
                      "phi3" = 2.6208)

# Net Exports
NX_real_test_trend <- ur.df(NX_real_ts,selectlags ='BIC', type='trend')
NX_real_test_drift <- ur.df(NX_real_ts,selectlags ='BIC', type='drift')
NX_real_test_none <- ur.df(NX_real_ts,selectlags ='BIC', type='none')

# Summarizes the first Dickey-Fuller Test to console
summary(NX_real_test_trend)

# Summarizes the second Dickey-Fuller Test to console
summary(NX_real_test_drift)

# Summarizes the final Dickey-Fuller Test to console
summary(NX_real_test_none)
# Records the test statistic values for later calculation
NX_real_tstats <- list("tau1" = 1.1326,
                      "tau2" = 0.227, 
                      "tau3" = -1.4084,
                      "phi1" = 1.2511, 
                      "phi2" = 2.0645, 
                      "phi3" = 1.8606)

# Perform Checks on Hypothesis for each real GDP component
# Initializes Data Frame for storing the results of test statistic analysis
C_real_results <- data.frame(col=c("phi","tau"))

# Iterates through all tests and compares them to the t-stats that were taken from
#   the various Dickey-Fuller tests that were performed
for(i in 1:length(tau_critical$Test)){
  # Gets a phi and tau test to begin with
  current_phi <- phi_critical$Test[i]
  current_tau <- tau_critical$Test[i]
  
  # Logic for phi comparisons. Compares the associated phiX test statistic from
  #   the Y_tstat list with the critical values for our given sample size. Stores
  #   the result of the comparison as a string, either confirming a rejection of
  #   the null or a failure to reject, and if the null was able to be rejected
  #   the level of confidence is also recorded. 
  if (C_real_tstats[current_phi] >= phi_critical$OnePCT[i]){
    df_phi_test_result <- "Reject (1%)"
  } else if (C_real_tstats[current_phi] >= phi_critical$FivePCT[i]){
    df_phi_test_result <- "Reject (5%)"
  } else if (C_real_tstats[current_phi] >= phi_critical$TenPCT[i]){
    df_phi_test_result <- "Reject (10%)"
  } else {
    df_phi_test_result <- "Fail"
  }
  
  # Performs above logic again, but for the tauX test
  if (C_real_tstats[current_tau] <= tau_critical$OnePCT[i]){
    df_tau_test_result <- "Reject (1%)"
  } else if (C_real_tstats[current_tau] <= tau_critical$FivePCT[i]){
    df_tau_test_result <- "Reject (5%)"
  } else if (C_real_tstats[current_tau] <= tau_critical$TenPCT[i]){
    df_tau_test_result <- "Reject (10%)"
  } else {
    df_tau_test_result <- "Fail"
  }
  
  # Records the results of the above logic in a new column of the Y_results data frame
  C_real_results <- cbind(C_real_results,c(df_phi_test_result,df_tau_test_result))
}
# Appropriately names the columns based on the phiX/tauX test performed
colnames(C_real_results) <- c("Test",1:length(phi_critical$Test))
# Names the rows to match the phi or tau tests, respectively

# Outputs the results to the console
C_real_results


# Initializes Data Frame for storing the results of test statistic analysis
I_real_results <- data.frame(col=c("phi","tau"))

# Iterates through all tests and compares them to the t-stats that were taken from
#   the various Dickey-Fuller tests that were performed
for(i in 1:length(tau_critical$Test)){
  # Gets a phi and tau test to begin with
  current_phi <- phi_critical$Test[i]
  current_tau <- tau_critical$Test[i]
  
  # Logic for phi comparisons. Compares the associated phiX test statistic from
  #   the Y_tstat list with the critical values for our given sample size. Stores
  #   the result of the comparison as a string, either confirming a rejection of
  #   the null or a failure to reject, and if the null was able to be rejected
  #   the level of confidence is also recorded. 
  if (I_real_tstats[current_phi] >= phi_critical$OnePCT[i]){
    df_phi_test_result <- "Reject (1%)"
  } else if (I_real_tstats[current_phi] >= phi_critical$FivePCT[i]){
    df_phi_test_result <- "Reject (5%)"
  } else if (I_real_tstats[current_phi] >= phi_critical$TenPCT[i]){
    df_phi_test_result <- "Reject (10%)"
  } else {
    df_phi_test_result <- "Fail"
  }
  
  # Performs above logic again, but for the tauX test
  if (I_real_tstats[current_tau] <= tau_critical$OnePCT[i]){
    df_tau_test_result <- "Reject (1%)"
  } else if (I_real_tstats[current_tau] <= tau_critical$FivePCT[i]){
    df_tau_test_result <- "Reject (5%)"
  } else if (I_real_tstats[current_tau] <= tau_critical$TenPCT[i]){
    df_tau_test_result <- "Reject (10%)"
  } else {
    df_tau_test_result <- "Fail"
  }
  
  # Records the results of the above logic in a new column of the Y_results data frame
  I_real_results <- cbind(I_real_results,c(df_phi_test_result,df_tau_test_result))
}
# Appropriately names the columns based on the phiX/tauX test performed
colnames(I_real_results) <- c("Test",1:length(phi_critical$Test))
# Names the rows to match the phi or tau tests, respectively

# Outputs the results to the console
I_real_results

# Initializes Data Frame for storing the results of test statistic analysis
G_real_results <- data.frame(col=c("phi","tau"))

# Iterates through all tests and compares them to the t-stats that were taken from
#   the various Dickey-Fuller tests that were performed
for(i in 1:length(tau_critical$Test)){
  # Gets a phi and tau test to begin with
  current_phi <- phi_critical$Test[i]
  current_tau <- tau_critical$Test[i]
  
  # Logic for phi comparisons. Compares the associated phiX test statistic from
  #   the Y_tstat list with the critical values for our given sample size. Stores
  #   the result of the comparison as a string, either confirming a rejection of
  #   the null or a failure to reject, and if the null was able to be rejected
  #   the level of confidence is also recorded. 
  if (G_real_tstats[current_phi] >= phi_critical$OnePCT[i]){
    df_phi_test_result <- "Reject (1%)"
  } else if (G_real_tstats[current_phi] >= phi_critical$FivePCT[i]){
    df_phi_test_result <- "Reject (5%)"
  } else if (G_real_tstats[current_phi] >= phi_critical$TenPCT[i]){
    df_phi_test_result <- "Reject (10%)"
  } else {
    df_phi_test_result <- "Fail"
  }
  
  # Performs above logic again, but for the tauX test
  if (G_real_tstats[current_tau] <= tau_critical$OnePCT[i]){
    df_tau_test_result <- "Reject (1%)"
  } else if (G_real_tstats[current_tau] <= tau_critical$FivePCT[i]){
    df_tau_test_result <- "Reject (5%)"
  } else if (G_real_tstats[current_tau] <= tau_critical$TenPCT[i]){
    df_tau_test_result <- "Reject (10%)"
  } else {
    df_tau_test_result <- "Fail"
  }
  
  # Records the results of the above logic in a new column of the Y_results data frame
  G_real_results <- cbind(G_real_results,c(df_phi_test_result,df_tau_test_result))
}
# Appropriately names the columns based on the phiX/tauX test performed
colnames(G_real_results) <- c("Test",1:length(phi_critical$Test))
# Names the rows to match the phi or tau tests, respectively

# Outputs the results to the console
G_real_results

# Initializes Data Frame for storing the results of test statistic analysis
NX_real_results <- data.frame(col=c("phi","tau"))

# Iterates through all tests and compares them to the t-stats that were taken from
#   the various Dickey-Fuller tests that were performed
for(i in 1:length(tau_critical$Test)){
  # Gets a phi and tau test to begin with
  current_phi <- phi_critical$Test[i]
  current_tau <- tau_critical$Test[i]
  
  # Logic for phi comparisons. Compares the associated phiX test statistic from
  #   the Y_tstat list with the critical values for our given sample size. Stores
  #   the result of the comparison as a string, either confirming a rejection of
  #   the null or a failure to reject, and if the null was able to be rejected
  #   the level of confidence is also recorded. 
  if (NX_real_tstats[current_phi] >= phi_critical$OnePCT[i]){
    df_phi_test_result <- "Reject (1%)"
  } else if (NX_real_tstats[current_phi] >= phi_critical$FivePCT[i]){
    df_phi_test_result <- "Reject (5%)"
  } else if (NX_real_tstats[current_phi] >= phi_critical$TenPCT[i]){
    df_phi_test_result <- "Reject (10%)"
  } else {
    df_phi_test_result <- "Fail"
  }
  
  # Performs above logic again, but for the tauX test
  if (NX_real_tstats[current_tau] <= tau_critical$OnePCT[i]){
    df_tau_test_result <- "Reject (1%)"
  } else if (NX_real_tstats[current_tau] <= tau_critical$FivePCT[i]){
    df_tau_test_result <- "Reject (5%)"
  } else if (NX_real_tstats[current_tau] <= tau_critical$TenPCT[i]){
    df_tau_test_result <- "Reject (10%)"
  } else {
    df_tau_test_result <- "Fail"
  }
  
  # Records the results of the above logic in a new column of the Y_results data frame
  NX_real_results <- cbind(NX_real_results,c(df_phi_test_result,df_tau_test_result))
}
# Appropriately names the columns based on the phiX/tauX test performed
colnames(NX_real_results) <- c("Test",1:length(phi_critical$Test))
# Names the rows to match the phi or tau tests, respectively

# Outputs the results to the console
NX_real_results

# Utilize the AutoARIMA function for C and the ACF/PACF for each component
auto.arima(C_real_ts, trace=F, stepwise = FALSE, ic="bic", seasonal=FALSE, stationary = F)
auto.arima(C_real_ts, trace=F, stepwise = FALSE, ic="aic", seasonal=FALSE, stationary = F)
auto.arima(C_real_ts, trace=F, stepwise = FALSE, ic="bic", seasonal=FALSE, stationary = T)
auto.arima(C_real_ts, trace=F, stepwise = FALSE, ic="aic", seasonal=FALSE, stationary = T)

acf2(C_real_ts, max.lag=30, plot=TRUE)
acf2(I_real_ts, max.lag=30, plot=TRUE)
acf2(G_real_ts, max.lag=30, plot=TRUE)
acf2(NX_real_ts, max.lag=30, plot=TRUE)

# Given the outcome for the ACF/PACF for each component seems to follow a similar
#   trend to the GDP measure, the results of the Dickey-Fuller tests, and the 
#   response of the AutoARIMA function, I will be moving forward with an ARIMA(0,2,2)
#   model for Y, C, I, G, and NX

acf2(chained_residual_ts, max.lag=30, plot=TRUE)

# Chained residuals are another matter entirely, but I think that if I am going to
#   use an ARIMA(0,2,2) model for almost everything else, then I should also do the
#   same on this unofficial GDP "component". I will have to think on this a bit though

################################################################################
# Dickey-Fuller on 2nd Differenced Components
# Creates the 2nd differenced time series of each component
## Consumption
C_d1_real_ts <- diff(C_real_ts)
C_d2_real_ts <- diff(C_d1_real_ts)
## Investment
I_d1_real_ts <- diff(I_real_ts)
I_d2_real_ts <- diff(I_d1_real_ts)
## Government
G_d1_real_ts <- diff(G_real_ts)
G_d2_real_ts <- diff(G_d1_real_ts)
## Net Exports
NX_d1_real_ts <- diff(NX_real_ts)
NX_d2_real_ts <- diff(NX_d1_real_ts)

# Consumption
C_d2_real_test_trend <- ur.df(C_d2_real_ts,selectlags ='BIC', type='trend')
C_d2_real_test_drift <- ur.df(C_d2_real_ts,selectlags ='BIC', type='drift')
C_d2_real_test_none <- ur.df(C_d2_real_ts,selectlags ='BIC', type='none')

# Summarizes the first Dickey-Fuller Test to console
summary(C_d2_real_test_trend)

# Summarizes the second Dickey-Fuller Test to console
summary(C_d2_real_test_drift)

# Summarizes the final Dickey-Fuller Test to console
summary(C_d2_real_test_none)
# Records the test statistic values for later calculation
C_d2_real_tstats <- list("tau1" = -23.6745,
                      "tau2" = -23.6359, 
                      "tau3" = -23.5966,
                      "phi1" = 279.3284, 
                      "phi2" = 185.6013, 
                      "phi3" = 278.4012)

# Investment
I_d2_real_test_trend <- ur.df(I_d2_real_ts,selectlags ='BIC', type='trend')
I_d2_real_test_drift <- ur.df(I_d2_real_ts,selectlags ='BIC', type='drift')
I_d2_real_test_none <- ur.df(I_d2_real_ts,selectlags ='BIC', type='none')

# Summarizes the first Dickey-Fuller Test to console
summary(I_d2_real_test_trend)

# Summarizes the second Dickey-Fuller Test to console
summary(I_d2_real_test_drift)

# Summarizes the final Dickey-Fuller Test to console
summary(I_d2_real_test_none)
# Records the test statistic values for later calculation
I_d2_real_tstats <- list("tau1" = -20.9169,
                      "tau2" = -20.8822,
                      "tau3" = -20.8478,
                      "phi1" = 218.0355, 
                      "phi2" = 144.8792, 
                      "phi3" = 217.3157)

# Government Expenditures
G_d2_real_test_trend <- ur.df(G_d2_real_ts,selectlags ='BIC', type='trend')
G_d2_real_test_drift <- ur.df(G_d2_real_ts,selectlags ='BIC', type='drift')
G_d2_real_test_none <- ur.df(G_d2_real_ts,selectlags ='BIC', type='none')

# Summarizes the first Dickey-Fuller Test to console
summary(G_d2_real_test_trend)

# Summarizes the second Dickey-Fuller Test to console
summary(G_d2_real_test_drift)

# Summarizes the final Dickey-Fuller Test to console
summary(G_d2_real_test_none)
# Records the test statistic values for later calculation
G_d2_real_tstats <- list("tau1" = -21.8995,
                      "tau2" = -21.8666, 
                      "tau3" = -21.8322,
                      "phi1" = 239.0749, 
                      "phi2" = 158.8824, 
                      "phi3" = 238.3234)

# Net Exports
NX_d2_real_test_trend <- ur.df(NX_d2_real_ts,selectlags ='BIC', type='trend')
NX_d2_real_test_drift <- ur.df(NX_d2_real_ts,selectlags ='BIC', type='drift')
NX_d2_real_test_none <- ur.df(NX_d2_real_ts,selectlags ='BIC', type='none')

# Summarizes the first Dickey-Fuller Test to console
summary(NX_d2_real_test_trend)

# Summarizes the second Dickey-Fuller Test to console
summary(NX_d2_real_test_drift)

# Summarizes the final Dickey-Fuller Test to console
summary(NX_d2_real_test_none)
# Records the test statistic values for later calculation
NX_d2_real_tstats <- list("tau1" = -21.3399,
                       "tau2" = -21.3047, 
                       "tau3" = -21.2694,
                       "phi1" = 226.9451, 
                       "phi2" = 150.796, 
                       "phi3" = 226.193)
C_d2_real_results <- data.frame(col=c("phi","tau"))

# Iterates through all tests and compares them to the t-stats that were taken from
#   the various Dickey-Fuller tests that were performed
for(i in 1:length(tau_critical$Test)){
  # Gets a phi and tau test to begin with
  current_phi <- phi_critical$Test[i]
  current_tau <- tau_critical$Test[i]
  
  # Logic for phi comparisons. Compares the associated phiX test statistic from
  #   the Y_tstat list with the critical values for our given sample size. Stores
  #   the result of the comparison as a string, either confirming a rejection of
  #   the null or a failure to reject, and if the null was able to be rejected
  #   the level of confidence is also recorded. 
  if (C_d2_real_tstats[current_phi] >= phi_critical$OnePCT[i]){
    df_phi_test_result <- "Reject (1%)"
  } else if (C_d2_real_tstats[current_phi] >= phi_critical$FivePCT[i]){
    df_phi_test_result <- "Reject (5%)"
  } else if (C_d2_real_tstats[current_phi] >= phi_critical$TenPCT[i]){
    df_phi_test_result <- "Reject (10%)"
  } else {
    df_phi_test_result <- "Fail"
  }
  
  # Performs above logic again, but for the tauX test
  if (C_d2_real_tstats[current_tau] <= tau_critical$OnePCT[i]){
    df_tau_test_result <- "Reject (1%)"
  } else if (C_d2_real_tstats[current_tau] <= tau_critical$FivePCT[i]){
    df_tau_test_result <- "Reject (5%)"
  } else if (C_d2_real_tstats[current_tau] <= tau_critical$TenPCT[i]){
    df_tau_test_result <- "Reject (10%)"
  } else {
    df_tau_test_result <- "Fail"
  }
  
  # Records the results of the above logic in a new column of the Y_results data frame
  C_d2_real_results <- cbind(C_d2_real_results,c(df_phi_test_result,df_tau_test_result))
}
# Appropriately names the columns based on the phiX/tauX test performed
colnames(C_d2_real_results) <- c("Test",1:length(phi_critical$Test))
# Names the rows to match the phi or tau tests, respectively

# Outputs the results to the console
C_d2_real_results


# Initializes Data Frame for storing the results of test statistic analysis
I_d2_real_results <- data.frame(col=c("phi","tau"))

# Iterates through all tests and compares them to the t-stats that were taken from
#   the various Dickey-Fuller tests that were performed
for(i in 1:length(tau_critical$Test)){
  # Gets a phi and tau test to begin with
  current_phi <- phi_critical$Test[i]
  current_tau <- tau_critical$Test[i]
  
  # Logic for phi comparisons. Compares the associated phiX test statistic from
  #   the Y_tstat list with the critical values for our given sample size. Stores
  #   the result of the comparison as a string, either confirming a rejection of
  #   the null or a failure to reject, and if the null was able to be rejected
  #   the level of confidence is also recorded. 
  if (I_d2_real_tstats[current_phi] >= phi_critical$OnePCT[i]){
    df_phi_test_result <- "Reject (1%)"
  } else if (I_d2_real_tstats[current_phi] >= phi_critical$FivePCT[i]){
    df_phi_test_result <- "Reject (5%)"
  } else if (I_d2_real_tstats[current_phi] >= phi_critical$TenPCT[i]){
    df_phi_test_result <- "Reject (10%)"
  } else {
    df_phi_test_result <- "Fail"
  }
  
  # Performs above logic again, but for the tauX test
  if (I_d2_real_tstats[current_tau] <= tau_critical$OnePCT[i]){
    df_tau_test_result <- "Reject (1%)"
  } else if (I_d2_real_tstats[current_tau] <= tau_critical$FivePCT[i]){
    df_tau_test_result <- "Reject (5%)"
  } else if (I_d2_real_tstats[current_tau] <= tau_critical$TenPCT[i]){
    df_tau_test_result <- "Reject (10%)"
  } else {
    df_tau_test_result <- "Fail"
  }
  
  # Records the results of the above logic in a new column of the Y_results data frame
  I_d2_real_results <- cbind(I_d2_real_results,c(df_phi_test_result,df_tau_test_result))
}
# Appropriately names the columns based on the phiX/tauX test performed
colnames(I_d2_real_results) <- c("Test",1:length(phi_critical$Test))
# Names the rows to match the phi or tau tests, respectively

# Outputs the results to the console
I_d2_real_results

# Initializes Data Frame for storing the results of test statistic analysis
G_d2_real_results <- data.frame(col=c("phi","tau"))

# Iterates through all tests and compares them to the t-stats that were taken from
#   the various Dickey-Fuller tests that were performed
for(i in 1:length(tau_critical$Test)){
  # Gets a phi and tau test to begin with
  current_phi <- phi_critical$Test[i]
  current_tau <- tau_critical$Test[i]
  
  # Logic for phi comparisons. Compares the associated phiX test statistic from
  #   the Y_tstat list with the critical values for our given sample size. Stores
  #   the result of the comparison as a string, either confirming a rejection of
  #   the null or a failure to reject, and if the null was able to be rejected
  #   the level of confidence is also recorded. 
  if (G_d2_real_tstats[current_phi] >= phi_critical$OnePCT[i]){
    df_phi_test_result <- "Reject (1%)"
  } else if (G_d2_real_tstats[current_phi] >= phi_critical$FivePCT[i]){
    df_phi_test_result <- "Reject (5%)"
  } else if (G_d2_real_tstats[current_phi] >= phi_critical$TenPCT[i]){
    df_phi_test_result <- "Reject (10%)"
  } else {
    df_phi_test_result <- "Fail"
  }
  
  # Performs above logic again, but for the tauX test
  if (G_d2_real_tstats[current_tau] <= tau_critical$OnePCT[i]){
    df_tau_test_result <- "Reject (1%)"
  } else if (G_d2_real_tstats[current_tau] <= tau_critical$FivePCT[i]){
    df_tau_test_result <- "Reject (5%)"
  } else if (G_d2_real_tstats[current_tau] <= tau_critical$TenPCT[i]){
    df_tau_test_result <- "Reject (10%)"
  } else {
    df_tau_test_result <- "Fail"
  }
  
  # Records the results of the above logic in a new column of the Y_results data frame
  G_d2_real_results <- cbind(G_d2_real_results,c(df_phi_test_result,df_tau_test_result))
}
# Appropriately names the columns based on the phiX/tauX test performed
colnames(G_d2_real_results) <- c("Test",1:length(phi_critical$Test))
# Names the rows to match the phi or tau tests, respectively

# Outputs the results to the console
G_d2_real_results

# Initializes Data Frame for storing the results of test statistic analysis
NX_d2_real_results <- data.frame(col=c("phi","tau"))

# Iterates through all tests and compares them to the t-stats that were taken from
#   the various Dickey-Fuller tests that were performed
for(i in 1:length(tau_critical$Test)){
  # Gets a phi and tau test to begin with
  current_phi <- phi_critical$Test[i]
  current_tau <- tau_critical$Test[i]
  
  # Logic for phi comparisons. Compares the associated phiX test statistic from
  #   the Y_tstat list with the critical values for our given sample size. Stores
  #   the result of the comparison as a string, either confirming a rejection of
  #   the null or a failure to reject, and if the null was able to be rejected
  #   the level of confidence is also recorded. 
  if (NX_d2_real_tstats[current_phi] >= phi_critical$OnePCT[i]){
    df_phi_test_result <- "Reject (1%)"
  } else if (NX_d2_real_tstats[current_phi] >= phi_critical$FivePCT[i]){
    df_phi_test_result <- "Reject (5%)"
  } else if (NX_d2_real_tstats[current_phi] >= phi_critical$TenPCT[i]){
    df_phi_test_result <- "Reject (10%)"
  } else {
    df_phi_test_result <- "Fail"
  }
  
  # Performs above logic again, but for the tauX test
  if (NX_d2_real_tstats[current_tau] <= tau_critical$OnePCT[i]){
    df_tau_test_result <- "Reject (1%)"
  } else if (NX_d2_real_tstats[current_tau] <= tau_critical$FivePCT[i]){
    df_tau_test_result <- "Reject (5%)"
  } else if (NX_d2_real_tstats[current_tau] <= tau_critical$TenPCT[i]){
    df_tau_test_result <- "Reject (10%)"
  } else {
    df_tau_test_result <- "Fail"
  }
  
  # Records the results of the above logic in a new column of the Y_results data frame
  NX_d2_real_results <- cbind(NX_d2_real_results,c(df_phi_test_result,df_tau_test_result))
}
# Appropriately names the columns based on the phiX/tauX test performed
colnames(NX_d2_real_results) <- c("Test",1:length(phi_critical$Test))
# Names the rows to match the phi or tau tests, respectively

# Outputs the results to the console
NX_d2_real_results
################################################################################
# Model Creation
################################################################################
# GDP ARIMA(0,2,2) Model w/o constant
Y_model <- Arima(Y_real_ts, order = c(0,2,2), method = c("ML"), include.constant = FALSE )


# Component Models
## Consumption ARIMA(0,2,2) Model w/o constant
C_model <- Arima(C_real_ts, order = c(0,2,2), method = c("ML"), include.constant = FALSE )
## Investment ARIMA(0,2,2) Model w/o constant
I_model <- Arima(I_real_ts, order = c(0,2,1), method = c("ML"), include.constant = FALSE )
## Government Expenditures ARIMA(0,2,2) Model w/o constant
G_model <- Arima(G_real_ts, order = c(0,2,1), method = c("ML"), include.constant = FALSE )
## Net Exports ARIMA(0,2,2) Model w/o constant
NX_model <- Arima(NX_real_ts, order = c(0,2,2), method = c("ML"), include.constant = FALSE )
## Chained Residual ARIMA(0,2,2) Model w/o constant
CR_model <- Arima(chained_residual_ts,order = c(0,2,2), method = c("ML"), include.constant = F)

all_models <- list("GDP" = Y_model,
                   "Consumption" = C_model, 
                   "Investment" = I_model,
                   "Government" = G_model,
                   "Net Exports" = NX_model,
                   "Chained Residulas" = CR_model)
modelsummary(all_models, stars = T)

