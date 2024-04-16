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
gdp_data = read_excel(file.choose())

# Attach the file for more ease of data calling for the next few lines
attach(gdp_data)
gdp_nobs <- nrow(gdp_data)
# Create time series objects for GDP (Y), Consumption (C), Investment (I), Government
#   Expenditures (G), and Net Exports (NX)
Y_ts <- ts(data = `Gross Domestic Product (Y)`,start = c(1947,1),end = c(2023,4),frequency = 4)
C_ts <- ts(data = `Personal Consumption Expenditures (C)`,start = c(1947,1),end = c(2023,4),frequency = 4)
I_ts <- ts(data = `Gross Private Domestic Investment (I)`,start = c(1947,1),end = c(2023,4),frequency = 4)
G_ts <- ts(data = `Government Consumption Expenditures and Gross Investment (G)`,start = c(1947,1),end = c(2023,4),frequency = 4)
NX_ts <- ts(data = `Net Exports (NX)`,start = c(1947,1),end = c(2023,4),frequency = 4)


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
composite_ts <- ts.union(C_ts, I_ts, G_ts, NX_ts)
plot.ts(composite_ts,main = "Nominal GDP Components (Quarterly, 1947-2023)")
ts.plot(composite_ts,main = "Nominal GDP Components (Quarterly, 1947-2023)")

plot.ts(Y_ts,main = "Nominal GDP (Quarterly, 1947-2023)")

# Plots for real GDP data
composite_real_ts <- ts.union(C_ts, I_ts, G_ts, NX_ts)
plot.ts(composite_real_ts,main = "Real GDP Components (Quarterly, 1947-2023)")
ts.plot(composite_real_ts,main = "Real GDP Components (Quarterly, 1947-2023)")

plot.ts(Y_real_ts,main="Real GDP (Quarterly, 1947-2023)")

plot.ts(chained_residual_ts,main="Real GDP Residual from Chain 2017 Conversion")

# Upon initial grataucal inspection, it appears that both real and nominal GDP
#   and their components have some type of time trend in the data. Given this,
#   I suspect that differencing would be appropriate. I will continue by performing
#   a Dickey-Fuller test on this data to confirm those suspicion
################################################################################
# Dickey-Fuller Tests

# Creates data frames that store the relevant tau and tau critical values based on 
# the number of observations in the data. Labels the row and column names accordingly
tau_critical <- data.frame()
tau_critical <- rbind(tau_critical,c("tau1",-2.58,-1.95,-1.62),c("tau2",-3.44,-2.87,-2.57),c("tau3",-3.98,-3.42,-3.13))
colnames(tau_critical) <- c("Test","OnePCT","FivePCT","TenPCT")


phi_critical <- data.frame()
phi_critical <- rbind(phi_critical,c("phi1",6.47,4.61,3.79),c("phi2",6.15,4.71,4.05),c("phi3",8.34,6.30,5.36))
colnames(phi_critical) <- c("Test","OnePCT","FivePCT","TenPCT")


# Performs Dickey-Fuller Tests on nominal GDP for each category of possible model
Y_test_trend <- ur.df(Y_ts,selectlags ='BIC', type='trend')
Y_test_drift <- ur.df(Y_ts,selectlags ='BIC', type='drift')
Y_test_none <- ur.df(Y_ts,selectlags ='BIC', type='none')

# Summarizes the first Dickey-Fuller Test to console
summary(Y_test_trend)

# Summarizes the second Dickey-Fuller Test to console
summary(Y_test_drift)

# Summarizes the final Dickey-Fuller Test to console
summary(Y_test_none)
# Records the test statistic values for later calculation
Y_tstats <- list("tau1" = 11.5295,"tau2" = 8.721, "tau3" = 3.4189,"phi1" = 66.6145, "phi2" = 44.2641, "phi3" = 37.9034)

# Initializes Data Frame for storing the results of test statistic analysis
Y_results <- data.frame(col=c("phi","tau"))

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
  if (Y_tstats[current_phi] >= phi_critical$OnePCT[i]){
    df_phi_test_result <- "Reject (1%)"
  } else if (Y_tstats[current_phi] >= phi_critical$FivePCT[i]){
    df_phi_test_result <- "Reject (5%)"
  } else if (Y_tstats[current_phi] >= phi_critical$TenPCT[i]){
    df_phi_test_result <- "Reject (10%)"
  } else {
    df_phi_test_result <- "Fail"
  }
  
  # Performs above logic again, but for the tauX test
  if (Y_tstats[current_tau] <= tau_critical$OnePCT[i]){
    df_tau_test_result <- "Reject (1%)"
  } else if (Y_tstats[current_tau] <= tau_critical$FivePCT[i]){
    df_tau_test_result <- "Reject (5%)"
  } else if (Y_tstats[current_tau] <= tau_critical$TenPCT[i]){
    df_tau_test_result <- "Reject (10%)"
  } else {
    df_tau_test_result <- "Fail"
  }
  
  # Records the results of the above logic in a new column of the Y_results data frame
  Y_results <- cbind(Y_results,c(df_phi_test_result,df_tau_test_result))
}
# Appropriately names the columns based on the phiX/tauX test performed
colnames(Y_results) <- c("Test",1:length(phi_critical$Test))
# Names the rows to match the phi or tau tests, respectively
 
# Outputs the results to the console
Y_results

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
Y_real_tstats <- list("tau1" = 8.6073,"tau2" = 3.4598, "tau3" = -0.9695,"phi1" = 38.7108, "phi2" = 26.8724, "phi3" = 7.3942)

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





