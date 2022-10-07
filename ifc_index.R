
#         Author: Pedro Salas-Rojo 
#         Date: 06/10/2022
#         Name of project: Adapt Marco Ranaldi's Stata code to R (IFC Index)

# Ranaldi, M. (2021) Income Composition Inequality. Review of Income and Wealth, 68 (1), 139-160.

library(tidyverse)

# Function IFC ----

ifc <- function(data, var1, var2, wts = NULL){
  
  data <- data
  
  if(is.null(wts)){
    data$wts <- rep(1, times = nrow(data))
    print("No weights provided. A vector of 1s is used as weights.")
  } else {
    data$wts <- data[[wts]]
  }
  
  # Generate the variables total income, capital and labor income
  
  data$SEI <- data[[var1]]
  data$EI <- data[[var2]]
  data$TEI <- data$SEI + data$EI
  
  # Generate the weighted cumulative shares of the concentration curves 
  # for capital and labor à la Kakwani (scaled up to 1)
  
  data <- data %>%
    arrange(TEI)
  
  data$rank <- cumsum(data$wts)/sum(data$wts)
  data$TEI_mean <- weighted.mean(data$TEI, data$wts)
  data$TEI_Tot <- data$TEI_mean*length(data$TEI)
  data$TEI_sh_cum <- cumsum(data$wts*data$TEI)/sum(data$wts*data$TEI)
  
  data$SEI_mean <- weighted.mean(data$SEI, data$wts)
  data$SEI_Tot <- data$SEI_mean*length(data$SEI)
  data$SEI_sh_cum <- cumsum(data$wts*data$SEI)/sum(data$wts*data$SEI)
  
  data$EI_mean <- weighted.mean(data$EI, data$wts)
  data$EI_Tot <- data$EI_mean*length(data$EI)
  data$EI_sh_cum <- cumsum(data$wts*data$EI)/sum(data$wts*data$EI)
  
  # Generate the weighted survey capital and labor shares
  
  data$SEI_Tot_sh = data$SEI_Tot/data$TEI_Tot
  data$EI_Tot_sh = data$EI_Tot/data$TEI_Tot
  
  # Generate the weighted cumulative shares of the zero concentration curves 
  # (scaled up to the factor shares)
  
  data$Zeroconc_SEI_cum = data$TEI_sh_cum*data$SEI_Tot_sh
  data$Zeroconc_EI_cum = data$TEI_sh_cum*data$EI_Tot_sh
  
  # Generate the weighted cumulative shares of the concentration curves for 
  # capital and labor (scaled up to the factor shares)	
  
  data$SEI_sh_cum_2 = data$SEI_sh_cum*data$SEI_Tot_sh
  data$EI_sh_cum_2 = data$EI_sh_cum*data$EI_Tot_sh
  
  # Create empty variables for the maximum concentration curves
  
  data$Maxconc_SEI = NA
  data$Maxconc_EI = NA
  
  # Rank according to "rank", get individual and individual share
  
  data <- data %>%
    arrange(rank)
  data = tibble::rowid_to_column(data, "ind")
  data$ind_sh = data$ind/length(data$ind)
  
  # Create lagged variables.
  
  data$lagZeroconc_SEI_cum = lag(data$Zeroconc_SEI_cum)
  data$lagZeroconc_EI_cum = lag(data$Zeroconc_EI_cum)				
  data$lagSEI_sh_cum_2 = lag(data$SEI_sh_cum_2)
  data$lagEI_sh_cum_2 = lag(data$EI_sh_cum_2)
  data$lagTEI_sh_cum = lag(data$TEI_sh_cum)
  data$lagEI_sh_cum = lag(data$EI_sh_cum)
  data$lagSEI_sh_cum = lag(data$SEI_sh_cum)	
  
  # Generate the difference with respect to the lag variable
  data$A = data$lagZeroconc_SEI_cum + data$Zeroconc_SEI_cum
  data$B = data$lagZeroconc_EI_cum + data$Zeroconc_EI_cum
  data$C = data$lagSEI_sh_cum_2 + data$SEI_sh_cum_2
  data$D = data$lagEI_sh_cum_2 + data$EI_sh_cum_2	
  data$E = data$lagTEI_sh_cum + data$TEI_sh_cum
  data$F = data$lagEI_sh_cum + data$EI_sh_cum
  data$G = data$lagSEI_sh_cum + data$SEI_sh_cum
  
  # Generate the sum of all h
  data$H = sum(data$A, na.rm = TRUE)
  data$I = sum(data$B, na.rm = TRUE)
  data$L = sum(data$C, na.rm = TRUE)
  data$M = sum(data$D, na.rm = TRUE)
  data$N = sum(data$E, na.rm = TRUE)
  data$O = sum(data$F, na.rm = TRUE)
  data$P = sum(data$G, na.rm = TRUE)		
  
  if(mean(data$H, na.rm = TRUE) > mean(data$L, na.rm = TRUE)){
    
    data$Maxconc_SEI <- ifelse(data$TEI_sh_cum<data$EI_Tot_sh, 
                               0, 
                               data$Maxconc_SEI)
    data$Maxconc_SEI <- ifelse(data$TEI_sh_cum>=data$EI_Tot_sh,
                               data$TEI_sh_cum-data$EI_Tot_sh,
                               data$Maxconc_SEI)
    
  } else {
    
    data$Maxconc_SEI <- ifelse(data$TEI_sh_cum<data$SEI_Tot_sh, 
                               data$TEI_sh_cum, 
                               data$Maxconc_SEI)
    data$Maxconc_SEI <- ifelse(data$TEI_sh_cum>=data$SEI_Tot_sh,
                               data$SEI_Tot_sh,
                               data$Maxconc_SEI)
    
  }
  
  data$lagMaxconc_SEI <- lag(data$Maxconc_SEI)	
  data$Q <- data$lagMaxconc_SEI + data$Maxconc_SEI
  data$R <- sum(data$Q, na.rm = TRUE)
  
  if(mean(data$I, na.rm = TRUE) > mean(data$M, na.rm = TRUE)){
    
    data$Maxconc_EI <- ifelse(data$TEI_sh_cum<data$SEI_Tot_sh, 
                              0, 
                              data$Maxconc_EI)
    data$Maxconc_EI <- ifelse(data$TEI_sh_cum>=data$SEI_Tot_sh,
                              data$TEI_sh_cum-data$SEI_Tot_sh,
                              data$Maxconc_EI)
    
  } else {
    
    data$Maxconc_EI <- ifelse(data$TEI_sh_cum<data$EI_Tot_sh, 
                              data$TEI_sh_cum, 
                              data$Maxconc_EI)
    data$Maxconc_EI <- ifelse(data$TEI_sh_cum>=data$EI_Tot_sh,
                              data$EI_Tot_sh,
                              data$Maxconc_EI)
     
  }
  
  data$lagMaxconc_EI <- lag(data$Maxconc_EI)	
  data$S <- data$lagMaxconc_EI + data$Maxconc_EI
  data$T <- sum(data$S, na.rm = TRUE)
  
  # Gini Coefficient
  
  gini <- mean(1- (1/length(data$ind))*data$N)
  
  # Income-Factor Concentration Index (ICF)
  
  ifc  <- mean(data$H - data$L)/mean(abs(data$H - data$R))
  
  # Area of the concentration curve for capital income
  
  mu_p <- mean(1/(length(data$ind)*2)*data$P)
  
  # Area of the concentration curve for labor income
  
  mu_w <- mean(1/(length(data$ind)*2)*data$O)
  
  
  return(list(`data` = data, `gini` = gini, `ifc` = ifc, 
              `mu_var1` = mu_p, `mu_var2` = mu_w))
  
}

# Get data ----

data <- read.csv("your_path/data_ifc.csv")

# Get results ----

# Weights should be analytics. If no weights are provided, use NULL.

results <- ifc(data, var1 = "cap_inc", var2 = "lab_inc", wts = "weights")

data_ifc <- results[["data"]]
print(results[["gini"]])
print(results[["ifc"]])
print(results[["mu_var1"]])
print(results[["mu_var2"]])

results <- ifc(data, var1 = "financial_assets", var2 = "nonfinanc_assets", wts = NULL)

data_ifc <- results[["data"]]
print(results[["gini"]])
print(results[["ifc"]])
print(results[["mu_var1"]])
print(results[["mu_var2"]])
