library(deSolve)
library(rstan)

# The present script calculates the original PBPK Index by Krishnan et al.1995

# Set the path containg the results of fit, the "Kreyling IV data.xlsx" and "Rat physiological parameters.xlsx" files
setwd("___")
load(file = "Fit TiO2 IV_results.RData")



###Prepare experimental data
data <- openxlsx::read.xlsx("Kreyling IV data.xlsx", sheet = 1, colNames = T, rowNames = T)
sd <- openxlsx::read.xlsx("Kreyling IV data.xlsx", sheet = 2, colNames = T, rowNames = T)
feces_exp <- openxlsx::read.xlsx("Kreyling IV data.xlsx", sheet = 3, colNames = T, rowNames = F)
urine_exp <- openxlsx::read.xlsx("Kreyling IV data.xlsx", sheet = 4, colNames = T, rowNames = F)
colnames(feces_exp) <- c("Time", "Feces") # Time in days, Feces in micro_g
colnames(urine_exp) <- c("Time", "Urine_excretion_rate", "Cumulative_ID") # Time in days
#Transform data
Transformed_data <- data
for (d in 1:length(doses)) {
  Transformed_data[d,] <- ((data[d,]/100)*doses[1])  #results give TiO2 in micro grams
  sd[d,] <- sd[d,]*doses[1]/100
}
#Drop Carcass column
Transformed_data <- subset(Transformed_data, select = -(Carcass))
Transformed_data <- Transformed_data[,c(4,2,1,3,5,6,7,9,10,8)]
rownames(Transformed_data) <- NULL  
Transformed_data <- as.matrix(Transformed_data)
sd <- subset(sd, select = -(Carcass))
feces_exp[,1] <- feces_exp[,1]*24 #transform time to hours
urine_exp[,1] <- urine_exp[,1]*24 #transform time to hours
feces_exp[,2] <- feces_exp[,2]*doses[1]/100
urine_exp <- subset(urine_exp, select = -(Urine_excretion_rate))
urine_exp[,2] <- urine_exp[,2]*doses[1]/100

urine_time <- urine_exp[,1] #in hours
feces_time <- feces_exp[,1] #in hours
exp_list <- list(Transformed_data, urine_exp[,2], feces_exp[,2])
names(exp_list) <- c("Transformed_data", "urine_exp", "feces_exp")
#################################################################################################################



stan_results <- extract(fit)
# Calculate the mean values of the posterior distributions
fitted_params <- c(uptake = mean(stan_results$theta[,1]),
                   uptake_spl = mean(stan_results$theta[,2]),
                   uptake_li = mean(stan_results$theta[,3]),
                   P_rest = mean(stan_results$theta[,4])
)

# Store in 1 vector both physiological and Substance-Specific parameters
params <- c(params, fitted_params)

# Define the dose
dose <- doses[1]
# Define the initial conditions of the ODEs system
inits <- c(Mcap_lu=0, Mcap_spl=0, Mcap_li=0, Mcap_ki=0, Mcap_ht=0, Mcap_br=0, Mcap_ut=0, Mcap_skel=0, Mcap_st=0,
           Mtis_lu=0, Mtis_spl=0, Mtis_li=0, Mtis_ki=0, Mtis_ht=0, Mtis_br=0, Mtis_ut=0, Mtis_skel=0, Mtis_st=0,
           Mm_lu=0, Mm_spl=0, Mm_li=0, Mm_ki=0, Mm_ht=0, Mm_br=0, Mm_ut=0, Mm_skel=0, Mm_st=0, 
           M_ven=dose, M_art=0, Mm_ven=0, Mm_art=0, M_feces=0, M_urine=0)


###############
# ODEs system #
###############
ode.func <- function(time, inits, params){
  with( as.list(c(inits,params)),{
    
    #Permeability coefficients (unitless)
    x_lu = xslow
    x_spl = xfast
    x_li = xfast
    x_ki = xslow
    x_ht = xslow
    x_br = xbr
    x_ut = xslow
    x_skel = xfast
    x_st = xslow
    #Partition Coefficients (unitless)
    P_lu = P_rest
    P_spl = P_rest
    P_li = P_li
    P_ki = P_rest
    P_ht = P_rest
    P_br = P_rest
    P_ut = P_rest
    P_skel = P_rest
    P_st = P_st
    #Capacity ug of NPs per g of PCs on organ
    uptake_spl = uptake_spl
    uptake_ki = uptake
    uptake_ht = uptake
    uptake_br = uptake
    uptake_ut = uptake
    uptake_skel = uptake
    uptake_st = uptake_st
    uptake_blood = uptake
    uptake_lu = uptake
    uptake_li = uptake_li
    
    #Capillary concentrations
    Ccap_lu <- Mcap_lu/Vcap_lu
    Ccap_spl <- Mcap_spl/Vcap_spl
    Ccap_li <- Mcap_li/Vcap_li
    Ccap_ki <- Mcap_ki/Vcap_ki
    Ccap_ht <- Mcap_ht/Vcap_ht
    Ccap_br <- Mcap_br/Vcap_br
    Ccap_ut <- Mcap_ut/Vcap_ut
    Ccap_skel <- Mcap_skel/Vcap_skel
    Ccap_st <- Mcap_st/Vcap_st
    
    #Tissue concentration
    Ctis_lu <- Mtis_lu/Vtis_lu
    Ctis_spl <- Mtis_spl/Vtis_spl
    Ctis_li <- Mtis_li/Vtis_li
    Ctis_ki <- Mtis_ki/Vtis_ki
    Ctis_ht <- Mtis_ht/Vtis_ht
    Ctis_br <- Mtis_br/Vtis_br
    Ctis_ut <- Mtis_ut/Vtis_ut
    Ctis_skel <- Mtis_skel/Vtis_skel
    Ctis_st <- Mtis_st/Vtis_st
    
    C_ven <- M_ven/Vven
    C_art <- M_art/Vart
    
    PA_lu<-x_lu*Q_total
    PA_spl<-x_spl*Q_spl
    PA_li<-x_li*Q_li
    PA_ki<-x_ki*Q_ki
    PA_ht<-x_ht*Q_ht
    PA_br<-x_br*Q_br
    PA_ut<-x_ut*Q_ut
    PA_skel<-x_skel*Q_skel
    PA_st<-x_st*Q_st
    
    P_lu <- P_lu 
    P_spl <- P_spl 
    P_li <- P_li 
    P_ki <- P_ki 
    P_ht <- P_ht 
    P_br <- P_br 
    P_ut <- P_ut 
    P_skel <- P_skel 
    P_st <- P_st
    
    Pup_lu<-Pup_max*(1-(Mm_lu/(Wm_lu*uptake_lu))) #1/h
    Pup_spl<-Pup_max*(1-(Mm_spl/(Wm_spl*uptake_spl))) #1/h
    Pup_li<-Pup_max*(1-(Mm_li/(Wm_li*uptake_li))) #1/h
    Pup_ki<-Pup_max*(1-(Mm_ki/(Wm_ki*uptake_ki))) #1/h
    Pup_ht<-Pup_max*(1-(Mm_ht/(Wm_ht*uptake_ht))) #1/h
    Pup_br<-Pup_max*(1-(Mm_br/(Wm_br*uptake_br))) #1/h
    Pup_ut<-Pup_max*(1-(Mm_ut/(Wm_ut*uptake_ut))) #1/h
    Pup_skel<-Pup_max*(1-(Mm_skel/(Wm_skel*uptake_skel))) #1/h
    Pup_st<-Pup_max*(1-(Mm_st/(Wm_st*uptake_st))) #1/h
    
    Pup_ven<-Pup_max*(1-(Mm_ven/(Wm_ven*uptake_blood))) #1/h
    Pup_art<-Pup_max*(1-(Mm_art/(Wm_art*uptake_blood))) #1/h
    
    #Lungs
    dMcap_lu <- Q_total*C_ven - Q_total*Ccap_lu - PA_lu*Ccap_lu + PA_lu*Ctis_lu/P_lu
    dMtis_lu <- PA_lu*Ccap_lu - PA_lu*Ctis_lu/P_lu - Pup_lu*Vtis_lu*Ctis_lu + k_de*Mm_lu
    dMm_lu   <- Pup_lu*Vtis_lu*Ctis_lu - k_de*Mm_lu
    
    #Spleen
    dMcap_spl <- Q_spl*C_art - Q_spl*Ccap_spl - PA_spl*Ccap_spl + PA_spl*Ctis_spl/P_spl
    dMtis_spl <- PA_spl*Ccap_spl - PA_spl*Ctis_spl/P_spl - Pup_spl*Vtis_spl*Ctis_spl + k_de*Mm_spl
    dMm_spl   <- Pup_spl*Vtis_spl*Ctis_spl - k_de*Mm_spl
    
    #Liver
    dMcap_li <- Q_li*C_art + Q_spl*Ccap_spl - (Q_li+Q_spl)*Ccap_li - PA_li*Ccap_li + PA_li*Ctis_li/P_li
    dMtis_li <- PA_li*Ccap_li - PA_li*Ctis_li/P_li - Pup_li*Vtis_li*Ctis_li + k_de*Mm_li - CLE_f*Mtis_li
    dMm_li   <- Pup_li*Vtis_li*Ctis_li - k_de*Mm_li
    dM_feces <- CLE_f*Mtis_li
    
    #Kidneys
    dMcap_ki <- Q_ki*C_art - Q_ki*Ccap_ki - PA_ki*Ccap_ki + PA_ki*Ctis_ki/P_ki - CLE_u*Mcap_ki
    dMtis_ki <- PA_ki*Ccap_ki - PA_ki*Ctis_ki/P_ki - Pup_ki*Vtis_ki*Ctis_ki + k_de*Mm_ki 
    dMm_ki   <- Pup_ki*Vtis_ki*Ctis_ki - k_de*Mm_ki
    dM_urine <- CLE_u*Mcap_ki
    
    #Heart
    dMcap_ht <- Q_ht*C_art - Q_ht*Ccap_ht - PA_ht*Ccap_ht + PA_ht*Ctis_ht/P_ht
    dMtis_ht <- PA_ht*Ccap_ht - PA_ht*Ctis_ht/P_ht - Pup_ht*Vtis_ht*Ctis_ht + k_de*Mm_ht
    dMm_ht   <- Pup_ht*Vtis_ht*Ctis_ht - k_de*Mm_ht
    
    #Brain
    dMcap_br <- Q_br*C_art - Q_br*Ccap_br - PA_br*Ccap_br + PA_br*Ctis_br/P_br
    dMtis_br <- PA_br*Ccap_br - PA_br*Ctis_br/P_br - Pup_br*Vtis_br*Ctis_br + k_de*Mm_br
    dMm_br   <- Pup_br*Vtis_br*Ctis_br - k_de*Mm_br
    
    #Uterus
    dMcap_ut <- Q_ut*C_art - Q_ut*Ccap_ut - PA_ut*Ccap_ut + PA_ut*Ctis_ut/P_ut
    dMtis_ut <- PA_ut*Ccap_ut - PA_ut*Ctis_ut/P_ut - Pup_ut*Vtis_ut*Ctis_ut + k_de*Mm_ut
    dMm_ut   <- Pup_ut*Vtis_ut*Ctis_ut - k_de*Mm_ut
    
    #Skeleton
    dMcap_skel <- Q_skel*C_art - Q_skel*Ccap_skel - PA_skel*Ccap_skel + PA_skel*Ctis_skel/P_skel
    dMtis_skel <- PA_skel*Ccap_skel - PA_skel*Ctis_skel/P_skel - Pup_skel*Vtis_skel*Ctis_skel + k_de*Mm_skel
    dMm_skel   <- Pup_skel*Vtis_skel*Ctis_skel - k_de*Mm_skel
    
    #Soft tissue
    dMcap_st <- Q_st*C_art - Q_st*Ccap_st - PA_st*Ccap_st + PA_st*Ctis_st/P_st
    dMtis_st <- PA_st*Ccap_st - PA_st*Ctis_st/P_st - Pup_st*Vtis_st*Ctis_st + k_de*Mm_st
    dMm_st   <- Pup_st*Vtis_st*Ctis_st - k_de*Mm_st
    
    
    #Veins
    dM_ven <- - Q_total*C_ven + (Q_li+Q_spl)*Ccap_li + Q_ki*Ccap_ki + Q_ht*Ccap_ht + Q_br*Ccap_br + Q_ut*Ccap_ut + Q_skel*Ccap_skel + 
      Q_st*Ccap_st  - Pup_ven*Vven*C_ven + k_de*Mm_ven
    dMm_ven <- Pup_ven*Vven*C_ven - k_de*Mm_ven
    
    #Arteries
    dM_art <- Q_total*Ccap_lu - Q_spl*C_art - Q_li*C_art - Q_ki*C_art - Q_ht*C_art - Q_br*C_art - Q_ut*C_art - Q_skel*C_art - Q_st*C_art - 
      Pup_art*Vart*C_art + k_de*Mm_art
    dMm_art <-Pup_art*Vart*C_art - k_de*Mm_art
    
    
    #Total amounts in each compartment
    #Lungs
    Lu_total <-  Mtis_lu + Mm_lu #Mcap_lu +
    
    #Spleen
    Spl_total <-  Mtis_spl + Mm_spl #Mcap_spl +
    
    #Liver
    Li_total <- Mtis_li + Mm_li # Mcap_li +
    
    #Kidneys
    Ki_total <-  Mtis_ki + Mm_ki #Mcap_ki +
    
    #Heart
    Ht_total <-  Mtis_ht + Mm_ht #Mcap_ht +
    
    #Brain
    Br_total <-  Mtis_br + Mm_br #Mcap_br +
    
    #Uterus
    Ut_total <-  Mtis_ut + Mm_ut #Mcap_ut +
    
    #Skeleton
    Skel_total <-  Mtis_skel + Mm_skel  #Mcap_skel +
    
    #Soft tissue
    St_total <-  Mtis_st + Mm_st #Mcap_st +
    
    #Blood
    Blood_total <- M_ven + Mm_ven + M_art + Mm_art
    
    Feces_total = M_feces
    Urine_total = M_urine
    
    list(c(dMcap_lu=dMcap_lu, dMcap_spl=dMcap_spl, dMcap_li=dMcap_li, dMcap_ki=dMcap_ki, dMcap_ht=dMcap_ht, dMcap_br=dMcap_br, dMcap_ut=dMcap_ut, dMcap_skel=dMcap_skel, dMcap_st=dMcap_st, 
           dMtis_lu=dMtis_lu, dMtis_spl=dMtis_spl, dMtis_li=dMtis_li, dMtis_ki=dMtis_ki, dMtis_ht=dMtis_ht, dMtis_br=dMtis_br, dMtis_ut=dMtis_ut, dMtis_skel=dMtis_skel, dMtis_st=dMtis_st,
           dMm_lu=dMm_lu, dMm_spl=dMm_spl, dMm_li=dMm_li, dMm_ki=dMm_ki, dMm_ht=dMm_ht, dMm_br=dMm_br, dMm_ut=dMm_ut, dMm_skel=dMm_skel, dMm_st=dMm_st, 
           dM_ven=dM_ven, dM_art=dM_art, dMm_ven=dMm_ven, dMm_art=dMm_art, dM_feces=dM_feces, dM_urine=dM_urine), 
         Lu_total=Lu_total, Spl_total=Spl_total, Li_total=Li_total, Ki_total=Ki_total, Ht_total=Ht_total, Br_total=Br_total, Ut_total=Ut_total, 
         Skel_total=Skel_total, St_total=St_total, Blood_total=Blood_total, Feces_total=Feces_total, Urine_total=Urine_total)
    
    
  })}
# Time points of integration
sample_time <- c(0, 10/60, 1, 1*24, 7*24, 28*24) #in hours
sample_time_urine <- c(0,urine_time)
sample_time_feces <- c(0, feces_time)
solution <- ode(times = sample_time, func = ode.func, y = inits, parms = params, method = "bdf")
urine_solution <- ode(times = sample_time_urine, func = ode.func, y = inits, parms = params, method = "bdf")
feces_solution <- ode(times = sample_time_feces, func = ode.func, y = inits, parms = params, method = "bdf")
# Total mass in each compartment
Total_amounts <- solution[2:length(sample_time),35:44]
Total_urine <- urine_solution[2:length(sample_time_urine),46]
Total_feces <- feces_solution[2:length(sample_time_feces),45]
sim_list <- list(Total_amounts, Total_urine, Total_feces)
names(sim_list) <- c("Total_amounts", "Total_urine", "Total_feces")

###########################################################################################################################################################################################################


############# Calculate PBPK indices #############
# pbpk.index a function the returns the compartment and consolidated (Total) discrepancy index
# of a PBPK model, given some experimental data. It follows the paper of Krishnan et al.1995.
# experimental: list of vectors containing the experimental data
# predictions: list of vectors containing the predicted data
# names of the compartments
pbpk.index <- function(observed, predicted, comp.names =NULL){
  # Check if the user provided the correct input format
  if (!is.list(observed) || !is.list(predicted)){
    stop(" The observations and predictions must be lists")
  }
  # Check if the user provided equal length lists
  if (length(observed) != length(predicted)){
    stop(" The observations and predictions must have the same compartments")
  }
  Ncomp <- length(observed) # Number of compartments
  I <- rep(NA, Ncomp) # Compartment discrepancy index
  N_obs <- rep(NA, Ncomp) #Number of observations per compartment
  #loop over the compartments
  for (i in 1:Ncomp){
    et <- 0
    Et <-0
    N <- length(observed[[i]]) # number of observations for compartment i
    # Check if observations and predictions have equal length
    if(N != length(predicted[[i]])){
      stop(paste0("Compartment ",i," had different length in the observations and predictions"))
    }
    N_obs[i] <- N # populate tne N_obs vector
    for (j in 1:N){
      # sum of absolute squared errors (error = observed - predicted)
      et <- et + (abs(observed[[i]][j] - predicted[[i]][j]))^2
      # Sum of squared observed values
      Et <- Et + (observed[[i]][j])^2
    }
    # root mean square of the absolute error
    RMet2 <-sqrt(et/N)
    # root mean of the square of observed values
    RMEt2 <- sqrt(Et/N)
    I[i] <- RMet2/RMEt2
  }
  # Total number of observations
  Ntot <- sum(N_obs)
  # Initialise the consolidated discrepancy index
  Ic <-0
  for (i in 1:Ncomp){
    Ic <- Ic +  I[i]* N_obs[i]/Ntot
  }
  # Name the list of compartment discrepancy indices
  if ( !is.null(comp.names)){
    names(I) <- comp.names
  }else if (!is.null(names(observed))){
    names(I) <- names(observed)
  } else if (!is.null(names(predicted)) && is.null(comp.names) ){
    names(I) <- names(predicted)
  }
  
  return(list(Total_index = Ic, Compartment_index= I))
}


#prepare input
observed <- list()
predicted <- list()
for (i in 1:dim(Total_amounts)[2]){
  observed[[i]] <- Transformed_data[,i]
  predicted[[i]] <- Total_amounts[,i]
}
observed[[i+1]] <- urine_exp[,2]
observed[[i+2]] <-feces_exp[,2]
predicted[[i+1]] <- Total_urine
predicted[[i+2]] <-Total_feces

names(predicted) <- c(colnames(Total_amounts), "urine", "feces")

discrepancy <- pbpk.index(observed, predicted)