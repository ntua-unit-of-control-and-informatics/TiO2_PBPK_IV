library(rstan)
library(ggplot2)
library(parallel)

options(max.print=1e05)
#Set the path of "Kreyling IV data.xlsx" and "Rat physiological parameters.xlsx" files
setwd("___")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

####################################################################
###################             DATA             ###################
####################################################################

doses <- c(18.15, 11.29, 16.53, 108.5, 46.92)

### Important!!! each compartment has a specific index vectors Tissue_fractions, Regional_flow_fractions, Capillary_fractions and cannot be changed
# The index of each compartment:
#Rest of Body (rob) --> 1
#Heart (ht) --> 2
#Kidneys (ki) --> 3
#Brain (br) --> 4
#Spleen (spl) --> 5
#Lungs (lu) --> 6
#Liver (li) --> 7
#Uterus (ut) --> 8
#Skeleton (skel) --> 9
#Adipose (ad) --> 10
#Skin (skin) --> 11
#Muscles (mu) --> 12


####################
### User's INPUT ###
####################
#### If any of these compartments don not exist in pbpk, just give it the value NA in compartments vector, example: "Heart" = NA and it will remove it 
#### from the equilibriums and the corresponding V_tis, V_cap, Q will be equal to NA.


compartments <- list( "St"="St","Heart"="Heart", "Kidneys"="Kidneys", "Brain"="Brain", "Spleen"="Spleen",
                      "Lungs"="Lungs", "Liver"="Liver", "Uterus"="Uterus", "Skeleton"="Skeleton", "Adipose"=NA, "Skin"=NA, "Muscles"=NA) #used as input in function, compartments that are used in pbpk

BW <- 263 # Total Body weight of rat in g

#####################################
### Function to create Parameters ###
#####################################

create.params <- function(comp_names, w){
  
  # List with names of all possible compartments
  all_comps <- list("St"="St","Heart"="Heart", "Kidneys"="Kidneys", "Brain"="Brain", "Spleen"="Spleen",
                    "Lungs"="Lungs", "Liver"="Liver", "Uterus"="Uterus", "Skeleton"="Skeleton", "Adipose"="Adipose", "Skin"="Skin", "Muscles"="Muscles") # List with names of all possible compartments
  
  ### Density of tissues/organs
  d_tissue <- 1 #g/ml
  d_skeleton <- 1.92 #g/ml
  d_adipose <- 0.940 #g/ml
  
  Q_total <- (1.54*w^0.75)*60 # Total Cardiac Output (ml/h)
  
  Total_Blood <- 0.06*w+0.77 # Total blood volume (ml)
  
  #Arterial blood volume
  #Vart <- 1.2905*w/100 #0.15*Total_Blood #(ml)
  
  #Veins blood volume
  #Vven <- 2.968*w/100 #0.64*Total_Blood #(ml)
  
  fr_ad <- 0.0199*w + 1.644 # w in g,  Brown et al.1997 p.420. This equation gives the  adipose % of body weight 
  
  #read data from excel
  fractions <- openxlsx::read.xlsx("Rat physiological parameters.xlsx", sheet = 1, colNames = T, rowNames = T)
  fractions <- as.matrix(sapply(fractions, as.numeric))
  
  #Tissue weight fraction 
  Tissue_fractions <- fractions[,1]/100 # % of BW. Na values refers to the volume of the rest organs(RoB)
  Tissue_fractions[10] <- fr_ad/100
  #Regional blood flow fraction
  Regional_flow_fractions <- fractions[,2]/100 # % of total cardiac output
  #Capillary volume fractions (fractions of tissue volume)
  Capillary_fractions <- fractions[,3] # of tissue volume
  #Macrophage content as fraction tissue volume for each tissue/organ
  Macrophage_fractions <- fractions[,4] 
  
  W_tis <- rep(0,length(comp_names))
  V_tis <- rep(0,length(comp_names))
  V_cap <- rep(0,length(comp_names))
  W_macro <- rep(0,length(comp_names))  #one more for blood compartment
  Q <- rep(0,length(comp_names))
  
  
  for (i in 1:length(comp_names)) {
    control <- comp_names[i]
    
    Tissue_fractions[i] <- ifelse(is.na(control), NA, Tissue_fractions[i])
    Regional_flow_fractions[i] <- ifelse(is.na(control), NA, Regional_flow_fractions[i])
    Capillary_fractions[i] <- ifelse(is.na(control), NA, Capillary_fractions[i])
    Macrophage_fractions[i] <- ifelse(is.na(control), NA, Macrophage_fractions[i])
    
    ### Calculation of tissue weights  
    W_tis[i] <- w*Tissue_fractions[i]
    
    
    ###Calculation of tissue volumes
    
    if (i==9){
      V_tis[i] <- W_tis[i]/d_skeleton
    } else if(i==10){
      V_tis[i] <- W_tis[i]/d_adipose
    } else{
      V_tis[i] <- W_tis[i]/d_tissue 
    }
    
    ###Calculation of capillary volumes
    V_cap[i] <- V_tis[i]*Capillary_fractions[i]
    
    ###Volume of macrophage contents
    W_macro[i] <- W_tis[i]*Macrophage_fractions[i]
    
    ###Calculation of regional blood flows
    Q[i] <- Q_total*Regional_flow_fractions[i]
  }
  
  #Vm_ven <- 0.01*Vven #macrophage content in veins
  #Vm_art <- 0.01*Vart #0.02*Vart #macrophage content in arteries
  
  ### Calculations for "Soft tissue" compartment
  W_tis[1] <- w - sum(W_tis[2:length(W_tis)], na.rm = TRUE)
  V_tis[1] <- W_tis[1]/d_adipose     
  Q[1] <- Q_total - sum(Q[2:length(Q)],na.rm = TRUE) + Q[6]
  V_cap[1] <- V_tis[1]*Capillary_fractions[1] #Total_Blood - Vven - Vart - sum(V_cap[2:length(V_cap)], na.rm = TRUE)
  W_macro[1] <- W_tis[1]*Macrophage_fractions[1]
  #Capillary_fractions[1] <- V_cap[1]/V_tis[1]
  
  
  parameters <- matrix(c(W_tis[],V_tis[],V_cap[],Q[],W_macro[]), ncol = 5)
  colnames(parameters) <- c("W_tis", "V_tis", "V_cap", "Q", "W_macro")
  rownames(parameters) <- all_comps
  
  Vven=0.64*Total_Blood
  Vart=0.15*Total_Blood
  Wm_ven=0.01*Vven
  Wm_art=0.01*Vart
  
  return(c(
    "Q_total"=Q_total, "V_blood"=Total_Blood, "Vven"=Vven, "Vart"=Vart, "Wm_ven"=Wm_ven, "Wm_art"=Wm_art,
    
    "W_st"=parameters[1,1], "W_ht"=parameters[2,1], "W_ki"=parameters[3,1], "W_br"=parameters[4,1], "W_spl"=parameters[5,1], "W_lu"=parameters[6,1], "W_li"=parameters[7,1], "W_ut"=parameters[8,1], "W_skel"=parameters[9,1],
    
    "Vtis_st"=parameters[1,2], "Vtis_ht"=parameters[2,2], "Vtis_ki"=parameters[3,2], "Vtis_br"=parameters[4,2], "Vtis_spl"=parameters[5,2], "Vtis_lu"=parameters[6,2], "Vtis_li"=parameters[7,2], "Vtis_ut"=parameters[8,2], "Vtis_skel"=parameters[9,2],
    
    "Vcap_st"=parameters[1,3], "Vcap_ht"=parameters[2,3], "Vcap_ki"=parameters[3,3], "Vcap_br"=parameters[4,3], "Vcap_spl"=parameters[5,3], "Vcap_lu"=parameters[6,3], "Vcap_li"=parameters[7,3], "Vcap_ut"=parameters[8,3], "Vcap_skel"=parameters[9,3],
    
    "Wm_st"=parameters[1,5], "Wm_ht"=parameters[2,5], "Wm_ki"=parameters[3,5], "Wm_br"=parameters[4,5], "Wm_spl"=parameters[5,5], "Wm_lu"=parameters[6,5], "Wm_li"=parameters[7,5], "Wm_ut"=parameters[8,5], "Wm_skel"=parameters[9,5],
    
    "Q_st"=parameters[1,4], "Q_ht"=parameters[2,4], "Q_ki"=parameters[3,4], "Q_br"=parameters[4,4], "Q_spl"=parameters[5,4], "Q_lu"=parameters[6,4], "Q_li"=parameters[7,4], "Q_ut"=parameters[8,4], "Q_skel"=parameters[9,4],
    
    xfast = 3,
    xslow = 1e-4,
    xbr = 1e-5,
    Pup_max = 20,#82,
    uptake_st = 1e-2,
    k_de = 4.9*1e-19, #1/h  desorption rate\P_li <<- 20,
    P_st = 6e-2,
    
    CLE_u = 2.4e-1,
    P_li = 20
  ))
}
# Create the parametric vector (Physiological parameters and Substance-Specific parameters
# considered as fixed values)  
params<-create.params(compartments,BW)


# Create a vector containig the mean values of the stohastic parameters --> Parameters to be estimated through the data 
eta_tr <- c(uptake <<- 1e-1, 
            uptake_spl <<- 3,
            uptake_li <<- 16,
            P_rest <<- 0.8,
            CLE_f  <<- 1e-3)

# Define the Coefficient of Variance (CV) of each stohastic parameter
cv <- c(2, 0.5, 0.5, 1, 2)
uptake_std <- uptake*cv[1]
uptake_spl_std <- uptake_spl*cv[2]
uptake_li_std <- uptake_li*cv[3]
P_rest_std <- P_rest * cv[4]
CLE_f_std <- CLE_f * cv[5]
eta_tr_std <- c( uptake_std, uptake_spl_std, uptake_li_std, P_rest_std, CLE_f_std)

# Times of experimental data 
#sampling time
time <- c(10/60, 1, 1*24, 7*24, 28*24) #in hours

#######################
###Experimental data###
#######################
data <- openxlsx::read.xlsx("Kreyling IV data.xlsx", sheet = 1, colNames = T, rowNames = T)
sd <- openxlsx::read.xlsx("Kreyling IV data.xlsx", sheet = 2, colNames = T, rowNames = T)
feces <- openxlsx::read.xlsx("Kreyling IV data.xlsx", sheet = 3, colNames = T, rowNames = F)
urine <- openxlsx::read.xlsx("Kreyling IV data.xlsx", sheet = 4, colNames = T, rowNames = F)
colnames(feces) <- c("Time", "Feces") # Time in days, Feces in micro_g
colnames(urine) <- c("Time", "Urine_excretion_rate", "Cumulative_ID") # Time in days
Transformed_data <- data
for (d in 1:length(doses)) {
  Transformed_data[d,] <- ((data[d,]/100)*doses[1])  #results give TiO2 in micro grams
  sd[d,] <- sd[d,]*doses[1]/100
}
#Drop Carcass column
Transformed_data <- subset(Transformed_data, select = -(Carcass))
sd <- subset(sd, select = -(Carcass))
feces[,2] <- feces[,2]*doses[1]/100
feces_time <- feces[,1]*24 #transform time to hours
feces <- feces[,2]
urine <- subset(urine, select = -(Urine_excretion_rate))
urine_time <- urine[,1]*24 #transform time to hours
urine <- urine[,2]*doses[1]/100


#############################################################################################
# Prepare the list of data that are needed in Stan 
DataList <- list(
  eta_tr = eta_tr, # Mean values of stohastic parameters vector
  eta_tr_std = eta_tr_std, # Standard deviation of stohastic parameters vector
  N_diff = 33, # Number of ODEs
  N_compart= 10, # Number of compartments
  N_param = length(eta_tr),
  params = params,
  time = time,
  feces_time = feces_time,
  urine_time = urine_time,
  mass = Transformed_data, #MassData,
  feces = feces,
  urine = urine,
  t_init = 0, # initial time point
  m0 = c(rep(0,27),doses[1], rep(0,5)), #micrograms #initial conditions of ODEs
  rel_tol = 5e-08, 
  abs_tol = 5e-08,
  max_num_steps = 1e05)

tic = proc.time()

fit <- stan(file = 'Fit TiO2 IV.stan', data = DataList, iter = 1000, warmup=500, chains=4, control = list(adapt_delta = 0.9))
options(max.print=5.5E5) 

#check_hmc_diagnostics(fit)
#check_divergences(fit)
#check_treedepth(fit)
#check_energy(fit)
#if I had set  control=list(max_treedepth=15) the the correct command would be check_treedepth(fit, 15) 

#library(shinystan)
#launch_shinystan(fit)
#print(fit)


#clock<-proc.time() - tic
#print(clock)
#pairs(fit, pars=c("par"))
#traceplot(fit,c("epsilon"));
#stan_dens(fit,("par[1]"),separate_chains = TRUE)
#exp_fit <- extract(fit)
#mean(exp_fit$theta[,1,2])
