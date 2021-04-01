library(deSolve)
library(sensitivity)
library(Bolstad)
library(pracma)
library(ggplot2)
library(reshape2)
options(max.print = 1e5)

#Set the path of "Kreyling IV data.xlsx" and "Rat physiological parameters.xlsx" files
setwd("___")

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
    
    "Q_st"=parameters[1,4], "Q_ht"=parameters[2,4], "Q_ki"=parameters[3,4], "Q_br"=parameters[4,4], "Q_spl"=parameters[5,4], "Q_lu"=parameters[6,4], "Q_li"=parameters[7,4], "Q_ut"=parameters[8,4], "Q_skel"=parameters[9,4]
    
    
    
  ))
}
params<-create.params(compartments,BW)
init_params<-params

# Declare the mean value of each Substance-Specific parameter 
groupped.prior <- c(
  xfast <<- 3,
  xslow <<- 1e-4,
  x_br <<- 1e-5,
  P_li <<- 20,
  P_st <<- 6e-2,
  P_rest <<- 0.8,
  Pup_max <<- 20,#82,
  CLE_f<<- 6e-4,# 1/h
  CLE_u<<-2.4e-1,# 1/h
  uptake <<- 1e-1, 
  uptake_spl <<- 3,
  uptake_li <<- 16,
  uptake_st <<- 1e-2,
  k_de <<- 4.9e-19
)

# "init_prior" is the vector of all Substance-Specific parameters, 
# considering that each compartment has its own independent parameters
init_prior <- c(
  x_lu = xslow,
  x_spl = xfast,
  x_li = xfast,
  x_ki = xslow,
  x_ht = xslow,
  x_br = x_br,
  x_ut = xslow,
  x_skel = xfast,
  x_st = xslow,
  P_lu = P_rest, 
  P_spl = P_rest, 
  P_li = P_li ,
  P_ki = P_rest, 
  P_ht = P_rest,
  P_br = P_rest ,
  P_ut = P_rest ,
  P_skel = P_rest, 
  P_st = P_st,
  uptake_spl = uptake_spl,
  uptake_ki = uptake,
  uptake_ht = uptake,
  uptake_br = uptake,
  uptake_ut = uptake,
  uptake_skel = uptake,
  uptake_st = uptake_st,
  uptake_blood = uptake,
  uptake_lu = uptake,
  uptake_li = uptake_li,
  Pup_max = Pup_max,
  CLE_f = CLE_f,
  CLE_u = CLE_u,
  k_de = k_de
)


###############
# ODEs system #
###############
ode.func <- function(time, inits, params){
  with( as.list(c(inits,params)),{
    
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
    
    #uptake_spl = uptake_spl
    #uptake_ki = uptake_ki
    #uptake_ht = uptake_ht
    #uptake_br = uptake_br
    #uptake_ut = uptake_ut
    #uptake_skel = uptake_skel
    #uptake_st = uptake_st
    #uptake_blood = uptake_blood
    #uptake_lu = uptake_lu
    #uptake_li = uptake_li
    
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
    Li_total <-  Mtis_li + Mm_li #Mcap_li +
    
    #Kidneys
    Ki_total <-  Mtis_ki + Mm_ki #Mcap_ki +
    
    #Heart
    Ht_total <-  Mtis_ht + Mm_ht  #Mcap_ht +
    
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

# Create the "PBPK_function". The input of PBPK_function is a data frame,
# that each row is a sample for the Substance_specific parameters. So the input "X" 
# is a data frame [N_samples x N_params] created by the GSA method
PBPK <- function(X){
  # Define the dose 
  dose <- doses[1]
  # Define the initial conditions 
  inits <- c(Mcap_lu=0, Mcap_spl=0, Mcap_li=0, Mcap_ki=0, Mcap_ht=0, Mcap_br=0, Mcap_ut=0, Mcap_skel=0, Mcap_st=0,
             Mtis_lu=0, Mtis_spl=0, Mtis_li=0, Mtis_ki=0, Mtis_ht=0, Mtis_br=0, Mtis_ut=0, Mtis_skel=0, Mtis_st=0,
             Mm_lu=0, Mm_spl=0, Mm_li=0, Mm_ki=0, Mm_ht=0, Mm_br=0, Mm_ut=0, Mm_skel=0, Mm_st=0, 
             M_ven=dose, M_art=0, Mm_ven=0, Mm_art=0, M_feces=0, M_urine=0)
  #sample_time <- c(0, 10/60, 1, 1*24, 7*24, 28*24) #in hours
  sample_time <- seq(0, 28*24, 1) #in hours
  
  # Create a matrix to store the results of the ODEs system for each parametric sample
  output <- matrix(0, nrow = nrow(X), ncol = 12)
  
  # Define the CV of each parameter
  cv <- rep(NA, length(init_prior))
  for (p in 1:length(init_prior)) {
    cv[p] <- 50/100
  }
  
  
  eta <- init_prior
  std <- cv*eta 
  #Tranform from Normal to lognormal
  eta_tr <- log((eta^2)/sqrt((eta^2)+std^2))
  std_tr <- sqrt(log(1 + (std^2)/(eta^2)))
  
  # Calculate the solution of the ODEs sytem, considering the initial mean values of the parameters.
  # This solution will replace the "NA" values which will occur from solution of the ODEs system 
  # for some parametric samples
  average_run <-  solution <- ode(times = sample_time, func = ode.func, y = inits,
                                  parms = c(eta,init_params), method = "bdf")
  error_counter = 0 # a counter to count how many times the integration failed
  success <- rep(NA, nrow(X))
  for (i in 1:nrow((X))) {
    print(paste("We are at iteration", i, sep = "  "))
    
    # Inverse transform the sample from uniform[0,1] to the lognormal distributions of the Substance-Specific parameters
    u <- as.numeric(X[i,])
    inv_cdf <- exp(eta_tr + sqrt(2)*std_tr*erfinv(2*u - 1)) #This is the inversed CDF function for lognormal
    names(inv_cdf) <-names(init_prior)
    params <- c(inv_cdf,init_params) #store in 1 parametric vector the physiological and the substance-specific parameters
    solution <- ode(times = sample_time, func = ode.func, y = inits, parms = params, method = "bdf")
    
    # Check if the current parametric sample gave a "normal" ODEs solution (not "NA" values).
    # If not, then replace this solution with the average solution (calculated above) in the output matrix
    if (dim(solution)[1] != dim(average_run)[1]){
      solution <- average_run
      error_counter = error_counter + 1
      success[i] = FALSE  
    }else{
      success[i] = TRUE
    }
    Total_amounts <- solution[,35:46]
    
    for (k in 1:ncol(Total_amounts)) {
      #Integrate the "Time" dimension by calculating the AUC of the "Total" curve of each compartment. 
      AUC <- sintegral(sample_time, Total_amounts[,k], n.pts = 256)
      output[i,k] <- as.numeric(AUC["value"])
    }
    success <<- success
  }
  #Return a matrix that each row contains the AUC results for all the compartments, of the corresponding parametric sample
  return(output)
}
tic = proc.time()
# Define the GSA method
# It is considered that the sampling distribution is uniform(0,1) and the sample is transformed to 
# lognormal inside the PBPK function. In addition, it is strongly recommended to set "scale=TRUE", in order to scale
# the output matrix
morris_results <- morris(model = PBPK, factors = length(init_prior), r = 50,
                         design = list(type = "oat", levels = 10, grid.jump = 5),
                         binf = 1e-10, bsup = 1-1e-10, scale = TRUE)
clock<-proc.time() - tic 
##########################################################################################################
# Export the "mu" values for each parameter
mu <- apply(morris_results$ee, 3, function(M){
  apply(M, 2, mean)
})
# Export the "mu*" values for each parameter
mu.star <- apply(abs(morris_results$ee), 3, function(M){
  apply(M, 2, mean)
})
# Export the "sigma" values for each parameter
sigma <- apply(morris_results$ee, 3, function(M){
  apply(M, 2, sd)
})
row.names(mu) <- c("x_lu","x_spl","x_li", "x_ki", "x_ht", "x_br", "x_ut", "x_skel", "x_st", "P_lu",  "P_spl", 
                   "P_li","P_ki", "P_ht", "P_br", "P_ut", "P_skel",  "P_st" , "uptake_spl", "uptake_ki", "uptake_ht", "uptake_br",
                   "uptake_ut", "uptake_skel", "uptake_st", "uptake_blood","uptake_lu", "uptake_li", "Pup_max", "CLE_f", "CLE_u",
                   "k_de")
row.names(mu.star) <- row.names(mu)
row.names(sigma) <- row.names(mu)
colnames(mu) <- c("Lungs", "Spleen", "Liver", "Kidneys", "Heart", "Brain", "Uterus", "Skeleton", "Soft tissue", "Blood", "Feces", "Urine")
colnames(mu.star) <- colnames(mu)
colnames(sigma) <- colnames(mu)

#Calculate Global Index (GI) metric.
# This metric combines both mu* and sigma of each parameter, in order to have one
# sensitivity index for each parameter 
GI <- sqrt(mu.star^2 + sigma^2)


# Integration to eliminate the "Compartments" parameters.
# The calculated weights are based on the AUC of mass collected in each compartment,
# according to the experimental data

#Calculate weights
###Prepare experimental data
#read data
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
sd <- subset(sd, select = -(Carcass))
feces_exp[,1] <- feces_exp[,1]*24 #transform time to hours
urine_exp[,1] <- urine_exp[,1]*24 #transform time to hours
feces_exp[,2] <- feces_exp[,2]*doses[1]/100
urine_exp <- subset(urine_exp, select = -(Urine_excretion_rate))
urine_exp[,2] <- urine_exp[,2]*doses[1]/100


# Time of sampling
time <-  c(1, 4, 24, 7*24, 28*24) #in hours
urine_time <- urine_exp[,1]*24 #in hours
feces_time <- feces_exp[,1]*24 #in hours


Liver_data <- as.data.frame(cbind(time,Transformed_data[,1], sd[,1]))
Spleen_data <- as.data.frame(cbind(time,Transformed_data[,2], sd[,2]))
Kidneys_data <- as.data.frame(cbind(time,Transformed_data[,3], sd[,3]))
Lungs_data <- as.data.frame(cbind(time,Transformed_data[,4], sd[,4]))
Heart_data <- as.data.frame(cbind(time,Transformed_data[,5], sd[,5]))
Brain_data <- as.data.frame(cbind(time,Transformed_data[,6], sd[,6]))
Uterus_data <- as.data.frame(cbind(time,Transformed_data[,7], sd[,7]))
Blood_data <- as.data.frame(cbind(time,Transformed_data[,8], sd[,8]))
Skeleton_data <- as.data.frame(cbind(time,Transformed_data[,9], sd[,9]))
Soft_tissue_data <- as.data.frame(cbind(time,Transformed_data[,10], sd[,10]))
Feces_data <- feces_exp
Urine_data <- urine_exp


colnames(Liver_data) <- c("Time", "Mass_data", "SD")
colnames(Spleen_data) <- colnames(Liver_data)
colnames(Kidneys_data) <- colnames(Liver_data)
colnames(Lungs_data) <- colnames(Liver_data)
colnames(Heart_data) <- colnames(Liver_data)
colnames(Brain_data) <- colnames(Liver_data)
colnames(Uterus_data) <- colnames(Liver_data)
colnames(Blood_data) <- colnames(Liver_data)
colnames(Skeleton_data) <- colnames(Liver_data)
colnames(Soft_tissue_data) <- colnames(Liver_data)
colnames(Feces_data) <- c("time", "Feces")
colnames(Urine_data) <- c("time", "Urine")

# A list with the experimental data of each compartment
observed_data <- list(Lungs_data, Spleen_data, Liver_data, Kidneys_data, Heart_data, Brain_data, Uterus_data,
                      Skeleton_data, Soft_tissue_data, Blood_data, Feces_data, Urine_data)

# Calculate the AUC of the observed mass (experimental data points) of TiO2 in each compartment
AUC_obs <- c(rep(0,12))
counter <- 1
for (data in observed_data) {
  xy_values <- data[,1:2] # The time and mass values of the model
  AUC_integration <- sintegral(xy_values[,1], xy_values[,2], n.pts = 256) #integrate the curve
  AUC_obs[counter] <- as.numeric(AUC_integration["value"]) #Taek only the AUC value
  counter <- counter+1
}
# Prepare the Vtis (volume of each compartment) values
params_list <- as.list(init_params) 
#calculate the total volume for urine and feces based on Bellamy et al.1970 rates for Sprague-Dawley rats (table 1)
u_rate <- (0.295+0.468)/2 #(ml/h) urinary rate
f_rate <- (0.133+0.247)/2 #(gr/h) fecal rate equal to ml/h considering density of feces close to 1
V_urine <- u_rate*Urine_data$time #ml total volume of urine @ moments that we have observes
V_urine <- sintegral(Urine_data$time, V_urine) #Integrate AUC of Urine - time
V_urine <- as.numeric(V_urine["value"])
V_feces <- f_rate*Feces_data$time #ml total volume of feces @ moments that we have observes
V_feces <- sintegral(Feces_data$time, V_feces) #Integrate AUC of feces - time
V_feces <- as.numeric(V_feces["value"])

# Vector with the Volume values
Volumes <- c(params_list$Vtis_lu, params_list$Vtis_spl, params_list$Vtis_li, params_list$Vtis_ki, params_list$Vtis_ht, params_list$Vtis_br, params_list$Vtis_ut,
             params_list$Vtis_skel, params_list$Vtis_st , params_list$V_blood, V_feces, V_urine)
weights <- AUC_obs/Volumes # Calculate the scaling factors
#Calculate weighted GI values
GI_weighted <- rep(0,length(init_prior))
for (i in 1:length(init_prior)) {
  GI_weighted[i] <- sum(weights*GI[i,])/sum(weights)
}
xlab_names <- rownames(mu)

df <- data.frame( xlab_names, GI_weighted)
Morris_plot <- ggplot(df, aes(xlab_names, GI_weighted))+ geom_col()+
  labs(title = "Morris Global Sensitivity", y = "value", x = "Parameters" )+
  theme(axis.text.x = element_text(size=15 ,angle=45),
        axis.text.y = element_text(size=15),
        plot.title = element_text(hjust = 0.5,size=26),
        axis.title.y =element_text(hjust = 0.5,size=18,face="bold"),
        axis.title.x =element_text(hjust = 0.5,size=18,face="bold"),
        #axis.text.x=element_text(size=18),
        legend.title=element_text(hjust = 0.5,size=18))

print(Morris_plot)
