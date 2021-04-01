library(deSolve) 
library(ggplot2)
library(rstan)
library(truncnorm)
options(max.print=1e05)

# Set the path containg the results of fit, the "Kreyling IV data.xlsx" and "Rat physiological parameters.xlsx" files
setwd("___")
load(file = "Fit TiO2 IV_results.RData")
stan_results <- extract(fit)

# Set the path where the plots should be saved
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

create.params <- function(comp_names, w, stohastic=TRUE){
  
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
  
  if (stohastic == TRUE){
    theta <- rep(NA, 5)
    theta_std <- theta
    for (i in 1:5) {
      theta[i] <- mean(stan_results$theta[,i])
      theta_std <- sd(stan_results$theta[,i])
    }
    # Transform from Normal to Lognormal distribution
    theta_tr <- log((theta^2)/sqrt((theta^2)+theta_std^2))
    theta_std_tr <- sqrt(log(1 + (theta_std^2)/(theta^2)))
    
    # Sample from the Lognormal distributions
    uptake <- rlnorm(1,theta_tr[1], theta_std_tr[1])
    uptake_spl <- rlnorm(1,theta_tr[2], theta_std_tr[2])
    uptake_li <- rlnorm(1,theta_tr[3], theta_std_tr[3])
    P_rest <- rlnorm(1,theta_tr[4], theta_std_tr[4])
    CLE_f <- rlnorm(1,theta_tr[5], theta_std_tr[5])
  }else{
    uptake <- mean(stan_results$theta[,1])
    uptake_spl <- mean(stan_results$theta[,2])
    uptake_li <- mean(stan_results$theta[,3])
    P_rest <- mean(stan_results$theta[,4])
    CLE_f <- mean(stan_results$theta[,5])
  }
  # Sample from the errors distributions
  e1 <- rtruncnorm(1, a=0, b=Inf, mean = mean(stan_results$sigma[,1]), sd = sd(stan_results$sigma[,1]))
  e2 <- rtruncnorm(1, a=0, b=Inf, mean = mean(stan_results$sigma[,2]), sd = sd(stan_results$sigma[,2]))
  e3 <- rtruncnorm(1, a=0, b=Inf, mean = mean(stan_results$sigma[,3]), sd = sd(stan_results$sigma[,3]))
  xfast = 3
  xslow = 1e-4
  x_br = 1e-5
  P_li = 20
  P_st = 6e-2
  Pup_max = 20
  CLE_u=2.4e-1# 1/h
  uptake_st = 1e-2
  k_de = 4.9e-19
  
  
  
  return(list(
    "Q_total"=Q_total, "V_blood"=Total_Blood, "Vven"=Vven, "Vart"=Vart, "Wm_ven"=Wm_ven, "Wm_art"=Wm_art,
    
    "W_st"=parameters[1,1], "W_ht"=parameters[2,1], "W_ki"=parameters[3,1], "W_br"=parameters[4,1], "W_spl"=parameters[5,1], "W_lu"=parameters[6,1], "W_li"=parameters[7,1], "W_ut"=parameters[8,1], "W_skel"=parameters[9,1],
    
    "Vtis_st"=parameters[1,2], "Vtis_ht"=parameters[2,2], "Vtis_ki"=parameters[3,2], "Vtis_br"=parameters[4,2], "Vtis_spl"=parameters[5,2], "Vtis_lu"=parameters[6,2], "Vtis_li"=parameters[7,2], "Vtis_ut"=parameters[8,2], "Vtis_skel"=parameters[9,2],
    
    "Vcap_st"=parameters[1,3], "Vcap_ht"=parameters[2,3], "Vcap_ki"=parameters[3,3], "Vcap_br"=parameters[4,3], "Vcap_spl"=parameters[5,3], "Vcap_lu"=parameters[6,3], "Vcap_li"=parameters[7,3], "Vcap_ut"=parameters[8,3], "Vcap_skel"=parameters[9,3],
    
    "Wm_st"=parameters[1,5], "Wm_ht"=parameters[2,5], "Wm_ki"=parameters[3,5], "Wm_br"=parameters[4,5], "Wm_spl"=parameters[5,5], "Wm_lu"=parameters[6,5], "Wm_li"=parameters[7,5], "Wm_ut"=parameters[8,5], "Wm_skel"=parameters[9,5],
    
    "Q_st"=parameters[1,4], "Q_ht"=parameters[2,4], "Q_ki"=parameters[3,4], "Q_br"=parameters[4,4], "Q_spl"=parameters[5,4], "Q_lu"=parameters[6,4], "Q_li"=parameters[7,4], "Q_ut"=parameters[8,4], "Q_skel"=parameters[9,4],
    
    "xfast"=xfast, "xslow"=xslow, "x_br"=x_br, "P_li"=P_li, "P_st"=P_st, "P_rest"=P_rest, "Pup_max"=Pup_max,
    "CLE_f"=CLE_f, "CLE_u"=CLE_u, "uptake_st"=uptake_st, "k_de"=k_de,
    
    
    "uptakes"=uptake, "uptake_spl"=uptake_spl, "uptake_li"=uptake_li, "error1" = e1, "error2" = e2, "error3" = e3
     
  ))
}

params<-create.params(compartments,BW)


###############
# ODEs system #
###############
dose <- doses[1]
inits <- c(Mcap_lu=0, Mcap_spl=0, Mcap_li=0, Mcap_ki=0, Mcap_ht=0, Mcap_br=0, Mcap_ut=0, Mcap_skel=0, Mcap_st=0,
           Mtis_lu=0, Mtis_spl=0, Mtis_li=0, Mtis_ki=0, Mtis_ht=0, Mtis_br=0, Mtis_ut=0, Mtis_skel=0, Mtis_st=0,
           Mm_lu=0, Mm_spl=0, Mm_li=0, Mm_ki=0, Mm_ht=0, Mm_br=0, Mm_ut=0, Mm_skel=0, Mm_st=0, 
           M_ven=dose, M_art=0, Mm_ven=0, Mm_art=0, M_feces=0, M_urine=0)

ode.func <- function(time, inits, params){
  with( as.list(c(inits,params)),{
    
    #Permeability coefficients (unitless)
    x_lu = xslow
    x_spl = xfast
    x_li = xfast
    x_ki = xslow
    x_ht = xslow
    x_br = x_br
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


# Produce one solution with mean values
params_determ <- create.params(compartments, BW, stohastic=FALSE)

# Here we create a first solution to get the dimensions because the total instances will be 
# equal to the dimension of sample_time 
sample_time <- c(0, 10/60, 1, 1*24, 7*24, 28*24) #in hours
solution <- ode(times = sample_time, func = ode.func, y = inits, parms = params_determ, method = "bdf")
Total_amounts <- solution[,c(1,35:46)]
colnames(Total_amounts) <- c("Time", "Lungs", "Spleen", "Liver", "Kidneys", "Heart", "Brain", "Uterus", "Skeleton", "Soft tissue",
                             "Blood", "Feces", "Urine")
deterministic.df <- as.data.frame(Total_amounts)



Nsim <- 1000
Nrat <- 25
ltime <- length(sample_time)
data_time <- sample_time
Ncomp <- dim(deterministic.df)[2]
pred_lu <- array(rep(NA, Nrat*ltime*Nsim), dim = c(Nrat,ltime,Nsim))
pred_spl <- array(rep(NA, Nrat*ltime*Nsim), dim = c(Nrat,ltime,Nsim))
pred_li <- array(rep(NA, Nrat*ltime*Nsim), dim = c(Nrat,ltime,Nsim))
pred_ki <- array(rep(NA, Nrat*ltime*Nsim), dim = c(Nrat,ltime,Nsim))
pred_ht <- array(rep(NA, Nrat*ltime*Nsim), dim = c(Nrat,ltime,Nsim))
pred_br <- array(rep(NA, Nrat*ltime*Nsim), dim = c(Nrat,ltime,Nsim))
pred_ut <- array(rep(NA, Nrat*ltime*Nsim), dim = c(Nrat,ltime,Nsim))
pred_skel <- array(rep(NA, Nrat*ltime*Nsim), dim = c(Nrat,ltime,Nsim))
pred_st <- array(rep(NA, Nrat*ltime*Nsim), dim = c(Nrat,ltime,Nsim))
pred_blood <- array(rep(NA, Nrat*ltime*Nsim), dim = c(Nrat,ltime,Nsim))
pred_feces <- array(rep(NA, Nrat*ltime*Nsim), dim = c(Nrat,ltime,Nsim))
pred_urine <- array(rep(NA, Nrat*ltime*Nsim), dim = c(Nrat,ltime,Nsim))

for (sim in 1:Nsim) {
  for (rat in 1:Nrat) {
    print(paste("rat   ",rat, "Nsim    ",sim))
    pred <- matrix(rep(NA, ltime*Ncomp), nrow = ltime)
    # parameter to inform about failed attempts
    failed_attempt <- 0
    # if params result in odes not solving, resample and resolve
    while( (dim(pred)[1] != ltime) | sum(is.na(pred)) ) {
      #params are selected stochastically because create.params samples from stan estimates
      params <- create.params(compartments, BW)
      sol <- ode(times = sample_time, func = ode.func, y = inits, parms = params, method = "bdf")
      pred <- sol[,c(1,35:46)]
      
    }
    for (t in 1:ltime) {
      pred_lu[rat, t, sim] <- rtruncnorm(1, a=0, b=Inf, mean = pred[t,2], sd=params$error2)
      pred_spl[rat, t, sim] <- rtruncnorm(1, a=0, b=Inf, mean = pred[t,3], sd=params$error1)
      pred_li[rat, t, sim] <- rtruncnorm(1, a=0, b=Inf, mean = pred[t,4], sd=params$error1)
      pred_ki[rat, t, sim] <- rtruncnorm(1, a=0, b=Inf, mean = pred[t,5], sd=params$error2)
      pred_ht[rat, t, sim] <- rtruncnorm(1, a=0, b=Inf, mean = pred[t,6], sd=params$error3)
      pred_br[rat, t, sim] <- rtruncnorm(1, a=0, b=Inf, mean = pred[t,7], sd=params$error3)
      pred_ut[rat, t, sim] <- rtruncnorm(1, a=0, b=Inf, mean = pred[t,8], sd=params$error3)
      pred_skel[rat, t, sim] <- rtruncnorm(1, a=0, b=Inf, mean = pred[t,9], sd=params$error1)
      pred_st[rat, t, sim] <- rtruncnorm(1, a=0, b=Inf, mean = pred[t,10], sd=params$error2)
      pred_blood[rat, t, sim] <- rtruncnorm(1, a=0, b=Inf, mean = pred[t,11], sd=params$error1)
      pred_feces[rat, t, sim] <- rtruncnorm(1, a=0, b=Inf, mean = pred[t,12], sd=params$error2)
      pred_urine[rat, t, sim] <- rtruncnorm(1, a=0, b=Inf, mean = pred[t,13], sd=params$error2)
    }
  }
}
bag_of_data <- list(pred_lu,pred_spl,pred_li,pred_ki,pred_ht,pred_br,pred_ut,
                    pred_skel,pred_st,pred_blood,pred_feces,pred_urine)
names(bag_of_data) <- c("Lungs", "Spleen", "Liver", "Kidneys", "Heart", "Brain", "Uterus", "Skeleton", "Soft tissue",
                        "Blood", "Feces", "Urine")
comp_names <- names(bag_of_data)

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
data.df <- Transformed_data[,c(4,2,1,3,5,6,7,9,10,8)]

sd <- subset(sd, select = -(Carcass))
sd.data.df <- sd[,c(4,2,1,3,5,6,7,9,10,8)]
feces_exp[,1] <- feces_exp[,1]*24 #transform time to hours
urine_exp[,1] <- urine_exp[,1]*24 #transform time to hours
feces_exp[,2] <- feces_exp[,2]*doses[1]/100
urine_exp <- subset(urine_exp, select = -(Urine_excretion_rate))
urine_exp[,2] <- urine_exp[,2]*doses[1]/100
colnames(urine_exp) <- c("Time", "Urine")

feces.df <- feces_exp
urine.df <- urine_exp
Time <- sample_time[2:length(sample_time)]
data.df <- cbind(Time,data.df)
sd.data.df <- cbind(Time,sd.data.df)

counter <-1
for (dat in bag_of_data) {
  
  # Get the compartment name to be plotted
  comp_name <- comp_names[counter]
  # Define the save name
  save_name <- paste0(comp_name, "_vpc.png", sep="")
  pred_con <- dat
  
  
  
  
  ##########################
  # Calculation of 95 % CI of the 50th quantile
  con<-matrix(rep(NA,ltime*Nsim), nrow = ltime)
  ymin<- rep(NA,ltime)
  ymax<- rep(NA,ltime)
  
  for (i in 1:ltime){
    for (j in 1:Nsim){
      con[i,j]<-quantile(pred_con[,i,j],probs=0.5)
    }
  }
  
  for (i in 1:ltime){
    con[i,]<-sort(con[i,])
    ymin[i]<-quantile(con[i,],probs=0.025)
    ymax[i]<-quantile(con[i,],probs=0.975)
  }
  
  q50<-data.frame(time=data_time,lo=ymin,hi=ymax)
  ############################  
  # Calculation of 95 % CI of the 5th quantile
  
  con<-matrix(rep(NA,ltime*Nsim), nrow = ltime)
  ymin<- rep(NA,ltime)
  ymax<- rep(NA,ltime)
  
  for (i in 1:ltime){
    for (j in 1:Nsim){
      con[i,j]<-quantile(pred_con[,i,j],probs=0.05)
    }
  }
  
  for (i in 1:ltime){
    con[i,]<-sort(con[i,])
    ymin[i]<-quantile(con[i,],probs=0.025)
    ymax[i]<-quantile(con[i,],probs=0.975)
  }
  
  q05<-data.frame(time=data_time,lo=ymin,hi=ymax)
  #################################  
  # Calculation of 95 % CI of the 95th quantile
  con<-matrix(rep(NA,ltime*Nsim), nrow = ltime)
  ymin<- rep(NA,ltime)
  ymax<- rep(NA,ltime)
  
  for (i in 1:ltime){
    for (j in 1:Nsim){
      con[i,j]<-quantile(pred_con[,i,j],probs=0.95)
    }
  }
  
  for (i in 1:ltime){
    con[i,]<-sort(con[i,])
    ymin[i]<-quantile(con[i,],probs=0.025)
    ymax[i]<-quantile(con[i,],probs=0.975)
  }
  
  q95<-data.frame(time=data_time,lo=ymin,hi=ymax)
  ############################
  # Mean values of predictions from deterministic.df model
  use_comp <- colnames(deterministic.df)[counter+1]
  mean_pred = data.frame(deterministic.df[,c("Time", use_comp)])
  colnames(mean_pred) = c("time", "value")
  
  if(!(comp_name %in% c("Feces", "Urine"))){
    use_comp <- colnames(data.df)[counter+1]
    df1 <- as.data.frame(cbind(data.df[,c("Time",use_comp)], sd.data.df[,use_comp]))
    colnames(df1) <- c("Time", "Mean", "Sd")
    my_plot <- ggplot() +   
      geom_line(data=mean_pred, aes(x=time, y=value, linetype = " Prediction mean"),
                size=1.2)+
      geom_ribbon(data=q50,aes(x=time, ymin = lo, ymax = hi, fill = "50th percentile"),
                  inherit.aes=FALSE, alpha = 0.2) + 
      geom_ribbon(data=q05,aes(x=time, ymin = lo, ymax = hi, fill = "5th percentile"),
                  inherit.aes=FALSE, alpha = 0.2) + 
      geom_ribbon(data=q95,aes(x=time, ymin = lo, ymax = hi, fill="95th percentile"),
                  inherit.aes=FALSE,alpha = 0.2) + 
      geom_point(data = df1, aes(x = Time, y=Mean),size=5)+
      geom_errorbar(data = df1, aes(x = Time, ymin=ifelse((Mean-Sd)>0,Mean-Sd,0), ymax=Mean+Sd),size=1)+
      labs(title = rlang::expr(!!comp_name), y = "TiO2 (ug)", x = "Time (in hours)") +
      scale_fill_manual("Prediction Ribbons", values = c( "50th percentile" = 1, 
                                                          "5th percentile" = 2,"95th percentile" = 3)) + 
      scale_linetype_manual("Line", values = c(" Prediction mean" = 1))+
      scale_colour_manual("Point", values = c("Biodistribtion data" = 1))+
      theme(plot.title =element_text(hjust = 0.5, size=30, face="bold"),
            axis.title.y =element_text(hjust = 0.5, size=20, face="bold"),
            axis.text.y=element_text(size=18),
            axis.title.x =element_text(hjust = 0.5, size=20, face="bold"),
            axis.text.x=element_text(size=18),
            legend.title=element_text(hjust = 0.01, size=20), 
            legend.text=element_text(size=18))
    png(rlang::expr(!!save_name), width = 15, height = 10, units = 'in', res= 100)
    print(my_plot)
  } else if(comp_name=="Feces"){
    observed <- feces.df[,c("Time",colnames(feces.df)[2])]
    colnames(observed) <- c("Time", "mean")
    my_plot <- ggplot(observed, aes(x=Time, y=mean, colour="Biodistribtion data"))+
      geom_point(shape=19, size=4) +   
      #geom_bar(aes(ymin=ifelse((mean-sd)>0,mean-sd,0), ymax=mean+sd), width=.1) +
      geom_line(data=mean_pred, aes(x=time, y=value, linetype = " Prediction mean"), 
                size=1.2)+
      geom_ribbon(data=q50,aes(x=time, ymin = lo, ymax = hi, fill = "50th percentile"),
                  inherit.aes=FALSE, alpha = 0.2) + 
      geom_ribbon(data=q05,aes(x=time, ymin = lo, ymax = hi, fill = "5th percentile"),
                  inherit.aes=FALSE, alpha = 0.2) + 
      geom_ribbon(data=q95,aes(x=time, ymin = lo, ymax = hi, fill="95th percentile"),
                  inherit.aes=FALSE,alpha = 0.2) + 
      labs(title = rlang::expr(!!comp_name), y = "TiO2 (ug)", x = "Time (in hours)") +
      scale_fill_manual("Prediction Ribbons", values = c( "50th percentile" = 1, 
                                                          "5th percentile" = 2,"95th percentile" = 3)) +
      scale_linetype_manual("Line", values = c(" Prediction mean" = 1))+
      scale_colour_manual("Point", values = c("Biodistribtion data" = 1))+
      theme(plot.title =element_text(hjust = 0.5, size=30, face="bold"),
            axis.title.y =element_text(hjust = 0.5, size=20, face="bold"),
            axis.text.y=element_text(size=18),
            axis.title.x =element_text(hjust = 0.5, size=20, face="bold"),
            axis.text.x=element_text(size=18),
            legend.title=element_text(hjust = 0.01, size=20), 
            legend.text=element_text(size=18))
    png(rlang::expr(!!save_name), width = 15, height = 10, units = 'in', res= 100)
    print(my_plot)
    
  }else{
    observed <- urine.df[,c("Time",colnames(urine.df)[2])]
    colnames(observed) <- c("Time", "mean")
    my_plot <- ggplot(observed, aes(x=Time, y=mean, colour="Biodistribtion data"))+
      geom_point(shape=19, size=4) +   
      #geom_bar(aes(ymin=ifelse((mean-sd)>0,mean-sd,0), ymax=mean+sd), width=.1) +
      geom_line(data=mean_pred, aes(x=time, y=value, linetype = " Prediction mean"), 
                size=1.2)+
      geom_ribbon(data=q50,aes(x=time, ymin = lo, ymax = hi, fill = "50th percentile"),
                  inherit.aes=FALSE, alpha = 0.2) + 
      geom_ribbon(data=q05,aes(x=time, ymin = lo, ymax = hi, fill = "5th percentile"),
                  inherit.aes=FALSE, alpha = 0.2) + 
      geom_ribbon(data=q95,aes(x=time, ymin = lo, ymax = hi, fill="95th percentile"),
                  inherit.aes=FALSE,alpha = 0.2) + 
      labs(title = rlang::expr(!!comp_name), y = "TiO2 (ug)", x = "Time (in hours)") +
      scale_fill_manual("Prediction Ribbons", values = c( "50th percentile" = 1, 
                                                          "5th percentile" = 2,"95th percentile" = 3)) +
      scale_linetype_manual("Line", values = c(" Prediction mean" = 1))+
      scale_colour_manual("Point", values = c("Biodistribtion data" = 1))+
      theme(plot.title =element_text(hjust = 0.5, size=30, face="bold"),
            axis.title.y =element_text(hjust = 0.5, size=20, face="bold"),
            axis.text.y=element_text(size=18),
            axis.title.x =element_text(hjust = 0.5, size=20, face="bold"),
            axis.text.x=element_text(size=18),
            legend.title=element_text(hjust = 0.01, size=20), 
            legend.text=element_text(size=18))
    png(rlang::expr(!!save_name), width = 15, height = 10, units = 'in', res= 100)
    print(my_plot)
  }
  dev.off()
  counter <- counter +1
}
