library(ggplot2)
library(rstan)

options(max.print=1e5)

# Set the path containg the results of fit
setwd("___")
load(file = "Fit TiO2 IV_results.RData")
stan_results <- extract(fit)

# Posterior distributions
posterior <- cbind(uptake =(stan_results$theta[,1]),
                   uptake_spl = (stan_results$theta[,2]),
                   uptake_li = (stan_results$theta[,3]),
                   P_rest = (stan_results$theta[,4]),
                   CLE_f = (stan_results$theta[,5])
)

# Prior distributions
eta_tr <- c( uptake = 1e-1, uptake_spl = 3, uptake_li = 16, P_rest = 0.8,
                   CLE_f  = 6e-4)


eta_tr_std <- c( uptake = 2, uptake_spl = 0.5, uptake_li = 0.5, P_rest = 1,
                       CLE_f  = 0.5)

prior_pars_mu <- rep(NA,length(eta_tr))
prior_pars_std <- rep(NA,length(eta_tr_std))
prior <- matrix(rep(NA, 5*dim(posterior)[1]), ncol = 5)
N_param = 5
for (i in 1:N_param){
  prior_pars_std[i] <- sqrt(log(((eta_tr_std[i]^2)/(eta_tr[i])^2)+1));
  prior_pars_mu[i] <- log(((eta_tr[i])^2)/sqrt((eta_tr_std[i]^2)+(eta_tr[i])^2));
  prior[,i] <- rlnorm(dim(posterior)[1], prior_pars_mu[i],  prior_pars_std[i])
}


theta1_wide <- data.frame(prior =prior[,1], posterior = posterior[,1])
theta2_wide <- data.frame(prior =prior[,2], posterior = posterior[,2])
theta3_wide <- data.frame(prior =prior[,3], posterior = posterior[,3])
theta4_wide <- data.frame(prior =prior[,4], posterior = posterior[,4])
theta5_wide <- data.frame(prior =prior[,5], posterior = posterior[,5])

theta1 <- tidyr::gather(theta1_wide, type, value, prior:posterior, factor_key=TRUE)
theta2 <- tidyr::gather(theta2_wide, type, value, prior:posterior, factor_key=TRUE)
theta3 <- tidyr::gather(theta3_wide, type, value, prior:posterior, factor_key=TRUE)
theta4 <- tidyr::gather(theta4_wide, type, value, prior:posterior, factor_key=TRUE)
theta5 <- tidyr::gather(theta5_wide, type, value, prior:posterior, factor_key=TRUE)

bag_of_data <- list(theta1, theta2, theta3, theta4, theta5)#, Urine)
names(bag_of_data) <- c("Generic Uptake", "Uptake Spleen", "Uptake Liver",
                        "Generic Partition Coefficient","Hepatobiliary Clearance")
x_axis <- c("ug of TiO2 per g of phagocytes","ug of TiO2 per g of phagocytes",
            "ug of TiO2 per g of phagocytes", "uniteless", "1/h")
comp_names <- names(bag_of_data)
counter <-1

for (dat in bag_of_data) {
  comp_name<-comp_names[counter]
  x_axis_title <- x_axis[counter]
  save_name<-paste0(comp_name, ".png", sep="")
  data_to_plot <- dat
  
 
  my_plot <-ggplot(data_to_plot, aes(x=value, fill=type)) +
    geom_density()+
   
    labs(title = rlang::expr(!!comp_name),  y = "Density", x = rlang::expr(!!x_axis_title)) +
    
    
    theme(plot.title = element_text(hjust = 0.5,size=30), 
          axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
          axis.text.y=element_text(size=22),
          axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
          axis.text.x=element_text(size=22),
          legend.title=element_text(hjust = 0.5,size=25), 
          legend.text=element_text(size=22))
  png(rlang::expr(!!save_name),width = 15, height = 10, units = 'in', res = 100)
  print(my_plot)
  
  dev.off()
  counter <- counter+1
}






