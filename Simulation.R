################################################################################################
#                                       Simulation study                                       #
#                                           Data                                               #
################################################################################################
#Last file: NA
#Author: Alma Revers
#Libraries needed
library(MASS)

data_sim_total <- list()
input_sim_total <- list()

#parameters to set
N_simulations <- 100
NB_dist <- FALSE

N_subj = 50 #Klein
#N_subj = 250 #Middel
#N_subj = 500 #Groot

N_otu = 10 #Klein
#N_otu = 100 #Middel
#N_otu = 200 #Groot

N_family <- 2
#total_count_min = 10 #300 #600
#total_count_max = 50 #500 #1000


for ( s in 1:N_simulations){
    index_family <- round(runif(n=N_otu, min=1, max=N_family))
    if ( s < 51 ){
      x = rnorm(n= N_subj, mean=0, sd=1)
    }
    else{
      x = scale(rnorm(n = N_subj, mean=0, sd=1)^2)
    }
   
    if(round(runif(n=1, min=0, max=1)) == 0){
      mean_slope_fam_1 = 0
      mean_slope_fam_2 = 1
      mean_intercept_fam_1 = 1
      mean_intercept_fam_2 = 1
    }else {
      mean_slope_fam_1 = 0
      mean_slope_fam_2 = -1
      mean_intercept_fam_1 = 1
      mean_intercept_fam_2 = 1
    }
    
    slope = c()
    intercept = c()
    fam_slope <- rnorm(n=2, mean=c(mean_slope_fam_1,mean_slope_fam_2) , sd=0.1)
    fam_intercept <- rnorm(n=2, mean=c(mean_intercept_fam_1,mean_intercept_fam_2) , sd=0.1) 
    for ( j in 1:N_otu){
      Sigma <- matrix(c(0.1^2, 0.1*0.1*(-0.7),  0.1*0.1*(-0.7), 0.1^2), 2,2)
      intercept_slope_mvn <- mvrnorm(n=1, mu= c(fam_intercept[index_family[j]], fam_slope[index_family[j]]), Sigma=Sigma)
      intercept <- c(intercept, intercept_slope_mvn[1])
      slope <- c(slope, intercept_slope_mvn[2])
    }
    dispersion <- rlnorm(n = N_otu, meanlog = exp(0.01/intercept), sdlog=0.1)
    
    y <- c()
    if(NB_dist) {
      y_subj <- c()
      for(i in 1:N_subj){
        y_subj <- c()
        y_subj <- c(y_subj, rnegbin(n = N_otu, mu = exp(intercept + slope*x[i]), theta = dispersion))
        y <- cbind(y, y_subj)
        }
    }
    
    if(!NB_dist){
      for(i in 1:N_subj){
        y_subj <- c()
        y_subj <- c(y_subj, rpois(n = N_otu, lambda = exp(intercept + slope*x[i])))
        y <- cbind(y, y_subj)
      }
      for(j in 1:N_otu){
        extreme_max <- order(y[j,])[(round(N_subj*0.8)):N_subj]
        y[j,extreme_max] <- round(mean(y[j,]))
      }
    }
    
    data_sim <- list(
      N_subj = N_subj,
      N_otu = N_otu,
      N_family = N_family,
      N_variables = 2,
      y = t(y),
      x= cbind(rep(1, N_subj), x),
      offset= rep(0, N_subj),
      index_family = index_family
    )
    data_sim_total[[s]] <- data_sim
    
    input_sim_s <- list(
      index_family = index_family,
      slope = slope,
      intercept = intercept,
      y = t(y),
      x = x
    )
    input_sim_total[[s]] <- input_sim_s
}
################################################################################################
#                               Save output for next step                                      #
################################################################################################
saveRDS(data_sim_total, file="./Simulation_data/True_data/data_NB_K_K.rds")
saveRDS(input_sim_total, file="./Simulation_data/True_data/input_NB_K_K.rds")


