################################################################################################
#                                           Tutorial                                           #
#                                                                                              #
################################################################################################
#Author: Alma Revers
#Libraries needed
library(rstan) 
library(phyloseq)
library(broom.mixed)
library(ggplot2)

#Settings for parallel computing
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#load data
#setwd("")
data <- readRDS("example_data.rds")

##########################Some basics about the data########################################

#The OTU counts
head(data@otu_table)
dim(data@otu_table)

#The sample data
dim(data@sam_data)
hist(data@sam_data[,1])

OTU_to_show <- 1
x_to_show <- 1
plot(unlist(data@sam_data[,x_to_show]), unlist(data@otu_table[,OTU_to_show]),
     xlab="X",
     ylab="OTU count")

#The phylogenetic structure
dim(data@tax_table)
head(data@tax_table)
data@tax_table[,5] #The family per OTU
table(data@tax_table[,5])

###########################Pre-process data##################################################
threshold_abundance <- 0.1
threshold_structure <- 5

N_samples_per_OTU<- apply(data@otu_table, 1, function(x) sum(x>0))
OTU_above_threshold <- N_samples_per_OTU > (length(N_samples_per_OTU)*threshold_abundance)
table(OTU_above_threshold)

families_to_include <- table(data@tax_table[,5]) > threshold_structure
OTUs_to_include_based_on_family <- data@tax_table[,5] %in% names(families_to_include)[families_to_include]
table(OTUs_to_include_based_on_family)

#Format data to expected input for stan
x <- 1
stan_data <- list(
  N_subj = length(unlist(data@sam_data[,x])),
  N_otu = dim(data@otu_table[OTU_above_threshold&OTUs_to_include_based_on_family,])[1],
  N_family = sum(families_to_include),
  N_variables = 2,
  y = t(data@otu_table[OTU_above_threshold&OTUs_to_include_based_on_family,]),
  x = cbind(1, unlist(data@sam_data[,x])),
  offset= log(apply(data@otu_table[OTU_above_threshold&OTUs_to_include_based_on_family,], 2, sum)),
  index_family= as.numeric(as.factor(data@tax_table[OTU_above_threshold&OTUs_to_include_based_on_family,5]))
)

############################Run stan for estimating prosirior#################################### 
model_obj <- stan_model(file="NB.stan")
fit <- sampling(model_obj, data=stan_data)

#####################################Check convergency###########################################
#traceplots of the hyperparameters
traceplot(fit, pars="mu_beta")
traceplot(fit, pars="sigma_beta")
traceplot(fit, pars="a1")
traceplot(fit, pars="alpha_0")

beta <- tidy(fit, pars="beta_otu", 
             estimate.method = "mean", 
             conf.int=TRUE, 
             conf.level = 0.95, 
             rhat=TRUE, 
             ess=TRUE)
#Rhat and effective sample size of slopes
hist(beta$rhat[seq(2, length(beta$rhat), by=2)])
hist(beta$ess[seq(2, length(beta$rhat), by=2)])

#####################################CPlot results##############################################
volcano_plot_data <- data.frame(x = beta$estimate[seq(2, length(beta$estimate), by=2)],
                                y = 1/(beta$std.error[seq(2, length(beta$estimate), by=2)]^2),
                              family = as.factor(data@tax_table[OTU_above_threshold&OTUs_to_include_based_on_family,5])
)
ggplot(volcano_plot_data, aes(x = x, y=abs(y), color=family)) +
  geom_point() +
  theme_bw()+
  geom_vline(xintercept = 0, linetype = "dashed", alpha=0.4) +
  labs(x = "Estimated mean association", y = expression("1/"~SD^2), color ="Phylogenetic Family")






