library(dplyr)
library(tidyr)
library(survival)
library(stats)
library(doParallel)
library(parallel)

N = 10e5 # sample size 

tau_end = 100; # maximum possible follow-up time for each subject, which we specify is the same across subjects.

tau = 5; # end time for the SPCE

B = 100; # number of bootstrap replicates desired 

###############################################
# Covariate distribution Parmeters 
###############################################

# U1 ~ N(muU1, sdU1)
muU1 = 0; sdU1 = 1 # U1 ~ N(muU1, sdU1) Individual-level covariate

# U2 ~ Bernoulli(probU2)
probU2 = 0.5 # U2 ~ Bernoulli(probU2) Individual-level covariate

########################################
# Model Coefficients 
########################################

# a) treatment assignment probability model 
probA_coeffs = c(-0.5, log(1.25), log(0.8))

# b) hazard model coefficients for potential survival times (primary outcome)
outcome_coeffs_0 = c(log(1.5), log(1.5)) # control group
outcome_coeffs_1 = c(log(1.25), log(1.25))  # treatment group

# c) baseline hazard functions for potential survival times (primary outcome)
c0_outcome = 0.3 # control group (constant baseline hazard)
c1_outcome = 0.25 # treatment group (constant baseline hazard)

# d) hazard model coefficients for potential censoring times 
censor_coeffs_0 = c(log(1.25), log(1.25)) # control group
censor_coeffs_1 = c(log(1.25), log(1.5))  # treatment group

# e) baseline hazard functions for the potential censoring times 
c0_censor = 0.05  # control group (constant baseline hazard)
c1_censor = 0.05 # treatment group (constant baseline hazard)

########################################
# Calculate the true SPCE 
########################################

f0 = function(outcome_coeffs_0, c0_outcome, tau, x) {

# u2 = 0,1 
# tau = end time for the SPCE 

y = 

exp(-(c0_outcome * tau) * exp(outcome_coeffs_0[1]*x[1])) * 
dnorm(x[1], mean = 0, sd=1) * (1-probU2) + 

exp(-(c0_outcome * tau) * exp(outcome_coeffs_0[1]*x[1] + outcome_coeffs_0[2])) * 
dnorm(x[1], mean = 0, sd=1) * (probU2)

return(y)

}

f1 = function(outcome_coeffs_1, c1_outcome, tau, x) {

# u2 = 0,1 
# tau = end time for the SPCE 

y = 

exp(-(c1_outcome * tau) * exp(outcome_coeffs_1[1]*x[1])) * 
dnorm(x[1], mean = 0, sd=1) * (1-probU2) + 

exp(-(c1_outcome * tau) * exp(outcome_coeffs_1[1]*x[1] + outcome_coeffs_1[2])) * 
dnorm(x[1], mean = 0, sd=1) * (probU2)

return(y)

}

# Numerical Integration 
EY0 = pcubature(f0, c(-Inf), c(Inf), outcome_coeffs_0 = outcome_coeffs_0, c0_outcome = c0_outcome, tau = tau, tol = 1e-05)$integral
EY1 = pcubature(f1, c(-Inf), c(Inf), outcome_coeffs_1 = outcome_coeffs_1, c1_outcome = c1_outcome, tau = tau, tol = 1e-05)$integral

SPCE_true = EY1 - EY0

print(SPCE_true)

###########################################
# Data + Potential Outcome Generation
###########################################

seed = 111; 

set.seed(seed); # set seed for reproducibility 

########################################################	
# Generate covariates (U1 and U2)
########################################################
U1 <- rnorm(N, muU1, sdU1)
U2 <- rbinom(N, 1, prob = probU2) # Scenario where U2 is independent of V and U1
    
########################################################	
# Generate Treatment Assignments
########################################################

# Generate Treatment Assignment A:
prob_A <- plogis(probA_coeffs[1] + U1*probA_coeffs[2] + U1*probA_coeffs[3])
A <- rbinom(N, 1, prob = prob_A)

########################################################	
# Generate potential survival times
########################################################

# Treatment group      
T1 <- (1/c1_outcome) * (-log(1-runif(N)) / exp(U1*outcome_coeffs_1[1] + U2*outcome_coeffs_1[2]))

# Control group 
T0 <- (1/c0_outcome) * (-log(1-runif(N)) / exp(U1*outcome_coeffs_0[1] + U2*outcome_coeffs_0[2]))

# "Observed" Event Time (even though it may not actually be observed)
T <- T1 * A + T0 * (1-A)

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

########################################################
# Generate potential censoring times 
########################################################

# Treatment Group
C1 <- (1/c1_censor) * (-log(1-runif(N)) / exp(U1*censor_coeffs_1[1] + U2*censor_coeffs_1[2]))

# Control Group
C0 <- (1/c0_censor) * (-log(1-runif(N)) / exp(U1*censor_coeffs_0[1] + U2*censor_coeffs_0[2]))

# "Observed" Censoring Time (even though it may not actually be observed)
C <- C1 * A + C0 * (1-A)

########################################################
# Create the data frame for analysis 
########################################################

# Create indicator variables 
index_delta_1 = which(T <= C & T <= tau_end) 
index_gamma_1 = which(C < T & C < tau_end) 
index_EOFU = which(tau_end < T & tau_end < C) # EOFU = end of follow-up 

delta = rep(0, N); delta[index_delta_1] = 1; # equals 1 if the subject experienced the primary outcome before censoring occurred or before time tau_end, 0 otherwise 
gamma = rep(0, N); gamma[index_gamma_1] = 1; # equals 1 if the subject was censored before the primary outcome occured or before time tau_end, 0 otherwise
zeta = rep(0, N); zeta[index_EOFU] = 1; # equals 1 if end of follow-up (EOFU) is reached before censoring or primary outcome occurence 

index_censored = which(gamma == 1 | zeta == 1)
censored = rep(0, N); censored[index_censored] = 1
event = 1-censored;

# This creates the data frame for analysis 
data = 
data.frame(
survtime = pmin(T, C, tau_end), # This is M 
A = A, # Treatment 
U1 = U1, # Covariate 1
U2 = U2, # Covariate 2
delta = delta, # 1 = event occured, 0 = o.w.
gamma = gamma, # 1 = censored from LTFU, 0 = o.w.
zeta = zeta,   # 1 = end of follow-up reached, 0 o.w.
censored = censored # 1 = censored (including administrative censoring)
)

###############################################
###############################################
# Perform Analysis to Estimate SPCE 
###############################################
###############################################

# Create function to perform G-computation and IPTCW estimation for each bootstrap replicate of the data 

# Begin function
func = function(b, data) {

set.seed(b);

index = sample(1:N, N, replace = TRUE)

data_b = data[index,];
rm("index"); gc();

######################### Estimate the Propensity scores #########################

fit_ps = glm(A ~ U1 + U2, family = binomial, data = data_b);

data_b$PS <- stats::predict.glm(fit_ps, type = "response") %>% unname()

rm("fit_ps"); gc();

# Separate data by treatment status 
data0 = 
data_b %>% 
filter(A == 0);

data0$wa = 1/(1-data0$PS)
data0$id = 1:nrow(data0)

data1 = 
data_b %>% 
filter(A == 1);

data1$wa = 1/(data1$PS)
data1$id = 1:nrow(data1)

######################### Estimate the Censoring weights #########################

fit0_censoring = coxph(Surv(survtime, censored) ~ U1 + U2, data = data0) # fit the cox model for censoring for A=0
fit1_censoring = coxph(Surv(survtime, censored) ~ U1 + U2, data = data1) # fit the cox model for censoring for A=1

newdata0 = data0; newdata0$survtime = tau; 
 
newdata0$Sprob_tau_cens = predict(fit0_censoring, newdata = newdata0,
                 type = "survival", se.fit = FALSE) # obtain the inverse-censoring weights for A=0

newdata1 = data1; newdata1$survtime = tau;

newdata1$Sprob_tau_cens <- predict(fit1_censoring, newdata = newdata1,
                 type = "survival", se.fit = FALSE) # obtain the inverse-censoring weights for A=1
 
newdata0$wc <- 1/newdata0$Sprob_tau_cens # obtain the censoring weights for A=0
newdata1$wc <- 1/newdata1$Sprob_tau_cens # obtain the censoring weights for A=1

######################### Estimate the SPCE using IPTCW #################################

qq = which(data0$survtime >= tau)
PT0_wt = sum(data0$wa[qq] * newdata0$wc[qq]) / nrow(data_b)

qq = which(data1$survtime >= tau)
PT1_wt = sum(data1$wa[qq] * newdata1$wc[qq]) / nrow(data_b)

SPCE_wt = PT1_wt - PT0_wt # This is the causal estimand estimate using IPTCW

######################### Estimate the SPCE using G-computation #########################

fit0_outcome = coxph(Surv(survtime, delta) ~ U1 + U2, data = data0) # fit the cox model for the outcome for A=0
fit1_outcome = coxph(Surv(survtime, delta) ~ U1 + U2, data = data1) # fit the cox model for the outcome for A=1

newdata0$Sprob_tau_outcome = predict(fit0_outcome, newdata = newdata0,
                 type = "survival", se.fit = FALSE) # obtain the predicted Survival probability at time tau under A=0

newdata1$Sprob_tau_outcome <- predict(fit1_outcome, newdata = newdata1,
                 type = "survival", se.fit = FALSE) # obtain the predicted Survival probability at time tau under A=1

PT0_gcomp = mean(newdata0$Sprob_tau_outcome);
PT1_gcomp = mean(newdata1$Sprob_tau_outcome);

SPCE_gcomp = PT1_gcomp - PT0_gcomp

# Clean-up
rm("qq", "fit0_censoring", "fit1_censoring", "fit0_outcome", "fit1_outcome", "newdata0", "newdata1", "data0", "data1");
gc();

# Export results 
return(c(PT1_wt, PT0_wt, SPCE_wt, PT1_gcomp, PT0_gcomp, SPCE_gcomp))

}
# End function

# Perform the bootstrap using do-parallel

n.cores <- parallel::detectCores() - 1

my.cluster = 
parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
  )

packlist = c("dplyr", "tidyr", "survival", "stats")

doParallel::registerDoParallel(cl = my.cluster)

Results <- foreach(
  b = 1:B, .combine = 'rbind',
  .packages = packlist
) %dopar% {

func(b, data)

}

parallel::stopCluster(cl = my.cluster)

LL_SPCE_wt = quantile(Results[,3], 0.025);
UL_SPCE_wt = quantile(Results[,3], 0.975);
SPCE_wt = mean(Results[,3]);

LL_SPCE_gcomp = quantile(Results[,6], 0.025);
UL_SPCE_gcomp = quantile(Results[,6], 0.975);
SPCE_gcomp = mean(Results[,6]);

data.frame(
Method = c("G-computation", "IPTCW"),
LL = c(LL_SPCE_wt, LL_SPCE_gcomp),
UL = c(UL_SPCE_wt, UL_SPCE_gcomp)
)

         # Method         LL         UL
# 1 G-computation 0.06779914 0.07132417
# 2         IPTCW 0.06794244 0.07140384
