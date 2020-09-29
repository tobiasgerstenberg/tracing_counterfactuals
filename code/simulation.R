library("TreeBUGS")
source("simfun.R")
library("tidybayes")    # for extracting parameters from Bayesian model
library("tidyverse") 

# load("cache/fit_latent-trait_all.rda")
# summary(fit_bayes_all)

# str(fit_bayes_all$runjags, 1)
# #summary(fit_bayes_all$runjags)
# 
# dimnames(fit_bayes_all$runjags$mcmc[[1]])[[2]]
# 
# str(fit_bayes_all$runjags$mcmc)

mus <- c(m_o = 2.0, m_t = 0.1, s = 0)
sigmas <- c(m_o = 1.2, m_t = 2, s = 1.5)
rho <- matrix(c(   1, -0.2,  0.2,
                -0.2,    1,  0.3,
                 0.2,   .3,    1), nrow=3)
colnames(rho) <- rownames(rho) <- c("m_o", "m_t", "s")


### check a few visually
sim1 <- one_sim(mu = mus, sigma = sigmas, rho = rho, 
                diff_s_prob = 0.2, show_diag = TRUE)
(sim2 <- one_sim(mu = mus, sigma = sigmas, rho = rho, 
                diff_s_prob = 0.6, show_diag = TRUE))

## run on server:
NSIM <- 500
diffs <- seq(from = 0.1, to = 0.9, by = 0.1)

# NSIM <- 2
# diffs <- c(0.1, 0.2)


all_diffs <- rep(diffs, each = NSIM)

library(parallel)

### needs to be potentially run in individual chunks.
### because fitting uses 4 cores, setting mc.cores = 6 requires a 24 core machine
allfit <- mclapply(all_diffs, one_sim, 
                   mu = mus, sigma = sigmas, rho = rho, mc.cores = 6L)
save(allfit, file = "cache/sim_results.rda")
