
one_sim <- function(
  diff_s_prob,
  mu, sigma, rho, 
  N = 24, numItems = c("items" = 6), 
  eqnfile = "model_sim.eqn", 
  restrictions = list("g1=g2=0.5"),
  n.iter = 50000, 
  n.thin = 5, 
  n.chains = 4,
  n.adapt = 5000, 
  n.burnin = 5000,
  show_diag = FALSE
) {
  
  mu1 <- mu2 <- mu
  mu1[["s"]] <- qnorm(pnorm(mu[["s"]]) - diff_s_prob/2)
  mu2[["s"]] <- qnorm(pnorm(mu[["s"]]) + diff_s_prob/2)
  
  dat1 <- genTraitMPT(
    N = N, 
    numItems = numItems, 
    eqnfile = eqnfile, 
    restrictions = restrictions, 
    mu =  mu1, 
    sigma = sigmas, 
    rho = rho
  )
  dat2 <- genTraitMPT(
    N = N, 
    numItems = numItems, 
    eqnfile = eqnfile, 
    restrictions = restrictions, 
    mu =  mu2, 
    sigma = sigmas, 
    rho = rho
  )
  #browser()
  
  dfit <- rbind(dat1$data, dat2$data)
  dc <- data.frame(age = factor(rep(c("4", "5"), each = N)))
  
  fit_both <- traitMPT(eqnfile, dfit, 
                  restrictions = restrictions, 
                  covData = dc,
                  predStructure = list("m_o m_t s ; age"), 
                  predType = "f",
                  n.iter = n.iter, n.thin = n.thin, n.chains = n.chains,
                  n.adapt = n.adapt, n.burnin = n.burnin)
  
  #show_diag <- TRUE
  if (show_diag) {
    par(ask=TRUE)
    plot(fit_both)
    plot(fit_both, "sigma")
    plot(fit_both, "rho")
    par(ask=FALSE)
    print(summary(fit_both))
  }
  
  # fit1$runjags
  # dimnames(fit_both$runjags$mcmc[[1]])[[2]]
  
  means <- gather_draws(fit_both$runjags, mu[i]) %>% 
    mutate(
      parameter = factor(i, levels = c(3, 1, 2), labels = c("s", "m_o", "m_t"))
    ) %>% 
    rename(mean = .value) %>% 
    ungroup() %>% 
    select(-i, -.variable)
  means %>% 
    filter(parameter == "s")
  
  means %>% 
    group_by(parameter) %>% 
    mean_qi()
  
  pars = gather_draws(fit_both$runjags, 
                     `factor_.*_age`[i], regex = TRUE) %>% 
  ungroup() %>% 
  mutate(age_group = factor(letters[i], level = letters[1:2]),
         #estimate = pnorm(.value),
         parameter = factor(
           str_replace(str_replace(.variable, "factor_", ""), "_age", ""), 
           levels = c("s", "m_o", "m_t"))) 
  pars <- left_join(pars, means) %>% 
    mutate(.value = pnorm(.value + mean))
  
  pars %>% 
    group_by(parameter, age_group) %>% 
    mean_qi(.value)
  
  post <- pars %>% 
    filter(parameter == "s") %>% 
    pivot_wider(id_cols = .draw, names_from = "age_group", values_from = ".value") %>% 
    mutate(diff = b - a)
  
  ci80t <- post %>% 
    mean_qi(diff, .width = .8) %>% 
    mutate(incl_0 = sign(.lower) == sign(.upper))
  ci95t <- post %>% 
    mean_qi(diff, .width = .95) %>% 
    mutate(incl_0 = sign(.lower) == sign(.upper))
  
  tt1 <- as_tibble(as.list(mu))
  colnames(tt1) <- paste0("mu.", colnames(tt1))
  tt2 <- as_tibble(as.list(sigma))
  colnames(tt2) <- paste0("sigma.", colnames(tt2))

  out <- tibble(
    sim_diff = diff_s_prob,
    obs_diff = ci80t$diff,
    ci80 = ci80t$incl_0,
    ci95 = ci95t$incl_0,
    width_80 = ci80t$.upper - ci80t$.lower,
    width_95 = ci95t$.upper - ci95t$.lower,
    sa = mean_qi(post, a)$a - pnorm(mu1[["s"]]),
    sb = pnorm(mu2[["s"]]) - mean_qi(post, b)$b
  )
  out <- bind_cols(
    out, tt1, tt2
  )
  out$rho <- list(rho)
  return(out)
  
}
