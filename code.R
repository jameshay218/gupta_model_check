## install.packages("deSolve")
## devtools::install_github("jameshay218/lazymcmc")
library(tidyverse)
library(lazymcmc)
library(doParallel)
library(patchwork)
setwd("~/Documents/GitHub/gupta_model_check")

rerun_fits <- TRUE

n_clusters <- 9
cl <- makeCluster(n_clusters)
registerDoParallel(cl)

mcmcPars1 <- c("iterations"=100000,"popt"=0.44,"opt_freq"=1000,
              "thin"=10,"adaptive_period"=50000,"save_block"=1000)
mcmcPars2 <- c("iterations"=200000,"popt"=0.234,"opt_freq"=1000,
               "thin"=10,"adaptive_period"=100000,"save_block"=1000)

get_index_par <- function(chain, index){
  par_names <- colnames(chain)[2:(ncol(chain)-1)]
  par <- as.numeric(chain[chain$sampno == index, 2:(ncol(chain)-1)])
  names(par) <- par_names
  return(par)
}

runnames <- c("prior_only_lourenco1","prior_only_lourenco2","prior_only_me","fitting_lourenco1","fitting_lourenco2","fitting_me")
n_runs <- length(runnames)
prior_controls <- c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE)
prior_use <- c(1,2,3,1,2,3)
n_chains <- 3

runnames <- rep(runnames, each=n_chains)
prior_controls <- rep(prior_controls,each=n_chains)
prior_use <- rep(prior_use, each=n_chains)
chain_nos <- rep(1:n_chains, n_runs)

## Define the set of ODEs for the model. Return such that we can solve with deSolve
SIR_odes <- function(t, x, params) {
  y <- x[1]
  z <- x[2]
  dead <- x[3]
  
  ## Extract model parameters
  R0 <- params[1]
  sigma <- params[2]
  rho <- params[3]
  theta <- params[4]
  N <- params[5]
  
  beta <- R0*sigma
  
  dY <- beta*y*(1-z) - sigma*y
  dZ <- beta*y*(1-z)
  
  ## Note need to shift this number forward by psi days post-hoc
  list(c(dY,dZ))
}

## ROUGH PARAMETERS FROM THE PAPER
sigma <- 1/4.5 ## 1/Infectious period
R0 <- 2.5 ## Basic reproductive number
rho <- 0.01 ## Proportion of population at risk of severe disease
theta <- 0.14 ## Probability of dying with severe disease
psi <- 17 ## Time between infection and death
t0 <- 31 + 28
psi <- 17

## Population of UK
N <- 66870000

## EXTRACT AND CLEAN DATA
## Number dead in UK to date
obs_deaths <- read_csv("~/Documents/GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv")
obs_deaths <- obs_deaths %>% filter(`Country/Region` == "United Kingdom") %>% select(-c("Province/State","Country/Region","Lat","Long"))
obs_deaths <- colSums(obs_deaths)
obs_deaths <- data.frame(time=names(obs_deaths),dead=obs_deaths)
obs_deaths$time <- as.Date(as.character(obs_deaths$time),origin="2019-12-01", format="%m/%d/%y")
obs_deaths <- obs_deaths %>% filter(time <= "2020-03-19")
obs_deaths <- obs_deaths %>% filter(time >= "2020-03-05")
data <- obs_deaths

## Take a look to see that model makes sense
## Times to solve model over
t <- seq(0,100,by=1)

## Note starting conditions for population size - we're working per capita here
y0 <- 1/N
results <- as.data.frame(deSolve::ode(y=c(y=y0,z=y0),
                                      times=t, func=SIR_odes,
                                      parms=c(R0,sigma, rho, theta, N)))


## Deaths for that time are actually reported psi days later
results$time_deaths <- results$time + psi + t0
results$time <- results$time + t0
results$susc <- 1 - results$z
results$dead <- N*rho*theta*results$z
results$time <- as.Date(results$time, origin="2019-12-01")
results$time_deaths <- as.Date(results$time_deaths, origin="2019-12-01")

par(mfrow=c(2,1))
plot(results[,c("time","susc")], type='l',col="green",xlab="Time (days)", ylab="Proportion of population",ylim=c(0,1))
lines(results[,c("time","z")],col="blue")
lines(results[,c("time","y")],col="red")
legend("topright", title="Key",
       legend=c("Proportion infectious","Proportion no longer susceptible","Proportion susceptible"),
       col=c("red","blue","green"),
       lty=c(1,1,1))
plot(results[,c("time_deaths","dead")],type='l')


##################################
## NOW SET UP CODE TO FIT MODEL
t_long <- seq(0,365,by=1)
pars <- c(R0, sigma, rho, theta, N, psi, t0)
names(pars) <- c("R0","sigma","rho","theta","N","psi","t0")
## Putting model solving code in a function for later use
solve_model <- function(pars, t){
  N <- pars["N"]
  ## Note starting conditions for population size - we're working per capita here
  results <- as.data.frame(deSolve::ode(y=c(y=1/N,z=1/N),
                                        times=t, func=SIR_odes,
                                        parms=pars))

  return(results)
}

## We need to put our likelihood function in a closure environment
## It's important to write the function in this form!
create_lik <- function(parTab, data, PRIOR_FUNC,t, PRIOR_ONLY=FALSE){
  par_names <- parTab$names
  
  ## Extract observed incidence
  obs_dead <- data$dead
  max_time <- max(data$time)
  min_time <- min(data$time)
  
  likelihood_func <- function(pars){
    names(pars) <- par_names
    ## Solve model
    if(!PRIOR_ONLY){
      results <- solve_model(pars, t)
      
      ## Deaths for that time are actually reported psi days later
      psi <- pars["psi"]
      t0 <- pars["t0"]
      theta <- pars["theta"]
      rho <- pars["rho"]
      
      ## Shift times and convert to dates
      results$time_deaths <- results$time + as.integer(psi + t0)
      results$time <- results$time + as.integer(t0)
      results$dead <- N*rho*theta*results$z
      results$time <- as.Date(results$time, origin="2019-12-01")
      results$time_deaths <- as.Date(results$time_deaths, origin="2019-12-01")
      
      ## Get deaths that match data
      predicted1 <-  results %>%
        filter(time_deaths <= max_time & time_deaths >= min_time) %>% pull(dead)
      lik <- sum(dpois(x=obs_dead,lambda=predicted1,log=TRUE))
    } else {
      lik <- 0
    }
    if(!is.null(PRIOR_FUNC)) lik <- lik + PRIOR_FUNC(pars)
    lik
  }
}

## Control parameters in MCMC
parTab <- data.frame(names=c("R0","sigma", "rho", "theta", "N", "psi", "t0"),
                     values=pars,
                     fixed=c(0,0,0,0,1,0,0),
                     steps=c(0.1,0.1,0.1,0.1,0.1,0.1,0.1),
                     lower_bound=c(1,0,0,0,0,0,0),
                     upper_bound=c(10,10,0.1,10,100000000,25,75))


## Test the model solves
f <- create_lik(parTab, data, NULL, t=t_long)
f(pars)

## Starting points, chosen pseudo at random
## seeding chains for SIR models is hard with deSolve,
## so I've chosen points near the true values.
startTab <- parTab

prior_func_lourenco1 <- function(pars){
  names(pars) <- parTab$names
  p1 <- dnorm(pars["R0"], 2.25, 0.025,log=TRUE)
  p2 <- dnorm(1/pars["sigma"], 4.5, 1, log=TRUE)
  p3 <- dnorm(pars["psi"], 17, 2, log=TRUE)
  p4 <- dgamma(pars["rho"], shape=5,rate=5/0.001,log=TRUE)
  p5 <- dnorm(pars["theta"],0.14,0.007,log=TRUE)
  return(sum(p1,p2,p3,p4,p5))
}
prior_func_lourenco2 <- function(pars){
  names(pars) <- parTab$names
  p1 <- dnorm(pars["R0"], 2.75, 0.025,log=TRUE)
  p2 <- dnorm(1/pars["sigma"], 4.5, 1, log=TRUE)
  p3 <- dnorm(pars["psi"], 17, 2, log=TRUE)
  p4 <- dgamma(pars["rho"], shape=5,rate=5/0.01,log=TRUE)
  p5 <- dnorm(pars["theta"],0.14,0.007,log=TRUE)
  return(sum(p1,p2,p3,p4,p5))
}

prior_func_me <- function(pars){
  names(pars) <- parTab$names
  p1 <- dnorm(pars["R0"], 2.5, 0.5,log=TRUE)
  p2 <- dnorm(1/pars["sigma"], 4.5, 1, log=TRUE)
  p3 <- dnorm(pars["psi"], 17, 2, log=TRUE)
  p4 <- 0
  p5 <- dnorm(pars["theta"],0.14,0.007,log=TRUE)
  return(sum(p1,p2,p3,p4,p5))
}
top_wd <- getwd()

if(rerun_fits){
  
res <- foreach(i=seq_along(runnames),.packages=c("lazymcmc","tidyverse")) %dopar% {
    setwd(top_wd)
    setwd("chains")
    
    if(!file.exists(runnames[i])) {
      dir.create(runnames[i])
    }
    setwd(runnames[i])
    
    filename_tmp <- paste0(runnames[i],"_",chain_nos[i])
    prior_only <- prior_controls[i]
    prior_use_tmp <- prior_use[i]
    
    if(prior_use_tmp == 1){
      prior_func <- prior_func_lourenco1
    } else if(prior_use_tmp == 2) {
      prior_func <- prior_func_lourenco2
    } else {
      prior_func <- prior_func_me
    }
    
    output <- run_MCMC(parTab=startTab, data=obs_deaths, mcmcPars=mcmcPars1, filename=filename_tmp,
                       CREATE_POSTERIOR_FUNC = create_lik, mvrPars = NULL, PRIOR_FUNC=prior_func, t=t_long,
                       PRIOR_ONLY=prior_only)
    chain <- read.csv(output$file)
    best_pars <- get_best_pars(chain)
    chain <- chain[chain$sampno >= mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]
    covMat <- cov(chain)
    mvrPars <- list(covMat,0.5,w=0.8)
    
    ## Start from best location of previous chain
    startTab$values <- best_pars
    
    ## Run second chain
    output <- run_MCMC(parTab=startTab, data=obs_deaths, mcmcPars=mcmcPars2, filename=filename_tmp,
                       CREATE_POSTERIOR_FUNC = create_lik,PRIOR_FUNC=prior_func, t=t_long,
                       PRIOR_ONLY=prior_only, mvrPars=mvrPars)
  }
}
  

## Read in the MCMC chains
setwd(top_wd)
setwd("chains")

all_runs <- list.files()
all_chains <- NULL
for(run in all_runs){
  chains <- load_mcmc_chains(paste0(top_wd, "/chains/",run),parTab,TRUE,1,200000)
  
  tmp_chain <- as.data.frame(chains[[2]])
  tmp_chain$sampno <- 1:nrow(tmp_chain)
  tmp_chain$run <- run
  tmp_chain <- reshape2::melt(tmp_chain, id.vars=c("sampno","run"))
  all_chains <- rbind(all_chains, tmp_chain)
  
  pdf(paste0(top_wd,"/",run,"_chains.pdf"))
  plot(chains[[1]])
  dev.off()
}
all_chains1 <- all_chains
all_chains <- all_chains %>% mutate(value=ifelse(variable=="sigma", 1/value, value))
all_chains$group <- "new"
all_chains <- all_chains %>% mutate(group = ifelse(run %in% c("fitting_lourenco1","prior_only_lourenco1"), "original_low", group),
                                    group = ifelse(run %in% c("fitting_lourenco2","prior_only_lourenco2"), "original_high",group))

var_key <- c("R0"="R[0]",
             "sigma"="sigma",
             "rho"="rho",
             "theta"="theta",
             "psi"="phi",
             "t0"="Seed_date",
             "lnlike"="Posterior_probability")

run_names_key <- c("fitting_lourenco1"="After fitting to data, original",
                   "prior_only_lourenco1"="Before fitting to data, original",
                   "fitting_lourenco2"="After fitting to data, original",
                   "prior_only_lourenco2"="Before fitting to data, original",
                   "fitting_me"="After fitting to data, my analysis",
                   "prior_only_me"="Before fitting to data, new")

all_chains$run <- run_names_key[all_chains$run]

all_chains$variable <- var_key[all_chains$variable]
all_chains$variable <- factor(all_chains$variable, levels=var_key)

colnames(all_chains)[2] <- "Version"

blank_limits <- data.frame(variable=c("R[0]","R[0]","rho","rho"),x=c(1,4,0,0.1))


## Plot the posteriors and priors
p1 <- ggplot(all_chains[all_chains$group == "original_low",]) + 
  geom_density(aes(x=value,fill=Version), alpha=0.5) +
  geom_blank(data=blank_limits,aes(x=x))+
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  facet_wrap(~variable,scales="free",ncol=3, labeller=label_parsed) + 
  xlab("Estimate") + ylab("Posterior density") +
  theme_bw() +
  theme(legend.position=c(0.7,0.2))
p2 <- ggplot(all_chains[all_chains$group == "original_high",]) + 
  geom_density(aes(x=value,fill=Version), alpha=0.5) +
  geom_blank(data=blank_limits,aes(x=x))+
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  facet_wrap(~variable,scales="free",ncol=3, labeller=label_parsed) + 
  xlab("Estimate") + ylab("Posterior density") +
  theme_bw() +
  theme(legend.position=c(0.7,0.2))
p3 <- ggplot(all_chains[all_chains$group == "new",]) + 
  geom_density(aes(x=value,fill=Version), alpha=0.5) + 
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  geom_blank(data=blank_limits,aes(x=x))+
  facet_wrap(~variable,scales="free",ncol=3, labeller=label_parsed) +
  xlab("Estimate") + ylab("Posterior density") +
  theme_bw() +
  theme(legend.position=c(0.7,0.2))

png(paste0(top_wd,"/plots/densities_original_low.png"),height=5,width=8,res=300,units="in")
p1
dev.off()

png(paste0(top_wd,"/plots/densities_original_high.png"),height=5,width=8,res=300,units="in")
p2
dev.off()

png(paste0(top_wd,"/plots/densities_new.png"),height=5,width=8,res=300,units="in")
p3
dev.off()

all_chains_subset <- all_chains %>% filter(Version %in% c("After fitting to data, original", "After fitting to data, my analysis") & variable %in% var_key[names(var_key) %in% c("R0","sigma","rho","t0")])
all_chains_subset$variable <- as.character(all_chains_subset$variable)
all_chains_subset <- all_chains_subset %>% pivot_wider(names_from=variable,values_from=value)

## Plot the 2D densities
p4 <- ggplot(all_chains_subset) + 
  geom_hex(aes(x=`sigma`,y=`R[0]`,fill=..density..),bins=100) +
  ylab("Basic reproductive number, R[0]") + xlab("Infectious period in days") +
  facet_wrap(~Version) +
  scale_fill_gradient2(low="#5E4FA2",mid="#FAFDB8",high="#9E0142",midpoint= 0.008) +
  theme_bw()

p5 <- ggplot(all_chains_subset) + 
  geom_hex(aes(y=`rho`,x=as.Date(`Seed_date`,origin="2019-12-01"),fill=..density..),bins=100) +
  xlab("Seed date") + ylab("Proportion of population at\n risk of severe disease, ρ") +
  scale_fill_gradient2(low="#5E4FA2",mid="#FAFDB8",high="#9E0142",midpoint= 0.008) +
  facet_wrap(~Version) + theme_bw()

png(paste0(top_wd,"/plots/correlations.png"),height=6,width=8,res=300,units="in")
p4/p5
dev.off()


## Calculate 95% CI on proportion infected
setwd(top_wd)
setwd("chains")
run <- "fitting_me/"
chains <- load_mcmc_chains(paste0(top_wd, "/chains/",run),parTab,FALSE,1,200000)
chain <- as.data.frame(chains[[2]])  
chain$sampno <- 1:nrow(chain)
samps <- sample(unique(chain$sampno), 30000)

prop_infected <- numeric(30000)
rhos <- numeric(30000)
t0s <- numeric(30000)
for(i in seq_along(samps)){
  pars <- get_index_par(chain, samps[i])
  results <- solve_model(pars, t_long)
  rhos[i] <- pars["rho"]
  t0s[i] <- pars["t0"]
  psi <- pars["psi"]
  t0 <- pars["t0"]
  theta <- pars["theta"]
  rho <- pars["rho"]
  results$time_deaths <- results$time + as.integer(psi + t0)
  results$time <- results$time + as.integer(t0)
  results$dead <- N*rho*theta*results$z
  results$time <- as.Date(results$time, origin="2019-12-01")
  results$time_deaths <- as.Date(results$time_deaths, origin="2019-12-01")
  
  prop_infected[i] <- results[results$time == "2020-03-19","z.N"]
}
print(paste0("Proportion immune: ", quantile(prop_infected,c(0.01,0.025,0.5, 0.975,0.99))*100))

final_dat <- data.frame(rho=rhos,prop_infected=prop_infected)

p7 <- ggplot(final_dat) + 
  geom_hex(aes(y=`rho`,x=prop_infected,fill=..density..),bins=100) +
  xlab("Proportion of UK population immune by 19/03/2020") + ylab("Proportion of population at\n risk of severe disease ρ") +
  scale_fill_gradient2(low="#5E4FA2",mid="#FAFDB8",high="#9E0142",midpoint= 0.0008) + theme_bw() + theme(legend.position=c(0.8,0.7))

png(paste0(top_wd,"/plots/final_new.png"),height=4,width=5,res=300,units="in")
p7
dev.off()

