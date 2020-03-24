## install.packages("deSolve")

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
  dDead <- N*rho*theta*z
  list(c(dY,dZ,dDead))
}


sigma <- 1/4.5 ## 1/Infectious period
R0 <- 2.25 ## Basic reproductive number
rho <- 0.001 ## Proportion of population at risk of severe disease
theta <- 0.14 ## Probability of dying with severe disease
psi <- 17 ## Time between infection and death

## Population of UK
N <- 66800000
## Population of Hubei
#N <- 58500000

## Number dead in UK to date
number_dead <- 422
## Number dead by 30th Jan in Hubei
## number_dead <- 162

## Approximate infection prevalence in Hubei on 30th Jan
## approx_prev <- 0.01

## Times to solve model over
t <- seq(0,100,by=1)

## Note starting conditions for population size - we're working per capita here
results <- as.data.frame(deSolve::ode(y=c(y=1/N,z=1/N,Dead=0),
                                      times=t, func=SIR_odes,
                                      parms=c(R0,sigma, rho, theta, N)))

## Deaths for that time are actually reported psi days later
results$time_deaths <- results$time + psi
results$susc <- 1 - results$z

## Number of days since start
days_to_test <- as.integer(as.POSIXct("2020-03-03") - as.POSIXct("2020-01-27"))

par(mfrow=c(2,1))
plot(results[results$time <= days_to_test,c("time","z")], type='l',col="green",xlab="Time (days)", ylab="Per capita")
lines(results[results$time <= days_to_test,c("time","y")], col="red")
lines(results[results$time <= days_to_test,c("time","susc")], col="blue")
abline(v=days_to_test, lty=2)
legend("topleft", title="Key",
       legend=c("Proportion infectious","Proportion no longer susceptible","Proportion susceptible"),
       col=c("red","green","blue"),
       lty=c(1,1,1))

plot(results[results$time_deaths <= days_to_test, c("time_deaths","Dead")],type='l', xlab="Time (days)", ylab="Cumulative number of deaths",xlim=c(0,days_to_test))


