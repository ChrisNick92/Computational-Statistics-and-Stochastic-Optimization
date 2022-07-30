
# Exercise 1 - Stochastic Simulation

# Importing libraries

library(data.table)
library(ggplot2)
library(latex2exp)

# 1.a) Rejection sampling

# Proposal function --> Cauchy

# Sample from N(0,1)

M <- sqrt(2*pi/exp(1)) # Optimal M

dt <- data.table(x=-5:5)

ggplot(dt,aes(x))+
  stat_function(fun = function(x) dnorm(x)/dcauchy(x),
                colour = "blue", size = 1)+ylab("h(x)")+
  labs(title = TeX("The graph of the function $h(x)=f(x)/g(x)$"))



samples <- 1000

rejection_simulation <- function(num_samples, M,seed = 42){
  set.seed(seed)
  sampled <- 0
  rejections <- 0
  trials <- 0
  samples = c()
  while (sampled < num_samples){
    accept = FALSE
    while(!accept){
      u <- runif(2)
      y <- tan(pi*(u[1]-1/2))
      trials <- trials +1
      if (M*dcauchy(y)*u[2] <= dnorm(y)){
        accept = TRUE
        samples <- append(samples,y)
        sampled <- sampled +1
      }
      else{rejections <- rejections + 1}
    }
  }
  info <- c("samples", "trials", "rejections")
  l <- list(samples = samples, trials = trials,
            rejections = rejections, info = info)
  return(l)
}

l <- rejection_simulation(num_samples = 1000, M = M)
print(paste("Probability of acceptance: ", 1-l$rejections/l$trials))
print(paste("Theoretical acceptance: ", 1/M))

# Visualize the results of sampling
h_opt <- 3.491/(1000)^(1/3) # optimal bin-width for Normal Distribution
dataset <- data.frame(X = l$samples)
m <- ggplot(dataset, aes(x=X))
m + geom_histogram(aes(y = ..density..),binwidth = h_opt,
                   color = "#24438C",fill = "#77AABB")+
  geom_density(alpha = 0.5, size = 1, aes(color = "Estimated Density"),
               key_glyph = draw_key_path)+
  stat_function(fun = dnorm,size = 1.2, aes(color = "Normal PDF"), key_glyph =draw_key_path)+
  labs(title = "Histogram of samples generated from N(0,1)",
       subtitle = "Method: Rejection sampling")+
  scale_colour_manual("Densities",values = c("blue","black"))+
  theme(legend.position = c(0.9, 0.9))
  
  

# Show the results - Summary of Sampling
print(paste("- Total trials:", l$trials))
print(paste("- Total rejections:", l$rejections))
print(paste("- Acceptance probability:", 1-l$rejections/l$trials))
print(paste("- Theoretical acceptance prob:", 1/M))
print(paste("- Estimated mean:", mean(l$samples)))
print(paste("- Estimated standard deviation", sd(l$samples)))


# 1.b)
conditional_sampling <- function(n=5000,r=1,p=0.5,seed = 42){
  set.seed(seed)
  lambdas <- rgamma(n=n,shape = r, scale = (1-p)/p) # Sample lambdas
  x <- rep(0,n)
  for(i in 1:n){
    x[i] <- rpois(1,lambdas[i])
  }
  info <- c("x", "lambdas")
  l <- list(info = info, x = x,
            lambdas = lambdas)
  return (l)
}

# Converge of delta_n and delta_n_star
max_sample <- 5000
num_samples <- seq(20,max_sample,10)
iterations <- length(num_samples)
mean_x <- rep(0,iterations)
mean_lambdas <- rep(0,iterations)

for(i in 1:iterations){
  samples <- conditional_sampling(n=num_samples[i])
  mean_x[i] <- mean(samples$x)
  mean_lambdas[i] <- mean(samples$lambdas)
}

d1 <- data.table(values = mean_lambdas, Mean = rep("lambdas", iterations),
                 Iterations = num_samples)
d2 <- data.table(values = mean_x, Mean = rep("X's", iterations),
                 Iterations = num_samples)
d3 <- data.table(values = rep(1,iterations), Mean = rep("True",iterations),
                 Iterations = num_samples)

dt <- rbindlist(list(d1,d2,d3))
ggplot(data = dt, aes(x= Iterations, y = values, color = Mean))+
  geom_line(size = 1.3)+
  labs(title = TeX("The converge of $\\delta_n$ and $\\delta_n^*$ to $E(X)=2$."),
       x = "Number of samples", y = "Value of mean")+
  scale_colour_manual("Mean",values = c("#4FA3AB","#709716","#F07B7B"))+
  theme(legend.position = c(0.9, 0.9))

# Comparison of variances
samples <- conditional_sampling(n=5000)
sprintf("- The variance of samples from X is %.3f",var(samples$x))
sprintf("- The variance of samples from L is %.3f",var(samples$lambdas))

# More experiments 

n <- 5000
p <- 0.5
r <- 1
sequence <- seq(20,3000,1)
iterations <- length(sequence)
delta <- rep(0,3000)
delta_star <- rep(0,3000)
var_delta <- rep(0,iterations)
var_delta_star <- rep(0, iterations)

for(i in 1:19){
  delta_star[i] <- mean(rgamma(n = n, shape = r, 
                               scale = (1-p)/p))
  delta[i] <- mean(rnbinom(n=n, size = r, 
                           prob = p))
}

for(i in 1:iterations){
  delta_star[sequence[i]] <- mean(rgamma(n = n, shape = r, 
                                    scale = (1-p)/p))
  delta[sequence[i]] <- mean(rnbinom(n=n, size = r, 
                                          prob = p))
  var_delta[i] <- var(delta[1:sequence[i]])
  var_delta_star[i] <- var(delta_star[1:sequence[i]])
}
d1 <- data.table(values = var_delta, 
                 Mean = rep("Simulated delta",iterations),
                 sample_size = sequence)
d2 <- data.table(values = rep((1/n)*r*(1-p)/p^2,iterations),
                 Mean = rep("Delta"),
                 sample_size = sequence)
d3 <- data.table(values = var_delta_star,
                 Mean = rep("Simulated delta star",
                            iterations),
                 sample_size = sequence)
d4 <- data.table(values = rep((1/n)*r*(1-p)^2/p^2,iterations),
                 Mean = rep("Delta star"),
                 sample_size = sequence)
dt <- rbindlist(list(d1,d2,d3,d4))
ggplot(data = dt, aes(x= sample_size, y = values, color = Mean))+
  geom_line(size = 0.8)+
  labs(title = TeX("Simulation: Variance of $\\delta_n$ and $\\delta_n^*$."),
                   subtitle = TeX("$n=5000$"),
       x = "Number of samples", y = "Variance")+
  scale_colour_manual("Variance of",values = c("#4FA3AB","#4FA3AB",
                                               "#F07B7B", "#7AB163"))+
  theme(legend.position = c(0.85, 0.88))+
  scale_fill_discrete(labels=c(TeX("Delta"), 'Delta Star',
                               "Simulated delta",
                               "Simulated delta star"))



# 1.c)

T_simulation <- function(n=80, num_samples=10000, seed = 42){
  set.seed(seed)
  samples <- c()
  for(i in 1:num_samples){
    samples <- append(samples, sum(runif(n))^2/n)
  }
  estimated_mean <- mean(samples)
  estimated_sd <- sd(samples)
  true_mean <- n/4-1/12
  info <- c("estimated mean", "estimated sd", "samples", "true mean")
  l <- list(estimated_mean = estimated_mean,
            estimated_sd = estimated_sd, samples = samples,
            info = info, true_mean = true_mean)
}
l <- T_simulation()
print(paste("- Sample mean:", l$estimated_mean))
print(paste("- True mean: ", l$true_mean))
print(paste("- Estimated standard deviation: ", l$estimated_sd))


# Histogram of samples
data <- data.table(samples = l$samples)
ggplot(data = data, aes(x = samples))+
  geom_histogram(aes(x = samples),bins=30,
                 color = "#24438C",fill = "#77AABB")+
  labs(title = TeX("The histogram for 10000 samples from $T_{80}$."),
       x=TeX("Value of $T_{80}$"))

# 1.d)

data <- readRDS("data1.rds")
print(paste("- The first 5 elemets of the dataset are: ", data[1:5]))

T_bootstrap <- function(samples, B=10000, seed = 42){
  set.seed(42)
  bootstrap_samples <- c() # Array to store bootstrap samples, size Bx80
  for(i in 1:B){ # Generating a bootstrap sample
    # and append new sample to array
    bootstrap_samples <- c(bootstrap_samples,
                           sample(samples,
                                  size=length(samples),replace = TRUE))
  }
  # Reshape array 1x(Bx80)->(Bx80)
  bootstrap_samples <- matrix(bootstrap_samples, nrow = B, byrow = T)
  # Calculate the mean row-wise
  T_on_bootstrap <- apply(bootstrap_samples, 1, sum)^2/length(samples)
  info <- c("samples", "Estimated sdt")
  l <- list(info = info, samples = T_on_bootstrap,
           estimated_sd = 
             sqrt((1/(B-1))*sum((T_on_bootstrap-mean(T_on_bootstrap))^2))) # step 6 of ALG-2
  return(l)
}

# Jackknife resampling
T_jackknife <- function(samples){
  n <- length(samples)
  theta_ommited <- rep(0,n)
  for(i in 1:n){
    theta_ommited[i] <- (1/(n-1))*sum(samples[-i])^2
  }
  return(sqrt(((n-1)/n)*sum((mean(theta_ommited)-theta_ommited)^2)))
}


bootstrap <- T_bootstrap(data)
print(paste("- Standard error using Bootstrap sampling: ",
            bootstrap$estimated_sd))
print(paste("- Standard error using Jackknife sampling: ", T_jackknife(data)))

# Compare histograms with simulation sampling & Bootstrap sampling

simul <- T_simulation()

dt_1 <- data.table(samples = simul$samples,
                   Method = rep("R - Simulation",length(simul$samples)))
dt_2 <- data.table(samples = bootstrap$samples,
                   Method = rep("Bootstrap - Simulation",
                                length(bootstrap$samples)))
dt <- rbindlist(list(dt_1,dt_2))
ggplot(dt,aes(x = samples, color = Method, fill = Method))+
  geom_histogram(position = "dodge",
                 bins =30)+
  labs(title = TeX("Histogram comparison for $\\hat{\\theta}$"),
       x = "Estimator Value",y = "Number of samples",
       subtitle = TeX("Total samples: $10000$, $n=80$"))+
  theme(legend.position = c(0.9, 0.9))











