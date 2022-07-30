
# Exercise 2 - Density Estimation

# Importing libraries & data
library(data.table)
library(ggplot2)
library(latex2exp)

sample <- faithful$eruptions
print(paste("- Mean:", mean(sample)))
print(paste("- Standard deviation:", sd(sample)))
print(sample[1:5])

# Utility functions
epanechnikov <- function(x){
  return(ifelse(abs(x)<= 1, (3/4)*(1-x^2), 0))
}

density_estimate <- function(x,h,sample){
  n <- length(sample)
  epan <- rep(0,n)
  for(i in 1:n){
    epan[i]<-epanechnikov((x-sample[i])/h)
  }
  return((1/n)*sum(epan)*1/h)
}

cv_log_likelihood <- function(h,sample){
  n <- length(sample)
  cv_log <- rep(0,n)
  for(i in 1:n){
    cv_log[i] <- log(density_estimate(sample[i],h,sample[-i]))
  }
  return(sum(cv_log))
}

logspace <- function(start, stop, num = 50, base = 10.0) {
  s <- seq(log(start, base), log(stop, base), length.out=num)
  return(base ^ s)
}


tic = Sys.time()
# 1st pass
x <- logspace(start = 0.17, stop =2, num =100)
values <- rep(0, length(x))
for(i in 1:length(x)){values[i] = cv_log_likelihood(x[i],sample)}
toc <- Sys.time()
print(sprintf("- 1st pass finished in %.3f sec(s) and %d were examined.",
              toc-tic,length(x)))
x1 <- which(values == max(values), arr.ind = TRUE)
print(sprintf("- Maximum value %.4f at point %.2f",
              max(values),x[x1]))




# 2nd pass
tic = Sys.time()
x <- seq(0.20,0.22, 0.0001)
values <- rep(0,length(x))
for(i in 1:length(x)){values[i] = cv_log_likelihood(x[i],sample)}
toc <- Sys.time()
print(sprintf("- 1st pass finished in %.3f sec(s) and %d were examined.",
              toc-tic,length(x)))
x1 <- which(values == max(values), arr.ind = TRUE)
print(sprintf("- Maximum value %.4f at point %.2f",
              max(values),x[x1]))
d <- data.frame(cbind(x,values))
ggplot(d, aes(x = x, y = values))+
  geom_path(col = "#E75177", size = 1.2)+xlab("h")+ylab("log-likelikehood")+
  labs(title = 
         TeX("Crossvalidated log-likelihood $\\sum_{i=1}^N\\log \\hat{f}_{h,i}(x_i)$"))


binary_search <- function(a,b,sample){
  mid_point <- (a+b)/2
  ifelse(cv_log_likelihood(a,sample)>cv_log_likelihood(b,sample),
         return(c(a,mid_point)), return(c(mid_point,b)))
}


search_for_maximum <- function(a,b,sample,tol=1e-8){
  while(abs(a-b)>tol){
    interval <- binary_search(a,b,sample)
    a <- interval[1]
    b <- interval[2]
  }
  points <- c(a, (a+b)/2, b)
  values <- c(cv_log_likelihood(a, sample),
              cv_log_likelihood((a+b)/2,sample),
              cv_log_likelihood(b,sample))
  index <- which(values == max(values),arr.ind = TRUE)
  return(points[index])
}

#Find optimal h
h_opt <- search_for_maximum(a=0.20,b=0.21,sample)
print(h_opt)


estimated_density <- density(sample, kernel = "epanechnikov", bw = h_opt,
                             n = 2048)

d <- data.frame(cbind(estimated_density$x,estimated_density$y))
ggplot(d, aes(x = estimated_density$x, y = estimated_density$y))+
  geom_point(col = "#E75177")+xlab("x")+ylab("Likelihood")+
  labs(title = TeX("Estimated Density for bandwidth $h=0.20866$"))+
  theme(text=element_text(size = 10),
        plot.subtitle = element_text(size = 10))




# 2.b) Compare densities

num_pts <- length(estimated_density$x)
custom_estimations <- rep(0,num_pts)
for(i in 1:num_pts){
  custom_estimations[i]<-density_estimate(estimated_density$x[i]
                                          ,h_opt,sample)}
dt_1 <- data.table(x = estimated_density$x,
                   y = custom_estimations, Estimation = "Custom")
dt_2 <- data.table(x = estimated_density$x,
                   y = estimated_density$y, Estimation ="Actual")
d <- rbindlist(list(dt_1,dt_2))

ggplot(d, aes(x = x, y = y), color = Estimation)+
  geom_path(aes(color = Estimation),size=1.2)+xlab("x")+ylab("Density")+
  labs(title = "Estimated density vs actual density")+
  scale_colour_manual("Density",values = c("#F07B7B","#4FA3AB"))+
  theme(legend.position = c(0.92, 0.9))

# 2.c) Calculate P(X>3.5)

final_density <- function(x, seed = 42){
  set.seed(seed)
  return(density_estimate(x,h_opt,sample))
}

int <- integrate(Vectorize(final_density),lower =3.5, upper = Inf)


# 2.d) Sample from the estimated f

# Step 1 - Pick at random a point X from sample
# Step 2 - Pick at random a point Y from kernel
# Step 3 - Then Z = X+hY is a sample from estimated f
set.seed(42)
n <- 250
y <- apply(matrix(runif(3*n,-h_opt,h_opt),3), 2, median)
x <- sample(sample,size=n,replace = TRUE)
z <- x+y

print(paste("- Estimated probability from sample", length(z[z>3.5])/length(z)))


hist(z)
