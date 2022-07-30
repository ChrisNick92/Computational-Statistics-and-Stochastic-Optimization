# Implementation of EM for Poisson Mixture
samples <- c(2,7,3,9)

EM_poisson_mixture <- function(samples,p,l,tol=1e-10){
  n <- length(samples)
  lambdas <- matrix(l,ncol = 2)
  reps <- 0
  repeat{
    # E - Step
    w <- matrix(rep(0,2*n), ncol = 2)
    for(i in 1:n){ # fill the matrix w
      p_total <- p[1]*dpois(samples[i], l[1]) + 
        p[2]*dpois(samples[i],l[2])
      w[i,1] <- p[1]*dpois(samples[i],l[1])/p_total
      w[i,2] <- p[2]*dpois(samples[i], l[2])/p_total}
    # M - Step
    p_new <- apply(w,2,mean)
    l_new <- apply(diag(samples) %*% w,2,sum)/apply(w,2,sum)
    lambdas <- rbind(lambdas, l_new)
    reps <-reps +1
    if(sum((l_new - l)^2) <= tol){
      break
    }else{
      p <- p_new
      l <- l_new
    }}
  return(list(p=p_new, l=l_new, reps = reps,
              lambdas = lambdas))}

set.seed(42)
p <- runif(1)
l <- runif(2, min =1, max = 10)

result <- EM_poisson_mixture(samples, p = c(p,1-p), l,
                             tol=1e-10)

print(sprintf("- EM converged after %d iterations",
              result$reps))
sprintf("- Estimate for (p1,p2)=(%.4f, %.4f)",
        result$p[1],result$p[2])
sprintf("- Estimate for (l1,l2)=(%.4f, %.4f)",
        result$l[1],result$l[2])












