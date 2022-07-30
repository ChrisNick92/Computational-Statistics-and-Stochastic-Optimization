
# Exercise 4 - Variable selection
library(lars)
library(data.table)
library(sys)

# 4.a - Find the best model by minimizing BIC score

# Full enumeration

data("diabetes")

N <- length(diabetes$y)

dt <- data.table(y = diabetes$y, age = diabetes$x[1:N,1],
                 sex = diabetes$x[1:N,2], bmi = diabetes$x[1:N,3],
                 map = diabetes$x[1:N,4], tc = diabetes$x[1:N,5],
                 ldl = diabetes$x[1:N,6], hdl = diabetes$x[1:N,7],
                 tch = diabetes$x[1:N,8],ltg = diabetes$x[1:N,9],
                 glu = diabetes$x[1:N,10])

Minimize_bic <- function(data, target_name = "y"){
  dt <- data
  variables <- names(data)
  target_pos <- which(names(data) == target_name)
  variables <- replace(variables, c(1,target_pos),
                       c(target_name,variables[1]))
  variables[1] <- "y"
  names(dt) <- variables
  num_variables <- length(variables)
  dt <- dt[,..variables]
  x <- replicate(length(variables)-1,c(0,1))
  binarization <- expand.grid(lapply(seq_len(ncol(x)),
                                     function(i) x[,i]))
  num_models <- length(binarization[,1])
  
  best_model <- lm(y~., data = dt[,"y"])
  min_bic <- BIC(best_model)
  best_variables <- c()
  
  tic <- Sys.time()
  for(i in 2:num_models){
    indices <- which(binarization[i,]>0)+1
    temp_variables <- c("y", variables[indices])
    temp_model <- lm(y~., data = dt[,..temp_variables])
    temp_bic <- BIC(temp_model)
    if(temp_bic<min_bic){
      best_model <- temp_model
      min_bic <-temp_bic
      best_variables <- temp_variables[-1]
    }
  }
  tac <- Sys.time()
  print(sprintf("- Exhaustive search finished in %.5f sec(s).",
                tac-tic))
  print(sprintf("- %d models were examined!", num_models))
  print(sprintf("- The best model contains the variables %s",
          paste(best_variables,collapse=",")))
  print(sprintf("- BIC value: %.4f",BIC(best_model)))
  return(best_model)
}
M1 <- Minimize_bic(dt)



# 4.b

# Importing the required package
library(glmnet)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)
library(latex2exp)


lasso <- glmnet(data.matrix(dt[,2:11]),dt$y)

betas = as.matrix(lasso$beta)
lambdas = lasso$lambda
names(lambdas) = colnames(betas)


as.data.frame(betas) %>% 
  tibble::rownames_to_column("variable") %>% 
  pivot_longer(-variable) %>% 
  mutate(lambda=lambdas[name]) %>% 
  ggplot(aes(x=lambda,y=value,col=variable)) + 
  geom_line() + 
  geom_label_repel(data=~subset(.x,lambda==min(lambda)),
                   aes(label=variable),nudge_x=-0.5) +
  scale_x_log10()+
  labs(x=TeX("Penalty parameter $\\lambda$"),
       y="Value of estimated parameter")+
  theme(legend.position = "none")


set.seed(413)
lasso2 <- cv.glmnet(data.matrix(dt[,2:11]), dt$y,
                    num_folds = 10)
plot(lasso2)

# M2

coefs1 <- coef(lasso2, s = "lambda.min")
coefs1 <- c("y",rownames(coefs1)[coefs1[,1]!=0][-1])
print(sprintf("- The variable corresponds to lmin are: %s",
                  paste(coefs1[-1],collapse=",")))
M2 <- lm(y~., data = dt[,..coefs1])

summary(M2)


# M3

coefs2 <- coef(lasso2)
coefs2 <- c("y",rownames(coefs2)[coefs2[,1]!=0][-1])
print(sprintf("- The variable corresponds to lse are: %s",
              paste(coefs2[-1],collapse=",")))

M3 <- lm(y~., data = dt[,..coefs2])

summary(M3)


# 4.c

rmse <- function(actual,preds){
  return(sqrt(mean((actual-preds)^2)))}

cross_val <- function(data, models, num_folds = 5, seed=42){
  set.seed(seed)
  dt_shuffled <- data[sample(nrow(data)),]
  folds <- cut(seq(1,nrow(dt_shuffled)),
               breaks = num_folds, labels = FALSE)
  RMSE <- rep(0,length(models))
  for(i in 1:length(models)){
    names(RMSE)[i] <- paste("M", as.character(i), sep ="")
  }
  
  for(i in 1:num_folds){
    indices <- which(folds == i, arr.ind = TRUE)
    test_data <- dt_shuffled[c(indices),]
    train_data <- dt_shuffled[-c(indices),]
    
    for(i in 1:length(models)){
      variables <- c(array(unlist(models[i])),"y")
      model <- lm(y~., data = train_data[,..variables])
      preds <- predict(model, test_data[,-"y"])
      RMSE[i] <- RMSE[i] + rmse(array(unlist(test_data[,"y"])),preds)
    }}
  return(RMSE/num_folds)}

M1_vars <- names(M1$coefficients)[-1]
M2_vars <- coefs1[-1]
M3_vars <- coefs2[-1]

models <- list(M1 = M1_vars, M2 = M2_vars,
               M3 = M3_vars)

RMSE <- cross_val(dt, models)
RMSE

