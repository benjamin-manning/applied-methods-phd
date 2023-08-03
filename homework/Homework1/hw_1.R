library(tidyverse)
library(mets)
library(combinat)
library(ggplot2)
library(RCAL)
setwd('/Users/benjaminmanning/Desktop/applied-methods-phd/homework/Homework1')
df = read.csv('lalonde_nsw.csv')

####### 1 ##########
#### a
exper_df = df |> 
  select(treat, re78)

treated = exper_df |> filter(treat == 1) 
untreated = exper_df |> filter(treat == 0) 
ate_meana = mean(treated$re78 , na.rm = TRUE) - mean(untreated$re78, na.rm = TRUE)
ate_meana


#### b It should be the same right? the treatment was completely randomized

#### c
# Define the list
treatment_status <- df$treat

# Define the number of permutations
n_permutations <- 1000

# Initialize a list to store the permutations
my_permutations <- vector("list", n_permutations)

# Generate and store the permutations
for (i in 1:n_permutations) {
  my_permutations[[i]] <- sample(treatment_status)
}

# Print the first permutation to check
print(my_permutations[[1]])

# Initialize a variable to store the number of duplicates
num_duplicates <- 0

# Loop over all permutations
for (i in 1:(n_permutations-1)) {
  for (j in (i+1):n_permutations) {
    # Check if the i-th and j-th permutations are identical
    if (all(my_permutations[[i]] == my_permutations[[j]])) {
      num_duplicates <- num_duplicates + 1
      break
    }
  }
}

# Print the number of duplicates
print(paste("Number of duplicate permutations:", num_duplicates))



# Create an empty vector to store the ate_mean values
ate_means <- vector()

# Loop over the permutations
for (data in my_permutations){
  
  df_rand <- data.frame(
    treat = data,
    re78 = df$re78
  )
  
  treated <- df_rand |> filter(treat == 1) 
  untreated <- df_rand |> filter(treat == 0) 
  ate_mean <- mean(treated$re78 , na.rm = TRUE) - mean(untreated$re78, na.rm = TRUE)
  
  # Append the ate_mean value to the vector
  ate_means <- c(ate_means, ate_mean)
}

# Print the ate_mean values
print(ate_means)
# Create a histogram of the ate_means variable
hist(ate_means, main = "Histogram of ATE Means", xlab = "ATE Mean")

# Add a vertical line at a specific point, let's say at x = 2
abline(v = ate_meana, col = "green", lwd = 2)
abline(v = quantile(ate_means, probs = 0.975), col = "red", lwd = 2)
abline(v = quantile(ate_means, probs = 0.025), col = "red", lwd = 2)

ecdf_list = ecdf(ate_means)
ecdf_list(ate_meana)
# Define the observed difference in means
observed_diff <- ate_meana  # replace with your value

# Compute the p-value
p_value <- sum(abs(ate_means) >= abs(observed_diff)) / length(ate_means)

# Print the p-value
print(p_value)
##pvalue" 0.007

##### d
fitd <- lm(re78 ~ treat, data = df)
summary(fitd, robust = T)

## pvalue: 0.00479


####### 2 ##########
#### a
df2 <-  read.csv('lalonde_psid.csv')
df_combo = rbind(df2, df |> filter(treat == 1))
#treated
model <- glm(treat ~ age + education + hispanic + black + married + nodegree + re74 + re75, data = df_combo , family = binomial)

# Get the predicted log-odds
predicted_logodds <- predict(model)

# Get the predicted probabilities
predicted_probabilities <- predict(model, type = "response")

# Print or view the predicted values
print(predicted_logodds)
print(predicted_probabilities)
df_combo$pred_probs <- predicted_probabilities


ggplot(df_combo, aes(x = predicted_probabilities, group = treat, fill = as.factor(treat))) +
  geom_density(alpha = 0.5) +
  scale_fill_discrete(name = "Group") +
  theme_minimal()

#### b 
#control mean:
#(Y_i (1-D_i))/ 1-pi(X_i))
E_y0 <- (df_combo$re78 * (1-df_combo$treat))/(1-df_combo$pred_probs)
mean(E_y0)
E_y1 <- (df_combo$re78 * (df_combo$treat))/(df_combo$pred_probs)
mean(E_y1)
ate_p <- mean(E_y1)- mean(E_y0)
ate_p 

## sipw
ate_sipw <- mean(E_y1) / mean(((1-df_combo$treat))/(1-df_combo$pred_probs)) - mean(E_y0) / mean((df_combo$treat)/(df_combo$pred_probs))
ate_sipw


# Assume your data frame is df
df_combo$untreat_prob <- ifelse(df_combo$treat == 0, df_combo$pred_probs, NA)
df_combo$treat_prob <- ifelse(df_combo$treat == 1, df_combo$pred_probs, NA)

ate.ipw(df_combo$re78, df_combo$treat, df_combo[, (ncol(df_combo)-1):ncol(df_combo)])
# Available in randomized settings and observational settings with unconfoundedness+overlap
data = df_combo
# Estimate the propensity score e(X) via logistic regression using splines
fmla <- as.formula(paste0("~", paste0("bs(", covariates, ", df=3)", collapse="+")))
W <- df_combo$treat
Y <- df_combo$re78
XX <- model.matrix(fmla, data)
logit <- cv.glmnet(x=XX, y=W, family="binomial")
e.hat <- predict(logit, XX, s = "lambda.min", type="response")

# Using the fact that
z <- Y * (W/df_combo$pred_probs - (1-W)/(1-df_combo$pred_probs))
ate.est <- mean(z)
ate.se <- sd(z) / sqrt(length(z))
ate.tstat <- ate.est / ate.se
ate.pvalue <- 2*(pnorm(1 - abs(ate.est/ate.se)))
ate.results <- c(estimate=ate.est, std.error=ate.se, t.stat=ate.tstat, pvalue=ate.pvalue)
print(ate.results)
##SAME AS MINE! : ate_p

#### c
fit2c <- lm(re78 ~ treat + age + education + hispanic + black + married + nodegree + re74 + re75, data = df_combo)
summary(fit2c, robust = T)

#### d
df_crump <- df_combo |> 
  filter(pred_probs >= 0.1 & pred_probs <= 0.9 )
Y_i <- df_crump$re78 
D_i <- df_crump$treat
prop_hat <- df_crump$pred_probs
E_y0 <- (Y_i * (1-D_i))/(1-prop_hat)
mean(E_y0)
E_y1 <- (Y_i * (D_i))/(prop_hat)
mean(E_y1)
ate_p <- mean(E_y1)- mean(E_y0)
ate_p 

## sipw
ate_sipw <- mean(E_y1) / mean(((1-D_i))/(1-prop_hat)) - mean(E_y0) / mean((D_i)/(prop_hat))
ate_sipw

##### black and nonblack only
#### d
#black
df_black <- df_combo |> 
  filter(pred_probs >= 0.1 & pred_probs <= 0.9 ) |> 
  filter(black == 1)
Y_i <- df_black$re78 
D_i <- df_black$treat
prop_hat <- df_black$pred_probs
E_y0 <- (Y_i * (1-D_i))/(1-prop_hat)
mean(E_y0)
E_y1 <- (Y_i * (D_i))/(prop_hat)
mean(E_y1)
ate_p <- mean(E_y1)- mean(E_y0)
ate_p 

## sipw
ate_sipw <- mean(E_y1) / mean(((1-D_i))/(1-prop_hat)) - mean(E_y0) / mean((D_i)/(prop_hat))
ate_sipw

## COMPARISON
exper_df = df |> 
  filter(black == 1) |> 
  select(treat, re78) 

treated = exper_df |> filter(treat == 1) 
untreated = exper_df |> filter(treat == 0) 
ate_meana = mean(treated$re78 , na.rm = TRUE) - mean(untreated$re78, na.rm = TRUE)
ate_meana


#nonblack
df_nonblack <- df_combo |> 
  filter(pred_probs >= 0.1 & pred_probs <= 0.9 ) |> 
  filter(black == 0)
Y_i <- df_nonblack$re78 
D_i <- df_nonblack$treat
prop_hat <- df_nonblack$pred_probs
E_y0 <- (Y_i * (1-D_i))/(1-prop_hat)
mean(E_y0)
E_y1 <- (Y_i * (D_i))/(prop_hat)
mean(E_y1)
ate_p <- mean(E_y1)- mean(E_y0)
ate_p 

## sipw
ate_sipw <- mean(E_y1) / mean(((1-D_i))/(1-prop_hat)) - mean(E_y0) / mean((D_i)/(prop_hat))
ate_sipw

## COMPARISON
exper_df = df |> 
  filter(black == 0) |> 
  select(treat, re78) 

treated = exper_df |> filter(treat == 1) 
untreated = exper_df |> filter(treat == 0) 
ate_meana = mean(treated$re78 , na.rm = TRUE) - mean(untreated$re78, na.rm = TRUE)
ate_meana




