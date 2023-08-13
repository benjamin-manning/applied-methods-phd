library(tidyverse)
library(mets)
library(combinat)
library(ggplot2)
library(RCAL)
library(lmtest)
library(plm)
library(sandwich)
setwd('/Users/benjaminmanning/Desktop/applied-methods-phd/homework/Homework3')
df = read.csv('networth_delta_elas.csv')

##### Problem 1 #####
# a QUITE ESOTERIC QUESTION - no need to answer

##### Problem 2 #####
# a
df =  read.csv('yearly_hp.csv')

fit2a = lm(hpi ~ elasticity + as.factor(year) +as.factor(county_fips), data = df)
# save coefficients with year in name
coeffs2a = coef(summary(fit2a))
coefs = coeffs2a[3:17] + coeffs2a[2]
plot(2002:2016, coefs, type = 'l', xlab = 'Year', ylab = 'Elasticity', main = 'Elasticity of Housing Prices by Year')

# Load the dplyr package

df <- df %>%
  mutate(decile = ntile(elasticity, 10))
table(df$decile)

# Save coefs per year
matrix <- matrix(NA, nrow = 10 * 16, ncol = 3)
counter = 1
for (i in 2002:2016){
    df2 = df[df$year == i,]
    fit2 = lm(hpi ~ as.factor(decile) + as.factor(county_fips), data = df2)
    coeffs2 = coef(summary(fit2))
    coefs = coeffs2[2:11, 1] # Assuming you want the estimates, not other summary statistics
    matrix[counter:(counter+9), 1] = rep(i, 10) # Assign the year to the first column
    matrix[counter:(counter+9), 2] = coefs      # Assign the coefficients to the second column
    matrix[counter:(counter+9), 3] = 1:10       # Assign the decile to the third column
    counter = counter + 10
}
#plot color by decile
plot(matrix[,1], matrix[,2], xlab = 'Year', ylab = 'Elasticity', main = 'Elasticity of Housing Prices by Year')
# Convert the matrix to a data frame
data_for_plot <- as.data.frame(matrix)
colnames(data_for_plot) <- c("Year", "Elasticity", "Decile")

# Load ggplot2
# Create the plot
ggplot(data_for_plot, aes(x = Year, y = Elasticity, color = as.factor(Decile))) +
  geom_point() +
    geom_line() +
  labs(x = 'Year', y = 'Elasticity', title = 'Elasticity of Housing Prices by Year') +
  scale_color_discrete(name = "Color Legend") +
  theme_minimal()
