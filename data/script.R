#ANOVA and Welch tests
library(readr)
library(car)
library(onewaytests)

#games howell test after Schlegel (2016) retrieved from https://gist.github.com/aschleg/ea7942efc6108aedfa9ec98aeb6c2096
games.howell <- function(grp, obs) {
  
  #Create combinations
  combs <- combn(unique(grp), 2)
  
  # Statistics that will be used throughout the calculations:
  # n = sample size of each group
  # groups = number of groups in data
  # Mean = means of each group sample
  # std = variance of each group sample
  n <- tapply(obs, grp, length)
  groups <- length(tapply(obs, grp, length))
  Mean <- tapply(obs, grp, mean)
  std <- tapply(obs, grp, var)
  
  statistics <- lapply(1:ncol(combs), function(x) {
    
    mean.diff <- Mean[combs[2,x]] - Mean[combs[1,x]]
    
    #t-values
    t <- abs(Mean[combs[1,x]] - Mean[combs[2,x]]) / sqrt((std[combs[1,x]] / n[combs[1,x]]) + (std[combs[2,x]] / n[combs[2,x]]))
    
    # Degrees of Freedom
    df <- (std[combs[1,x]] / n[combs[1,x]] + std[combs[2,x]] / n[combs[2,x]])^2 / # Numerator Degrees of Freedom
      ((std[combs[1,x]] / n[combs[1,x]])^2 / (n[combs[1,x]] - 1) + # Part 1 of Denominator Degrees of Freedom 
         (std[combs[2,x]] / n[combs[2,x]])^2 / (n[combs[2,x]] - 1)) # Part 2 of Denominator Degrees of Freedom
    
    #p-values
    p <- ptukey(t * sqrt(2), groups, df, lower.tail = FALSE)
    
    # Sigma standard error
    se <- sqrt(0.5 * (std[combs[1,x]] / n[combs[1,x]] + std[combs[2,x]] / n[combs[2,x]]))
    
    # Upper Confidence Limit
    upper.conf <- lapply(1:ncol(combs), function(x) {
      mean.diff + qtukey(p = 0.95, nmeans = groups, df = df) * se
    })[[1]]
    
    # Lower Confidence Limit
    lower.conf <- lapply(1:ncol(combs), function(x) {
      mean.diff - qtukey(p = 0.95, nmeans = groups, df = df) * se
    })[[1]]
    
    # Group Combinations
    grp.comb <- paste(combs[1,x], ':', combs[2,x])
    
    # Collect all statistics into list
    stats <- list(grp.comb, mean.diff, se, t, df, p, upper.conf, lower.conf)
  })
  
  # Unlist statistics collected earlier
  stats.unlisted <- lapply(statistics, function(x) {
    unlist(x)
  })
  
  # Create dataframe from flattened list
  results <- data.frame(matrix(unlist(stats.unlisted), nrow = length(stats.unlisted), byrow=TRUE))
  
  # Select columns set as factors that should be numeric and change with as.numeric
  results[c(2, 3:ncol(results))] <- round(as.numeric(as.matrix(results[c(2, 3:ncol(results))])), digits = 3)
  
  # Rename data frame columns
  colnames(results) <- c('groups', 'Mean Difference', 'Standard Error', 't', 'df', 'p', 'upper limit', 'lower limit')
  
  return(results)
}

#ariaDNE test for normality of values for principal carnassials in 5 groups
normality.test <- read.csv("ariaDNE_anova_groups.csv", header = TRUE)

shapiro.test(normality.test$X01_carnivora_basal)
shapiro.test(normality.test$X02_carnivora_derived)
shapiro.test(normality.test$X03_dasy_basal)
shapiro.test(normality.test$X04_dasy_derived)
shapiro.test(normality.test$X06_hy_derived)

normality.test.2 <- read.csv("ariaDNE_anova.csv", header = TRUE)
leveneTest(normality.test.2$ariaDNE, normality.test.2$function_all, center=mean)

aov.results <- aov(normality.test.2$ariaDNE ~ normality.test.2$function_all)
summary(aov.results)
TukeyHSD(aov(normality.test.2$ariaDNE ~ normality.test.2$function_all))

welch.df <- normality.test.2
welch.test(ariaDNE ~ function_all, data=welch.df)

gh.output <- games.howell(welch.df$function_all, welch.df$ariaDNE)
gh.output


#DNE test for normality of values for principal carnassials in 5 groups
normality.test <- read.csv("DNE_anova_groups.csv", header = TRUE)

shapiro.test(normality.test$X01_carnivora_basal)
shapiro.test(normality.test$X02_carnivora_derived)
shapiro.test(normality.test$X03_dasy_basal)
shapiro.test(normality.test$X04_dasy_derived)
shapiro.test(normality.test$X05_hy_derived)

normality.test.2 <- read.csv("DNE_anova.csv", header = TRUE)
leveneTest(normality.test.2$DNE, normality.test.2$function_all, center=mean)

aov.results <- aov(normality.test.2$DNE ~ normality.test.2$function_all)
summary(aov.results)
TukeyHSD(aov(normality.test.2$DNE ~ normality.test.2$function_all))

welch.df <- normality.test.2
welch.test(DNE ~ function_all, data=welch.df)

gh.output <- games.howell(welch.df$function_all, welch.df$DNE)
gh.output


# linear regression with alpha angle
lm.test <- read.csv("alpha_angle_regression.csv", header = TRUE)

## OLS regression ariaDNE
lm.results <- lm(ariaDNE ~ angle, data = lm.test)
summary(lm.results)

## WLS regression ariaDNE, script after Foley (2019) retreived from https://rpubs.com/mpfoley73/500818
cal.weights <- 1 / lm(abs(lm.results$residuals) ~ lm.results$fitted.values)$fitted.values^2
cal.lmw <- lm(ariaDNE ~ angle, data = lm.test, weights = cal.weights)
summary(cal.lmw)

## OLS regression DNE
lm.results <- lm(DNE ~ angle, data = lm.test)
summary(lm.results)

# LDA
library(MASS)
## LDA of DNE values
normality.test.lda <- read.csv("DNE_lda_principals_groups.csv", header = TRUE)
shapiro.test(normality.test.lda$basal)
shapiro.test(normality.test.lda$derived)
lda.test <- read.csv("DNE_lda_principals.csv", header = TRUE)
dne.lda.principals <- lda(function_all ~ DNE, data = lda.test)
dne.lda.p <- predict(dne.lda.principals, newdata=lda.test[1])$x
table(dne.lda.p, lda.test[,2])
table.predicted <- table(dne.lda.p, lda.test[,2])

#### calcutate posterior probabilities using CV
dne.cv <- lda(function_all ~ DNE, data = lda.test, CV = TRUE)
dne.cv[["posterior"]]

#### predict group for Proviverra principal carnassial
proviverra.dne<- read.csv("DNE_proviverra_principal.csv")
proviverra.predict <- predict(dne.lda.principals, newdata = proviverra.dne)
proviverra.predict[["posterior"]]

### all carnassials
normality.test.lda <- read.csv("DNE_lda_all_groups.csv", header = TRUE)
shapiro.test(normality.test.lda$basal)
shapiro.test(normality.test.lda$derived)
lda.test <- read.csv("DNE_lda_all.csv", header = TRUE)
dne.lda.principals <- lda(function_all ~ DNE, data = lda.test)
dne.lda.p <- predict(dne.lda.principals, newdata=lda.test[1])$x
table(dne.lda.p, lda.test[,2])
table.predicted <- table(dne.lda.p, lda.test[,2])

#### calcutate posterior probabilities using CV
dne.cv <- lda(function_all ~ DNE, data = lda.test, CV = TRUE)
dne.cv[["posterior"]]

#### predict group for Proviverra carnassials
proviverra.dne<- read.csv("DNE_proviverra.csv")
proviverra.predict <- predict(dne.lda.principals, newdata = proviverra.dne)
proviverra.predict[["posterior"]]


# LDA of ariaDNE values

### confusion matrix and bootstrap scripts after Maindonald (2008) retrieved from https://maths-people.anu.edu.au/~johnm/courses/dm/math3346/2008/pdf/
confusion <- function(actual, predicted, names = NULL, printit = TRUE, prior = NULL) {
  if (is.null(names))
    names <- levels(actual)
  tab <- table(actual, predicted)
  acctab <- t(apply(tab, 1, function(x) x/sum(x)))
  dimnames(acctab) <- list(Actual = names, "Predicted (cv)" = names)
  if (is.null(prior)){
    relnum <- table(actual)
    prior <- relnum/sum(relnum)
    acc <- sum(tab[row(tab) == col(tab)])/sum(tab)
  }
  else {
    acc <- sum(prior * diag(acctab))
    names(prior) <- names
  }
  if (printit)
    print(round(c("Overall accuracy" = acc, "Prior accuracy" = prior), 4))
  if (printit) {
    cat("\nConfusion matrix", "\n")
    print(round(acctab, 4))
  }
  invisible(acctab)
}

boot.function <- function(data) {
  prior <- table(data$function_all)
  prior <- prior/sum(prior)
  filename <- sample(1:1000000, 1)
  index <- sample(1:dim(lda.test)[1], replace = TRUE)
  boot.lda <- lda(function_all ~ ariaDNE, data = data[index, ], CV = TRUE)
  cmat <- confusion(data[index, "function_all"], boot.lda$class, printit = FALSE)
  print(c(acc = round(prior[1] * cmat[1, 1] + prior[2] * cmat[2, 2], 4)))
  boot.export <- (c(acc = round(prior[1] * cmat[1, 1] + prior[2] * cmat[2, 2], 4)))
  #write.csv(boot.export, file = filename, row)}
  return(boot.export)}

### principal carnassials
normality.test.lda <- read.csv("lda_principals_groups.csv", header = TRUE)
shapiro.test(normality.test.lda$basal)
shapiro.test(normality.test.lda$derived)
lda.test <- read.csv("lda_principals.csv", header = TRUE)
ariadne.lda.principals <- lda(function_all ~ ariaDNE, data = lda.test)
ariadne.lda.p <- predict(ariadne.lda.principals, newdata=lda.test[1])$x
table(ariadne.lda.p, lda.test[,2])
table.predicted <- table(ariadne.lda.p, lda.test[,2])

#### calcutate posterior probabilities using CV
ariadne.cv <- lda(function_all ~ ariaDNE, data = lda.test, CV = TRUE)
ariadne.cv[["posterior"]]

#### predict group for Proviverra principal carnassial
proviverra.ariadne<- read.csv("lda_proviverra_principal.csv")
proviverra.predict <- predict(ariadne.lda.principals, newdata = proviverra.ariadne)
proviverra.predict[["posterior"]]

#### bootstrap CV lda
output <- replicate(1000, boot.function(lda.test))
output

### all carnassials
normality.test.lda <- read.csv("lda_all_groups.csv", header = TRUE)
shapiro.test(normality.test.lda$basal)
shapiro.test(normality.test.lda$derived)
lda.test <- read.csv("RStudio/Statistik/statistik_final/lda.csv", header = TRUE)
ariadne.lda.principals <- lda(function_all ~ ariaDNE, data = lda.test)
ariadne.lda.p <- predict(ariadne.lda.principals, newdata=lda.test[1])$x
table(ariadne.lda.p, lda.test[,2])
table.predicted <- table(ariadne.lda.p, lda.test[,2])

#### calcutate posterior probabilities using CV
ariadne.cv <- lda(function_all ~ ariaDNE, data = lda.test, CV = TRUE)
ariadne.cv[["posterior"]]

#### predict group for Proviverra carnassials
proviverra.ariadne<- read.csv("lda_proviverra.csv")
proviverra.predict <- predict(ariadne.lda.principals, newdata = proviverra.ariadne)
proviverra.predict[["posterior"]]

#### bootstrap CV lda
output <- replicate(1000, boot.function(lda.test))
output

