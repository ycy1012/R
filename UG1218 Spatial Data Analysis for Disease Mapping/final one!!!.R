
library(spdep)
library(CARBayes)
library(dplyr)
library(coda)
library(sf)


# start with data clean
data <- read.csv("~/Desktop/lung cancer/my_data.csv")
shape <- read_sf(dsn = "~/Desktop/pub_las/pub_las.shp")
shape$local_auth <- gsub("Eilean Siar", "Na h-Eileanan Siar", shape$local_auth)

# imputataion
imputed_data <- mice(data, m=40, method='norm', maxit=10, seed=500)
plot(imputed_data)
stripplot(imputed_data) # Create strip plots
reg.fit.mi <- with(imputed_data, lm(Lung_cancer_count~ Manufacturing + Construction
                                    + Emissions + Earnings + Smoking_Rate)) # Fit a linear regression model
summary(reg.fit.mi$analyses[[3]])
pool.fit <- pool(reg.fit.mi) # Pool the results of the fitted linear models from all imputed datasets to provide
summary(pool.fit)

# Extract the completed data

completed_data <- complete(imputed_data, action = 1) # or use "long" or "all" for different analysis strategies
merged_data <- merge(completed_data, shape, by.x = "Council.Areas", by.y = "local_auth")

# add ratio
merged_data$lung_cancer_ratio <- merged_data$Lung_cancer_count / merged_data$Pop_count

# mean ratio for do neighbours
ratio_sf <- read.csv('~/Desktop/lung cancer/lung_cancer_ratios_2008_2017.csv')
# mean ratio和shape
merged_data_sf <- merge(shape, ratio_sf, by.x = "local_auth", by.y = "council.areas")
head(merged_data_sf)

sf_data <- merged_data_sf %>%
  dplyr::select('local_auth', 'code', 'hectares', 'mean_ratio', 'geometry')
head(sf_data)
class(sf_data) #sf

# finish data clean --------------------------------------------------------------------------------

# neighbour
neighbors <- poly2nb(sf_data)
summary(neighbors)

# Find areas with no links
isolated_area_names <- sf_data$local_auth[c(20, 23, 27)]
isolated_area_names

neighbors[[23]] <- c(neighbors[[23]], 1, 16, 27)
neighbors[[27]] <- c(neighbors[[27]], 23, 1)
neighbors[[1]] <- c(neighbors[[1]], 23, 27)
neighbors[[16]] <- c(neighbors[[16]], 20, 23)
neighbors[[20]] <- c(neighbors[[20]], 4, 16)
neighbors[[4]] <- c(neighbors[[4]], 20)

neighbors[[23]] <- neighbors[[23]][neighbors[[23]] != 0]
neighbors[[27]] <- neighbors[[27]][neighbors[[27]] != 0]
neighbors[[20]] <- neighbors[[20]][neighbors[[20]] != 0]

# Check neighbors
print(neighbors[[23]])
print(neighbors[[27]])
print(neighbors[[20]])

neighbors <- lapply(neighbors, as.integer)
class(neighbors) <- "nb"
summary(neighbors)

# neighbor part over ----------------------------------------------------------------------------------

# weight matrix
weights <- nb2listw(neighbors, style="B", zero.policy=TRUE)
neighbors_list <- as(weights, "listw")$neighbours
W_matrix <- nb2mat(neighbors_list, style="B", zero.policy=TRUE)

# MCMC parameters
M.burnin <- 10000       # Number of burn-in iterations (discarded)
M <- 100000               # Number of iterations retained

# use 16 model predict 17

# car model for 2016 since need to predict for 2017
combined_data_2016 <- filter(completed_data, Year == "X2016")

# BYM
set.seed(444)          # For reproducability
MCMC_bym <- S.CARbym(
  formula = Lung_cancer_count ~ offset(log(Pop_count)), 
  data = combined_data_2016,
  family = "poisson",
  W = W_matrix,
  burnin = M.burnin,
  n.sample = M.burnin + M,    # Total iterations
  n.chains = 3,
  n.cores = 3,
  verbose = FALSE
)

# Broad summary of results
print(MCMC_bym$summary.results)

# Additional results on model fit criteria
MCMC_bym$modelfit

# BYM with all
set.seed(444)          # For reproducability
MCMC_bym_all <- S.CARbym(
  formula = Lung_cancer_count ~ offset(log(Pop_count)) + Smoking_Rate + Manufacturing + 
    Construction + Emissions + Earnings, 
  data = combined_data_2016,
  family = "poisson",
  W = W_matrix,
  burnin = M.burnin,
  n.sample = M.burnin + M,    # Total iterations
  n.chains = 3,
  n.cores = 3,
  verbose = FALSE
)

# Broad summary of results
print(MCMC_bym_all$summary.results)

# Additional results on model fit criteria
MCMC_bym_all$modelfit

# model check 
beta_bym_all <- MCMC_bym_all$samples$beta
gelman_results_bym_all <- gelman.diag(beta_bym_all)
print(gelman_results_bym_all) #1 是1 所以收敛 表明链之间的变异与链内变异相比较小，即模型已经收敛

plot(beta_bym_all)
autocorr.plot(beta_bym_all)
gelman.plot(beta_bym_all)

# check posterior predict model
psi <- MCMC_bym_all$samples$psi

samples_beta_bym_all <- do.call(rbind, lapply(MCMC_bym_all$samples$beta, function(x) x))
samples_psi_bym_all <- do.call(rbind, lapply(MCMC_bym_all$samples$psi, function(x) x))

Nsamples = 1000
postpred_Y <- array(NA,dim=c(32,Nsamples))
postpred_mean <- array(NA,dim=c(32,Nsamples))

for (i_s in 1:Nsamples) {
  # 从合并后的样本中随机选择一组参数样本
  s = 10*i_s
  intercept_s <- samples_beta_bym_all[s, 1]
  beta_smoking_s <- samples_beta_bym_all[s, 2]
  beta_construction_s <- samples_beta_bym_all[s, 3]
  beta_earnings_s <- samples_beta_bym_all[s, 4]
  beta_manufacturing_s <- samples_beta_bym_all[s, 5]
  beta_emissions_s <- samples_beta_bym_all[s,6]
  psi_means_s <- samples_psi_bym_all[s, ]
  
  # 计算预测的均值
  postpred_mean[, i_s] <- with(combined_data_2016, 
                               exp(intercept_s + 
                                     beta_smoking_s * Smoking_Rate + 
                                     beta_construction_s * Construction +
                                     beta_manufacturing_s * Manufacturing +
                                     beta_earnings_s * Earnings + 
                                     beta_emissions_s * Emissions + psi_means_s + 
                                     offset(log(Pop_count))
                               ))
  # 生成预测计数
  postpred_Y[, i_s] <- rpois(32, lambda = postpred_mean[, i_s])
}

# 比较预测数据和实际观测数据
hist(colMeans(postpred_Y))
abline(v = mean(combined_data_2016$Lung_cancer_count), col = "red")

hist(apply(postpred_Y, 2, sd))
abline(v = sd(combined_data_2016$Lung_cancer_count), col = "red")

# predict
car_2017 <- filter(completed_data, Year == "X2017")

# Extract posterior means 
intercept_bym <- MCMC_bym_all$summary.results["(Intercept)", "Mean"]
beta_smoking_bym <- MCMC_bym_all$summary.results["Smoking_Rate", "Mean"]
beta_construction_bym <- MCMC_bym_all$summary.results["Construction", "Mean"]
beta_earnings_bym <- MCMC_bym_all$summary.results["Earnings", "Mean"]
beta_manufacturing_bym <- MCMC_bym_all$summary.results["Manufacturing", "Mean"]
beta_emissions_bym <- MCMC_bym_all$summary.results["Emissions", "Mean"]
psi_means_bym = colMeans(samples_psi_bym_all)

# Predicting lung cancer cases for 2017
car_2017$predicted_lung_cancer_bym <- with(car_2017, 
                                           exp(intercept_bym + 
                                                 beta_smoking_bym * Smoking_Rate + 
                                                 beta_construction_bym * Construction +
                                                 beta_manufacturing_bym * Manufacturing +
                                                 beta_earnings_bym * Earnings + 
                                                 beta_emissions_bym * Emissions + psi_means_bym + 
                                                 offset(log(Pop_count))
                                           )
)


# Add a column to compare the predicted cases to actual cases
car_2017$comparison_b <- with(car_2017, predicted_lung_cancer_bym - Lung_cancer_count)
mean(car_2017$comparison_b^2) # 489.4742

# Create a scatter plot 
ggplot(car_2017, aes(x = Lung_cancer_count, y = predicted_lung_cancer_bym)) +
  geom_point() +  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Actual Lung Cancer Count", y = "Predicted Lung Cancer Count", title = "Actual vs Predicted Lung Cancer Cases") +
  theme_minimal()  


# LEROUX
set.seed(444)          # For reproducability
MCMC_leroux <- S.CARleroux(
  formula = Lung_cancer_count ~ offset(log(Pop_count)), 
  data = combined_data_2016,
  family = "poisson",
  W = W_matrix,
  burnin = M.burnin,
  n.sample = M.burnin + M,    # Total iterations
  n.chains = 3,
  n.cores = 3,
  verbose = FALSE
)

# Broad summary of results
print(MCMC_leroux$summary.results)

# Additional results on model fit criteria
MCMC_leroux$modelfit

# LEROUX 加上相关变量
set.seed(444)          # For reproducability
MCMC_leroux_all <- S.CARleroux(
  formula = Lung_cancer_count ~ offset(log(Pop_count)) + Smoking_Rate + Manufacturing + 
    Construction + Emissions + Earnings, 
  data = combined_data_2016,
  family = "poisson",
  W = W_matrix,
  burnin = M.burnin,
  n.sample = M.burnin + M,    # Total iterations
  n.chains = 3,
  n.cores = 3,
  verbose = FALSE
)

# Broad summary of results
print(MCMC_leroux_all$summary.results)

# Additional results on model fit criteria
MCMC_leroux_all$modelfit

# model check 
beta_le_all <- MCMC_leroux_all$samples$beta
gelman_results_le_all <- gelman.diag(beta_le_all)
print(gelman_results_le_all) #1 是1 所以收敛 表明链之间的变异与链内变异相比较小，即模型已经收敛

plot(beta_le_all)
autocorr.plot(beta_le_all)
gelman.plot(beta_le_all)

# 检查 posterior predict model
phi <- MCMC_leroux_all$samples$phi

samples_beta_le_all <- do.call(rbind, lapply(MCMC_leroux_all$samples$beta, function(x) x))
samples_phi_le_all <- do.call(rbind, lapply(MCMC_leroux_all$samples$phi, function(x) x))

Nsamples = 1000
postpred_Y <- array(NA,dim=c(32,Nsamples))
postpred_mean <- array(NA,dim=c(32,Nsamples))

for (i_s in 1:Nsamples) {
  # 从合并后的样本中随机选择一组参数样本
  s = 10*i_s
  intercept_s_le <- samples_beta_le_all[s, 1]
  beta_smoking_s_le <- samples_beta_le_all[s, 2]
  beta_construction_s_le <- samples_beta_le_all[s, 3]
  beta_earnings_s_le <- samples_beta_le_all[s, 4]
  beta_manufacturing_s_le <- samples_beta_le_all[s, 5]
  beta_emissions_s_le <- samples_beta_le_all[s, 6]
  phi_means_s_le <- samples_phi_le_all[s, ]
  
  # 计算预测的均值
  postpred_mean[, i_s] <- with(combined_data_2016, 
                               exp(intercept_s_le + 
                                     beta_smoking_s_le * Smoking_Rate + 
                                     beta_construction_s_le * Construction +
                                     beta_manufacturing_s_le * Manufacturing +
                                     beta_earnings_s_le * Earnings + 
                                     beta_emissions_s_le * Emissions + phi_means_s_le + 
                                     offset(log(Pop_count))
                               ))
  # 生成预测计数
  postpred_Y[, i_s] <- rpois(32, lambda = postpred_mean[, i_s])
}

# 比较预测数据和实际观测数据
hist(colMeans(postpred_Y))
abline(v = mean(combined_data_2016$Lung_cancer_count), col = "red")

hist(apply(postpred_Y, 2, sd))
abline(v = sd(combined_data_2016$Lung_cancer_count), col = "red")

# predict for 2017
# Extract posterior means 
intercept_le <- MCMC_leroux_all$summary.results["(Intercept)", "Mean"]
beta_smoking_le <- MCMC_leroux_all$summary.results["Smoking_Rate", "Mean"]
beta_construction_le <- MCMC_leroux_all$summary.results["Construction", "Mean"]
beta_earnings_le <- MCMC_leroux_all$summary.results["Earnings", "Mean"]
beta_manufacturing_le <- MCMC_leroux_all$summary.results["Manufacturing", "Mean"]
beta_emissions_le <- MCMC_leroux_all$summary.results["Emissions", "Mean"]
phi_means_le <- colMeans(samples_phi_le_all)

# Predicting lung cancer cases for 2017
car_2017$predicted_lung_cancer_le <- with(car_2017, 
                                          exp(intercept_le + 
                                                beta_smoking_le * Smoking_Rate + 
                                                beta_construction_le * Construction +
                                                beta_manufacturing_le * Manufacturing +
                                                beta_earnings_le * Earnings + 
                                                beta_emissions_le * Emissions + phi_means_le + 
                                                offset(log(Pop_count))
                                          )
)


# Add a column to compare the predicted cases to actual cases
car_2017$comparison_le <- with(car_2017, predicted_lung_cancer_le - Lung_cancer_count)
mean(car_2017$comparison_le^2) # 498.3557

# Create a scatter plot 
ggplot(car_2017, aes(x = Lung_cancer_count, y = predicted_lung_cancer_le)) +
  geom_point() +  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Actual Lung Cancer Count", y = "Predicted Lung Cancer Count", title = "Actual vs Predicted Lung Cancer Cases") +
  theme_minimal()  

# 9年去预测2017年

# 先需要搞一个超级大的w matrix
# make a 288x288 matrix filled with 0
W_big_2 <- matrix(0, nrow = 288, ncol = 288)

# Replace diagonal with W_matrix
for (i in 0:8) {
  rows <- (1:32) + i * 32
  cols <- (1:32) + i * 32
  W_big_2[rows, cols] <- W_matrix
}

completed_data_car <-completed_data[order(completed_data$Year), ]
filtered_data_car <- completed_data_car[completed_data_car$Year != "X2017", ]

M.burnin <- 10000       # Number of burn-in iterations (discarded)
M <- 100000               # Number of iterations retained

# BYM
set.seed(444)          # For reproducability
MCMC_big_bym <- S.CARbym(
  formula = Lung_cancer_count ~ offset(log(Pop_count)), 
  data = filtered_data_car,
  family = "poisson",
  W = W_big_2,
  burnin = M.burnin,
  n.sample = M.burnin + M,    # Total iterations
  n.chains = 3,
  n.cores = 3,
  verbose = FALSE
)

# Broad summary of results
print(MCMC_big_bym$summary.results)

# Additional results on model fit criteria
MCMC_big_bym$modelfit

# BYM 包括所有的相关数据
set.seed(444)          # For reproducability
MCMC_big_bym_all <- S.CARbym(
  formula = Lung_cancer_count ~ offset(log(Pop_count)) + Smoking_Rate + Manufacturing + 
    Construction + Emissions + Earnings, 
  data = filtered_data_car,
  family = "poisson",
  W = W_big_2,
  burnin = M.burnin,
  n.sample = M.burnin + M,    # Total iterations
  n.chains = 3,
  n.cores = 3,
  verbose = FALSE
)

# Broad summary of results
print(MCMC_big_bym_all$summary.results)

# Additional results on model fit criteria
MCMC_big_bym_all$modelfit

# model check 
beta_bym_all_big <- MCMC_big_bym_all$samples$beta
gelman_results_bym_all_big <- gelman.diag(beta_bym_all_big)
print(gelman_results_bym_all_big) # 1.01 所以收敛 表明链之间的变异与链内变异相比较小，即模型已经收敛

plot(beta_bym_all_big)
autocorr.plot(beta_bym_all_big)
gelman.plot(beta_bym_all_big)

# 检查 posterior predict model
psi <- MCMC_big_bym_all$samples$psi

samples_beta_bym_all_big <- do.call(rbind, lapply(MCMC_big_bym_all$samples$beta, function(x) x))
samples_psi_bym_all_big <- do.call(rbind, lapply(MCMC_big_bym_all$samples$psi, function(x) x))

Nsamples = 1000
postpred_Y <- array(NA,dim=c(288,Nsamples))
postpred_mean <- array(NA,dim=c(288,Nsamples))

for (i_s in 1:Nsamples) {
  # 从合并后的样本中随机选择一组参数样本
  s = 10*i_s
  intercept_s <- samples_beta_bym_all_big[s, 1]
  beta_smoking_s <- samples_beta_bym_all_big[s, 2]
  beta_construction_s <- samples_beta_bym_all_big[s, 3]
  beta_earnings_s <- samples_beta_bym_all_big[s, 4]
  beta_manufacturing_s <- samples_beta_bym_all_big[s, 5]
  beta_emissions_s <- samples_beta_bym_all_big[s,6]
  psi_means_s <- samples_psi_bym_all_big[s, ]
  
  # 计算预测的均值
  postpred_mean[, i_s] <- with(filtered_data_car, 
                               exp(intercept_s + 
                                     beta_smoking_s * Smoking_Rate + 
                                     beta_construction_s * Construction +
                                     beta_manufacturing_s * Manufacturing +
                                     beta_earnings_s * Earnings + 
                                     beta_emissions_s * Emissions + psi_means_s + 
                                     offset(log(Pop_count))
                               ))
  # 生成预测计数
  postpred_Y[, i_s] <- rpois(288, lambda = postpred_mean[, i_s])
}

# 比较预测数据和实际观测数据
hist(colMeans(postpred_Y)) 
abline(v = mean(filtered_data_car$Lung_cancer_count), col = "red")

hist(apply(postpred_Y, 2, sd))
abline(v = sd(filtered_data_car$Lung_cancer_count), col = "red")

# predict!
car_big_2017 <- filter(completed_data, Year == "X2017")

# Extract posterior means 
intercept_big <- MCMC_big_bym_all$summary.results["(Intercept)", "Mean"]
beta_smoking_big <- MCMC_big_bym_all$summary.results["Smoking_Rate", "Mean"]
beta_construction_big <- MCMC_big_bym_all$summary.results["Construction", "Mean"]
beta_earnings_big <- MCMC_big_bym_all$summary.results["Earnings", "Mean"]
beta_manufacturing_big <- MCMC_big_bym_all$summary.results["Manufacturing", "Mean"]
beta_emissions_big <- MCMC_big_bym_all$summary.results["Emissions", "Mean"]

psi_all_chains <- do.call(rbind, MCMC_big_bym_all$samples$psi)
psi_means <- colMeans(psi_all_chains)
psi_means_big_chain <- matrix(psi_means, nrow = 32, ncol = 9, byrow = FALSE)
psi_row_averages <- rowMeans(psi_means_big_chain)

# Predicting lung cancer cases for 2017
car_big_2017$predicted_lung_cancer_bym_big <- with(car_big_2017, 
                                                   exp(intercept_big + 
                                                         beta_smoking_big * Smoking_Rate + 
                                                         beta_construction_big * Construction +
                                                         beta_manufacturing_big * Manufacturing +
                                                         beta_earnings_big * Earnings + 
                                                         beta_emissions_big * Emissions + psi_row_averages + 
                                                         offset(log(Pop_count))
                                                   )
)


# Add a column to compare the predicted cases to actual cases
car_big_2017$comparison_big_B <- with(car_big_2017, predicted_lung_cancer_bym_big - Lung_cancer_count)
mean(car_big_2017$comparison_big_B^2) # 671.8069

# Create a scatter plot 
ggplot(car_big_2017, aes(x = Lung_cancer_count, y = predicted_lung_cancer_bym_big)) +
  geom_point() +  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Actual Lung Cancer Count", y = "Predicted Lung Cancer Count", title = "Actual vs Predicted Lung Cancer Cases") +
  theme_minimal()  

# LEROUX
set.seed(444)          # For reproducability
MCMC_big_L <- S.CARleroux(
  formula = Lung_cancer_count ~ offset(log(Pop_count)), 
  data = filtered_data_car,
  family = "poisson",
  W = W_big_2,
  burnin = M.burnin,
  n.sample = M.burnin + M,    # Total iterations
  n.chains = 3,
  n.cores = 3,
  verbose = FALSE
)

# Broad summary of results
print(MCMC_big_L$summary.results)

# Additional results on model fit criteria
MCMC_big_L$modelfit

# LEROUX 加上所有变量数据
set.seed(444)          # For reproducability
MCMC_big_le_all <- S.CARleroux(
  formula = Lung_cancer_count ~ offset(log(Pop_count)) + Smoking_Rate + Manufacturing + 
    Construction + Emissions + Earnings, 
  data = filtered_data_car,
  family = "poisson",
  W = W_big_2,
  burnin = M.burnin,
  n.sample = M.burnin + M,    # Total iterations
  n.chains = 3,
  n.cores = 3,
  verbose = FALSE
)

# Broad summary of results
print(MCMC_big_le_all$summary.results)

# Additional results on model fit criteria
MCMC_big_le_all$modelfit

# model check 
beta_le_all_big <- MCMC_big_le_all$samples$beta
gelman_results_le_all_big <- gelman.diag(beta_le_all_big)
print(gelman_results_le_all_big) # 1.01 所以收敛 表明链之间的变异与链内变异相比较小，即模型已经收敛

plot(beta_le_all_big)
autocorr.plot(beta_le_all_big)
gelman.plot(beta_le_all_big)

# 检查 posterior predict model
phi <- MCMC_big_le_all$samples$phi

samples_beta_le_all_big <- do.call(rbind, lapply(MCMC_big_le_all$samples$beta, function(x) x))
samples_phi_bym_all_big <- do.call(rbind, lapply(MCMC_big_le_all$samples$phi, function(x) x))

Nsamples = 1000
postpred_Y <- array(NA,dim=c(288,Nsamples))
postpred_mean <- array(NA,dim=c(288,Nsamples))

for (i_s in 1:Nsamples) {
  # 从合并后的样本中随机选择一组参数样本
  s = 10*i_s
  intercept_s <- samples_beta_le_all_big[s, 1]
  beta_smoking_s <- samples_beta_le_all_big[s, 2]
  beta_construction_s <- samples_beta_le_all_big[s, 3]
  beta_earnings_s <- samples_beta_le_all_big[s, 4]
  beta_manufacturing_s <- samples_beta_le_all_big[s, 5]
  beta_emissions_s <- samples_beta_le_all_big[s,6]
  phi_means_s <- samples_phi_bym_all_big[s, ]
  
  # 计算预测的均值
  postpred_mean[, i_s] <- with(filtered_data_car, 
                               exp(intercept_s + 
                                     beta_smoking_s * Smoking_Rate + 
                                     beta_construction_s * Construction +
                                     beta_manufacturing_s * Manufacturing +
                                     beta_earnings_s * Earnings + 
                                     beta_emissions_s * Emissions + phi_means_s + 
                                     offset(log(Pop_count))
                               ))
  # 生成预测计数
  postpred_Y[, i_s] <- rpois(288, lambda = postpred_mean[, i_s])
}

# 比较预测数据和实际观测数据
hist(colMeans(postpred_Y)) 
abline(v = mean(filtered_data_car$Lung_cancer_count), col = "red")

hist(apply(postpred_Y, 2, sd))
abline(v = sd(filtered_data_car$Lung_cancer_count), col = "red")

# predict!
# Extract posterior means 
intercept_big <- MCMC_big_le_all$summary.results["(Intercept)", "Mean"]
beta_smoking_big <- MCMC_big_le_all$summary.results["Smoking_Rate", "Mean"]
beta_construction_big <- MCMC_big_le_all$summary.results["Construction", "Mean"]
beta_earnings_big <- MCMC_big_le_all$summary.results["Earnings", "Mean"]
beta_manufacturing_big <- MCMC_big_le_all$summary.results["Manufacturing", "Mean"]
beta_emissions_big <- MCMC_big_le_all$summary.results["Emissions", "Mean"]

phi_all_chains <- do.call(rbind, MCMC_big_le_all$samples$phi)
phi_means <- colMeans(phi_all_chains)
phi_means_big_chain <- matrix(phi_means, nrow = 32, ncol = 9, byrow = FALSE)
phi_row_averages <- rowMeans(phi_means_big_chain)

# Predicting lung cancer cases for 2017
car_big_2017$predicted_lung_cancer_le_big <- with(car_big_2017, 
                                                   exp(intercept_big + 
                                                         beta_smoking_big * Smoking_Rate + 
                                                         beta_construction_big * Construction +
                                                         beta_manufacturing_big * Manufacturing +
                                                         beta_earnings_big * Earnings + 
                                                         beta_emissions_big * Emissions + phi_row_averages + 
                                                         offset(log(Pop_count))
                                                   )
)


# Add a column to compare the predicted cases to actual cases
car_big_2017$comparison_big_L <- with(car_big_2017, predicted_lung_cancer_le_big - Lung_cancer_count)
mean(car_big_2017$comparison_big_L^2) # 666.2226

# Create a scatter plot 
ggplot(car_big_2017, aes(x = Lung_cancer_count, y = predicted_lung_cancer_le_big)) +
  geom_point() +  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Actual Lung Cancer Count", y = "Predicted Lung Cancer Count", title = "Actual vs Predicted Lung Cancer Cases") +
  theme_minimal()  

#####################
poisson_model_4 <- glm(Lung_cancer_count ~ Manufacturing + Construction
                       + Emissions + Earnings + Smoking_Rate, offset=offset(log(Pop_count)),
                       data = completed_data,
                       family = "poisson")
summary(poisson_model_4)
plot(poisson_model_4)

# prepare the model for 2008-2016 since we need to predict for 2017
filtered_data <- completed_data %>%
  filter(completed_data$Year >= "X2008" & completed_data$Year <= "X2016")

# final model for glm in poisson
poisson_model_4 <- glm(Lung_cancer_count ~ Manufacturing + Construction
                       + Emissions + Earnings + Smoking_Rate, offset=offset(log(Pop_count)),
                       data = filtered_data,
                       family = "poisson")
summary(poisson_model_4)
plot(poisson_model_4)

# predict for 2017
y_2017 <- completed_data %>% filter(Year == "X2017")

# Predicting lung cancer cases in 2017
predicted_counts <- predict(poisson_model_4, y_2017, type = "response")

# Add the predicted counts to your y_2017 dataset
y_2017$predicted_lung_cancer = predicted_counts

# Compute absolute value
y_2017$absolute_error <- abs(y_2017$Lung_cancer_count - y_2017$predicted_lung_cancer)

# View the results
View(y_2017)

# calculate the MSE
squared_errors <- (y_2017$Lung_cancer_count - y_2017$predicted_lung_cancer)^2
mse <- mean(squared_errors)
print(mse) ## 1311.372

########### poission predict end ---------------------------------------------------------------

library(knitr)

# Create a data frame
models_aic <- data.frame(
  Model = c("Leroux(multi-year)", "BYM(multi-year)", "Leroux(1 year)","BYM(1 year)","Poisson"),
  MSE = c(mean(car_big_2017$comparison_big_L^2), mean(car_big_2017$comparison_big_B^2) , 
          mean(car_2017$comparison_le^2), mean(car_2017$comparison_b^2), mse)
)
# Create a table using knitr::kable()
kable(models_aic, caption = "AIC Values for Models", format = "markdown")

############


