# p = 5
# 
# q_true_skew = 0
# q_true_kurt = 2
# t_df = 10
# 
# skew_or_kurt_or_both = 1
# 
# n = 200
# 
# alpha_value = 5

Omega_1 = matrix(0.5, nrow = q_true_skew, ncol = q_true_skew) + diag(0.5, nrow = q_true_skew)
Omega_2 = matrix(0.5, nrow = (p - q_true_skew), ncol = (p - q_true_skew)) + diag(0.5, nrow = (p - q_true_skew))

alpha_vector = rep(alpha_value, q_true_skew)

# ####
# signscale_vector = c(1, -0.5, 1, -0.5, 1, -0.5, 1, -0.5, 1, -0.5)
# alpha_vector = alpha_vector * signscale_vector[1:length(alpha_vector)]
# ####

Omega_3 = matrix(0.5, nrow = q_true_kurt, ncol = q_true_kurt) + diag(0.5, nrow = q_true_kurt)
Omega_4 = matrix(0.5, nrow = (p - q_true_kurt), ncol = (p - q_true_kurt)) + diag(0.5, nrow = (p - q_true_kurt))

allsamplenum = 2^p - 1
allsampleindices = matrix(nrow = 0, ncol = p)
for (q in 1:p){
  sampleindices_q = t(apply(combn(p, q), 2, function(x) (1:p %in% x)))
  allsampleindices = rbind(allsampleindices, sampleindices_q)
}

q_vector = rowSums(allsampleindices)

if (skew_or_kurt_or_both %in% c(1, 3)){
  if (q_true_skew >= 0 && q_true_skew <= p){
    sampleindices_true = c(rep(TRUE, q_true_skew), rep(FALSE, (p - q_true_skew)))
  }else{
    stop('ERROR in q_true_skew!!')
  }
}else if (skew_or_kurt_or_both == 2){
  if (q_true_kurt >= 0 && q_true_kurt <= p){
    sampleindices_true = c(rep(TRUE, q_true_kurt), rep(FALSE, (p - q_true_kurt)))
  }else{
    stop('ERROR in q_true_kurt!!')
  }
}else{
  stop('ERROR in skew_or_kurt_or_both!!')
}

index_true = 0
count_index = 0
for (i in 1:allsamplenum){
  if (all(sampleindices_true == allsampleindices[i,])){
    index_true = i
    count_index = count_index + 1
  }
}
if (count_index > 1)
  stop('ERROR!!!!')

num_repl = 1000

set.seed(seed = NULL)
Seed_vector = sample(10000000, num_repl)

require(foreach)
require(doParallel)

packagelist = c('MASS', 'sn')

ncores = detectCores()
cl = makeCluster(min(ncores, 120))
registerDoParallel(cl)

parallel_outputs = foreach(index_replicate = 1:num_repl, .combine = rbind, .packages = packagelist) %dopar%
  {
    seed = Seed_vector[index_replicate]
    
    if(!is.integer(seed))
      stop(paste('ERROR!!!!', as.character(seed)))
    set.seed(seed, kind = 'default', normal.kind = 'default')
    
    if (skew_or_kurt_or_both == 1){
      if (q_true_skew == 0){
        symmetric_components = MASS::mvrnorm(n = n, mu = rep(0, p - q_true_skew), Sigma = Omega_2)
        Data = symmetric_components
      }else if (q_true_skew == 1){
        skewed_components = sn::rsn(n = n, xi = rep(0,length(alpha_vector)), omega = Omega_1, alpha = alpha_vector)
        symmetric_components = MASS::mvrnorm(n = n, mu = rep(0, p - q_true_skew), Sigma = Omega_2)
        Data = cbind(skewed_components, symmetric_components)
      }else if (q_true_skew > 1 && q_true_skew < p){
        skewed_components = sn::rmsn(n = n, xi = rep(0,length(alpha_vector)), Omega = Omega_1, alpha = alpha_vector)
        symmetric_components = MASS::mvrnorm(n = n, mu = rep(0, p - q_true_skew), Sigma = Omega_2)
        Data = cbind(skewed_components, symmetric_components)
      }else if (q_true_skew == p){
        skewed_components = sn::rmsn(n = n, xi = rep(0,length(alpha_vector)), Omega = Omega_1, alpha = alpha_vector)
        Data = skewed_components
      }else{
        stop('ERROR in q_true_skew!!')
      }
    }else if (skew_or_kurt_or_both == 2){
      if (q_true_kurt == 0){
        normal_components = MASS::mvrnorm(n = n, mu = rep(0, p - q_true_kurt), Sigma = Omega_4)
        Data = normal_components
      }else if (q_true_kurt == 1){
        heavytail_components = rnorm(n, mean = 0, sd = 1) / sqrt(rchisq(n, t_df) / t_df)
        normal_components = MASS::mvrnorm(n = n, mu = rep(0, p - q_true_kurt), Sigma = Omega_4)
        Data = cbind(heavytail_components, normal_components)
      }else if (q_true_kurt > 1 && q_true_kurt < p){
        heavytail_components = MASS::mvrnorm(n = n, mu = rep(0, q_true_kurt), Sigma = Omega_3) /
          matrix(sqrt(rchisq(n, t_df) / t_df), nrow = n, ncol = q_true_kurt, byrow = FALSE)
        normal_components = MASS::mvrnorm(n = n, mu = rep(0, p - q_true_kurt), Sigma = Omega_4)
        Data = cbind(heavytail_components, normal_components)
      }else if (q_true_kurt == p){
        heavytail_components = MASS::mvrnorm(n = n, mu = rep(0, q_true_kurt), Sigma = Omega_3) /
          matrix(sqrt(rchisq(n, t_df) / t_df), nrow = n, ncol = q_true_kurt, byrow = FALSE)
        Data = heavytail_components
      }else{
        stop('ERROR in q_true_kurt!!')
      }
    }else if (skew_or_kurt_or_both == 3){
      if (q_true_skew == 0){
        symmetric_components = MASS::mvrnorm(n = n, mu = rep(0, p - q_true_skew), Sigma = Omega_2)
        Data_temp = symmetric_components
      }else if (q_true_skew == 1){
        skewed_components = sn::rsn(n = n, xi = rep(0,length(alpha_vector)), omega = Omega_1, alpha = alpha_vector)
        symmetric_components = MASS::mvrnorm(n = n, mu = rep(0, p - q_true_skew), Sigma = Omega_2)
        Data_temp = cbind(skewed_components, symmetric_components)
      }else if (q_true_skew > 1 && q_true_skew < p){
        skewed_components = sn::rmsn(n = n, xi = rep(0,length(alpha_vector)), Omega = Omega_1, alpha = alpha_vector)
        symmetric_components = MASS::mvrnorm(n = n, mu = rep(0, p - q_true_skew), Sigma = Omega_2)
        Data_temp = cbind(skewed_components, symmetric_components)
      }else if (q_true_skew == p){
        skewed_components = sn::rmsn(n = n, xi = rep(0,length(alpha_vector)), Omega = Omega_1, alpha = alpha_vector)
        Data_temp = skewed_components
      }else{
        stop('ERROR in q_true_skew!!')
      }
      
      if (q_true_kurt == 0){
        normal_components = Data_temp
        Data = normal_components
      }else if (q_true_kurt == 1){
        heavytail_components = Data_temp[, 1] / sqrt(rchisq(n, t_df) / t_df)
        normal_components = Data_temp[, (q_true_kurt + 1):p]
        Data = cbind(heavytail_components, normal_components)
      }else if (q_true_kurt > 1 && q_true_kurt < p){
        heavytail_components = Data_temp[, 1:q_true_kurt] /
          matrix(sqrt(rchisq(n, t_df) / t_df), nrow = n, ncol = q_true_kurt, byrow = FALSE)
        normal_components = Data_temp[, (q_true_kurt + 1):p]
        Data = cbind(heavytail_components, normal_components)
      }else if (q_true_kurt == p){
        heavytail_components = Data_temp[, 1:q_true_kurt] /
          matrix(sqrt(rchisq(n, t_df) / t_df), nrow = n, ncol = q_true_kurt, byrow = FALSE)
        Data = heavytail_components
      }else{
        stop('ERROR in q_true_kurt!!')
      }
    }else{
      stop('ERROR in skew_or_kurt_or_both!!')
    }
    
    m1 = mat.or.vec(allsamplenum, 1)
    mean_m1 = mat.or.vec(allsamplenum, 1)
    sd_m1 = mat.or.vec(allsamplenum, 1)
    for (i in 1:allsamplenum){
      Data_i = Data[, allsampleindices[i,]]
      if (is.vector(Data_i))
        Data_i = matrix(Data_i, nrow = length(Data_i), ncol = 1)
      
      q = q_vector[i]
      
      Data_i_centered = Data_i - matrix(colMeans(Data_i), nrow = n, ncol = q, byrow = TRUE)
      
      Sigma_hat_i_inverse = solve(cov(Data_i))
      
      temp_matrix = Data_i_centered %*% solve(cov(Data_i), t(Data_i_centered))
      
      b1 = sum((temp_matrix)^3) / (n^2)
      
      m1[i] = b1
      
      m4 = mean((diag(temp_matrix))^2)
      m6 = mean((diag(temp_matrix))^3)
      
      alpha_1 = (3/q) * ((m6 / (q + 2)) - (2 * m4) + (q * (q + 2)))
      alpha_2 = (6 * m6) / (q * (q + 2) * (q + 4))
      
      multiplicity_1 = q
      multiplicity_2 = (q * (q - 1) * (q + 4)) / 6
      
      Var_b1 = ((2 * multiplicity_1 * alpha_1^2) + (2 * multiplicity_2 * alpha_2^2)) / (n^2)
      
      mean_m1[i] = ((multiplicity_1 * alpha_1) + (multiplicity_2 * alpha_2)) / n
      
      sd_m1[i] = sqrt(Var_b1)
    }
    
    centered_scaled_m1 = abs((m1 - mean_m1) / sd_m1)
    
    max_centered_scaled_m1 = max(centered_scaled_m1)
    
    indexmax_centered_scaled_m1 = which(centered_scaled_m1 == max_centered_scaled_m1)
    
    if (length(indexmax_centered_scaled_m1) > 1)
      stop('Error!!!!')
    
    c(indexmax_centered_scaled_m1)
  }

stopCluster(cl)

indicesmax_centered_scaled_m1 = as.vector(parallel_outputs)

allsample_tally_centered_scaled_m1 = mat.or.vec(allsamplenum, 1)
for (i in 1:allsamplenum){
  allsample_tally_centered_scaled_m1[i] = mean(i == indicesmax_centered_scaled_m1)
}

print(cbind(allsample_tally_centered_scaled_m1, q_vector))

barplot(allsample_tally_centered_scaled_m1, names.arg = 1:allsamplenum, ylim = c(0,1),
        xlab = 'sample indices', ylab = 'proportion',
        main = paste('true index = ', index_true, ', n = ', n, ', \nalpha = (',
                     paste(alpha_vector, collapse = ', '), ')', sep = ''))
