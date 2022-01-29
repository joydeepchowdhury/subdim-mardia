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

####
signscale_vector = c(1, -0.5, 1, -0.5, 1, -0.5, 1, -0.5, 1, -0.5)
alpha_vector = alpha_vector * signscale_vector[1:length(alpha_vector)]
####

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
      