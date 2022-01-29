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
