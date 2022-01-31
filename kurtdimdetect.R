# p = 5
# 
# q_true_skew = 2
# q_true_kurt = 2
# t_df = 10
# 
# skew_or_kurt_or_both = 2
# 
# n = 200
# 
# alpha_value = 5
# 
# Omega_1 = matrix(0.5, nrow = q_true_skew, ncol = q_true_skew) + diag(0.5, nrow = q_true_skew)
# Omega_2 = matrix(0.5, nrow = (p - q_true_skew), ncol = (p - q_true_skew)) + diag(0.5, nrow = (p - q_true_skew))
# 
# alpha_vector = rep(alpha_value, q_true_skew)

Omega_3 = matrix(0.5, nrow = q_true_kurt, ncol = q_true_kurt) + diag(0.5, nrow = q_true_kurt)
Omega_4 = matrix(0.5, nrow = (p - q_true_kurt), ncol = (p - q_true_kurt)) + diag(0.5, nrow = (p - q_true_kurt))

allsamplenum = 2^p - 1
allsampleindices = matrix(nrow = 0, ncol = p)
for (q in 1:p){
  sampleindices_q = t(apply(combn(p, q), 2, function(x) (1:p %in% x)))
  allsampleindices = rbind(allsampleindices, sampleindices_q)
}

q_vector = rowSums(allsampleindices)
