# remove(list = ls())
# 
# p = 5
# 
# t_df = 10
# 
# skew_or_kurt_or_both = 2
# 
# n = 200
# 
# par(mfrow = c(3, 3))
# for (q_true_kurt in 0:2){
#   source("~/kurtdimdetect.R")
# }
# par(mfrow = c(3, 3))
# for (q_true_kurt in 3:p){
#   source("~/kurtdimdetect.R")
# }

remove(list = ls())

p = 5

skew_or_kurt_or_both = 2

n = 500

par(mfrow = c(3, 3))
for (q_true_kurt in 1:3){
  for (t_df in c(1, 5, 10)){
    source("~/kurtdimdetect1.R")
  }
}