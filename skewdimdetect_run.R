# remove(list = ls())
# 
# p = 5
# 
# q_true_kurt = 2
# t_df = 10
# 
# skew_or_kurt_or_both = 1
# 
# n = 200
# 
# alpha_value = 5
# 
# par(mfrow = c(3, 3))
# for (q_true_skew in 0:2){
#   source("~/skewdimdetect.R")
# }
# par(mfrow = c(3, 3))
# for (q_true_skew in 3:p){
#   source("~/skewdimdetect.R")
# }

remove(list = ls())

p = 5

q_true_kurt = 2
t_df = 10

skew_or_kurt_or_both = 1

n = 500

par(mfrow = c(3, 3))
for (q_true_skew in 1:3){
  for (alpha_value in c(1, 3, 5)){
    source("~/skewdimdetect1.R")
  }
}
par(mfrow = c(3, 3))
for (q_true_skew in 4:p){
  for (alpha_value in c(1, 3, 5)){
    source("~/skewdimdetect1.R")
  }
}
q_true_skew = 0
source("~/skewdimdetect1.R")

