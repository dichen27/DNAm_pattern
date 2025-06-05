rm(list=ls())
setwd('/Users/dichen/Documents/Methylation_Running/')

library(R.matlab) 
library(cocor)

##################### 
############## 
## Z-test
r1 <- 0.8623234  # 第一个相关系数
r2 <- 0.7879541  # 第二个相关系数
n1 <- 175  # 第一个样本大小
n2 <- 190  # 第二个样本大小

# Fisher z-transformation
z1 <- 0.5 * log((1 + r1) / (1 - r1))
z2 <- 0.5 * log((1 + r2) / (1 - r2))

# 计算 z 值之差和标准误差
difference <- z1 - z2
SE <- sqrt((1/(n1-3)) + (1/(n2-3)))

# 计算 z 统计量
z_statistic <- difference / SE

# 计算 p 值
p_value <- 2 * (1 - pnorm(abs(z_statistic)))

# 输出结果
cat("Z statistic:", z_statistic, "\nP-value:", p_value)

####



##################### 
# finished
################ 

