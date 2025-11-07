library(data.table)
library(dplyr)
library(stringr)

# 随机数种子
n <- sample(1:1000000, 1)
set.seed(n)

# 同态加密文件
args <- commandArgs(trailingOnly = TRUE)

# 解析命令行参数
options <- list()
for (arg in args) {
  kvp <- strsplit(arg, "=")[[1]]
  key <- sub("--", "", kvp[1])
  value <- kvp[2]
  options[[key]] <- value
}
print(args)

# 获取信息
work_path <- options[['workPath']]   # 获取工作路径
result_file <- options[['resultPath']]
matrix_file <- options[['matrixPath']]
samplesize <- options[['SampleSize']]

# 设置工作路径
setwd(work_path)

# 读取加密结果矩阵
orthogonal_result <- tryCatch({
  fread(result_file)
}, error = function(e) {
  stop(e)
})
colname_result <- colnames(orthogonal_result)
orthogonal_result <- orthogonal_result[orthogonal_result$V1 == 3 & orthogonal_result$V4 == 1 ,]
strat = nrow(orthogonal_result) - as.numeric(samplesize) + 1
orthogonal_result <- orthogonal_result[c(strat:nrow(orthogonal_result)),]

# 读取正交函数矩阵
orthogonal_matrix <- tryCatch({
  readRDS(matrix_file)
}, error = function(e) {
  stop(e)
})

# 解密
tryCatch({
  decryption_result <- t(orthogonal_matrix) %*% orthogonal_result$V8
  orthogonal_result$V8 <- decryption_result
  write.table(orthogonal_result,"decryptionResult.txt",row.names = F,quote = F,col.names = T)
}, error = function(e) {
  stop(e)
})
print("well done!")

# 输出结果信息到控制台
cat("The run has successfully finished. Results saved in decryptionResult.txt.\n")
