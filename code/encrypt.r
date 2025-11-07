# 加载工具
library(data.table)
library(Rcpp)
library(parallel)
library(bigstatsr)
library(bigsnpr)
library(bigreadr)
library(parallelly)
library(purrr)
library(dplyr)
library(foreach)
library(doParallel)

# 同态加密文件
args <- commandArgs(trailingOnly = TRUE)

# 随机数种子
num <- sample(1:1000000, 1)
set.seed(num)

# 解析参数
options <- list()
for (arg in args) {
  kvp <- strsplit(arg, "=")[[1]]
  key <- sub("--", "", kvp[1])
  value <- kvp[2]
  options[[key]] <- value
}


# 获取信息
work_path <- options[['workPath']]  # 获取工作路径
phePath <- options[['phenotypePath']]
genoPath <- options[['genotypePath']]
ped_file <- options[['pedigreePath']]
phneoPos <- options[['selectedTraits']]
dcovar <- options[['selectedFixedFactors']]
qcovar <- options[['selectedCovariates']]
RandomOtherEffects <- options[['selectedRandomOtherEffects']]
RandomGeneticEffects <- options[["selectedRandomGeneticEffects"]]
plink <- options[['plink']]
gcta <- options[['gcta']]
hiblup <- options[['hiblup']]
# 设置工作路径
setwd(work_path)

# 功能：将给定方阵 M 分解为 B，使得 M ≈ B %*% t(B)
decompose_to_B <- function(M, ridge = 0, ev_tol = 1e-12) {
  # 轻度对称化（避免数值非对称）
  M <- (M + t(M)) / 2
  
  # 可选对角加脊（改善数值正定性；对非常接近半正定的矩阵很有用）
  if (ridge > 0) {
    M <- M + diag(ridge, nrow(M))
  }
  
  # 1) 首选 Cholesky（允许主元）
  ch_try <- tryCatch(chol(M, pivot = TRUE), error = function(e) NULL)
  if (!is.null(ch_try)) {
    R   <- ch_try                     # A[piv, piv] = Rᵀ R
    piv <- attr(R, "pivot")
    Bp  <- t(R)                       # 先在置换后的顺序下得到 Bp
    # 还原到原始顺序：在对应行放回
    B <- matrix(0, nrow(M), ncol(M))
    B[piv, ] <- Bp
    return(B)
  }
  
  # 2) 回退：特征分解并裁剪负特征值为 0（最近 PSD）
  ee   <- eigen(M, symmetric = TRUE)
  vals <- pmax(ee$values, 0)          # 将微小负值截断为 0
  V    <- ee$vectors
  B    <- V %*% diag(sqrt(vals), nrow = length(vals))
  return(B)
}

# 初始化 Success 变量
Success <- 1

phe <- tryCatch({
  fread(phePath)
}, error = function(e) {
  Success <<- 2
  stop(e)
})
phe <- as.data.frame(phe)

# 提取表头内容并在末尾添加null，如果没有选中的表头则为null
header <- names(phe)
header <- c(header, NA)

system(paste0(plink," --vcf ",genoPath," --make-bed --out genotype"))
# 按照fam文件调整表型文件
fam <- fread("genotype.fam")
fam_ids <- fam$V2

# 保留在fam文件中的个体
phe_filtered <- phe[phe[,as.numeric(RandomGeneticEffects)] %in% fam_ids, ]

# 按照fam文件中的顺序进行排序
phe <- phe_filtered[match(fam_ids, phe_filtered[,as.numeric(RandomGeneticEffects)]), ]

# 检查是否有重复的表型值
duplicated_values <- any(duplicated(phe[,as.numeric(RandomGeneticEffects)]))

# 根据是否有重复的值决定是否添加 --rand 1 参数
rand_param <- ifelse(duplicated_values, "--rand 1", "")

# 检查 selectedFixedFactors 和 selectedCovariates 是否为空
dcovar_param <- ifelse(!is.null(dcovar) && dcovar != "", paste0(" --dcovar ", dcovar), "")
qcovar_param <- ifelse(!is.null(qcovar) && qcovar != "", paste0(" --qcovar ", qcovar), "")

# HIBLUP 命令
hiblup_command <- paste0(
  hiblup, " --single-trait",
  " --pheno ", phePath,
  " --pheno-pos 2",  # 表型位置（根据您的说明）
  dcovar_param,      # 根据条件添加 --dcovar
  qcovar_param,      # 根据条件添加 --qcovar
  rand_param,        # 根据是否有重复值添加 --rand 1
  " --bfile genotype",
  " --threads 32",
  " --out genotype"
)

# 执行 HIBLUP 命令
system(hiblup_command)
hiblup_result <- fread("genotype.rand")
hiblup_result <- hiblup_result[!is.na(hiblup_result$residuals)]
hiblup_result$ADphe <- hiblup_result$GA + hiblup_result$residuals


# 生成加密器
n <- dim(phe)[1]
A <- matrix(rnorm(n^2), nrow = n)  # 生成随机矩阵
qr_decomp <- qr(A)                 # QR分解
encrypter <- qr.Q(qr_decomp)

HE_phe <- encrypter %*% as.matrix(hiblup_result$ADphe)

phe[[colnames(phe)[as.numeric(phneoPos)]]] <- HE_phe

# 保存加密数据
write.table(phe,paste0("HE_phenotype.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE)
saveRDS(encrypter, file = paste0("encrypter.rds"))

# 检查是否提供了 pedigree 文件
if (file.exists(ped_file)) {
  ids <- fread("genotype.fam")
  ids <- unique(ids$V1)
  # P阵
  system(paste0(hiblup," --make-xrm",
                " --pedigree ",ped_file,
                " --add",
                " --write-txt", 
                " --out ped"))
  
  Pmat <- data.table::fread("ped.txt")
  id <- fread("ped.PA.id",header = F)
  Pmat <- as.matrix(Pmat)
  dimnames(Pmat) <- list(id$V1, id$V1)
  Pmat <- Pmat[ids, ids, drop = FALSE]
  
  
  # H阵
  system(paste0(hiblup," --make-xrm",
                " --pedigree ",ped_file,
                " --bfile genotype",
                " --add",
                " --write-txt", 
                " --out Hmatrix"))
  
  Hmat <- data.table::fread("Hmatrix.txt")
  id <- fread("Hmatrix.HA.id",header = F)
  Hmat <- as.matrix(Hmat)
  dimnames(Hmat) <- list(id$V1, id$V1)
  Hmat <- Hmat[ids, ids, drop = FALSE]
  
  # G阵
  system(paste0(hiblup," --make-xrm",
                " --bfile genotype",
                " --add",
                " --write-txt", 
                " --out Gmatrix"))
  
  Gmat <- data.table::fread("Gmatrix.txt")
  id <- fread("Gmatrix.GA.id",header = F)
  Gmat <- as.matrix(Gmat)
  dimnames(Gmat) <- list(id$V1, id$V1)
  Gmat <- Gmat[ids, ids, drop = FALSE]
  
  # Step 1: 读取 .bed 文件
  snp_readBed("genotype.bed")  # 生成 rds 文件
  bigSNP <- snp_attach("genotype.rds")  # 加载 RDS 文件
  G <- bigSNP$genotypes  # 提取基因型矩阵 (FBM)
  
  # Step 2: 创建标准化矩阵
  # 创建一个空的标准化矩阵，其维度与 G 相同
  G_std <- FBM(nrow = nrow(G), ncol = ncol(G))  # 维度 2000 × 409512
  
  print("G_std completed")
  
  # Step 3: 使用 big_apply 对矩阵进行标准化
  big_apply(G, a.FUN = function(G, ind, G_std) {
    # 对每列（SNP）进行标准化
    p <- colMeans(G[, ind], na.rm = TRUE) / 2  # 计算等位基因频率 p
    sd <- sqrt(2 * p * (1 - p))  # 计算标准差
    G_std[, ind] <- scale(G[, ind], center = 2 * p, scale = sd)  # 标准化
    NULL  # big_apply 要求返回 NULL
  }, a.combine = "c", G_std = G_std, ind = cols_along(G), ncores = 1)
  
  print("big_apply completed")
  
  block_matrix_multiply_parallel <- function(A, B, block_size, cores = 1) {
    n <- nrow(A)
    m <- ncol(A)
    p <- ncol(B)
    
    cl <- parallel::makeCluster(cores)
    
    # 仅导出必需的变量（避免导出整个A）
    parallel::clusterExport(cl, varlist = c("B", "block_size", "m", "p"), envir = environment())
    
    # 准备任务列表（避免主进程内存溢出）
    start_rows <- seq(1, n, by = block_size)
    tasks <- lapply(start_rows, function(i) {
      list(
        i = i,
        A_block = A[i:min(i + block_size - 1, n), , drop = FALSE]
      )
    })
    
    results <- parallel::parLapply(cl, tasks, function(task) {
      i <- task$i
      A_block <- task$A_block
      len_row <- nrow(A_block)
      block_C <- matrix(0, len_row, p)  # 仅分配当前块的内存
      
      # 三重循环计算局部块
      for (j in seq(1, p, by = block_size)) {
        col_end <- min(j + block_size - 1, p)
        col_range <- j:col_end
        
        for (k in seq(1, m, by = block_size)) {
          mid_end <- min(k + block_size - 1, m)
          mid_range <- k:mid_end
          
          # 矩阵乘法累加
          block_C[, col_range] <- block_C[, col_range] +
            A_block[, mid_range, drop = FALSE] %*%
            B[mid_range, col_range, drop = FALSE]
        }
      }
      list(i = i, len = len_row, block_C = block_C)
    })
    
    parallel::stopCluster(cl)
    
    # 预分配结果矩阵
    C <- matrix(0, n, p)
    
    # 组合结果
    for (res in results) {
      row_range <- res$i:(res$i + res$len - 1)
      C[row_range, ] <- res$block_C
    }
    
    return(C)
  }
  
  # 调用示例
  block_size <- 500
  gstd <- block_matrix_multiply_parallel(encrypter, G_std, block_size)
  saveRDS(gstd, file = paste0("scaled_genomic_matrix.rds"))
  
  # 分解矩阵
  B_P <- decompose_to_B(Pmat)
  B_H <- decompose_to_B(Hmat)
  B_P_HE <- encrypter %*% B_P
  B_H_HE <- encrypter %*% B_H
  # 保存H阵A阵
  saveRDS(B_P_HE, file = paste0("P_matrix.rds"))
  saveRDS(B_H_HE, file = paste0("H_matrix.rds"))
  
  B_G_HE <- encrypter %*% Gmat %*% t(encrypter)
  ginv <- solve(B_G_HE)
  # 构建GINV
  max_num <- dim(ginv)[1]
  # 第一列：生成 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, ...
  column1 <- rep(1:max_num, times = 1:max_num)
  # 第二列：根据第一列生成递增的子序列，跳过0
  column2 <- unlist(lapply(1:max_num, function(x) 1:x))
  # 组合成数据框
  ginv_dat<-data.frame(c(0,column1),c(0,column2),
                       c(determinant(ginv)$modulus,ginv[upper.tri(ginv,diag = T)]))
  write.table(ginv_dat,file="GINV",row.names = F,col.names = F,quote=F)
  system("bgzip GINV")  
}else{
  ids <- fread("genotype.fam")
  ids <- unique(ids$V1)
  # G阵
  system(paste0(hiblup," --make-xrm",
                " --bfile genotype",
                " --add",
                " --write-txt", 
                " --out Gmatrix"))
  
  Gmat <- data.table::fread("Gmatrix.txt")
  id <- fread("Gmatrix.GA.id",header = F)
  Gmat <- as.matrix(Gmat)
  dimnames(Gmat) <- list(id$V1, id$V1)
  Gmat <- Gmat[ids, ids, drop = FALSE]
  
  # 生成加密器
  n <- length(ids)
  A <- matrix(rnorm(n^2), nrow = n)  # 生成随机矩阵
  qr_decomp <- qr(A)                 # QR分解
  encrypter <- qr.Q(qr_decomp)
  
  # Step 1: 读取 .bed 文件
  snp_readBed("genotype.bed")  # 生成 rds 文件
  bigSNP <- snp_attach("genotype.rds")  # 加载 RDS 文件
  G <- bigSNP$genotypes  # 提取基因型矩阵 (FBM)
  
  # Step 2: 创建标准化矩阵
  # 创建一个空的标准化矩阵，其维度与 G 相同
  G_std <- FBM(nrow = nrow(G), ncol = ncol(G))  # 维度 2000 × 409512
  
  print("G_std completed")
  
  # Step 3: 使用 big_apply 对矩阵进行标准化
  big_apply(G, a.FUN = function(G, ind, G_std) {
    # 对每列（SNP）进行标准化
    p <- colMeans(G[, ind], na.rm = TRUE) / 2  # 计算等位基因频率 p
    sd <- sqrt(2 * p * (1 - p))  # 计算标准差
    G_std[, ind] <- scale(G[, ind], center = 2 * p, scale = sd)  # 标准化
    NULL  # big_apply 要求返回 NULL
  }, a.combine = "c", G_std = G_std, ind = cols_along(G), ncores = 1)
  
  print("big_apply completed")
  
  block_matrix_multiply_parallel <- function(A, B, block_size, cores = 1) {
    n <- nrow(A)
    m <- ncol(A)
    p <- ncol(B)
    
    cl <- parallel::makeCluster(cores)
    
    # 仅导出必需的变量（避免导出整个A）
    parallel::clusterExport(cl, varlist = c("B", "block_size", "m", "p"), envir = environment())
    
    # 准备任务列表（避免主进程内存溢出）
    start_rows <- seq(1, n, by = block_size)
    tasks <- lapply(start_rows, function(i) {
      list(
        i = i,
        A_block = A[i:min(i + block_size - 1, n), , drop = FALSE]
      )
    })
    
    results <- parallel::parLapply(cl, tasks, function(task) {
      i <- task$i
      A_block <- task$A_block
      len_row <- nrow(A_block)
      block_C <- matrix(0, len_row, p)  # 仅分配当前块的内存
      
      # 三重循环计算局部块
      for (j in seq(1, p, by = block_size)) {
        col_end <- min(j + block_size - 1, p)
        col_range <- j:col_end
        
        for (k in seq(1, m, by = block_size)) {
          mid_end <- min(k + block_size - 1, m)
          mid_range <- k:mid_end
          
          # 矩阵乘法累加
          block_C[, col_range] <- block_C[, col_range] +
            A_block[, mid_range, drop = FALSE] %*%
            B[mid_range, col_range, drop = FALSE]
        }
      }
      list(i = i, len = len_row, block_C = block_C)
    })
    
    parallel::stopCluster(cl)
    
    # 预分配结果矩阵
    C <- matrix(0, n, p)
    
    # 组合结果
    for (res in results) {
      row_range <- res$i:(res$i + res$len - 1)
      C[row_range, ] <- res$block_C
    }
    
    return(C)
  }
  
  # 调用示例
  block_size <- 500
  gstd <- block_matrix_multiply_parallel(encrypter, G_std, block_size)
  saveRDS(gstd, file = paste0("scaled_genomic_matrix.rds"))
  Gmat <- Gmat + diag(0.01,nrow(Gmat),nrow(Gmat))
  B_G_HE <- encrypter %*% Gmat %*% t(encrypter)
  ginv <- solve(B_G_HE)
  
  fam <- fread("genotype.fam")
  ord <- 1:nrow(fam)
  rnd_gen_name <- data.frame(Name=fam[[1]],Code=ord)
  fam_id <- rnd_gen_name[match(fam$V1,rnd_gen_name$Name),'Code']
  ginv_dat <- data.frame(c(0,fam_id[grm$V1]),c(0,fam_id[grm$V2]),
                         c(determinant(ginv)$modulus,ginv[upper.tri(ginv,diag = T)]))
  write.table(ginv_dat,file="GINV",row.names = F,col.names = F,quote=F)
  system("bgzip GINV")  
}