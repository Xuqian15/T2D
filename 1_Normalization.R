#!/usr/bin/env Rscript

if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T)
}

if (TRUE){
  option_list = list(
    make_option(c("-W", "--workplace"), type="character", default="/data3/Group5/xuqian/TD-fg/multi-kingdom/xMarkerFinder/data/input_files/",
                help="Input workplace [default %default]"),
    make_option(c("-p", "--profile"), type="character", default="T2DM_merged_Bacteria.csv",
                help="feature abundance profile[default %default]"),
    make_option(c("-m", "--method"), type="character", default="REL",
                help="normalization method[default %default]"),
    make_option(c("-o", "--output"), type="character", default="Bacteria",
                help="output file prefix [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))

  # 显示输入输出确认是否正确
  print(paste("The workplace is ", opts$workplace,  sep = ""))
  print(paste("The microbial profile is ", opts$profile,  sep = ""))
  print(paste("The normalization method is ", opts$method, sep = ""))
}


package_list = c("dplyr","compositions","edgeR")
###CHECK PACKAGES
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

#import data and preprocess
feat_abd <- read.csv(file = paste(opts$workplace,opts$profile,sep=''),sep = ',',header =  TRUE, stringsAsFactors = FALSE, check.names = FALSE)
# 删除第一列（Classification列）中值为空的行
feat_abd <- as.data.frame(feat_abd[!is.na(feat_abd$Classification), ])
# 将第一列作为新的行名
rownames(feat_abd) <- feat_abd$Classification
# 去除第一列
feat_abd <- feat_abd[, -1]
feat_abd[is.na(feat_abd)]<-0
feat_abd <- t(feat_abd)
write.csv(feat_abd,file = paste(opts$workplace,"t_T2DM_merged_Bacteria.csv",sep=''))
 
norm_method <- opts$method
#normalization
if(norm_method == "REL"){
  data_norm <- t(apply(feat_abd,1, function(x){x/sum(x)}))
} else if (norm_method =="AST"){
  data_norm <- t(apply(feat_abd,1, function(x){x/sum(x)}))
  data_norm <- t(apply(data_norm, 1, function(x) asin(sqrt(x))))
} else if (norm_method =="CLR"){
  data_norm <- t(clr(t(feat_abd)))
} else if (norm_method =="TMM"){
  data_norm = t(feat_abd)
  dge <- DGEList(counts = data_norm)
  dge <- calcNormFactors(dge,lib.size=NULL,method = "TMM")
  data_norm <- t(data_norm)/(dge$samples$norm.factors)
}


write.csv(data_norm,file = paste(opts$workplace,opts$output,"_normalized_abundance.csv",sep=''))

print("FINISH")

