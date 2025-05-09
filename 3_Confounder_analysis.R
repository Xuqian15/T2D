#!/usr/bin/env Rscript
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T)
}

if (TRUE){
  option_list = list(
    make_option(c("-W", "--workplace"), type="character", default="/data3/Group5/xuqian/TD-fg/multi-kingdom/xMarkerFinder/data/input_files/",
                help="Input workplace [default %default]"),
    make_option(c("-m", "--metadata"), type="character", default="merged-metadata.csv",
                help="metadata file [default %default]"),
    make_option(c("-p", "--profile"), type="character", default="Bacteria_filtered_abundance.csv",
                help="filtered relative abundance profile[default %default]"),
    make_option(c("-d", "--distance"), type="character", default="bc",
                help="the distance metric [default %default]"),
    make_option(c("-c", "--count"), type="character", default="999",
                help="the number of permutations [default %default]"),
    make_option(c("-g", "--group"), type="character", default="Group",
                help="the column name of experiment interest(group) in metadata [default %default]"),
    make_option(c("-o", "--output"), type="character", default="Bacteria",
                help="output file prefix [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))

  print(paste("The workplace is ", opts$workplace,  sep = ""))
  print(paste("The metadata file is ", opts$metadata,  sep = ""))
  print(paste("The profile of features is ", opts$profile,  sep = ""))
  print(paste("The distance metric is ", opts$distance,  sep = ""))
  print(paste("The permutation count is ", opts$count,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}


package_list = c("vegan","dplyr","labdsv","ggplot2","parallel","foreach","doParallel")

site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

metadata <- read.csv(file = paste(opts$workplace,opts$metadata,sep=''),sep = ',',header =  TRUE,row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
feat_abd <- read.csv(file = paste(opts$workplace,opts$profile,sep=''),sep = ',',header =  TRUE,row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
feat_abd <-feat_abd[rownames(feat_abd) %in% rownames(metadata), ]
##这里metadata有很多缺失值，字符向量手动填补，数值向量用均值填补
num_cols <- sapply(metadata, is.numeric)
# 对数值型列进行填补
metadata[num_cols] <- lapply(metadata[num_cols], function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
})
###adonis需要特征表的行列转换作为输入
#feat_abd = t(feat_abd)
group <-opts$group
distance <- opts$distance
count <- as.numeric(opts$count)


#设置并行计算的核心数
#numCores <- detectCores() - 1  # 使用全部核心数减去1
#cl <- makeCluster(numCores)
#registerDoParallel(cl)
#多线程计算
#distances <- foreach(method = c('bray', 'euclidean'), .packages = 'vegan', .combine = 'list') %dopar% {
#    if (method == 'bray') {
#        vegdist(feat_abd, method = 'bray', correction = 'lingoes', na.rm = TRUE)
#    } else if (method == 'euclidean') {
#        vegdist(feat_abd, method = 'euclidean', correction = 'lingoes', na.rm = TRUE)
#    }
#}

#stopCluster(cl)  # 关闭并行计算的集群

# 根据 distance 参数选择相应的距离矩阵
#if (distance == "bc") {
#    dist <- distances[1]$dist  # 'bray' 距离矩阵
#} else if (distance == "euclidean") {
#    dist <- distances[[2]]  # 'euclidean' 距离矩阵
#}


if (distance == "bc"){
    dist = vegdist(feat_abd, method = 'bray', correction = 'lingoes', na.rm = TRUE)
} else if (distance == "euclidean"){
    dist = vegdist(feat_abd, method = 'euclidean', correction = 'lingoes', na.rm = TRUE)
}

colnames

phenotype_microbiota_corr<- function(dist,microbiota){
    metadata2 <- metadata[,1:10]
    df_result <- matrix(nrow = 2, ncol = ncol(metadata2),0)
    row.names(df_result) <- c('R2','Pvalue')
    colnames(df_result) <- colnames(metadata2)
    for (i in 1:ncol(metadata2)){
        tmp = metadata2[,i]
        if (distance == "bc"){
            permanova = adonis(dist~tmp, permutations=count, method = "bray")
        } else if (distance == "euclidean"){
            permanova = adonis(dist~tmp, permutations=count, method = "euclidean")
        }
        df_result['R2',colnames(metadata2)[i]] <- as.data.frame(permanova[1])['tmp','aov.tab.R2']
        df_result['Pvalue',colnames(metadata2)[i]] <- as.data.frame(permanova[1])['tmp','aov.tab.Pr..F.']
    }
    return(df_result)    
}

#permanova test
dist2 <- as.matrix(dist)
result <- phenotype_microbiota_corr(dist2,feat_abd)
write.csv(result,file = paste(opts$workplace,opts$output,"_metadata_microbiota.csv",sep=''))

#plot
dist2[is.na(dist2)]<-0
pco.results=pco(dist,k=2)
axis.1.title <- paste('PCoA1 ', 
                      round((pco.results$eig[1]/sum(pco.results$eig))*100,1),
                      '%', sep='')
axis.2.title <- paste('PCoA2 ', 
                      round((pco.results$eig[2]/sum(pco.results$eig))*100,1),
                      '%', sep='')

df.plot <- data.frame(Axis1 = -1*pco.results$points[,1],
                  Axis2 = pco.results$points[,2],
                  Sample_ID = rownames(pco.results$points),
                  Group=unlist(metadata[group]))

pcoa_plot<-ggplot(df.plot,aes(x=Axis1,y=Axis2,colour=Group))+
                    geom_point(alpha =.7, size=4)+
                    scale_size_area()+scale_colour_brewer(type = "div", palette = "Dark2")+
                    xlab(axis.1.title)+
                    ylab(axis.2.title)+
                    stat_ellipse(level=0.95,linetype=2,type="norm")+
                    theme_classic()+geom_vline(xintercept = 0, color = 'black', size = 0.4, linetype = 4)+
                    geom_hline(yintercept = 0, color = 'black', size = 0.4, linetype = 4)+
                    theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
                        panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
                        panel.background = element_rect(color = 'black', fill = 'transparent'),
                        legend.title=element_blank(),aspect.ratio=1)+
                    ggtitle(paste0("PCoA plot based on the ",distance," distance"))
pdf(file = paste(opts$workplace,opts$output,"_pcoa_plot.pdf",sep=''),width=10,height=10)
pcoa_plot
while (!is.null(dev.list()))  dev.off()

svg(file = paste(opts$workplace,opts$output,"_pcoa_plot.svg",sep=''),width=10,height=10)
pcoa_plot
while (!is.null(dev.list()))  dev.off()

print("FINISH")