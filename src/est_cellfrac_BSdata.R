# lib, R-4.4.1
library(EpiDISH)
 
## 输入和输出
args <- commandArgs(TRUE)
ref.file  <- args[1] # 参考数据集文件
bulk.file <- args[2] # 带评估的输入样本文件
out.file  <- args[3] # 输出评估结果，组分评估表格
pdf.file  <- args[4] # 输出评估结果，组分分布箱型图

#ref.file  <- "ref/ref.file"
#bulk.file <- "example/input.tsv"  
#out.file  <- "example/input.est_cellfrac_res.xls"
#pdf.file  <- "example/input.est_cellfrac_res.boxplot.pdf"

## 读入 ref marker    
ref.dat = read.table(ref.file, header =T)
ref.dt <- ref.dat[,-c(1:3)]
  
celltypes <- colnames(ref.dt) <- read.table(text =colnames(ref.dt), sep = ".")[,2]
      
ref <- matrix(NA, ncol = length(unique(celltypes)), nrow = nrow(ref.dt))
colnames(ref) <- unique(celltypes)
for(i in colnames(ref)){
    ref[,i] <- apply(ref.dt[, celltypes == i], 1, mean)
}


## 读入待评估的数据   
df.dat = read.table(bulk.file, header =T, check.names = F)
df <- df.dat[,-c(1:4)]


## 组分评估  
rownames(df) <-   rownames(ref)  <- paste("R",c(1:nrow(df)), sep="_")
predicted_cell_proportions <- epidish(as.matrix(df), as.matrix(ref), method='RPC')$estF


## 输出保存   
out <- data.frame(sample.ids =colnames(df), predicted_cell_proportions)
write.table(out, out.file, row.names = F, sep = "\t", quote = F)


## 其他分析
### 宽边长
library(reshape2)
data_long <- melt(out,
             id.vars= c("sample.ids"),  #需要保留不参与聚合的变量,
             measure.vars= 2:14,        #用于聚合的变量, 也可以使用列标
             variable.name= 'cell_type', 
             value.name= 'cell_ratio')

### 绘制组间的boxplot
library(ggplot2)
# data_long$group <- factor(data_long$group, levels=c("HC","NPA","AD","CRC"))
p <- ggplot(data = data_long, aes(x= cell_type, y = cell_ratio)) +
  geom_boxplot(outlier.size = 0, outlier.color = NA) +
  geom_point(position = position_jitterdodge(), alpha=0.3, size=0.3) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  scale_y_continuous(limits = c(0, 1)) 
#  xlab("") +
#  ylab("") +

pdf(file=pdf.file, height=4, width=8)
print(p)
dev.off()

#END
