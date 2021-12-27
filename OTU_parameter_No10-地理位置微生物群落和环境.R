

#---library R  package
library(fs)# 文件夹处理
library(phyloseq)
library(tidyverse)

library(ggClusterNet)
library(EasyStat)
library(EasyMicrobiome)

library(ggthemes)# 主题
library(RColorBrewer)#调色板调用包

# # 建立第二级目录
result_path <- paste("./amplicon1/","/No_11_zhuyuxi",sep = "")
dir_create(result_path)

# 设置主题
mytheme1 = theme_bw() + theme(
  panel.background=element_blank(),
  panel.grid=element_blank(),
  legend.position="right",
  
  legend.title = element_blank(),
  legend.background=element_blank(),
  legend.key=element_blank(),
  # legend.text= element_text(size=7),
  # text=element_text(),
  # axis.text.x=element_text(angle=45,vjust=1, hjust=1)
  plot.title = element_text(vjust = -8.5,hjust = 0.1),
  axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
  axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
  axis.text = element_text(size = 20,face = "bold"),
  axis.text.x = element_text(colour = "black",size = 14,angle = 90),
  axis.text.y = element_text(colour = "black",size = 14),
  legend.text = element_text(size = 15,face = "bold")
)
# 设置主题
mytheme2 = theme_bw() + theme(
  panel.background=element_blank(),
  panel.grid=element_blank(),
  legend.position="right",
  
  legend.title = element_blank(),
  legend.background=element_blank(),
  legend.key=element_blank(),
  # legend.text= element_text(size=7),
  # text=element_text(),
  # axis.text.x=element_text(angle=45,vjust=1, hjust=1)
  plot.title = element_text(vjust = -8.5,hjust = 0.1),
  axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
  axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
  axis.text = element_text(size = 20,face = "bold"),
  axis.text.x = element_text(colour = "black",size = 14,angle = 90),
  axis.text.y = element_text(colour = "black",size = 14),
  legend.text = element_text(size = 15,face = "bold")
)


#----set color or fill
#调用所有这个包中的调色板
RColorBrewer::display.brewer.all()
#提取特定个数的调色板颜色，会出图显示
RColorBrewer::display.brewer.pal(9,"Set1")
# colset1 <- c(brewer.pal(12,"Paired"),brewer.pal(9,"Pastel1"))
colset1 <- brewer.pal(12,"Set3")
colset2 <- brewer.pal(12,"Paired")
colset3 <- c(brewer.pal(12,"Paired"),brewer.pal(9,"Pastel1"))
colset4 <- c(brewer.pal(12,"Paired"),brewer.pal(12,"Set3"),brewer.pal(9,"Pastel1"),brewer.pal(9,"Set1"))




#AM真菌#----------
Bac = T
Fun = F
AM = F



#------构建phyloseq对象

library(phyloseq)
library(ggClusterNet)
library(tidyverse)

ps0 = readRDS("./data_all/No_11_zhuyuxi//ps_insect.rds")

otu = t(vegan_otu(ps0))
head(otu)

map = sample_data(ps0)
head(map)

map$Group = gsub(", ","_",map$Location)


table(map$Group)
map$Group = as.factor(map$Group)

sample_data(ps0)  =map

tax = vegan_tax(ps0)
head(tax)
# tree  = read_tree("./data_all/No_9_jiangshangtao/AM/iqtree/otus.contree")
# library(Biostrings)
# rep = readDNAStringSet("./data_all/No_9_jiangshangtao/AM/reptree.fa")

# # 导入phyloseq(ps)对象
# library(phyloseq)
# 
# ps = phyloseq(sample_data(map),
#               otu_table(as.matrix(otu), taxa_are_rows=TRUE),
#               tax_table(as.matrix(tax)), phy_tree(tree)
#               
# )
# ps
# saveRDS(ps,"./data_all/No_9_jiangshangtao/AM/ps.rds")



if (Bac) {
  # ps0 = readRDS("./data_all/No_10_beigeshidi/16s///ps.rds")
  ps0 <- ps0 %>%
    subset_taxa(
      Kingdom == "d__Bacteria"
    )
  ps0 = filter_taxa(ps0, function(x) sum(x ) > 0 , TRUE);ps0
  
  TFdir.b <- as.data.frame(table(tax_table(ps0)[,1]))[,2][as.data.frame(table(tax_table(ps0)[,1]))[,1] %in% c("Bacteria","K__Bacteria","k__Bacteria")] > 10
  if (length(TFdir.b) != 0) {
    print("16s")
    res1path <- paste(result_path,"/OTU_base_diversity_16s",sep = "")
  }
  dir.create(res1path)
}

#--样本数量调整
ps0 <- prune_samples(sample_sums(ps0) >=2000,ps0);ps0
ps = ps0

head(map)


#--提取有多少个分组
map = sample_data(ps0)
gnum <- unique(map$Group) %>% length()

table(map$Group)
# 分组标签为Group


# 设定排序顺序
axis_order = unique(map$Group)


# 堆叠柱状图展示TOp前是个门,j 展示的水平为门水平
Top = 10

rank.names(ps0)
j = "Phylum"
# j = "Genus"
# -类似韦恩图的二分网络设置过滤阈值 序列数量
bionum = 0

# 差异分析 edger设定分组
# group1 = c("Gro1","Gro2")
# b= data.frame(group1)
# 如果每两个组之间都做差异，那就指定
b = NULL



# 如果多个组和CK比较，则可以使用这个功能Plot.CompareWithCK
# 那么上面edger的b必须设定为都和这个设定的CK比较
map
CK = "Laboratory"



# 热图和气泡图对丰度最高的前多少个OTU做
heatnum　=　50



#--R 语言做lefse法分析-过滤一下需要序列
# ps_sub <- ps0 %>%
#   subset_taxa(
#     # Kingdom == "Fungi"
#     Kingdom == "Bacteria"
#     )
ps_sub = filter_taxa(ps0, function(x) sum(x ) > 1000, TRUE);ps_sub
ps_sub <- prune_samples(sample_sums(ps_sub) >=200,ps_sub);ps_sub


#---机器学习部分
# ROC是三种机器学习的ROC曲线，但是只能跑两个分组，如果两组，可以选择为T。
ROC = T
# 是否做交叉检验
rfcv = F
# 选择多少个重要变量
optimal = 30




# 网络
# 过滤多少丰度的做网络
N = 0.001
zipi = F

# ## 环境因子
# data = read_xlsx("./data_all/No_9_jiangshangtao/AM/env.xlsx") %>% as.data.frame()
# head(data)
# 
# ps_env =  ps0
# 
# data$ID %in%sample_data(ps0)$ID
# 
# data <- data[data$ID %in%sample_data(ps0)$ID,]
# 
# # #---细菌和真菌网络
# # library(igraph)
# # library(sna)
# # library(ggrepel)
# # value1 = 0
