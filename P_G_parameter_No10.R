
library(fs)
#---library R  package-------
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(EasyStat)
library(EasyMicrobiome)

# 设置主题
library(ggthemes)
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
  axis.text.x = element_text(colour = "black",size = 14),
  axis.text.y = element_text(colour = "black",size = 14),
  legend.text = element_text(size = 15,face = "bold")
)


#----set color or fill#------------
library(RColorBrewer)#调色板调用包
#调用所有这个包中的调色板
display.brewer.all()
#提取特定个数的调色板颜色，会出图显示
display.brewer.pal(9,"Set1")
# colset1 <- c(brewer.pal(12,"Paired"),brewer.pal(9,"Pastel1"))
colset1 <- c(brewer.pal(12,"Paired"))
colset2 <- brewer.pal(12,"Paired")
colset3 <- c(brewer.pal(12,"Paired"),brewer.pal(9,"Pastel1"))





#--提取有多少个分组
map = sample_data(ps)
gnum <- unique(map$Group) %>% length()

# 分组标签为Group


# 设定排序顺序
axis_order = unique(map$Group)

# 堆叠柱状图展示TOp前是个门,j 展示的水平为门水平
Top = 10
# -类似韦恩图的二分网络设置过滤阈值 序列数量
# bionum = 200

# 差异分析 edger设定分组

# group1 = c("Gro1","Gro2")
# b= data.frame(group1)
# 如果每两个组之间都做差异，那就指定
b = NULL

# # 如果多个组和CK比较，则可以使用这个功能Plot.CompareWithCK
# # 那么上面edger的b必须设定为都和这个设定的CK比较
# CK = "Gro1"


# 热图和气泡图对丰度最高的前多少个OTU做
heatnum　=　20

