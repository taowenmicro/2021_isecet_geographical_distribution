#-分析共四个部分，这是第二个部分，代号2，第一部分是扩增子原始序列处理
library(EasyMicrobiome)
#---result1 base_diversity analyses-------------
otupath = paste(res1path,"/OTU_sourcetrack/",sep = "")
dir.create(otupath)

#------------虫子微生物群落来源#-----------
source("G:\\Shared_Folder\\Function_local\\R_function\\micro\\FEAST.R",encoding = "UTF-8")

head(map)


alppath = paste(otupath,"/Feast/",sep = "")
dir.create(alppath)


unique(map$Group)

ps_sub <- phyloseq::subset_samples(ps, !ID %in% c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12","E13"));ps_sub
map = sample_data(ps_sub)
table(map$Group)



result = FEAST(ps = ps_sub,group = "Group",sinkG = "L.sativae",sourceG = c("Root","Soil","Leaf"))
result

# -例子
p1 <- Plot_FEAST(data = result)

p2 <- MuiPlot_FEAST(data = result)


FileName <- paste(alppath,"Feast_all", ".pdf", sep = "")
ggsave(FileName, p1, width = 6, height =6,limitsize = FALSE)

FileName <- paste(alppath,"Feast_sample", ".pdf", sep = "")
ggsave(FileName, p2, width = 6, height =6,limitsize = FALSE)



FileName <- paste(alppath,"Feast_all", ".jpg", sep = "")
ggsave(FileName, p1, width = 6, height =6,limitsize = FALSE)

FileName <- paste(alppath,"Feast_sample", ".jpg", sep = "")
ggsave(FileName, p2, width = 6, height =6,limitsize = FALSE)


