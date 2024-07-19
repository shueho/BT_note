
#---00.初始化与加载包---

options(stringsAsFactors = F)
Sys.setenv(LANGUAGE = "en")
library(WGCNA)
library(FactoMineR) #PCA分析
library(factoextra)  
library(tidyverse) # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
### 启用WGCNA多核计算
enableWGCNAThreads(nThreads = 0.75*parallel::detectCores()) 

#---01.文件指定---

#指定基因表达矩阵和表现型矩阵（可为空）
geneExPath = 'DEG_fpkm.xls' #基因表达量路径
phenoPaths = c('Trait1.xls','Trait2.xls','Trait0.xls') #特征表文件名
# WGCNA原理
# 1.构建基因关系网络：
# ①基因间相关性系数；
# 【软阈值设定】，②幂指数接邻函数|Sij|^β；
# ③确定参数β；
# 【间接关系】，④拓扑矩阵（TOM）。

# ---02.基因相关性---

# 加载数据
dataExpr <- read.delim(geneExPath,header=T,row.names = 1,sep="\t",na.strings = "-")
dim(dataExpr)
dataExpr <- na.omit(dataExpr)
#dataExpr <- dataExpr[0:1000,]#测试时选择较少的基因
dim(dataExpr)

# 数据筛选
## 筛选中位绝对偏差前75%的基因，至少MAD大于0.01
## 筛选后会降低运算量，也会失去部分信息
## 也可不做筛选，使MAD大于0即可
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > 
                                max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]

## 转换为样品在行，基因在列的矩阵
dataExpr <- as.data.frame(t(dataExprVar))

## 检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3)


if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

dim(dataExpr)
head(dataExpr)[,1:8]

# 相关系数计算参数设置

# 官方推荐 "signed" 或 "signed hybrid"；
# 分析过程中unsigned无论是正负相关的基因，都取绝对值。
# signed会使得负相关得到的值很小，正相关的值很大，得到模块中的基因则都为正相关。
# 如果是单细胞的数据，那么就是因为一群基因正负调控得到的结果，这时候unsigned的比较合适。
type = "unsigned"
#type = "signed"
# 相关性计算，官方推荐 biweight mid-correlation & bicor，corType: pearson or bicor，为与原文档一致，故未修改，
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，或基因表达严重依赖于疾病状态时，需设置下面参数
#maxPOutliers = ifelse(corType=="pearson",1,0.05)
# 关联样品性状的二元变量时，设置
#robustY = ifelse(corType=="pearson",T,F)

# 软阈值筛选

#软阈值的筛选原则是使构建的网络更符合无标度网络特征。
## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")
pdf("01-Sample-clustering-to-detect-outliers.pdf",width = 6, height = 4)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()

powers = c(c(1:10), seq(from = 12, to=30, by=2)) #一般选取<16，推荐20
#需等待，和基因数目有关！！！
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5) #RsquaredCut = 0.85

pdf("02-chooing-power.pdf",width = 6, height = 6)
par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=0.85,col="red")  #一般选取>0.80

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")

dev.off()
power = sft$powerEstimate
power


#经验power (无满足条件的power时选用)
# 无向网络在power小于15或有向网络power小于30内，没有一个power值可以使
# 无标度网络图谱结构R^2达到0.8，平均连接度较高如在100以上，可能是由于
# 部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对
# 表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
# 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))       
                 )
  )
}


k = softConnectivity(dataExpr, power = power) -1
scaleFreePlot(k, main = paste("data set I, power=", power), truncated = F)

#---02.网络构建---

##一步法网络构建：One-step network construction and module detection##
# power: 上一步计算的软阈值
# maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)；
#  4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可
#  以处理3万个
#  计算资源允许的情况下最好放在一个block里面。
# corType: pearson or bicor
# numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
# saveTOMs：最耗费时间的计算，存储起来，供后续使用
# mergeCutHeight: 合并模块的阈值，越大模块越少


net = blockwiseModules(dataExpr, 
                       power = power, 
                       maxBlockSize = nGenes, #最大block大小，将所有基因放在一个block中
                       TOMType = type, 
                       #deepSplit = 2, 
                       minModuleSize = 30, #剪切树参数，deepSplit取值0-4
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.25, #模块合并参数，越大模块越少
                       numericLabels = TRUE, # T返回数字，F返回颜色
                       pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       #maxPOutliers=maxPOutliers, 
                       loadTOM =TRUE,
                       saveTOMFileBase = paste0(geneExPath, ".tom"),
                       verbose = 3)


table(net$colors)


# 层级聚类树展示各个模块
## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间

pdf("03-Cluster-Dendrogram.pdf",width = 12, height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#绘制模块之间相关性图

# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs

### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME.", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
pdf("04-Eigengene-adjacency-heatmap.pdf",width = 10, height = 6)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()

## 如果有表型数据，也可以跟ME数据放一起，一起出图
#MEs_colpheno = orderMEs(cbind(MEs_col, traitData))
#pdf("4-Cluster-Dendrogram.pdf",width = 6, height = 6)
#plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap", 
#                      marDendro = c(3,3,2,4),
#                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
#                      xLabelsAngle = 90)
#dev.off()


#---03.层次聚类（非常耗时）---

# 如果采用分步计算，或设置的blocksize>=总基因数，直接load计算好的TOM结果
# 否则需要再计算一遍，比较耗费时间
# TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
load(net$TOMFiles[1], verbose=T)
## Loading objects:
##   TOM
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong 
# connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function

if(F){
# ！可以先跳过，这一部分特别耗时，行列同时做层级聚类，最后做
pdf("05-Network-heatmap-plot.pdf",width = 6, height = 6)
TOMplot(plotTOM, net$dendrograms, moduleColors, 
        main = "Network heatmap plot, all genes")
dev.off()
}
#---04.导出网络用于Cytoscape---

probes = colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)

# Export the network into edge and node list files Cytoscape can read
# threshold 默认为0.5, 可以根据自己的需要调整，也可以都导出后在
# cytoscape中再调整
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste("06-",geneExPath, ".edges.txt", sep=""),
                               nodeFile = paste("07-",geneExPath, ".nodes.txt", sep=""),
                               weighted = TRUE, threshold = 0,  #阈值可以设置很小在软件上也可以设置
                               nodeNames = probes, nodeAttr = moduleColors)



#---05.自定义模块-表型关联图---

## 比如MEblue模块与Trait1相关

## 模块内基因与表型数据关联

# 性状跟模块虽然求出了相关性，可以挑选最相关的那些模块来分析，
# 但是模块本身仍然包含非常多的基因，还需进一步的寻找最重要的基因。
# 所有的模块都可以跟基因算出相关系数，所有的连续型性状也可以跟基因的表达
# 值算出相关系数。
# 如果跟性状显著相关基因也跟某个模块显著相关，那么这些基因可能就非常重要
# 。

### 计算模块与基因的相关性矩阵

if (corType=="pearsoon") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, 
                                         #robustY=robustY
                                         )
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}



if(F){
  
  # 计算性状与基因的相关性矩阵
  
  ## 只有连续型性状才能进行计算，如果是离散变量，在构建样品表时就转为0-1矩阵。
  
  if (corType=="pearsoon") {
    geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
    geneTraitP = as.data.frame(corPvalueStudent(
      as.matrix(geneTraitCor), nSamples))
  } else {
    geneTraitCorA = bicorAndPvalue(dataExpr, traitData, 
                                   # robustY=robustY
    )
    geneTraitCor = as.data.frame(geneTraitCorA$bicor)
    geneTraitP   = as.data.frame(geneTraitCorA$p)
  }
  
  ## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
  ## Pearson correlation was used for individual columns with zero (or missing)
  ## MAD.
# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
module = "blue"#去掉ME
pheno = "Trait1"
modNames = substring(colnames(MEs_col), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))
# 获取模块内的基因
moduleGenes = moduleColors == module

pdf("99-tt.pdf",width = 7,height = 7)
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
# 与性状高度相关的基因，也是与性状相关的模型的关键基因
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
}
#---07.其他图片---

#===============模块绘图==================
#将表达矩阵转换为一个颜色矩阵，使用log10（FPKM+1）
expColor=t(numbers2colors(log10(dataExpr+1),colors=blueWhiteRed(100),naColor="grey"))
colnames(expColor)=rownames(dataExpr)
#绘制基因的树形图，模块图，以及每个样品的表达图
pdf("08-wgcna.dendroColors.pdf",height = 7,width = 12)
plotDendroAndColors(net$dendrograms[[1]], 
                    colors=cbind(moduleColors[net$blockGenes[[1]]],expColor),
                    c("Module",colnames(expColor)),
                    dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05,
                    cex.rowText=0.5)
dev.off()

#绘制两两模块间的邻接矩阵
pdf("09-wgcna.adjacency.heatmap.pdf",height = 6,width = 9)
plotEigengeneNetworks(MEs_col, # MEs输出的是文本
                      "Eigengene adjacency heatmap",plotDendrograms = F,
                      marDendro = c(4,4,2,4))
dev.off()

#绘制所有模块的表达值热图与特征值条形图
flag = 0
for(module in substring(colnames(MEs_col),4)){ #注意substring的第二个参数是颜色前的前缀字数比如ME填写3，ME.填写4
  if(module == "grey") next
  #module="red"
  flag = flag+1
  ME=MEs_col[,paste("ME.",module,sep="")]
  pdf(paste("10-",flag,"-wgcna.", module, ".express.barplot.pdf", sep=""),height = 7,width = 9)
  par(mfrow=c(2,1),mar=c(0.3,5.5,3,2))
  plotMat(t(scale(dataExpr[,moduleColors==module])),
          rlabels=F,main=module,cex.main=2,clabels=F)
  
  par(mar=c(5,4.2,0,0.7))
  barplot(ME,col=module,main="",cex.main=2,ylab="eigengene expression",xlab="sample")
  dev.off()
}

#重新计算TOM矩阵

#TOM = TOMsimilarityFromExpr(dataExpr, power =sft$powerEstimate,TOMType = "unsigned"); 
flag = 0
for(module in substring(colnames(MEs_col),4)){  #注意substring的第二个参数是颜色前的前缀字数比如ME填写3，ME.填写4
  if(module == "grey") next
  flag = flag+1
  probes = colnames(dataExpr)
  inModule = is.finite(match(moduleColors, module))
  modProbes = probes[inModule]
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("11-",flag,"-CytoscapeInput-edges-", module, ".txt", sep=""),
                                 nodeFile = paste("12-",flag,"-CytoscapeInput-nodes-", module, ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0,   #TOMcutoff
                                 nodeNames = modProbes,
                                 nodeAttr = moduleColors[inModule])
}


#---关联表型---

#####################如果更改特征表，只需要重新执行下面的分析！！！
#计算ME与表型的相关系数，并计算p值
flag = 0
for(phenoPath in phenoPaths){
  #phenoPath = phenoPaths[1]
  flag = flag+1
  allTraits<-read.delim(phenoPath,header=T,sep="\t")  #读取表型值
  dataTraits = allTraits[, -1]
  tem = MEs_col[as.list(allTraits[1])$SampleID, ] #SampleID是对应特征表的样本列列名！
  moduleTraitCor = cor(tem, dataTraits, use = "p")# use=p代表去掉缺失计算
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(tem))
  #整理要显示在图中的数字,第一行为相关系数，第二行为p值
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  
  #绘制关联热图
  pdf(paste("13-",flag,"-wgcna.Module-trait.heatmap.pdf"), width = 10, height = 10)
  par(mar = c(6, 8.8, 3, 2.2))  
  labeledHeatmap(Matrix = moduleTraitCor,   #相关系数
                 xLabels = colnames(dataTraits), #x轴为表型
                 yLabels = names(MEs_col),#y轴为模块
                 ySymbols = names(MEs_col),
                 colorLabels = FALSE, 
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix, #每个单元格的内容
                 setStdMargins = FALSE,
                 cex.text = 0.8,
                 zlim = c(-1,1),
                 main = "Module-trait relationships")
  dev.off()
  
  pdf(paste("13-",flag,".2-wgcna.Module-trait.heatmap.pdf"), width = 15, height = 8)
  par(mar = c(6, 8.8, 3, 2.2))  
  labeledHeatmap(Matrix = moduleTraitCor,   #相关系数
                 xLabels = colnames(dataTraits), #x轴为表型
                 yLabels = names(MEs_col),#y轴为模块
                 ySymbols = names(MEs_col),
                 colorLabels = FALSE, 
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix, #每个单元格的内容
                 setStdMargins = FALSE,
                 cex.text = 0.8,
                 zlim = c(-1,1),
                 main = "Module-trait relationships")
  dev.off()
  write.table(moduleTraitCor,file = paste("14-",flag,"-Module-trait-Correlation-coefficient.xls"),sep = "\t")
  write.table(moduleTraitCor,file = paste("15-",flag,"-Module-trait-Correlation-Pvalue.xls"),sep = "\t")
  #绘制相关性系数图不需要重新运行
}





#---08.Hubgene提取---
# https://www.jianshu.com/p/8aa0c2b1851b
#1. 基于连通度connectivity，等同于常规网络中的degree。
#第一种： 先计算邻接矩阵，再计算连通度，推荐此种,邻接矩阵(Adjacency matrix)指基因和基因之间的加权相关性值取power次方即无尺度化之后构成的矩阵。
adjacency = abs(cor(dataExpr,use="p"))^power #datExpr为表达矩阵，power请与一步法中的软阈值保持一致
kIM <- intramodularConnectivity(adjacency, moduleColors)
write.table(kIM,file = "16-KIM.xls",sep = "\t")

#第二种：直接输入表达矩阵计算
#kIM <- intramodularConnectivity.fromExpr(dataExpr, colors, power = power)

#2. WGCNA自带一个函数，可以提取每个模块中连通度最高的基因
hub = chooseTopHubInEachModule(dataExpr,moduleColors, omitColors = "grey", power = power)
write.table(hub,file = "17-Hubgene.xls",sep = "\t")

#3. 有的文献中认为，hubgene 应该在某一模块中，与核心性状关联且与此模块也同样关联,比如：|GS|>.2&|MM|>.8 。
#GS 为 gene significance 缩写,反映某基因表达量与对应表型数值的相关性，调用 cor 函数计算相关性，绝对值越大，相关性越大，越趋近0，越不相关。
#MM 为 module membership 缩写，反映某基因表达量与模块特征值(module eigengenes)的相关。同样调用 cor 函数计算相关性，绝对值越大，相关性越大，越趋近0，越不相关。
#ps: 模块特征值 (module eigengenes) 是指基因和样本矩阵进行PCA分析后，每个模块的第一主成分即PC1,是此模块特性的高度代表,可以理解为一个模块就是一个超级代表基因。
#3.1.首先计算模块特征值(module eigengenes)
MEs2 = moduleEigengenes(dataExpr, moduleColors)$eigengenes

#3.2.计算module membership即MM
datKME =signedKME(dataExpr, MEs2, outputColumnName="MM.")
write.table(datKME,file = "18-MM.xls",sep = "\t")
#等价于：
#modNames = substring(names(MEs_col), 4)
#MM = as.data.frame(cor(dataExpr, MEs_col, use = "p"))


#3.3.gene significance即GS
flag = 0
for(phenoPath in phenoPaths){
  flag = flag+1
  allTraits<-read.delim(phenoPath,header=T,sep="\t")  #读取表型值
  dataTraits = allTraits[, -1]
  tem = dataExpr[as.list(allTraits[1])$SampleID, ] #SampleID是对应特征表的样本列列名！
  
  GS = as.data.frame(cor(tem, dataTraits, use = "p")) #datTraits为目标性状的矩阵文件
  write.table(GS,file=paste("19-",flag,"-GS.xls"),sep="\t")
}

## 模块内基因连接度

Alldegrees =intramodularConnectivity(adjacency, net$colors)
write.table(Alldegrees,file = "20-intramodularConnectivity.xls",sep = "\t")
#### kTotal:基因在整个网络中的连接度
#### kWithin: 基因在所属模块中的连接度，即Intramodular connectivity
#### kOut: kTotal-kWithin
#### kDiff: kIn-kOut

#也可以绘制一个模块基因的Intramodular connectivity与Gene significance的散点图
#3.4.循环提取每个模块的hubgene
MMname = colnames(datKME)
#for (mm in MMname){
  mm = MMname[1]
  FilterGenes =abs(GS)> .2 & abs(datKME$mm)>.8
  hubgenes = dimnames(data.frame(dataExpr))[[2]][FilterGenes]  
#}
#需要注意的是这里并没有添加显著检验，正确的做法应该是同时针对pvalue进行控制，比如：|GS|>.2 & |MM|>.8 & pvalue<0.05

#4. 同时，WGCNA官方建议结合 GS 与 kWithin 进行筛选
#比如: kWithin > 30 & |GS| >0.5





if(F){
  #5. 实际操作中，由于导出到 cytoscape 等可视化工具时边的数量格外多，就需要额外进行一个过滤步骤。此时常用的是针对 TOM 值做一个阈值的筛选，例如只保留 TOM>0.8 的基因关系对，此时计算模块内基因的 degree ，值越高，则越可能为hubgene 。
  #5.1. 载入TOM值
  # 如果采用分步计算，或设置的blocksize>=总基因数，直接load计算好的TOM结果
  # 否则需要再计算一遍，比较耗费时间
  #TOM = TOMsimilarityFromExpr(datExpr, power = power)
  
  #假如之前一步法中设定一个block计算，则可以直接载入节省时间
  # load(net_power$TOMFiles[1], verbose=T)
  # TOM <- as.matrix(TOM)
  
  #5.2. 设定阈值并输出文件
  geneid_allnet <- names(dataExpr)
  MEs = moduleEigengenes(dataExpr, moduleColors)$eigengenes
  modNames <- substring(names(MEs),4) 
  TOMcutoff <- 0.8 #实际操作中我建议设定一个小的值，到cytoscape中再进行个性化调整
  
  #循环输出所有模块
  for (mod in 1:nrow(table(moduleColors))){
    modules = names(table(moduleColors))[mod]
    # Select module probes
    probes = names(dataExpr)
    inModule = (moduleColors == modules)
    modProbes = probes[inModule]
    modGenes = modProbes
    # Select the corresponding Topological Overlap
    modTOM = TOM[inModule, inModule]
    
    dimnames(modTOM) = list(modProbes, modProbes)
    # Export the network into edge and node list files Cytoscape can read
    cyt = exportNetworkToCytoscape(modTOM,
                                   edgeFile = paste("CytoscapeInput-edges-", modules , ".txt", sep=""),
                                   nodeFile = paste("CytoscapeInput-nodes-", modules, ".txt", sep=""),
                                   weighted = TRUE,
                                   threshold = TOMcutoff,
                                   nodeNames = modProbes,
                                   altNodeNames = modGenes,
                                   nodeAttr = moduleColors[inModule])
  }
  
  #5.3，接下来针对CytoscapeInput-edges-*txt文件，进行统计，看哪个基因degree（边的数目总和）高，则此基因为hubgene
  
  #6. 除开以上连通度，度的衡量和过滤，hubgene的筛选还要采取灵活的方法，结合自己的研究目的及更多的数据进行。
  

#分步法展示每一步都做了什么
### 计算邻接矩阵
adjacency = adjacency(dataExpr, power = power)

### 把邻接矩阵转换为拓扑重叠矩阵，以降低噪音和假相关，获得距离矩阵。
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

### 层级聚类计算基因之间的距离树 
geneTree = hclust(as.dist(dissTOM), method = "average")

### 模块合并
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)

### 通过计算模块的代表性模式和模块之间的定量相似性评估，合并表达图谱相似的模块
MEList = moduleEigengenes(dataExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
MEDissThres = 0.25

# Call an automatic merging function
merge = mergeCloseModules(dataExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged

## 分步法完结
}

