if("sas7bdat" %in% rownames(installed.packages()) == FALSE) {source("http://bioconductor.org/biocLite.R");biocLite('sas7bdat')}
setwd("C:\\Users\\Zongliang\\Desktop\\new files\\AC study")
library("sas7bdat")
#AllInOne=read.sas7bdat("all_in_one.sas7bdat")
AllInOne=read.table("clearall_v2.txt",header=T,sep="\t")
head(AllInOne)
rowName=AllInOne[,1]
rowName
row.names(AllInOne)=make.names(rowName,unique=TRUE)

Data=AllInOne[,-1]
which(AllInOne=="A4032")
AllInOne=AllInOne[-1941,]
na_count <-sapply(AllInOne, function(AllInOne) sum(length(which(is.na(AllInOne)))))
StatOfNonNAN=-(na_count-nrow(AllInOne))
StatOfNonNAN
blank_count <-sapply(AllInOne, function(AllInOne) sum(length(which(AllInOne==""))))
na_count <- data.frame(na_count)
length(which(AllInOne==NA))
length(which(AllInOne==""))

###blank value ####
AllInOne=as.matrix(AllInOne)
head(AllInOne)
AllInOne=replace(AllInOne,AllInOne=="",NA)

function(AllInOne) sum(length(which(is.na(AllInOne))))




table(AllInOne[,"CCP_status"])
CCP_status=AllInOne[,"CCP_status"]

transfer=function(CCP_status){
  
  CCP_status=as.matrix(CCP_status)
  CCP_status[CCP_status=='POS']=as.numeric(1)
  CCP_status[CCP_status=='NEG']=as.numeric(-1)
  CCP_status[CCP_status=='UNK']=as.numeric(0)
  CCP_status=apply(CCP_status,c(1,2),as.numeric)
  return(CCP_status)
  
}
transfer(CCP_status)
table(AllInOne[,"RFIgM_status"])
RFIgM_status=AllInOne[,"RFIgM_status"]
transfer(RFIgM_status)
CCPAndRFIgM=transfer(CCP_status)+transfer(RFIgM_status)
colName=colnames(AllInOne)
type="AgeRA"
type="Gender"
res=matrix(nrow=length(colName),ncol=3)
ncount=0
NumericAllInOne=AllInOne
for(type in colName){
  ncount=ncount+1
  x1=AllInOne[,type]
  print(paste(ncount,type,sep=">"))
  
  NonNANPos=which(x1!="NaN")
  length(NonNANPos)
  if(typeof(x1[NonNANPos][1])=='character'){
    Num=as.numeric(factor(x1[NonNANPos]))
    x1[NonNANPos]=as.matrix(Num)
    x1[-NonNANPos]
    pvalue=cor(Num,CCPAndRFIgM[NonNANPos],method="spearman")
  }else{
    pvalue=cor(x1[NonNANPos],CCPAndRFIgM[NonNANPos],method="spearman")
  }
  
  res[ncount,]=c(type,length(NonNANPos),pvalue)
  NumericAllInOne[,type]=x1
}
head(NumericAllInOne)
res
PvalueMatrix=as.matrix(as.numeric(res[,3]))
rownames(PvalueMatrix)=res[,1]
### distribution of the spearman correlation in the parameters ###
#,xaxt = 'n'
hist(PvalueMatrix,breaks=seq(-1,1,0.05), freq = TRUE,xlim=c(-1,1),main="Histogram of the spearman correlation",xlab = "the spearman correlation")

hist(PvalueMatrix)
par(mar=c(5.1,4.1,4.1,2.1))

### threshold of the spearman correlation ###
plot(PvalueMatrix[abs(PvalueMatrix)>0.4],xaxt="n",cex=0.8)
axis(1, at=1:10, labels=as.vector(colName[which(abs(PvalueMatrix)>0.4)]),las = 2,cex.axis=0.5)


##### colnum cluster #####
#SampleCluster=read.table("Data_column_coding.txt",head=T,sep="\t")
SampleCluster=read.table("Column_Heading_categorization_v2.txt",head=T,sep=" ", fill = TRUE)
SampleCluster
SampleGroup=SampleCluster[-1]
rownames(SampleGroup)=SampleCluster[,1]
table(SampleGroup)
SampleGroup=as.matrix(SampleGroup)
row_na_count <-sapply(SampleGroup, function(AllInOne) sum(length(which(is.na(AllInOne)))))


category=unique(SampleGroup)
i=10
NumericAllInOne[,"VAS_Disease_Act_3"]

### test
array=Value
head(array)
TreeClust=function(array,xi){
  rown=nrow(array)
  rowNum=1
  array[1,]
  #row_na_count <-apply(array,1,function(AllInOne) sum(length(which(is.na(AllInOne)))))
  RemoveNA=function(x){if(is.na(x)){0}else{x}}
  ValueMatrix=apply(array,c(1,2),RemoveNA)
  #### tree ####
  sampleTree = hclust(dist(ValueMatrix,method = "euclidean"), method = "average")
  pdf(paste0("biTree",xi,".pdf"))
  plot(sampleTree, main = "hcluster", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2,)
  dev.off()
  groups <- cutree(sampleTree, k=2)
  tree1Name=names(groups)[which(groups==1)]
  tree2Name=names(groups)[which(groups==2)]
  Group1Name=rownames(CCPAndRFIgM)[which(CCPAndRFIgM==-2)]
  Group2Name=rownames(CCPAndRFIgM)[which(CCPAndRFIgM==2)]
  correct1Case1=length(which(tree1Name%in%Group1Name))
  correct2Case1=length(which(tree2Name%in%Group2Name))
  case1=correct1Case1+ correct2Case1
  correct1Case2=length(which(tree1Name%in%Group2Name))
  correct2Case2=length(which(tree2Name%in%Group1Name))
  case2=correct1Case2+ correct2Case2
  total=length(groups)
  rateOfCorrect=max(case1,case2)/total
  
  return(rateOfCorrect)
}
Group1Name=rownames(CCPAndRFIgM)[which(CCPAndRFIgM==-2)]
Group2Name=rownames(CCPAndRFIgM)[which(CCPAndRFIgM==2)]

resCorrect=matrix(ncol=2,nrow=length(category))
n=0
for(i in category){
  n=n+1
  GroupName=rownames(SampleGroup)[which(SampleGroup==i)]
  Value=NumericAllInOne[,which(colnames(NumericAllInOne)%in%GroupName)]
  Value=as.matrix(Value)
  rateOfCorrect=TreeClust(Value,i)
  resCorrect[n,]=c(i,rateOfCorrect)
}


baseline=min(length(Group1Name),length(Group2Name))/length(CCPAndRFIgM)
rownames(resCorrect)=resCorrect[,1]
resCorrect=resCorrect[,-1]
plot(resCorrect,xaxt="n",cex=0.8)
axis(1, at=1:12, labels=as.vector(names(resCorrect)),las = 2,cex.axis=0.5)
