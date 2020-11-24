WindowInfo=function(target.df, RefFaFileName="chr20.fa", WindowSize=100){
library(exomeCopy)
library(IRanges)

colnames(target.df)=c("seqnames","start","end")
RefFaFile =scanFa(RefFaFileName)
CalIRange=function(start,end,WindowSize){
  NofW= ceiling((end-start+1)/ WindowSize)
  Istart= start+(c(1:(NofW))-1)* WindowSize    
  IRangesSep = IRanges(Istart, width= WindowSize)
end(IRangesSep[length(IRangesSep)])= end
  return(c(IRangesSep= IRangesSep, NofW = NofW))
}

IRangesSepF=CalIRange(target.df$start[1], target.df$end[1],WindowSize)
IRangesSep=IRangesSepF$IRangesSep
IRNumber= IRangesSepF$ NofW
if(nrow(target.df)>1){
 for(i in 2:nrow(target.df)){
   IRangesSepF=CalIRange(target.df$start[i],target.df$end[i],WindowSize)   
IRangesSepT=IRangesSepF $ IRangesSep
   IRangesSep =c(IRangesSep, IRangesSepT) 
IRNumberT=IRangesSepF $ NofW
IRNumber=c(IRNumber, IRNumberT)
 }}
target <- GRanges(seqname= rep(target.df$seqname, IRNumber), ranges=IRangesSep)

rdata <- RangedData(space=seqnames(target),ranges=ranges(target))

IRangeSE=cbind(start(rdata), end(rdata))
colnames(IRangeSE)=c("start","end")
IRangeSE=as.data.frame(IRangeSE)

nonAmb=function(IRangeSE, start,end, RefFaFile){
return(letterFrequency (narrow(RefFaFile, start= IRangeSE[start], end= IRangeSE [end]),  "ATCG", as.prob=T))
} 
#rdata[[nonAmb]] =apply(IRangeSE, MARGIN=1, nonAmb, start= "start", end="end", RefFaFile = RefFaFile)
nonAmb=apply(IRangeSE, MARGIN=1, nonAmb, start= "start", end="end", RefFaFile = RefFaFile)

GCcontent=function(IRangeSE ,start, end, Fa){
return(letterFrequency (narrow(RefFaFile, start= IRangeSE[start], end= IRangeSE [end]),  "GC", as.prob=T))
} 
GC= apply(IRangeSE, MARGIN=1, GCcontent, start= "start", end="end", Fa= RefFaFile)

CONY.TempRegion=cbind(rep(noquote(names(rdata)),nrow(rdata)),start(rdata), end(rdata),width(rdata), nonAmb,GC)
return(CONY.TempRegion)
write.table(CONY.TempRegion, "CONY.TempRegion.txt",quote = F, sep = " ",row.names = F,col.names =c("seq","start","end","width","nonAmb","GC"))
}


################################################################################################################################################################################################################

CalRD=function(TempRegion, CRDMethod=c("PointR","SumUp"), SampleBamFileName= "aln1101.final.sorted.bam", MPileCountFileName="aln1101.Used.txt",SampleName="Sample1101",TargetChr="chr20", WindowSize=100){

library(exomeCopy)
library(IRanges)

if(CRDMethod=="PointR"){
IRangesSep =IRanges(start=TempRegion$start, end=TempRegion$end)
target <- GRanges(seqname= rep(TargetChr, nrow(TempRegion)), ranges=IRangesSep)
rdata <- RangedData(space=seqnames(target),ranges=ranges(target))
scanBamHeader(SampleBamFileName)[[1]]$targets
levels(seqnames(target))
ARD=countBamInGRanges(SampleBamFileName,target)
write.table(cbind(TempRegion,ARD), paste("CONY.1-TempRegion.",TargetChr,".",SampleName,".PointR.RD.txt", sep=""),quote = F, sep = " ",row.names = F,col.names =T) 
return(cbind(TempRegion,ARD))

}



if(CRDMethod=="SumUp"){
MPileCount=read.table(MPileCountFileName, header=F)
colnames(MPileCount)=c("location","count") 
MPileCount=as.data.frame(MPileCount)

RDVec=rep(0,TempRegion$end[nrow(TempRegion)])
RDVec[MPileCount[,1]]= MPileCount[,2]
RDVec=c(RDVec, rep(0,(100-length(RDVec)%%100)) )
RDMatrix=matrix(RDVec,100,ceiling(length(RDVec)/100))
SRD=colSums(RDMatrix)
#ARD= SRD* WindowSize /(TempRegion$nonAmb * TempRegion$width)
ARD= SRD /(TempRegion$nonAmb * TempRegion$width)

write.table(cbind(TempRegion,ARD), paste("CONY.1-TempRegion.",TargetChr,".",SampleName,".SumUp.RD.txt", sep=""),quote = F, sep = " ",row.names = F,col.names =T) 
return(cbind(TempRegion,ARD))

}
}


##################################################################################################################################################################################################################

AdjRD=function(CRDMethod=c("PointR","SumUp"), TargetChr, SampleName){
if(CRDMethod=="PointR"){
TempRegion=read.table(paste("CONY.1-TempRegion.",TargetChr,".",SampleName,".PointR.RD.txt", sep=""),header=T)
}
if(CRDMethod=="SumUp"){
TempRegion=read.table(paste("CONY.1-TempRegion.",TargetChr,".",SampleName,".SumUp.RD.txt", sep=""), header=T) 
}

TempInx=which( (TempRegion$GC>=0.2) & (TempRegion$GC<=0.8) & (TempRegion$ARD!="NaN") &(TempRegion$ARD!=0))
LOWS=lowess(TempRegion$ARD [TempInx]~TempRegion$GC[TempInx])
MLS=mean(LOWS$y)
AdjRD= TempRegion$ARD 
T.GC=c(20:80)/100
T.RD=rep(0,length(T.GC))
for( i in 1:length(T.GC)){
      T.RD[i]= LOWS$y[which(LOWS$x== T.GC [i])] [1]
}
AdjRD[TempInx]= TempRegion$ARD[TempInx]* MLS /T.RD[ match(round(TempRegion$GC[TempInx],2) ,T.GC)]
TempRegion=cbind(TempRegion, AdjRD)
write.table(TempRegion, paste("CONY.2-TempRegion.",TargetChr,".",SampleName, ".", CRDMethod ,".AdjRD.txt", sep="")  ,quote = F, sep = " ",row.names = F,col.names =T) 
}

####################################################################################################################################################################################################################
UsedRD=function(CRDMethod=c("PointR","SumUp"), AnaMethod=c("Single","Paired"), TargetChr, SampleName,ControlName){

  if(AnaMethod=="Single"){
    TempRegion= read.table(paste("CONY.2-TempRegion.",TargetChr,".",SampleName, ".", CRDMethod , ".AdjRD.txt", sep=""),header=T)
    NonInf=which((is.na(TempRegion$AdjRD)==T) | (TempRegion$nonAmb<0.5))
    NonInfRegion=TempRegion[NonInf,]
    write.table(NonInfRegion, paste("CONY.3-NonInfRegion.",TargetChr,".",SampleName, ".", CRDMethod , ".", AnaMethod, ".txt", sep=""),quote = F, sep = " ",row.names = F,col.names =T)

    # only discard the NonInf regions if there are any
    if (length(NonInf)>0) { TempRegion = TempRegion[-NonInf,] }

    # get the CN0 region 
    CN0=which(TempRegion$AdjRD==0)
    CN0Region=TempRegion[CN0,]
    write.table(CN0Region, paste("CONY.3-CN0Region.",TargetChr,".",SampleName, ".", CRDMethod , ".", AnaMethod,".txt", sep=""),quote = F, sep = " ",row.names = F,col.names =T)
    
    # only discard the CN0 region if there is any 
    if (length(CN0)>0) { TempRegion = TempRegion[-CN0,] }
      
    RD=TempRegion$AdjRD
    target=rep(1,length(RD))
    write.table(cbind(TempRegion,RD, target), paste("CONY.3-TempRegion.",TargetChr,".",SampleName, ".", CRDMethod , ".", AnaMethod ,".UsedRD.txt", sep=""),quote = F, sep = " ",row.names = F,col.names =T)
  }

  if(AnaMethod=="Paired"){

    # leave an error message
    print("AnaMethod Paired should be debugged taking into consideration that there might not be any 0 or Tmp regions ")
    quit(status=1)

    TempRegion1= read.table(paste("CONY.2-TempRegion.",TargetChr,".",SampleName, ".", CRDMethod ,".AdjRD.txt", sep=""),header=T) 
    TempRegion2= read.table(paste("CONY.2-TempRegion.",TargetChr,".",ControlName, ".", CRDMethod ,".AdjRD.txt", sep=""),header=T) 

    NonInf=which((is.na(TempRegion1$AdjRD)==T)| (TempRegion1$nonAmb<0.5))
    NonInfRegion=TempRegion1[NonInf,]
    write.table(NonInfRegion, paste("CONY.3-NonInfRegion.",TargetChr,".",SampleName, ".", CRDMethod , ".", AnaMethod, ".txt", sep=""),quote = F, sep = " ",row.names = F,col.names =T)

    TempRegion1= TempRegion1[-NonInf,]
    TempRegion2= TempRegion2[-NonInf,]
    Both0=which((TempRegion1$AdjRD==0)& (TempRegion2$AdjRD==0))
    CN0=which((TempRegion1$AdjRD==0)& (TempRegion2$AdjRD!=0))
    CNG=which((TempRegion1$AdjRD!=0)& (TempRegion2$AdjRD==0))
    Both0Region= TempRegion1[Both0,]
    CN0Region= TempRegion1[CN0,]
    CNGRegion= TempRegion1[CNG,]

    write.table(Both0Region, paste("CONY.3-Both0Region.",TargetChr,".",SampleName, ".", CRDMethod , ".", AnaMethod,".txt", sep=""),quote = F, sep = " ",row.names = F,col.names =T)
    write.table(CN0Region, paste("CONY.3-CN0Region.",TargetChr,".",SampleName, ".", CRDMethod , ".", AnaMethod,".txt", sep=""),quote = F, sep = " ",row.names = F,col.names =T)
    write.table(CNGRegion, paste("CONY.3-CNGRegion.",TargetChr,".",SampleName, ".", CRDMethod , ".", AnaMethod,".txt", sep=""),quote = F, sep = " ",row.names = F,col.names =T)
    Del=sort(c(Both0, CN0, CNG))
    TempRegion1= TempRegion1[-Del,]
    TempRegion2= TempRegion2[-Del,]
    RD=TempRegion1$AdjRD/ TempRegion2$AdjRD
    target=rep(1,length(RD))
    TempRegion=cbind(TempRegion1,RD,target)
    write.table(TempRegion, paste("CONY.3-TempRegion.",TargetChr,".",SampleName, ".", CRDMethod , ".", AnaMethod ,".UsedRD.txt", sep=""),quote = F, sep = " ",row.names = F,col.names =T)
  }
  }

#######################################################################################################################################################################################################################
EstPar=
function(CRDMethod=c("PointR","SumUp"), AnaMethod=c("Single","Paired"), TargetChr, SampleName, NCN) {

TempRegion=read.table(paste("CONY.3-TempRegion.",TargetChr,".",SampleName, ".", CRDMethod , ".", AnaMethod ,".UsedRD.txt", sep=""),header=T)

# mean of CN2
Med0=median(na.omit(log(TempRegion$RD)))
# var of CN2
logV=((sort(log(TempRegion$RD),decreasing = T)[round(nrow(TempRegion) *.025)]- Med0)/1.96)^2
#estimate the mean and variance for each CN state
Mean2=exp(Med0+logV/2)
Var2=(exp(logV) -1)*exp(2*Med0+logV)

OMean=Mean2*(c(1: NCN)/2)
OVar= Var2*(c(1: NCN)/2)

GroupMean =log(OMean^2/((OVar+OMean^2)^0.5))
GroupVar =log(1+(OVar/((OMean)^2)))

#estimate of proportion of each CN state
AllP=matrix(rep(0,NCN),NCN,length(TempRegion$RD),byrow=T) 
for(i in 1:NCN){
AllP[i,]= dlnorm((TempRegion$RD),GroupMean[i],(GroupVar[i])^0.5)
}
AsP=as.numeric(apply(t(abs(AllP)), 1, which.max))
GroupPro=rep(0,NCN)
for(i in 1:NCN){
GroupPro[i]=length(which(AsP==i))/length(TempRegion$RD)
}

GroupCNV=c(1:NCN)
(NCN);(GroupMean);(GroupVar);(GroupPro);(GroupCNV)


GroupSumm=cbind(GroupCNV,GroupMean,GroupVar,GroupPro)
colnames(GroupSumm)=c("GroupCNV","GroupMean","GroupVar","GroupPro")
write.table(GroupSumm, paste("CONY.3-GroupSumm.",TargetChr,".",SampleName, ".", CRDMethod , ".", AnaMethod ,".txt", sep=""),quote = F, sep = " ",row.names =F,col.names =T)
}
#######################################################################################################################################################################################################################

RunCONY= function(CRDMethod=c("PointR","SumUp"), AnaMethod=c("Single","Paired"), TargetChr, SampleName, RunTime = 300000,BurnN = 5000,RTN = 1000,BCPoint = 20, FragLength=500000){




TempRegion=read.table(paste("CONY.3-TempRegion.",TargetChr,".",SampleName, ".", CRDMethod , ".", AnaMethod ,".UsedRD.txt", sep=""),header=T)
GroupSumm=read.table(paste("CONY.3-GroupSumm.",TargetChr,".",SampleName, ".", CRDMethod , ".", AnaMethod ,".txt", sep=""),header =T)


RunCONYParallel=function(Frag,AFrag,TempRegion, GroupSumm, RunTime = 300000,BurnN = 5000,RTN = 1000,BCPoint = 20, TargetChr, SampleName) {

TempRegion= TempRegion[(ceiling((Frag-1)*nrow(TempRegion)/AFrag)+1): (ceiling(Frag*nrow(TempRegion)/AFrag)),]

SetIni=function(TempRegion,GroupSumm){
### Setting fixed parameters 
print("setting fixed parameters")
RD=log(TempRegion$RD)
if(length(which(RD==Inf))>0){RD[(which(RD==Inf))]=NA }
if(length(which(RD==(-Inf)))>0){RD[(which(RD==(-Inf)))]=NA }
Inx_AD=TempRegion$start
RDsq=RD^2
target=TempRegion$target
MeanF= GroupSumm$GroupMean
GroupN= nrow(GroupSumm)
AlphaF=rep(5, GroupN)
GroupPro= GroupSumm$ GroupPro
GroupVar= GroupSumm$GroupVar
GroupCNV= GroupSumm$ GroupCNV
KappaF=GroupPro*length(RD) 
KappaF[which(GroupPro==0)]=1
BetaF=((AlphaF-1)*GroupVar*KappaF)/( KappaF+1)
LambdaF=0.5
WeightIni=GroupPro
WeightUp=matrix(0,length(GroupPro),length(GroupPro))
for(i in 1: length(GroupPro)){
WeightUp[i,]=GroupPro
for(j in 1: length(GroupPro)){
WeightUp[i,j]=GroupPro[j]/(1-GroupPro[i])
}
WeightUp[i,i]=0
}

GroupPartLog=rep(0, GroupN)
for(k in 1:GroupN){
GroupPartLog[k]=log((BetaF[k]^AlphaF[k])* (KappaF[k]^0.5) / gamma(AlphaF[k]))
}
#### Breakpoint initial
### setting breakpoint gap breakpoint & last one
print("setting initial")

if(TempRegion$target[nrow(TempRegion)]==0){
 BrFix=sort(c(which(TempRegion$target==0), (which(TempRegion$target==0))-1) ) 
if(TempRegion$target[1]==0){ 
  BrFix= BrFix[-1]
   }
}

if(TempRegion$target[nrow(TempRegion)]!=0){
 BrFix=sort(c(which(TempRegion$target==0), (which(TempRegion$target==0))-1)  )
if(TempRegion$target[1]==0){ 
  BrFix= BrFix[-1]
   }
BrFix=c(BrFix, nrow(TempRegion))
}

if(length(BrFix)==1){CNRUpBr=matrix(c(0,BrFix),1,2)}
if(length(BrFix)>1){
if(TempRegion$target[1]==0){
     CNRUpBr=matrix(0, floor(length(BrFix) /2),2)
for(i in 1: (floor(length(BrFix) /2))){
CNRUpBr[i,]=c(BrFix[2*i-1],BrFix[2*i])
}
   }
if(TempRegion$target[1]==1){
     CNRUpBr=matrix(0, floor(length(BrFix) /2)+1,2)
CNRUpBr[1,]=c(0,BrFix[1])
for(i in 2: (floor(length(BrFix) /2)+1)){
CNRUpBr[i,]=c(BrFix[2*i-2],BrFix[2*i-1])
}
   }
}

#### breakpoint

SRM=smooth.spline(na.omit(RD))
SRD=rep(NA,length(RD))

SRMX =which(is.na(RD)==F)
SRD[SRMX]= SRM$y

SRDT=SRM$y

SRD1=SRDT[1: (length(SRDT)-1)]
SRD2=SRDT[2: (length(SRDT))]

LowB=0.02
UppB=0.98

BrInxT=sort(which(   ((SRD1<quantile(SRDT,LowB))&(SRD2>=quantile(SRDT,LowB)))
|
((SRD1>quantile(SRDT,LowB))&(SRD2<=quantile(SRDT,LowB)))
|
((SRD1<quantile(SRDT,UppB))&(SRD2>=quantile(SRDT,UppB)))
|
((SRD1>quantile(SRDT,UppB))&(SRD2<=quantile(SRDT,UppB)))
))

BrInxO= SRMX [BrInxT]

if(length(which(is.na(match(BrInxO,BrFix))==F))!=0){
BrInxO= BrInxO[-(which(is.na(match(BrInxO,BrFix))==F))]
}
BrInx=sort(c(BrInxO,BrFix))

if(is.na(match(length(RD),BrInx))==T){ BrInx =c(BrInx,length(RD))}
Br=rep(0,length(RD)) 
Br[BrInx]=1
## CN initial
CNVTargetInx=rep(0,length(BrInx))
if(BrFix[1]==1){
for(i in 2:length(BrInx)){
CNVTargetInx[i]=TempRegion$target[BrInx[i]]
}
}
if(BrFix[1]!=1){
for(i in 1:length(BrInx)){
CNVTargetInx[i]=TempRegion$target[BrInx[i]]
}}

CNR=rep(0,length(BrInx))
CNR[which(CNVTargetInx==1)]= rep(c(1:GroupN), ceiling (length(which(CNVTargetInx==1))/GroupN)) [1: (length(which(CNVTargetInx==1)))]
return(list("TempRegion"=TempRegion,"BrFix"=BrFix, "BrInx"=BrInx, "Br"= Br , "CNRUpBr"= CNRUpBr , "CNR"=CNR ,"RD"=RD ,"RDsq"=RDsq ,"GroupN"=GroupN ,"MeanF"=MeanF ,"AlphaF"=AlphaF ,"KappaF"=KappaF ,"BetaF"=BetaF ,"LambdaF"=LambdaF ,"WeightIni"=WeightIni, "WeightUp"=WeightUp, "GroupPartLog"= GroupPartLog, "GroupCNV"= GroupCNV))
}
DP3=SetIni(TempRegion,GroupSumm)



TempRegion=DP3$TempRegion
RD=DP3$RD
RDsq=DP3$ RDsq
GroupN=DP3$ GroupN
MeanF=DP3$ MeanF
BetaF=DP3$ BetaF
AlphaF=DP3$ AlphaF
KappaF=DP3$ KappaF
LambdaF=DP3$ LambdaF
WeightIni=DP3$ WeightIni
WeightUp=DP3$ WeightUp
GroupPartLog=DP3$ GroupPartLog
CNR=DP3$ CNR
BrInx=DP3$ BrInx
GroupCNV=DP3$ GroupCNV
BrFix=DP3$ BrFix
CNRUpBr=DP3$ CNRUpBr
target=TempRegion$target

#source("DataPrep5N.r")
source("CONY.R")

tic<- Sys.time()
CNRUp=UpC(RD,RDsq,GroupN,MeanF,BetaF,AlphaF,KappaF,LambdaF,WeightIni,WeightUp ,GroupPartLog,CNR,BrInx,GroupCNV,BrFix,CNRUpBr)
toc<- Sys.time()
comp.time<- toc - tic
print(comp.time)

CNMode=matrix(0,GroupN,length(RD))
BFInxTemp=rep(0,length(RD))

BFDiff=matrix(NA,1,2)
colnames(BFDiff)=c("times","diff")
tor=0/length(RD)

for(i in 2:BurnN){

tp=sample(c(1,2,3,4),1,prob=c(   ( (GroupN-1)/(GroupN*3)),(1/3),(1/(GroupN*3)), (1/3)   )    )
print(c(i,tp))
BB=switch(tp,Spt(RD, RDsq, GroupN,MeanF, BetaF,AlphaF,KappaF,LambdaF, WeightIni,WeightUp, GroupPartLog,CNRUp,BrInx,GroupCNV,BrFix,CNRUpBr), Meg(RD, RDsq, GroupN,MeanF, BetaF,AlphaF,KappaF,LambdaF, WeightIni,WeightUp, GroupPartLog,CNRUp,BrInx,GroupCNV,BrFix,CNRUpBr),DBSpt(RD, RDsq, GroupN,MeanF, BetaF,AlphaF,KappaF,LambdaF, WeightIni,WeightUp, GroupPartLog,CNRUp,BrInx,GroupCNV,BrFix,CNRUpBr),Cbound(RD, RDsq, GroupN,MeanF, BetaF,AlphaF,KappaF,LambdaF, WeightIni,WeightUp, GroupPartLog,CNRUp,BrInx,GroupCNV,BrFix,CNRUpBr) )

CNRUp=BB$CNRUpN
BrInx=BB$BrInxN
#ChangInx=BB$ChangInx

print(c (BB$Type,BB$Result))
print(c("Number of Breakpoint",length(BB$BrInxN)))
if(BB$Result=="Success"     ){
tic<- Sys.time()
CNRUp= UpC (RD,RDsq,GroupN,MeanF,BetaF,AlphaF,KappaF,LambdaF,WeightIni,WeightUp ,GroupPartLog,CNR= BB$CNRUpN,BrInx= BB$BrInxN,GroupCNV,BrFix,CNRUpBr)
toc<- Sys.time()
comp.time<- toc - tic
print(comp.time)
}
gc()
}

SuccessPro=0
SuccessProVec=0

for(i in (BurnN+1):RunTime){

tp=sample(c(1,2,3,4),1,prob=c(( (GroupN-1)/(GroupN*3)),(1/3),(1/(GroupN*3)), (1/3)))
print(c(i,tp))
BB=switch(tp,Spt(RD, RDsq, GroupN,MeanF, BetaF,AlphaF,KappaF,LambdaF, WeightIni,WeightUp, GroupPartLog,CNRUp,BrInx,GroupCNV,BrFix,CNRUpBr), Meg(RD, RDsq, GroupN,MeanF, BetaF,AlphaF,KappaF,LambdaF, WeightIni,WeightUp, GroupPartLog,CNRUp,BrInx,GroupCNV,BrFix,CNRUpBr),DBSpt(RD, RDsq, GroupN,MeanF, BetaF,AlphaF,KappaF,LambdaF, WeightIni,WeightUp, GroupPartLog,CNRUp,BrInx,GroupCNV,BrFix,CNRUpBr),Cbound(RD, RDsq, GroupN,MeanF, BetaF,AlphaF,KappaF,LambdaF, WeightIni,WeightUp, GroupPartLog,CNRUp,BrInx,GroupCNV,BrFix,CNRUpBr) )

CNRUp=BB$CNRUpN
BrInx=BB$BrInxN

if(BB$Result=="Success"){
SuccessPro=SuccessPro+1
CNRUp= UpC (RD,RDsq,GroupN,MeanF,BetaF,AlphaF,KappaF,LambdaF,WeightIni,WeightUp ,GroupPartLog,CNR= BB$CNRUpN,BrInx= BB$BrInxN,GroupCNV,BrFix,CNRUpBr)

}

CNModeTemp= CNRWin (RD,CNRUp,BrInx)
CNModeTS=matrix(0,GroupN,length(CNModeTemp))
for(j in 1:GroupN){
CNModeTS[j,which(CNModeTemp==j)]=1
}
CNMode= CNMode+ CNModeTS
rm(CNModeTS)

if( (i%%RTN)==0){
BFDec=matrix(0,GroupN,length(RD))
BFInx=rep(0,length(RD))

SuccessProVec=c(SuccessProVec,SuccessPro/(i-BurnN))


for(m in 1:GroupN){
BFDec[m,]=(CNMode[m,]/CNMode[2,])/( WeightIni[m]/ WeightIni[2])
}
for(j in 1:length(RD)){
if(target[j]==0){BFInx[j]=0}
if(target[j]!=0) {
if( max(na.omit(BFDec[,j])) >BCPoint ){ BFInx[j]=which.max(BFDec[,j]) }
if( max(na.omit(BFDec[,j])) <=BCPoint ){ BFInx[j]=2 }
  }
}
rm(BFDec)
BFDiffT= length(which(BFInx!=BFInxTemp))/length(RD)
BFDiffTM=matrix(c(i, BFDiffT),1,2)
BFDiff=rbind(BFDiff,BFDiffTM)
BFInxTemp=BFInx
print(BFDiff)
gc()

if(BFDiffT<=tor){break}

}
gc()
}


return(BFInxTemp)
}

#source("http://bioconductor.org/biocLite.R")
#biocLite("snow")


library(snow)
AFrag=floor(sum(TempRegion$width)/FragLength)
cl <- makeCluster(AFrag)
AA=matrix(c(1:AFrag), AFrag, 1)
source("CONY.R")
CN =parRapply(cl, AA, RunCONYParallel, AFrag,TempRegion, GroupSumm, RunTime,BurnN,RTN,BCPoint, TargetChr, SampleName)
stopCluster(cl)

write.table(cbind(TempRegion,CN),paste("CONY.4-Result.",TargetChr,".",SampleName, ".", CRDMethod , ".", AnaMethod,".txt", sep=""),quote = F, sep = " ",row.names = F,col.names =T)

}


#######################################################################################################################################################################################################################



ComResult=function(CRDMethod="PointR", AnaMethod="Single", TargetChr="chr20", SampleName="Sample1101"){

CONYResultT=read.table(paste("CONY.4-Result.",TargetChr,".",SampleName, ".", CRDMethod , ".", AnaMethod,".txt", sep=""),header=T)
TempRegionFull=read.table(paste("CONY.2-TempRegion.",TargetChr,".",SampleName,".", CRDMethod ,".AdjRD.txt", sep=""),header=T)

CN=rep(NA,nrow(TempRegionFull))
TempRegionFull=cbind(TempRegionFull,CN)
ResultInx=match(CONYResultT$start, TempRegionFull$start)
TempRegionFull$CN[ResultInx]= CONYResultT$CN

if(AnaMethod=="Single"){
CN0Region=read.table(paste("CONY.3-CN0Region.",TargetChr,".",SampleName, ".", CRDMethod , ".", AnaMethod,".txt", sep=""),header=T)
CN0Inx= match(CN0Region$start, TempRegionFull$start)
TempRegionFull$CN[CN0Inx]= 0
}




if(AnaMethod=="Paired"){
Both0Region= read.table(paste("CONY.3-Both0Region.",TargetChr,".",SampleName, ".", CRDMethod , ".", AnaMethod,".txt", sep=""),header=T)
CN0Region= read.table(paste("CONY.3-CN0Region.",TargetChr,".",SampleName[1], ".", CRDMethod , ".", AnaMethod,".txt", sep=""),header=T)
CNGRegion= read.table(paste("CONY.3-CNGRegion.",TargetChr,".",SampleName[1], ".", CRDMethod , ".", AnaMethod,".txt", sep=""),header=T)
Both0Inx= match(Both0Region$start, TempRegionFull$start)
CN0Inx= match(CN0Region$start, TempRegionFull$start)
CNGInx= match(CNGRegion$start, TempRegionFull$start)
TempRegionFull$CN[Both0Inx]=2
TempRegionFull$CN[CN0Inx]=0
TempRegionFull$CN[CNGInx]=9
}

write.table(TempRegionFull, paste("CONY.Result.",TargetChr,".",SampleName,".", CRDMethod, ".", AnaMethod ,".Window.txt", sep=""),quote = F, sep = " ",row.names = F,col.names =T)

TempRegionFullCNTemp=TempRegionFull$CN
TempRegionFullCNTemp [which(is.na(TempRegionFullCNTemp)==T)]=-1

TempRegionFullCN1= TempRegionFullCNTemp [-1]
TempRegionFullCN2= TempRegionFullCNTemp [-nrow(TempRegionFull)]
TempRegionFullCND=TempRegionFullCN1- TempRegionFullCN2

Breakpoint=c(0,which(TempRegionFullCND!=0), nrow(TempRegionFull))

CONY.CNVResult=matrix(,length(Breakpoint)-1,3)
colnames(CONY.CNVResult)=c("start","end","CN")
CONY.CNVResult=as.data.frame(CONY.CNVResult)

for(i in 1: (length(Breakpoint)-1)){
   CONY.CNVResult$start[i]= TempRegionFull$start[(Breakpoint[i]+1) ]
   CONY.CNVResult$end[i]= TempRegionFull$end[Breakpoint[(i+1)]]
   CONY.CNVResult$CN[i]= TempRegionFull$CN[(Breakpoint[i]+1)]
}
CONY.CNV= CONY.CNVResult [which((is.na(CONY.CNVResult$CN)==F)& (CONY.CNVResult$CN!=2)),]

write.table(CONY.CNVResult, paste("CONY.Result.",TargetChr,".",SampleName,".", CRDMethod, ".", AnaMethod ,".CNRegionAll.txt", sep=""),quote = F, sep = " ",row.names = F,col.names =T)
write.table(CONY.CNV, paste("CONY.Result.",TargetChr,".",SampleName,".", CRDMethod, ".", AnaMethod ,".CNV.txt", sep=""),quote = F, sep = " ",row.names = F,col.names =T)

}



########################################################################################################################################################################################################################

### Update C step by step 1010 2013 delete front and back

UpC<-function(RD,RDsq,GroupN,MeanF,BetaF,AlphaF,KappaF,LambdaF,WeightIni,WeightUp ,GroupPartLog,CNR,BrInx,GroupCNV,BrFix,CNRUpBr){

BrInxT=c(0,BrInx)
for(i in 1:length(BrInx)){
   if(CNR[i]!=0){
CNVComb=matrix(CNR,GroupN,length(CNR),byrow=T)
TS=rep(0,GroupN)

for(j in 1:GroupN){ 
CNVComb[j,i]=GroupCNV[j]
GS=sum(GroupPartLog[CNVComb[j,]])
WS=WeightPartL (WeightIni,WeightUp, CNR= CNVComb[j,],GroupN,CountMatrix= CountC (CNR= CNVComb[j,],GroupN))  
RS=RegionPartL (RD,RDsq, Wstart= BrInxT[i]+1 ,Wend= BrInxT[i+1],CNR= CNVComb[j,i],MeanF,BetaF,AlphaF,KappaF)   

TS[j]=GS+WS+RS
#print(c(GS,WS,RS,TS))
}

inx=order(TS,decreasing=T)


if( (i!=1) & (i!= length(BrInx))){ 
print("1")

  if((CNR[i-1]!=0) & ( CNR[i+1]!=0)){
print("1.1")
   CND=c(CNR[i-1], CNR[i+1])
CNR=CNVComb [inx[which(is.na(match(inx,CND)!=T))[1]],]
}
  if((CNR[i-1]!=0) & ( CNR[i+1]==0)){
print("1.2")
   CND=CNR[i-1]
CNR=CNVComb [inx[which(is.na(match(inx,CND)!=T))[1]],]
}
  if((CNR[i-1]==0) & ( CNR[i+1]!=0)){
print("1.3")
   CND=CNR[i+1]
CNR=CNVComb [inx[which(is.na(match(inx,CND)!=T))[1]],]
}
  if((CNR[i-1]==0) & ( CNR[i+1]==0)){
print("1.4")
CNR=CNVComb [inx[1],]
}
}


if( (i==1) & (i!= length(BrInx))){ 
print("2")
  if(CNR[i+1]!=0){
print("2.1")
   CND=CNR[i+1]
CNR=CNVComb [inx[which(is.na(match(inx,CND)!=T))[1]],]
}
  if( CNR[i+1]==0){
print("2.2")
CNR=CNVComb [inx[1],]
}
}

if( (i!=1) & (i== length(BrInx))){ 
print("3")
  if(CNR[i-1]!=0){
print("3.1")
   CND=CNR[i-1]
CNR=CNVComb [inx[which(is.na(match(inx,CND)!=T))[1]],]
}
  if( CNR[i-1]==0){
print("3.2")
CNR=CNVComb [inx[1],]
}
}

if( (i==1) & (i== length(BrInx))){ 
print("4")
CNR=CNVComb [inx[1],]
}

}
   if(CNR[i]==0){CNR=CNR}

}
return(CNR)
}



### Step 4.1 Merge 

Meg<-function(RD, RDsq, GroupN,MeanF, BetaF,AlphaF,KappaF,LambdaF, WeightIni,WeightUp, GroupPartLog,CNRUp,BrInx,GroupCNV,BrFix,CNRUpBr){

BrInxO=sort(BrInx[-(match(BrFix,BrInx))]) 

DecPro=0
if(length(BrInxO)==0){DecPro =0}
if(length(BrInxO)!=0){

if(length(BrInxO)==1){CPoint=BrInxO}
if(length(BrInxO)>1){CPoint=sample(BrInxO,1)}

ChangInx=which((CNRUpBr[,1]<CPoint) &(CNRUpBr[,2]>CPoint))

#CN Region
UPCNInx=which(BrInx==CPoint)
CNRTC=sample(c(1,2),1)
CNRT=switch(CNRTC, CNRUp[UPCNInx], CNRUp[UPCNInx+1])


















#General Form
if(  (  (UPCNInx==1)& (is.na(match(CNRT,CNRUp[UPCNInx+2]))==T)  ) |  (  (UPCNInx>1)& (is.na(match(CNRT,c(CNRUp[UPCNInx-1],CNRUp[UPCNInx+2])))==T)   )    ){
print("GeneralFormMerge")
CNRUpN=CNRUp;CNRUpN[UPCNInx]=CNRT;CNRUpN=CNRUpN[-(UPCNInx+1)]
BrInxN=BrInx[-(UPCNInx)]

WSN=WeightPartL(WeightIni,WeightUp, CNRUpN,GroupN, CountMatrix= as.matrix(CountC(CNRUpN,GroupN)))

WSO=WeightPartL(WeightIni,WeightUp, CNRUp,GroupN, CountMatrix= as.matrix(CountC(CNRUp,GroupN)))

GSN= GroupPartLog[CNRT]
GSO=sum(GroupPartLog[c(CNRUp[UPCNInx], CNRUp[UPCNInx+1])])

if(UPCNInx==1){WstartN=1;WstartO1=1}
if(UPCNInx>1){WstartN=BrInx[UPCNInx-1]+1; WstartO1= BrInx[UPCNInx-1]+1}
WendN= BrInx[UPCNInx+1]
WstartO2= BrInx[UPCNInx]+1
WendO1= BrInx[UPCNInx]; WendO2= BrInx[UPCNInx+1];

RSN=RegionPartL (RD,RDsq, WstartN , WendN ,CNR= CNRT,MeanF,BetaF,AlphaF,KappaF)


RSO=RegionPartL (RD,RDsq, WstartO1 , WendO1 ,CNR= CNRUp[UPCNInx],MeanF,BetaF,AlphaF,KappaF)+RegionPartL (RD,RDsq, WstartO2 , WendO2 ,CNR= CNRUp[UPCNInx+1],MeanF,BetaF,AlphaF,KappaF)


if( (CNRTC==1)&&( is.na(match(0,CNRUp[UPCNInx+2]))==T)){TranP= -log(GroupN-2)}
if( (CNRTC==2)&&( is.na(match(0,CNRUp[UPCNInx-1]))==T)){TranP= -log(GroupN-2)}
if( (CNRTC==1)&&( is.na(match(0,CNRUp[UPCNInx+2]))==F)){TranP= -log(GroupN-1)}
if( (CNRTC==2)&&( is.na(match(0,CNRUp[UPCNInx-1]))==F)){TranP= -log(GroupN-1)}

DecPro=exp(GSN+WSN+RSN- GSO-WSO-RSO+log(1- LambdaF)-log(LambdaF)+TranP)
DecPro
}

































#Double Merge

if(   ((UPCNInx==1)& (is.na(match(CNRT,CNRUp[UPCNInx+2]))==F))   |  ((UPCNInx>1)& (is.na(match(CNRT,c(CNRUp[UPCNInx-1],CNRUp[UPCNInx+2])))==F))    ) {


# Form1
if(CNRTC==1){
print("DoubleMerge Form1")
   CNRUpN=CNRUp;CNRUpN[UPCNInx]=CNRT;CNRUpN=CNRUpN[-((UPCNInx+1): (UPCNInx+2))]
BrInxN=BrInx[-( (UPCNInx) : (UPCNInx+1))]

WSN=WeightPartL(WeightIni,WeightUp, CNRUpN,GroupN, CountMatrix= as.matrix(CountC(CNRUpN,GroupN)))

WSO=WeightPartL(WeightIni,WeightUp, CNRUp,GroupN, CountMatrix= as.matrix(CountC(CNRUp,GroupN)))

GSN= GroupPartLog[CNRT]
GSO=sum(GroupPartLog[c(CNRUp[UPCNInx], CNRUp[UPCNInx+1], CNRUp[UPCNInx+2])])

if(UPCNInx==1){WstartN=1;WstartO1=1}
if(UPCNInx>1){WstartN=BrInx[UPCNInx-1]+1; WstartO1= BrInx[UPCNInx-1]+1}
WendN= BrInx[UPCNInx+2]
WstartO2= BrInx[UPCNInx]+1; WstartO3= BrInx[UPCNInx+1]+1;
WendO1= BrInx[UPCNInx]; WendO2= BrInx[UPCNInx+1]; WendO3= BrInx[UPCNInx+2]

RSN=RegionPartL (RD,RDsq, WstartN , WendN ,CNR= CNRT,MeanF,BetaF,AlphaF,KappaF)


RSO=RegionPartL (RD,RDsq, WstartO1 , WendO1 ,CNR= CNRUp[UPCNInx],MeanF,BetaF,AlphaF,KappaF)+RegionPartL (RD,RDsq, WstartO2 , WendO2 ,CNR= CNRUp[UPCNInx+1],MeanF,BetaF,AlphaF,KappaF) +RegionPartL (RD,RDsq, WstartO3 , WendO3 ,CNR=CNRUp[UPCNInx+2],MeanF,BetaF,AlphaF,KappaF)


if(is.na(match(0,CNRUp[length(CNRUp)]))==F){TL=1}
if(is.na(match(0,CNRUp[length(CNRUp)]))==T){TL=0}
if(is.na(match(0,CNRUp[1]))==F){T1=1}
if(is.na(match(0,CNRUp[1]))==T){T1=0}

TranP =(length(RD)-2-rowSums(CountC(CNRUpN,GroupN))[GroupN+1]-length(which(is.na(RD)==T))+T1+TL)*4/(GroupN-1)

DecPro=exp(GSN+WSN+RSN- GSO-WSO-RSO+2*log(1- LambdaF)-2*log(LambdaF)+log(TranP))

}






















# Form2
if(CNRTC==2){
print("DoubleMerge Form2")
   CNRUpN=CNRUp;CNRUpN[UPCNInx-1]=CNRT;CNRUpN=CNRUpN[-((UPCNInx): (UPCNInx+1))]
BrInxN=BrInx[-((UPCNInx-1): (UPCNInx))]

WSN=WeightPartL(WeightIni,WeightUp, CNRUpN,GroupN, CountMatrix= as.matrix(CountC(CNRUpN,GroupN)))

WSO=WeightPartL(WeightIni,WeightUp, CNRUp,GroupN, CountMatrix= as.matrix(CountC(CNRUp,GroupN)))

GSN= GroupPartLog[CNRT]
GSO=sum(GroupPartLog[c(CNRUp[UPCNInx-1], CNRUp[UPCNInx], CNRUp[UPCNInx+1])])

if((UPCNInx-1)==1){WstartN=1;WstartO1=1}
if((UPCNInx-1)>1){WstartN=BrInx[UPCNInx-2]+1; WstartO1= BrInx[UPCNInx-2]+1}
WendN= BrInx[UPCNInx+1]
WstartO2= BrInx[UPCNInx-1]+1; WstartO3= BrInx[UPCNInx]+1;
WendO1= BrInx[UPCNInx-1]; WendO2= BrInx[UPCNInx]; WendO3= BrInx[UPCNInx+1]

RSN=RegionPartL (RD,RDsq, WstartN , WendN ,CNR= CNRT,MeanF,BetaF,AlphaF,KappaF)


RSO=RegionPartL (RD,RDsq, WstartO1 , WendO1 ,CNR= CNRUp[UPCNInx-1],MeanF,BetaF,AlphaF,KappaF)+RegionPartL (RD,RDsq, WstartO2 , WendO2 ,CNR= CNRUp[UPCNInx],MeanF,BetaF,AlphaF,KappaF) +RegionPartL (RD,RDsq, WstartO3 , WendO3 ,CNR=CNRUp[UPCNInx+1],MeanF,BetaF,AlphaF,KappaF)


if(is.na(match(0,CNRUp[length(CNRUp)]))==F){TL=1}
if(is.na(match(0,CNRUp[length(CNRUp)]))==T){TL=0}
if(is.na(match(0,CNRUp[1]))==F){T1=1}
if(is.na(match(0,CNRUp[1]))==T){T1=0}
TranP =4/((GroupN-1)*(length(RD)-2-rowSums(CountC(CNRUpN,GroupN))[GroupN+1]-length(which(is.na(RD)==T))+T1+TL))

DecPro=exp(GSN+WSN+RSN- GSO-WSO-RSO+2*log(1- LambdaF) -2*log(LambdaF)+ log(TranP))

}
}
}

#if(is.na(DecPro)==T){DecPro=0}
MAccp=rbinom(1,1,min(1, DecPro))

if (MAccp==1){
return(list("Type"= "Merge","Result"="Success", "CNRUpN"=CNRUpN,"BrInxN" =BrInxN, "CNRO "= CNRUp, "BrInxO" =BrInx))
}
if (MAccp==0){
return(list("Type"="Merge" ,"Result"="Fail" , "CNRUpN"=CNRUp,"BrInxN" =BrInx, "CNRO"= CNRUp,"BrInxO" =BrInx))
}



}











# split
Spt<-function(RD, RDsq, GroupN,MeanF, BetaF,AlphaF,KappaF,LambdaF, WeightIni,WeightUp, GroupPartLog,CNRUp,BrInx,GroupCNV,BrFix,CNRUpBr){

BrInxO=BrInx[-(match(BrFix,BrInx))] 
CNDiff=rep(0,(length(CNRUp)))
if(length(CNRUp)>1){
for (i in 2: (length(CNRUp))){
CNDiff[i]=BrInx[i]-BrInx[i-1] }
}
CNDiff[1]= BrInx[1]
DecPro =0
if(length(which((CNDiff>1)&(CNRUp!=0)))==0 ){DecPro=0}
if(length(which((CNDiff>1)&(CNRUp!=0)))>0 ){

if(length(which((CNDiff>1)&(CNRUp!=0)))==1 ){SptReg= which((CNDiff>1)&(CNRUp!=0))}
if(length(which((CNDiff>1)&(CNRUp!=0)))>1 ){
SptReg=sample(which((CNDiff>1)&(CNRUp!=0)) ,1)}

if(SptReg==1){CPoint=sample(c(1: (BrInx[1]-1)),1)}
if(SptReg!=1){
   if( (BrInx[SptReg-1]+1)== (BrInx[SptReg]-1)  ){CPoint= (BrInx[SptReg-1]+1)}
   if( (BrInx[SptReg-1]+1)!= (BrInx[SptReg]-1)  ){
CPoint=sample(c( (BrInx[SptReg-1]+1): (BrInx[SptReg]-1) ),1)
}
}
ChangInx=which((CNRUpBr[,1]<CPoint) &(CNRUpBr[,2]>CPoint))

BrInxON=sort(c(BrInxO,CPoint))
  BrInxN= sort(c(BrInx,CPoint))

#Region
UPCNInx=which(BrInxN==CPoint)
CNRUpN=rep(NA,length(CNRUp)+1)
if(UPCNInx>1){CNRUpN[1: (UPCNInx-1)]=CNRUp[1: (UPCNInx-1)]}
if(UPCNInx<(length(CNRUpN)-1)){CNRUpN[(UPCNInx+2): (length(CNRUpN))]=CNRUp[(UPCNInx+1): (length(CNRUp))]}
CNRTC=sample(c(1,2),1)
CNRT=CNRUp[UPCNInx]

if(CNRTC==1){ 
CNRUpN[UPCNInx]= CNRUp [UPCNInx]

if(UPCNInx<(length(CNRUpN)-1)){
CNRUpN[UPCNInx+1]=sample(GroupCNV[-(which  ( (GroupCNV==CNRUpN[UPCNInx+2]) | (GroupCNV==CNRUpN[UPCNInx])))] ,1)
}

if(UPCNInx==(length(CNRUpN)-1)){
CNRUpN[UPCNInx+1]=sample(GroupCNV[-(which ( (GroupCNV==CNRUpN[UPCNInx])))] ,1)
}

}
if(CNRTC==2){
    CNRUpN[UPCNInx+1]= CNRUp [UPCNInx]

if(UPCNInx>1){
CNRUpN[UPCNInx]=sample(GroupCNV[-(which  ( (GroupCNV==CNRUpN[UPCNInx+1]) | (GroupCNV==CNRUpN[UPCNInx-1])))] ,1)
}
if(UPCNInx==1){
CNRUpN[UPCNInx]=sample(GroupCNV[-(which ( (GroupCNV==CNRUp[1])))] ,1)
}
}

WSN=WeightPartL(WeightIni,WeightUp, CNRUpN,GroupN, CountMatrix= as.matrix(CountC(CNRUpN,GroupN)))

WSO=WeightPartL(WeightIni,WeightUp, CNRUp,GroupN, CountMatrix= as.matrix(CountC(CNRUp,GroupN)))

GSO= GroupPartLog[CNRT]
GSN=sum(GroupPartLog[c(CNRUpN[UPCNInx], CNRUpN[UPCNInx+1])])

if(UPCNInx==1){WstartN1=1;WstartO=1}
if(UPCNInx>1){WstartN1=BrInxN[UPCNInx-1]+1; WstartO= BrInx[UPCNInx-1]+1}
WstartN2= BrInxN[UPCNInx]+1

WendN1= BrInxN[UPCNInx]
WendN2= BrInxN[UPCNInx+1]
WendO= BrInx[UPCNInx]



RSO=RegionPartL (RD,RDsq, WstartO , WendO ,CNR= CNRT,MeanF,BetaF,AlphaF,KappaF)

RSN=RegionPartL (RD,RDsq, WstartN1 , WendN1 ,CNR= CNRUpN[UPCNInx],MeanF,BetaF,AlphaF,KappaF)+RegionPartL (RD,RDsq, WstartN2 , WendN2 ,CNR= CNRUpN[UPCNInx+1],MeanF,BetaF,AlphaF,KappaF)

if( (CNRTC==1)&&( is.na(match(0,CNRUp[UPCNInx+2]))==T)){TranP= log(GroupN-2)}
if( (CNRTC==2)&&( is.na(match(0,CNRUp[UPCNInx-1]))==T)){TranP= log(GroupN-2)}
if( (CNRTC==1)&&( is.na(match(0,CNRUp[UPCNInx+2]))==F)){TranP= log(GroupN-1)}
if( (CNRTC==2)&&( is.na(match(0,CNRUp[UPCNInx-1]))==F)){TranP= log(GroupN-1)}
DecPro=exp(GSN+WSN+RSN- GSO-WSO-RSO-log(1- LambdaF)+log(LambdaF)+TranP)
}
if(is.na(DecPro)==T){DecPro=0}
MAccp=rbinom(1,1,min(1, DecPro))

if (MAccp==1){
return(list("Type"= "Split","Result"="Success", "CNRUpN"=CNRUpN,"BrInxN" =BrInxN, "CNRO "= CNRUp, "BrInxO" =BrInx))
}
if (MAccp==0){
return(list("Type"="Split" ,"Result"="Fail" , "CNRUpN"=CNRUp,"BrInxN" =BrInx, "CNRO"= CNRUp,"BrInxO" =BrInx))
}
}
### Step 4.3Double Split 
DBSpt<-function(RD, RDsq, GroupN,MeanF, BetaF,AlphaF,KappaF,LambdaF, WeightIni,WeightUp, GroupPartLog,CNRUp,BrInx,GroupCNV,BrFix,CNRUpBr){


BrInxO=BrInx[-(match(BrFix,BrInx))] 

CNDiff=rep(0,(length(CNRUp)))
if(length(CNRUp)>1){
for (i in 2: (length(CNRUp))){
CNDiff[i]=BrInx[i]-BrInx[i-1] }
}
CNDiff[1]= BrInx[1]
DecPro =0
if(  length(which((CNDiff>2)&(CNRUp!=0)))  ==0    ){ DecPro =0}

if(  length(which((CNDiff>2)&(CNRUp!=0)))  >0    ){
if(  length(which((CNDiff>2)&(CNRUp!=0)))  ==1    ){ SptReg= which((CNDiff>2)&(CNRUp!=0))}
if(  length(which((CNDiff>2)&(CNRUp!=0)))  >1    ){
SptReg=sample(which((CNDiff>2)&(CNRUp!=0)) ,1) }

if(SptReg==1){CPoint=sort(sample(c(1: (BrInx[1]-1)),2))}
if(SptReg!=1) {CPoint=sort(sample(c( (BrInx[SptReg-1]+1): (BrInx[SptReg]-1) ),2))}
ChangInx=which((CNRUpBr[,1]<CPoint) &(CNRUpBr[,2]>CPoint))

BrInxON=sort(c(BrInxO,CPoint))
  BrInxN= sort(c(BrInx,CPoint))
UPCNInx= match(CPoint,BrInxN)

CNRUpN=rep(NA,length(CNRUp)+2)
if(UPCNInx[1]>1){CNRUpN[1: (UPCNInx[1]-1)]=CNRUp[1: (UPCNInx[1]-1)]}
if(UPCNInx[1]<(length(CNRUpN)-2)){CNRUpN[(UPCNInx[1]+3): (length(CNRUpN))]=CNRUp[(UPCNInx[1]+1): (length(CNRUp))]}

CNRUpN[UPCNInx[1]]=CNRUp[UPCNInx[1]]
CNRUpN[UPCNInx[1]+2]=CNRUp[UPCNInx[1]]
CNRUpN[UPCNInx[1]+1]=sample(GroupCNV[-(which ( (GroupCNV== CNRUp[UPCNInx[1]])))] ,1)
WSN=WeightPartL(WeightIni,WeightUp, CNRUpN,GroupN, CountMatrix= as.matrix(CountC(CNRUpN,GroupN)))

WSO=WeightPartL(WeightIni,WeightUp, CNRUp,GroupN, CountMatrix= as.matrix(CountC(CNRUp,GroupN)))

GSN=sum( GroupPartLog[c(CNRUpN[UPCNInx[1]], CNRUpN[UPCNInx[1]+1], CNRUpN[UPCNInx[1]+2] ) ])
GSO=(GroupPartLog[CNRUpN[UPCNInx[1]]])

if((UPCNInx[1])==1){WstartN1=1;WstartO=1}
if((UPCNInx[1])>1){WstartN1=BrInx[UPCNInx[1]-1]+1; WstartO= BrInxN[UPCNInx[1]-1]+1}
WstartN2= BrInxN[UPCNInx[1]]+1; WstartN3= BrInxN[UPCNInx[1]+1]+1;
WendO= BrInx[UPCNInx[1]]
WendN1= BrInxN[UPCNInx[1]]; WendN2= BrInxN[UPCNInx[1]+1]; WendN3= BrInxN[UPCNInx[1]+2]

RSO=RegionPartL (RD,RDsq, WstartO , WendO ,CNR= CNRUp[UPCNInx[1]],MeanF,BetaF,AlphaF,KappaF)

RSN=RegionPartL (RD,RDsq, WstartN1 , WendN1 ,CNR= CNRUpN[UPCNInx[1]],MeanF,BetaF,AlphaF,KappaF) +RegionPartL (RD,RDsq, WstartN2 , WendN2 ,CNR= CNRUpN[UPCNInx[1]+1],MeanF,BetaF,AlphaF,KappaF) +RegionPartL (RD,RDsq, WstartN3 , WendN3 ,CNR= CNRUpN[UPCNInx[1]+2],MeanF,BetaF,AlphaF,KappaF)


if(is.na(match(0,CNRUp[length(CNRUp)]))==F){TL=1}
if(is.na(match(0,CNRUp[length(CNRUp)]))==T){TL=0}
if(is.na(match(0,CNRUp[1]))==F){T1=1}
if(is.na(match(0,CNRUp[1]))==T){T1=0}

TranP =(length(RD)-2-rowSums(CountC(CNRUpN,GroupN))[GroupN+1]-length(which(is.na(RD)==T))+T1+TL)* (GroupN-1)/4
DecPro=exp(GSN+WSN+RSN- GSO-WSO-RSO-2*log(1- LambdaF)+2*log(LambdaF)+log(TranP))
}
if(is.na(DecPro)==T){DecPro=0}
MAccp=rbinom(1,1,min(1, DecPro))

if (MAccp==1){
return(list("Type"= "Trifid","Result"="Success", "CNRUpN"=CNRUpN,"BrInxN" =BrInxN, "CNRO "= CNRUp, "BrInxO" =BrInx))
}
if (MAccp==0){
return(list("Type"="Trifid" ,"Result"="Fail" , "CNRUpN"=CNRUp,"BrInxN" =BrInx, "CNRO"= CNRUp,"BrInxO" =BrInx))
}

}





















Cbound<-function(RD, RDsq, GroupN,MeanF, BetaF,AlphaF,KappaF,LambdaF, WeightIni,WeightUp, GroupPartLog,CNRUp,BrInx,GroupCNV,BrFix,CNRUpBr){

BrInxO=BrInx[-(match(BrFix,BrInx))] 

CNRUpN=CNRUp
if(length(BrInxO)==0){ DecPro =0}
if(length(BrInxO)!=0){

if( length(BrInxO)==1){CBD=BrInxO}
if(length(BrInxO)>1){CBD=sample(BrInxO,1)}

CP=match(CBD,BrInx)

if(CBD==1){LR=1}
if(CBD==length(RD)){LR=(-1)}
if(  (CBD!=1)& (CBD!=length(RD)) ){LR=sample(c(-1,1),1)}

ChangInx=which((CNRUpBr[,1]<CBD) &(CNRUpBr[,2]>CBD))


if (is.na(match(CBD+LR,BrInx))==F){DecPro =0}
if (is.na(match(CBD+LR,BrInx))==T){

CBN=CBD+LR
BrInxN=BrInx; BrInxN[CP]=CBN
if(CP==1){WstartN1=1;WstartO1=1}
if(CP>1){WstartN1=BrInx[CP-1]+1; WstartO1=BrInx[CP-1]+1}
WendN1=CBN
WendO1=CBD
WstartN2=CBN+1
WstartO2=CBD+1
WendN2=BrInx[CP+1]
WendO2=BrInx[CP+1]


RSN=RegionPartL (RD,RDsq, WstartN1 , WendN1 ,CNR= CNRUpN[CP],MeanF,BetaF,AlphaF,KappaF) +RegionPartL (RD,RDsq, WstartN2 , WendN2 ,CNR= CNRUpN[CP+1],MeanF,BetaF,AlphaF,KappaF)

RSO=RegionPartL (RD,RDsq, WstartO1 , WendO1 ,CNR= CNRUpN[CP],MeanF,BetaF,AlphaF,KappaF) +RegionPartL (RD,RDsq, WstartO2 , WendO2 ,CNR= CNRUpN[CP+1],MeanF,BetaF,AlphaF,KappaF)
DecPro=exp(RSN-RSO)
 DecPro
}

}
#if(is.na(DecPro)==T){DecPro=0}
MAccp=rbinom(1,1,min(1, DecPro))

if (MAccp==1){
return(list("Type"= "BoundaryC","Result"="Success", "CNRUpN"=CNRUpN,"BrInxN" =BrInxN, "CNRO"= CNRUp, "BrInxO" =BrInx))
}
if (MAccp==0){
return(list("Type"="BoundaryC" ,"Result"="Fail" , "CNRUpN"=CNRUp,"BrInxN" =BrInx, "CNRO"= CNRUp,"BrInxO" =BrInx))
}


}








CountC=function(CNR,GroupN=5){

CountMatrix=matrix(0,GroupN+1,GroupN)
CountInx=which(CNR!=0)
 
if(is.na(match(1, CountInx))==T) {
for(j in 1:GroupN){          
       for(i in 1:GroupN){
CountMatrix[i,j]=
length(which((CNR[CountInx-1]==i)&(CNR[CountInx]==j)))
        }
CountMatrix[GroupN+1,j]= length(which((CNR[CountInx-1]==0)&(CNR[CountInx]==j)))
}
}
if(is.na(match(1, CountInx))==F) {
   CNR=c(0,CNR)
CountInx=which(CNR!=0)
for(j in 1:GroupN){          
       for(i in 1:GroupN){
CountMatrix[i,j]=
length(which((CNR[CountInx-1]==i)&(CNR[CountInx]==j)))
        }
CountMatrix[GroupN+1,j]= length(which((CNR[CountInx-1]==0)&(CNR[CountInx]==j)))
}
}

diag(CountMatrix)=0
return(CountMatrix)


}




WeightPartL=function(WeightIni,WeightUp, CNR,GroupN, CountMatrix){

WeightUp=rbind(WeightUp, WeightIni)
ToUp=WeightUp+CountMatrix
SumV=rep(0,GroupN+1)

for(i in 1: (GroupN)){
   WM=WeightUp[i,]; WM=WM[-i]
   TM=ToUp[i,] ; TM=TM[-i]
   SumV[i]= sum(ApLogGamma(TM))- sum(ApLogGamma(WM))
}
   WML=WeightUp[(GroupN+1),]
   TML=ToUp[(GroupN+1),] 
   SumV[(GroupN+1)]= sum(ApLogGamma(TML))- sum(ApLogGamma(WML))

WeightPartLog=sum(SumV)+sum(-ApLogGamma((rowSums(CountMatrix)+1 )))


return(WeightPartLog)

}











RegionPartL=function(RD,RDsq, Wstart,Wend,CNR,MeanF,BetaF,AlphaF,KappaF)
{

if((Wend>=Wstart)&(Wstart>0)&(Wend>0)){
Leng=Wend-Wstart+1
SRD=sum(RD[Wstart: Wend], na.rm=T)
SRDsq=sum(RDsq[Wstart: Wend],na.rm=T)
}
if((Wend<Wstart)| (Wstart<1)|(Wend<1) ){
Leng=0
SRD=0
SRDsq=0
}
RegionPartLog=ApLogGamma(AlphaF[CNR]+ 0.5*Leng)-(0.5)*log(KappaF[CNR]+ Leng)-(AlphaF[CNR]+0.5*Leng)*log(BetaF[CNR]+(KappaF[CNR]*MeanF[CNR]*MeanF[CNR]/2)+(SRDsq/2)-((KappaF[CNR]*MeanF[CNR]+SRD)^2)/(2*(KappaF[CNR]+ Leng)))

return(RegionPartLog)
}





















ApLogGamma=function(x){

if(length(x)==1){
if(x>1){
LogGamma=0.5*log(2*pi*(x-1))+(x-1)*((log(x-1))-1)
}
if((x>0)&&(x<=1)){
LogGamma=log(gamma(x))
}
if(x==0){
LogGamma=0
}
return(LogGamma)
}
if(length(x)>1){
LogGamma=rep(0, length(x))
for(i in 1:length(x)){
if(x[i]>1){
LogGamma[i]=0.5*log(2*pi*(x[i]-1))+(x[i]-1)*((log(x[i]-1))-1)
}
if((x[i]>0)&&(x[i]<=1)){
LogGamma[i]=log(gamma(x[i]))
}
if(x[i]==0){
LogGamma[i]=0
}
}
return(LogGamma)
}

}


CNRWin=function(RD,CNR,BrInx){

CNRW=length(RD)
if(BrInx[1]!=1){
CNRW[1: (BrInx[1])]=CNR[1]
}
if(BrInx[1]<=1){
CNRW[1]=CNR[1]
print(BrInx)
}

if(length(CNR)>1){
for( i in 2: (length(CNR))){
CNRW[(BrInx[i-1]+1): (BrInx[i])]=CNR[i]
}}
return(CNRW)

}























WeiPartLogCountC=function(CNR,GroupN=5, WeightUp,WeightIni){

CountMatrix=matrix(0,GroupN+1,GroupN)
CountInx=which(CNR!=0)
 
if(is.na(match(1, CountInx))==T) {
for(j in 1:GroupN){          
       for(i in 1:GroupN){
CountMatrix[i,j]=
length(which((CNR[CountInx-1]==i)&(CNR[CountInx]==j)))
        }
CountMatrix[GroupN+1,j]= length(which((CNR[CountInx-1]==0)&(CNR[CountInx]==j)))
}
}
if(is.na(match(1, CountInx))==F) {
   CNR=c(0,CNR)
CountInx=which(CNR!=0)
for(j in 1:GroupN){          
       for(i in 1:GroupN){
CountMatrix[i,j]=
length(which((CNR[CountInx-1]==i)&(CNR[CountInx]==j)))
        }
CountMatrix[GroupN+1,j]= length(which((CNR[CountInx-1]==0)&(CNR[CountInx]==j)))
}
}

WeightUp=rbind(WeightUp,WeightIni)
ToUp=WeightUp+CountMatrix
SumV=rep(0,GroupN+1)
for(i in 1: (GroupN)){
   WM=WeightUp[i,]; WM=WM[-i]
   TM=ToUp[i,] ; TM=TM[-i]
   SumV[i]= sum(ApLogGamma(TM))- sum(ApLogGamma(WM))
}
   WML=WeightUp[(GroupN+1),]
   TML=ToUp[(GroupN+1),] 
   SumV[(GroupN+1)]= sum(ApLogGamma(TML))- sum(ApLogGamma(WML))

WeightPartLog=sum(SumV)+sum(-ApLogGamma((rowSums(CountMatrix)+1 )))

return(WeightPartLog)

}






GroupPartLogS=function(CNR, GroupPartLog){

GroupPartLogSum=sum(GroupPartLog[CNR])
return(GroupPartLogSum)
}








