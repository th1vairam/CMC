#script to grab SVA CMC data

sczOutSva <- as.character(commandArgs(TRUE)[[1]])
controlOutSva <- as.character(commandArgs(TRUE)[[2]])
sczOut <- as.character(commandArgs(TRUE)[[3]])
controlOut <- as.character(commandArgs(TRUE)[[4]])

require(synapseClient)
synapseLogin()


makeDataMatrices <- function(geneSyn,tfSyn,sczFile,contFile,bpFile){
  SVAGeneExpression <- synGet('syn2757147')
  cmcSVAGeneExpression <- read.delim(SVAGeneExpression@filePath)
  
  SVATFExpression <- synGet('syn2757149')
  cmcSVATFExpression <- read.delim(SVATFExpression@filePath)
  
  #combine gene expression data
  cmcExpression <- rbind(cmcSVATFExpression,cmcSVAGeneExpression)
  
  #download metadata
  syncmcMetaData <- synGet('syn2299154')
  metadata <- read.delim(syncmcMetaData@filePath,sep=',',stringsAsFactors=FALSE)
  
  metadata_freeze = metadata[-which(metadata$DLPFC_RNA_report..Exclude. == 1), ]
  
  cmcExpressionFreeze <- cmcExpression[,which(colnames(cmcExpression)%in%as.character(metadata_freeze$DLPFC_RNA_isolation..Sample.RNA.ID))]
  
  scz_samples <- metadata_freeze$DLPFC_RNA_isolation..Sample.RNA.ID[which(metadata_freeze$Dx == "SCZ")]
  
  control_samples <- metadata_freeze$DLPFC_RNA_isolation..Sample.RNA.ID[which(metadata_freeze$Dx == "Control")]
  
  #extract schizophrenia cases only
  cmcExpressionFreezeSCZ <- scale(t(data.matrix(cmcExpressionFreeze[,colnames(cmcExpressionFreeze)%in%scz_samples])))
  
  write.csv(cmcExpressionFreezeSCZ,file=file,quote=FALSE)
}