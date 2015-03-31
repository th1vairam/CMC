#script to grab SVA CMC data

require(synapseClient)
synapseLogin()

SVAGeneExpressionSummary <- synGet('syn2757151')
cmcSVAGeneSummary <- read.delim(SVAGeneExpressionSummary@filePath)

SVATFExpressionSummary <- synGet('syn2757153')
cmcSVATFSummary <- read.delim(SVATFExpressionSummary@filePath)

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

cmcExpressionFreezeSCZ <- scale(data.matrix(cmcExpressionFreeze[,colnames(cmcExpressionFreeze)%in%scz_samples]))

#extract schizophrenia cases only
#