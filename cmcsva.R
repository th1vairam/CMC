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

#