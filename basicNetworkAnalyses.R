require(synapseClient)
synapseLogin()

res <- synQuery('select name,id,method from file where projectId==\'syn3455058\' and matrixType==\'sparse\'')

loadData <- function(syn){
  synMeta<-synGet(syn)
  load(synMeta@filePath)
  return(network)
}

SVAGeneExpressionSummary <- synGet('syn2757151')
cmcSVAGeneSummary <- read.delim(SVAGeneExpressionSummary@filePath,stringsAsFactors = FALSE)

SVATFExpressionSummary <- synGet('syn2757153')
cmcSVATFSummary <- read.delim(SVATFExpressionSummary@filePath, stringsAsFactors = FALSE)



require(metaNet)
networks <- sapply(res$file.id,loadData)
names(networks) <- res$file.method
#netStat
netStats <- simplify2array(sapply(networks,networkStatistics))
netStats <- cbind(rownames(netStats),c(cmcSVATFSummary$MAPPED_genes, cmcSVAGeneSummary$MAPPED_genes),netStats)
colnames(netStats) <- c('EnsemblId','GeneId',res$file.method)
netStats <- data.frame(netStats,stringsAsFactors = FALSE)
netStats[,3:9]<-apply(netStats[,3:9],2,as.numeric)
netStats[,3:9] <- apply(netStats[,3:9],2,scale)

tcresult <- as.tableColumns(netStats)
cols <- tcresult$tableColumns
cols[[2]]$maximumSize <- as.integer(80)
fileHandleId <- tcresult$fileHandleId

projectId <- "syn3455058"

schema<-TableSchema(name="SVA Schizophrenia Network Statistics", parent=projectId, columns=cols)
table<-Table(schema, fileHandleId)
table<-synStore(table, retrieveData=TRUE)


geneRanks <- apply(netStats[,3:9],2,order,decreasing=T)
geneRanksName <- apply(geneRanks,2,function(x,g) return(g[x]),netStats$GeneId)

enrichRes <- apply(geneRanksName,2,enrichmentPath,schizophreniaHits)
enrichResPval <- sapply(enrichRes,function(x) return(-log10(as.numeric(x$pval))))


matplot(1:5000,enrichResPval[1:5000,],'l',lwd=3)
legend('topleft',colnames(enrichResPval),col=c(1:6,1),lty=1:7,lwd=2,cex=.5)


#get differentially expressed genes

deSyn <- synGet('syn3493933')
deSynFile <- read.delim(deSyn@filePath,stringsAsFactors=F)

targetListDE <- as.character(deSynFile$MAPPED_genes[1:488])
targetListDEens <- as.character(deSynFile$genes[1:488])
enrichResDe <- apply(geneRanksName,2,enrichmentPath,targetListDE)
enrichResPvalDE <- sapply(enrichResDe,function(x) return(-log10(as.numeric(x$pval))))

matplot(1:16423,enrichResPvalDE,'l',lwd=3)
legend('topleft',colnames(enrichResPvalDE),col=c(1:6,1),lty=1:7,lwd=2,cex=.5)


def <- synTableQuery('SELECT * FROM syn3546459')
def@values$isGWAS <- netStats$GeneId%in%schizophreniaHits
def@values$isDE <- netStats$EnsemblId%in%targetListDEens
def <- synStore(def)


def <- synTableQuery('SELECT * FROM syn3546459')
colnames(def@values)
defMethod <- def@values[,c('sparrow','wgcna','aracne','genie3')]
a3 <- rowSums(defMethod)
plot(a3)
def@values <- cbind(def@values,a3)
