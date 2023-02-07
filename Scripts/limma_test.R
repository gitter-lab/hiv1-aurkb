#!/usr/bin/env Rscript

#Got help from here:
#http://www.biostat.jhsph.edu/~kkammers/software/eupa/R_guide.html 
#Guide above "is freely available as Free Software under the terms of the Free Software Foundation's GNU General Public License in source code form."

#Installs packages if not present and loads them
install.packages("librarian", quiet = TRUE, repos = "http://cran.us.r-project.org")
librarian::shelf(qvalue, limma, readr,quiet=TRUE)

#Load Data
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  stop("At least an input file needs to be provided.")
} else if (length(args)==1){
  args[2] = "output"
}
fName = args[1]

data <- read_csv(fName,col_names = TRUE,show_col_types=FALSE)

#Limma expects log-expression data
data[,] <- log(data[,])

cNames <- colnames(data)
design <- c()

for (i in 1:length(cNames)){
  #tmp <- unlist(strsplit(cNames[i],"_",fixed=TRUE))
  tmp <- substr(cNames[i], 1, nchar(cNames[i])-2)
  design[i] <- tmp[1]
}
fDesign = model.matrix(~factor(design))

fit <- lmFit(data,fDesign)
fit2 <- eBayes(fit)

outData <- data.frame(fit2$p.value[,2])
colnames(outData) <- c("pval")
outData$qval <- qvalue(outData$pval)$qvalues

outF <- unlist(strsplit(fName,"_",fixed=TRUE))
outF <- append(outF,"tmpLimma",after=length(outF))
outF <- paste(outF,sep="_",collapse="_")
print(outF)
write.csv(outData,file=outF,quote=FALSE,row.names=FALSE)
