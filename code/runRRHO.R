# run RRHO comparing to brainSpan data

options(stringsAsFactors = FALSE)
args = commandArgs(trailingOnly = TRUE)

#load RRHO function
source("../code/TransitionMapping_function.R")

#load your_data data
input_dir <- args[[1]]
if (!dir.exists(input_dir)) {
  stop("input directory doesn't exist")
}
input_files <- list.files(input_dir, pattern = "csv", full.names = TRUE)
if (length(input_files) == 0) {
  stop("input file is empty")
}
test_data <- lapply(input_files, function(f) {
  read.csv(f, header = T)
})
names(test_data) <- gsub(".csv$","", basename(input_files))

output_dir <- args[[2]]
dir.create(output_dir, showWarnings = F,recursive = T)

# load brainSPan data
load("../data/brainSpan_pariedVoom_results.rdata")

stepsize <-  200

# run RRHO
for (j in c(1:length(brainSpanVoom))) {

  ref <- as.data.frame(brainSpanVoom[[j]])
  ref$geneID <- rownames(ref)

  for (i in c(1:length(test_data))) {
    name_i <- names(test_data)[i]
    name_j <- names(brainSpanVoom)[j]
    print(paste(name_i,"......",name_j))

    test <- as.data.frame(test_data[[i]])
    gene_to_keep <- intersect(ref$geneID, test$geneID)
    ref <- ref[ref$geneID %in% gene_to_keep,c("geneID","logFC")]
    test <- test[test$geneID %in% gene_to_keep,c("geneID","logFC")]

    if (j == 1 & i == 1) {
      hypermat.all = array(data = NA,
                           dim = c(length(test_data),
                                   length(brainSpanVoom),
                                   length(seq(1,length(ref[,1]),stepsize)),
                                   length(seq(1,length(test[,1]),stepsize))))
    }

    tmp <- RRHO(test,ref,stepsize = stepsize,labels = c(paste0("iPSC_DE_",name_i),paste0("brianSpan_DE_",name_j)),plots = F)
    hypermat.all[i,j,,] = tmp$hypermat

  }
}
save(hypermat.all,file = paste0(output_dir,"/hypermatAll.rdata"))

#load(file=paste0(output_dir,"hypermatAll.rdata"))

### Plot the output

comparisons <- names(brainSpanVoom)
nr <- length(unique(gsub("-.*","",comparisons)))
nc <- length(unique(gsub(".*-","",comparisons)))

jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
colormap = jet.colors(100)

minhypermat <-  min(hypermat.all,na.rm = TRUE)
maxhypermat <-  max(hypermat.all,na.rm = TRUE)
# for one comparison in the CIRMA data
for (k in c(1:length(test_data))) {
  transition <- names(test_data)[k]
  png(paste0(output_dir,"/",transition,"_brainSpan_RRHOMap.png"),width = 10, height = 10, units = "in", res = 300)
  par(mfrow = c(nr + 1, nc + 1),
    mar = c(0, 0, 1, 1),
    oma = c(0, 0, 1, 0))
  count = 1
  for (j in 1:(nr + 1)) {
    for (i in 1:(nc + 1)) {
      if (j < i) {
        ## Scale the output to the max across all hypermats so that the colorbar is the same for all hypermats
        image(hypermat.all[k, count, , ],
              xlab = '',
              ylab = '',
              axes = FALSE,
              col = colormap,
              zlim = c(minhypermat,maxhypermat))
        mtext(names(brainSpanVoom)[count],3,cex = 0.5, line = 0)
        count = count + 1
      } else {
        plot(0,type = "n", axes = F, ylab = '', xlab = '')
      }
    }
  }
  dev.off()
}
## Function to plot color bar
## Modified from http://www.colbyimaging.com/wiki/statistics/color-bars
png(paste0(output_dir,"/brainSpan_RRHOMap_color_bar.png",sep = ""), width = 2, height = 6, units = "in", res = 300)
color.bar(jet.colors(100),min = minhypermat, max = maxhypermat, nticks =  6, title = "-log10(Nominal P-value)")
dev.off()
