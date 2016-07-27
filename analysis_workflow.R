# data analysis workflow for MCF
# written by Felix J. Hartmann


# set working directory (wd) to where the files are located
# files should only contain events of interest (eg live, single, CD45+)
setwd("/Users/hartmann/Desktop/data/")



# load (or install) all neccessary packages
source("http://bioconductor.org/biocLite.R")
if (!require(flowCore)) {biocLite("flowCore")} 
if (!require(gplots)) {biocLite("gplots")} 
if (!require(ggplot2)) {biocLite("ggplot2")} 
if (!require(Rtsne)) {biocLite("Rtsne")}
if (!require(FlowSOM)) {biocLite("FlowSOM")}
if (!require(RColorBrewer)) {biocLite("RColorBrewer")}
if (!require(gdata)) {biocLite("gdata")}
if (!require(reshape2)) {biocLite("reshape2")}
if (!require(matlab)) {biocLite("matlab")}
if (!require(grid)) {biocLite("grid")}
if (!require(Rmisc)) {biocLite("Rmisc")}



# load the fcs-files into a flowSet (fs)
fs <- read.flowSet(path = getwd(), pattern = "*.fcs", alter.names = T)
fs



# combine the flowset into a single flowframe with an additional gate.source channel
dim <- fsApply(fs, dim)
dim <- as.numeric(dim[,1])
dim
gate.source <- as.vector(x = NULL)
for(i in 1:length(dim)) {temp.source <- rep(i, dim[i])
  gate.source <- c(gate.source, temp.source)}



# combine data into a matrix
data <- fsApply(fs, exprs)
data <- cbind(data, gate.source)
head(data)
dim(data)



# rename channels and select clustering channels
# maybe here implement substring
(channel.names <- make.names(as.vector(parameters(fs[[1]])@data$desc), allow_ = F))
(clustering.channels <- channel.names[c(8,11,14,15,19,21,25,26,30,31:34,37)])
colnames(data) <- c(channel.names, "gate.source")



# do arcsinh transformation only for the clustering.channels
asinh_scale <- 5
data[,clustering.channels] <- asinh(data[,clustering.channels] / asinh_scale)



# percentile normalize the data (optional)
quantile_value <- 0.9999
percentile.vector <- apply(data[,clustering.channels], 2, function(x) quantile(x, quantile_value, names = F))
percentile.vector
data[,clustering.channels] <- t(t(data[,clustering.channels]) / as.numeric(percentile.vector))
head(data)
dim(data)



#write.FCS(x = flowFrame(exprs = data), filename = "transformed.fcs")
data.df <- data.frame(data)
dl <- split(data.df, data.df$gate.source)
sampleNames(fs)
dl1 <- dl$`13`
data <- data.matrix(dl1)
head(data)
dim(data)



# subsample data (optional) for t-SNE
n_sub <- 20000
n <- nrow(data)
set.seed(123)
ix <- sample(1:n, n_sub) #if you use subsampled data
#ix <- 1:nrow(data) #if you dont subsample the data
#data.sub <- data[ix,]
data_rtsne <- data.matrix(data[ix,clustering.channels])



# run FlowSOM (with set.seed for reproducibility)
set.seed(123)
out_fSOM <- FlowSOM::ReadInput(flowFrame(exprs = data[ix,], desc = list(FIL = 1)), transform = FALSE, scale = FALSE, compensate = FALSE)
out_fSOM <- FlowSOM::BuildSOM(out_fSOM, colsToUse = clustering.channels)
out_fSOM <- FlowSOM::BuildMST(out_fSOM)
labels <- out_fSOM$map$mapping[,1]


# define the range of k.values for metaclusters
min.cluster <- 3
max.cluster  <- 20


# try the suggested automatic metaclustering method for a hint for k
auto_meta <- MetaClustering(out_fSOM$map$codes, method = "metaClustering_consensus", max = max.cluster)
max(auto_meta)


# do a manual metaclustering for all values up to max
meta_results <- data.frame(as.factor(data[ix,"gate.source"]))
for (i in min.cluster:max.cluster) {
  set.seed(123)
  out_meta <- FlowSOM::metaClustering_consensus(out_fSOM$map$codes, k = i)
  meta_results <- cbind(meta_results, as.factor(out_meta[labels]))}

meta_results <- meta_results[,2:ncol(meta_results)]
colnames(meta_results) <- paste("k.", min.cluster:max.cluster, sep = "")
head(meta_results)



# run bh SNE
set.seed(123)
out_rtsne <- Rtsne(data_rtsne, dims = 2, perplexity = 50, theta = 0.5, max_iter = 1000, verbose = T, pca = F, check_duplicates = F)
#save(out_rtsne, file = "out_rtsne.robject", compress=F)
#save.image(file = "exported/image", compress = F)


# load and prepare metadata
#md2 <- read.xls("md2.xlsx", sheet = 1, header = T, verbose = F)
#md2
gate.df <- as.data.frame(data[ix,"gate.source"])
colnames(gate.df) <- "gate.source"
head(gate.df)
gate.df$cell.id <-  as.factor(ix) #as.factor(row.names(gate.df))
#gate.df <- merge(gate.df, md2, by = "gate.source")


# prepare the tSNE data
tsne <- as.data.frame(out_rtsne$Y)
colnames(tsne) <- c("tSNE1", "tSNE2")
tsne$cell.id <- as.factor(ix)
head(tsne)


# prepare the expression data
data.ix.df <- data.frame(data[ix,])
data.ix.df$cell.id <- as.factor(ix)
data.melt <- melt(data.ix.df, variable.name = "antigen", value.name = "expression", id.vars = c("cell.id", "gate.source"))
joined.expr <- merge(data.melt, tsne, by = "cell.id")
joined.expr <- merge(joined.expr, gate.df, by = "cell.id")


# prepare the metaclustering data
meta.sub <- meta_results
meta.sub$cell.id <- as.factor(ix)
meta.melt <- melt(meta.sub, id.vars = "cell.id", variable.name = "k.value", value.name = "cluster.assigment")
meta.melt$cluster.assigment <- factor(meta.melt$cluster.assigment, levels = c(1:length(meta.melt$cluster.assigment)), ordered = T)
joined.meta <- merge(meta.melt, tsne, by = "cell.id")
joined.meta <- merge(joined.meta, gate.df, by = "cell.id")
#joined.meta$genotype <- factor(joined.meta$genotype, levels(joined.meta$genotype)[c(2,1)])


# define the plotting stuff
theme_tsne <-  theme(panel.margin = unit(1.3, "lines"), 
                     strip.text = element_text(size = rel(0.8)), 
                     axis.ticks.length = unit(0.3, "lines"),
                     axis.text = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.background = element_rect(color="white"),
                     strip.background = element_blank(), 
                     panel.background = element_blank(),
                     axis.ticks = element_blank(),
                     aspect.ratio = 1)
db1 <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"))
red.blue <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(9))


# plot tSNE black
t1 <- ggplot(tsne, aes(x = tSNE1, y = tSNE2)) +
  geom_point(size = 1) +
  coord_fixed(ratio = 1) +
  ggtitle("tSNE map black") +
  theme_tsne
t1



# plot tSNEs with expression overlayed
t2 <- ggplot(droplevels(subset(joined.expr, antigen %in% clustering.channels)), aes(x = tSNE1, y = tSNE2, color = expression)) +
  geom_point(size = 0.025) +
  coord_fixed(ratio = 1) +
  scale_colour_gradientn(colours = jet.colors(100), limits = c(0,1)) +
  facet_wrap(~ antigen, ncol = 5, scales = "free") +
  ggtitle("tSNE map with expression values") +
  theme_tsne + theme(strip.text = element_text(size = rel(0.8)))
t2



# plot tSNEs with k.values overlayed
t3 <- ggplot(joined.meta, aes(x = tSNE1, y = tSNE2, color = cluster.assigment)) +
  geom_point(size = 0.025) +
  coord_fixed(ratio = 1) +
  scale_colour_manual(name = NULL,values = db1) +
  facet_wrap(~ k.value, ncol = 6, scales = "free") + 
  ggtitle("tSNE map with different k.values") +
  theme_tsne
t3



# save
save.name <- "irina_2nd"
ggsave(filename = paste(save.name, "tnse_black.png", sep=""), plot = t1, scale = 1) #, width = 6, height = 6, units = c("in")) 
ggsave(filename = paste(save.name, "tnse_exprs.png", sep=""), plot = t2, scale = 1) #, width = 8, height = 5, units = c("in")) 
ggsave(filename = paste(save.name, "tnse_kvalues.png", sep=""), plot = t3, scale = 1) #, width = 8, height = 6, units = c("in")) 



# make combined expression matrix
ce <- data.matrix(cbind(data[ix,], meta_results))
head(ce)
dim(ce)



# go through all k.values and make mean matrices
plot.list <- list()
for(ia in min.cluster:max.cluster) {
    cc <- paste("k.", ia, sep = "")
    cluster_num <- ia
    return_mat <- matrix(, nrow = cluster_num, ncol = length(clustering.channels))
    colnames(return_mat) <- clustering.channels
    for(ib in 1:cluster_num) {
        temp_mat <- ce[which(ce[,cc] == ib),clustering.channels]
        return_mat[ib,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))}
    rownames(return_mat) <- paste("cluster", 1:cluster_num, sep="")
    plot.list[[ia-min.cluster+1]] <- return_mat
}    


  
# make a list of ggplot2 heatmap items
pdf(paste(save.name, "heatmaps.pdf", sep=""), paper = "a4r", useDingbats = T, width = 8, height = 5)
temp <- lapply(1:length(plot.list), function(i) {
      heatmap.2(plot.list[[i]][,colnames(return_mat) %in% clustering.channels], 
      scale = "none",
      main = names(plot.list)[i],
      Colv = F,
      Rowv = T,
      dendrogram =  "row",
      trace = "none",
      col = red.blue,
      breaks = seq(0, 1, by = 0.1111111111))
})
dev.off()


