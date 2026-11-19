\name{taxa.binphy.big.cm}
\alias{taxa.binphy.big.cm}
\title{
Phylogenetic binning based on phylogenetic tree for multi-metacommunity data
}
\description{
Phylogenetic binning for iCAMP analysis given multiple metacommunities, to ensure the taxa number of each bin and the detected taxa number of each none-zero bin in each metacommunity reach the minimal requirements. To handle large phylogenetic tree, phylogenetic distance matrix should be calculated and saved using the package 'bigmemory' in advance. 
}
\usage{
taxa.binphy.big.cm(tree, comm = NULL, meta.group = NULL, meta.com = NULL,
                   meta.frequency = NULL, meta.ab = NULL,
                   pd.desc, pd.spname, pd.wd,
                   outgroup.tip = NA, outgroup.rm = TRUE,
                   d.cut = NULL, ds = 0.2, bin.size.limit = 24, 
                   bin.size.limit.eachmeta = 6, nworker = 4,
                   d.cut.method = c("maxpd", "maxdroot"))
}
\arguments{
  \item{tree}{phylogenetic tree, an object of class "phylo".}
  \item{comm}{matrix or data.frame, community data, each row is a sample or site, each colname is a taxon (a species or OTU or ASV), thus rownames should be sample IDs, colnames should be taxa IDs.}
  \item{meta.group}{matrix or data.frame, a one-column (n x 1) matrix indicating which metacommunity each sample belongs to. Rownames are sample IDs. The first column is metacommunity names. Such that different samples can belong to different metacommunities. If input a n x m matrix, only the first column is used. NULL means all samples belong to the same metacommunity. Default is NULL, means all samples from the same metacommunity.}
  \item{meta.com}{a list object, each element is a matrix or data.frame to define abundance (or relative abundance) of taxa in a metacommunity (regional pool). The element names indicate metacommunity names, which should be consistent with the metacommunity names defined in meta.group. If there is only one metacommunity, meta.com can be a matrix or data.frame to define taxa abundance (or relative abundance) in the metacommunity. Default is NULL, means to calculate metacommunity structure from comm according to metacommunities defined in meta.group.}
  \item{meta.frequency}{a list object, each element is a matrix or data.frame to define abundance (or relative abundance) of taxa in a metacommunity (regional pool). The element names indicate metacommunity names, which should be consistent with the metacommunity names defined in meta.group. If there is only one metacommunity, meta.com can be a matrix or data.frame to define taxa abundance (or relative abundance) in the metacommunity. Default is NULL, means to calculate metacommunity structure from comm according to metacommunities defined in meta.group, or calculate from meta.com.}
  \item{meta.ab}{matrix or data.frame, each column represents a taxon, each row represents a metacommunity (regional pool), to define the aubndance (or relative abundance) of each taxon in each metacommunity. The rownames indicate metacommunity names, which should be the same as the metacommunity names in meta.group. Default setting is NULL, means to calculate meta.ab as average relative abundance of each taxon in comm across the samples within each metacommunity defined by meta.group, or calculate from meta.com.}
  \item{pd.desc}{the name of the file to hold the backingfile description of the phylogenetic distance matrix, it is usually "pd.desc" if using default setting in pdist.big function.}
  \item{pd.spname}{character vector, taxa id in the same rank as the big matrix of phylogenetic distances.}
  \item{pd.wd}{folder path, where the bigmemmory file of the phylogenetic distance matrix are saved.}
  \item{outgroup.tip}{a vector of tip names (i.e. OTU IDs) which is in totally different lineage from all other tips, thus can be used as outgroup to root the tree. For example, Archaeal OTUs may be set as outgroup tips when analyzing Bacterial OTUs. Default is NA, means no need to set outgroup tip.}
  \item{outgroup.rm}{logic, whether to remove the outgroup.tip after the tree is rooted. Default is TRUE.}
  \item{d.cut}{numeric, the distance from root to the truncating point of the tree.}
  \item{ds}{numeric, the general threshold of phylogenetic distance within which the phylogenetic signal is significant. default is 0.2.}
  \item{bin.size.limit}{integer, the minimal requirement of bin size (taxa numer in a bin). Default setting is 24.}
  \item{bin.size.limit.eachmeta}{integer, the minimal requirement of none-zero bin size (taxa numer in a bin) in each metacommunity. Default setting is 6.}
  \item{nworker}{integer, for parallel computing. Either a character vector of host names on which to run the worker copies of R, or a positive integer (in which case that number of copies is run on localhost). default is 4, means 4 threads will be run.}
  \item{d.cut.method}{character, to specify the method to calculate d.cut from ds. 'maxpd' means based on maximum phylogenetic distance, d.cut = (maxpd - ds)/2. 'maxdroot' means based on maximum distance to root, d.cut = maxdroot - (ds/2), which is preferred if the tree only has one edge from the root.}
}
\details{
The phylogenetic tree is truncated at a certain phylogenetic distance (as short as necessary) to the root (d.cut), by which all the rest connections between tips (taxa) are lower than a threshold. Within the threshold, phylogenetic signal is generally significant. The taxa derived from the same ancestor after the truncating point are grouped to the same strict bin. Then, each small bin is merged into the bin with the nearest relatives. This procedure is repeated until all merged bins have enough taxa (>= bin.size.limit) in each metacommunity. Bigmemory (Kane et al 2013) is used to deal with large datasets.
}
\value{
Output is a list.
\item{sp.bin}{matrix, rownames are taxa IDs; the first column is strict bin IDs; the second column indicates which strict bin the taxon is merged into; the third column is the final bin IDs.}
\item{bin.united.sp}{list, each element is a vector of taxa IDs, indicating the taxa in a final bin (after small bins are merged into nearest large bins).}
\item{bin.strict.sp}{list, each element is a vector of taxa ID(s), indicating the taxa in a strict bin (before small bins are merged into large bins).}
\item{state.strict}{matrix, status of each strict bin. bin.strict.id, the strict bin ID; bin.strict.taxa.num, taxa number in each strict bin; bin.pd.max, bin.pd.mean, and bin.pd.sd, the maximum, mean, and standard deviation of the pairwise phylogenetic distances in each strict bin; the last column bin.strict.taxa.num.metamin shows the minimal taxa number of the strict bin across different metacommunities.}
\item{state.united}{matrix, status of each final bin. bin.united.id.old, the ID of the largest strict bin in each final bin; bin.united.tax.num, taxa number in each final bin; bin.pd.max, bin.pd.mean, and bin.pd.sd, the maximum, mean, and standard deviation of the pairwise phylogenetic distances in each strict bin; the last column bin.strict.taxa.num.metamin shows the minimal taxa number of the final bin across different metacommunities.}
}
\references{
Ning, D., Yuan, M., Wu, L., Zhang, Y., Guo, X., Zhou, X. et al. (2020). A quantitative framework reveals ecological drivers of grassland microbial community assembly in response to warming. Nature Communications, 11, 4717.
}
\author{
Daliang Ning
}
\note{
Version 2: 2024.12.26, add bin.size.limit.eachmeta.
Version 1: 2024.12.22.
}
\seealso{
\code{\link{taxa.binphy.big}}, \code{\link{icamp.cm}}
}
\examples{
data("example.data")
comm=example.data$comm
tree=example.data$tree

# in this example, 10 samples from one metacommunity,
# the other 10 samples from another metacommunity.
meta.group=data.frame(meta.com=c(rep("meta1",10),rep("meta2",10)))
rownames(meta.group)=rownames(comm)

# since pd.big need to specify a certain folder,
# the following code is set as 'not test'.
# but you may test the code on your computer
# after change the folder path for 'save.wd'.
\donttest{
  wd0=getwd()
  save.wd=paste0(tempdir(),"/pdbig.taxa.binphy")
  # please change to the folder you want to save the big niche difference matrix.
  
  nworker=2 # parallel computing thread number
  pd.big=pdist.big(tree = tree, wd=save.wd, nworker = nworker)
  
  ds = 0.2 # setting can be changed to explore the best choice
  
  bin.size.limit = 5 # setting can be changed to explore the best choice.
  bin.size.limit.eachmeta = 3 # setting can be changed to explore the best choice.
  # here set as 5 just for the small example dataset.
  # For real data, usually try 12 to 48.
  
  pbin=taxa.binphy.big.cm(tree = tree, comm=comm, meta.group=meta.group,
                          pd.desc = pd.big$pd.file, pd.spname = pd.big$tip.label,
                          pd.wd = pd.big$pd.wd, ds = ds,
                          bin.size.limit = bin.size.limit,
                          bin.size.limit.eachmeta = bin.size.limit.eachmeta,
                          nworker = nworker)
  setwd(wd0)
}
}
\keyword{main}
