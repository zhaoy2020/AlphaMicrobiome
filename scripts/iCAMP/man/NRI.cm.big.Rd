\name{NRI.cm.big}
\alias{NRI.cm.big}
\title{
Calculate net relatedness index (NRI) under multiple metacommunities for big data
}
\description{
Calculate net relatedness index (NRI) or other index of null model significance test based on mean pairwise distance (MPD) by parallel computing, for big dataset. This function can deal with local communities under different metacommunities (regional pools).
}
\usage{
NRI.cm.big(comm, pd.desc, pd.spname, pd.wd,
           meta.group = NULL, meta.spool = NULL,
           nworker = 4, memo.size.GB = 50,
           weighted = c(TRUE, FALSE), check.name = TRUE,
           rand = 1000, output.MPD = c(FALSE, TRUE),
           silent = FALSE, sig.index = c("SES", "NRI", "Confidence", "RC", "all"),
           sp.limit = 10000)
}
\arguments{
  \item{comm}{community data matrix. rownames are sample names. colnames are species names.}
  \item{pd.desc}{character, the name to describe bigmemory file of phylogenetic distance matrix, default is "pd.desc".}
  \item{pd.spname}{vector, the OTU ids (species names) in exactly the same order as the phylogenetic matrix rows or columns}
  \item{pd.wd}{the path of the folder saving the phylogenetic distance matrix.}
  \item{meta.group}{matrix or data.frame, a one-column (n x 1) matrix indicating which metacommunity each sample belongs to. rownames are sample IDs. first column is metacommunity names. Such that different samples can belong to different metacommunities. If input a n x m matrix, only the first column is used. NULL means all samples belong to the same metacommunity. Default is NULL, means all samples from the same metacommunity.}
  \item{meta.spool}{a list object, each element is a character vector listing all taxa IDs in a metacommunity. The names of the elements indicate metacommunity names, which should be the same as the metacommunity names in meta.group. Default is NULL, means to use the observed taxa in comm across samples within the same metacommunity that is defined by meta.group.}
  \item{nworker}{for parallel computing. Either a character vector of host names on which to run the worker copies of R, or a positive integer (in which case that number of copies is run on localhost). default is 4, means 4 threads will be run.}
  \item{memo.size.GB}{numeric, to set the memory size as you need, so that calculation of large tree will not be limited by physical memory. unit is Gb. default is 50Gb.}
  \item{weighted}{Logic, consider abundances or not (just presence/absence). default is TRUE.}
  \item{check.name}{Logic, whether to check the species names in comm and dis. default is TRUE.}
  \item{rand}{integer, randomization times. default is 1000.}
  \item{output.MPD}{Logic, whether to output observed MPD, so that you do not need to calculate observed MNTD alone. default is FALSE.}
  \item{silent}{Logic, if FALSE, some messages will be showed during calculation. Default is FALSE.}
  \item{sig.index}{character, the index for null model significance test. SES or NRI, standard effect size, i.e. net relatedness index (NRI); Confidence, percentage of null values less extreme than the observed value, i.e. non-parametric one-side confidence level; RC, modified Raup-Crick index (RC) based on MPD, i.e. count the number of null MPD lower than observed MPD plus a half of the number of null MPD equal to observed MPD, to get alpha, then calculate MPD-based RC as (2 x alpha - 1); all, output all the three indexes. default is SES. If input a vector, only the first element will be used.}
  \item{sp.limit}{numeric, the limitation of species number in each round of the looping inside this funciton, default is 10000, generally fine for most computer, if returns error, set lower value.}
}
\details{
This function is particularly designed for samples from different metacommunities. The null model "taxa shuffle" will be done under different metacommunities, separately (and independently). All other details are the same as the function NRI.p.
}
\value{
Output can be a data.frame with each row representing a sample and only one column of index values, or a list of several data.frame objects.
\item{SES}{output if sig.index is Confidence or all, a data.frame with NRI value for each sample.}
\item{Confidence}{output if sig.index is SES or all, a data.frame showing confidence level based on MPD for each sample.}
\item{RC}{output if sig.index is RC or all, a data.frame showing RC based on MPD for each sample.}
\item{MPD.obs}{output if output.MPD is TRUE, a data.frame showing observed MPD for each sample.}
\item{MPD.rand}{output if output.MPD is TRUE, a matrix showing all null MPD values.}
}
\references{
Webb CO, Ackerly DD, and Kembel SW. 2008. Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Bioinformatics 18:2098-2100

Kembel, S.W. (2009). Disentangling niche and neutral influences on community assembly: assessing the performance of community phylogenetic structure tests. Ecol Lett, 12, 949-960.

Chase, J.M., Kraft, N.J.B., Smith, K.G., Vellend, M. & Inouye, B.D. (2011). Using null models to disentangle variation in community dissimilarity from variation in alpha-diversity. Ecosphere, 2, 1-11.

Ning, D., Yuan, M., Wu, L., Zhang, Y., Guo, X., Zhou, X. et al. (2020). A quantitative framework reveals ecological drivers of grassland microbial community assembly in response to warming. Nature Communications, 11, 4717.
}
\author{
Daliang Ning
}
\note{
Version 2: 2023.12.23, fix a bug when input comm is data.frame.
Version 1: 2023.9.4
}
\seealso{
\code{\link{NRI.p}}, \code{\link{NRI.cm}}
}
\examples{
data("example.data")
comm=example.data$comm
pd=example.data$pd

# in this example, 10 samples from one metacommunity,
# the other 10 samples from another metacommunity.
meta.group=data.frame(meta.com=c(rep("meta1",10),rep("meta2",10)))
rownames(meta.group)=rownames(comm)

# since it needs to save some file to a certain folder,
# the following code is set as 'not test'.
# but you may test the code on your computer
# after change the folder path for 'save.wd'.
\donttest{
wd0=getwd()
save.wd=paste0(tempdir(),"/pdbig.nricm.big")
nworker=2 # parallel computing thread number.
pd.big=pdist.big(tree = tree, wd=save.wd, nworker = nworker)
rand.time=20 # usually use 1000 for real data.
sigmpd=NRI.cm.big(comm=comm, pd.desc = pd.big$pd.file,
               pd.spname = pd.big$tip.label, pd.wd = pd.big$pd.wd,
               meta.group=meta.group, nworker=nworker,
               weighted=TRUE, rand=rand.time,
               sig.index="all")
NRI=sigmpd$SES
CMPD=sigmpd$Confidence
RCMPD=sigmpd$RC
}
}
\keyword{ phylogenetic }
