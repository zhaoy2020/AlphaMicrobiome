\name{bmpd.big}
\alias{bmpd.big}
\title{
Beta mean pairwise distance (betaMPD) for big data
}
\description{
Calculates mean pairwise distance separating taxa in each pair of communities in a given community matrix, for big data.
}
\usage{
bmpd.big(comm, pd.desc = "pd.desc", pd.spname, pd.wd,
         spname.check = FALSE, abundance.weighted = TRUE,
         na.zero = TRUE, unit.sum = NULL)
}
\arguments{
  \item{comm}{matrix or data.frame, community data matrix, rownames are sample names, colnames are taxa ids.}
  \item{pd.desc}{character, the name to describe bigmemory file of phylogenetic distance matrix, default is "pd.desc".}
  \item{pd.spname}{vector, the OTU ids (species names) in exactly the same order as the phylogenetic matrix rows or columns}
  \item{pd.wd}{the path of the folder saving the phylogenetic distance matrix.}
  \item{spname.check}{logic, whether to check the OTU ids (species names) in community matrix and phylogenetic distance matrix are the same.}
  \item{abundance.weighted}{logic, whether weighted by species abundance, default is TRUE, means weighted.}
  \item{na.zero}{logic. when the sum of a row (a sample) is zero in community data matrix, the relative abundance will be NAN. Sometimes, to avoid some problem in following calculation, this kind of NAN value need be set as zero. Defalt is TRUE.}
  \item{unit.sum}{When a beta diversity index is calculated for a bin, the taxa abundances will be divided by unit.sum to calculate the relative abundances. usually, unit.sum can be set as the sequencing depth in each sample. Default setting is NULL, means not to do this special transformation.}
}
\details{
beta mean pairwise distance.
}
\value{
Output is a distance object of pairwise betaMPD between samples.
}
\references{
Webb CO, Ackerly DD, and Kembel SW. 2008. Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Bioinformatics 18:2098-2100
}
\author{
Daliang Ning
}
\note{
Version 2: 2023.9.4.
Version 1: 2017.3.13.
}
\seealso{
\code{\link{bNRI.big.cm}}
}
\examples{
data("example.data")
comm=example.data$comm
tree=example.data$tree

# since it needs to save some file to a certain folder,
# the following code is set as 'not test'.
# but you may test the code on your computer
# after change the folder path for 'save.wd'.
\donttest{
wd0=getwd()
save.wd=paste0(tempdir(),"/pdbig.bmpd.big")
# you may change save.wd to the folder you want to save the pd.big output.
nworker=2 # parallel computing thread number
pd.big=pdist.big(tree = tree, wd=save.wd, nworker = nworker)
bmpd.wt=bmpd.big(comm=comm, pd.desc = pd.big$pd.file,
                 pd.spname = pd.big$tip.label, pd.wd = pd.big$pd.wd,
                 abundance.weighted = TRUE)
setwd(wd0)
}
}
\keyword{phylogenetic}
