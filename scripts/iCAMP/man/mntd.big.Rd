\name{mntd.big}
\alias{mntd.big}
\title{
Mean nearest taxon distance (MNTD) for big data
}
\description{
Calculate mean nearest taxon distance (MNTD) in each community in a given community matrix, using 'bigmemory' to deal with too large dataset.
}
\usage{
mntd.big(comm, pd.desc = "pd.desc", pd.spname = NULL,
         pd.wd = getwd(), spname.check = FALSE,
         abundance.weighted = TRUE, silent = TRUE,
         memory.G = 50)
}
\arguments{
  \item{comm}{matrix or data.frame, community data matrix, rownames are sample names, colnames are OTU ids.}
  \item{pd.desc}{character, the name to describe bigmemory file of phylogenetic distance matrix, default is "pd.desc".}
  \item{pd.spname}{vector, the OTU ids (species names) in exactly the same order as the phylogenetic matrix rows or columns}
  \item{pd.wd}{the path of the folder saving the phylogenetic distance matrix.}
  \item{spname.check}{logic, whether to check the OTU ids (species names) in community matrix and phylogenetic distance matrix are the same.}
  \item{abundance.weighted}{logic, whether weighted by species abundance, default is TRUE, means weighted.}
  \item{silent}{Logic, if FALSE, some messages will be showed during calculation. Default is FALSE.}
  \item{memory.G}{numeric, to set the memory size as you need, so that calculation of large tree will not be limited by physical memory. unit is Gb. default is 50Gb.}
}
\details{
mean nearest taxon distance (MNTD) in each community, using the same algrithm as the function 'mntd' in package 'picante'.
}
\value{
result is a numeric vector with sample names
}
\references{
Webb CO, Ackerly DD, and Kembel SW. (2008). Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Bioinformatics 18:2098-2100

Kembel, S.W., Cowan, P.D., Helmus, M.R., Cornwell, W.K., Morlon, H., Ackerly, D.D. et al. (2010). Picante: R tools for integrating phylogenies and ecology. Bioinformatics, 26, 1463-1464.
}
\author{
Daliang Ning
}
\note{
Version 2: 2023.9.4.
Version 1: 2017.3.13.
}
\seealso{
\code{\link{NTI.cm.big}}
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
save.wd=paste0(tempdir(),"/pdbig.mntd.big")
# you may change save.wd to the folder you want to save the pd.big output.
nworker=2 # parallel computing thread number
pd.big=pdist.big(tree = tree, wd=save.wd, nworker = nworker)
mntd.wt=mntd.big(comm=comm, pd.desc = pd.big$pd.file,
               pd.spname = pd.big$tip.label, pd.wd = pd.big$pd.wd,
               spname.check = FALSE, abundance.weighted = TRUE)
setwd(wd0)
}
}
\keyword{ phylogenetic }
