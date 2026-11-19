\name{mpd.big}
\alias{mpd.big}
\title{
Mean pairwise distance (MPD) from big data.
}
\description{
Calculate mean pairwise distance (MNTD) in each community in a given community matrix, using 'bigmemory' to deal with too large dataset.
}
\usage{
mpd.big(comm, pd.desc, pd.spname = NULL, pd.wd,
        sp.limit = 10000, spname.check = FALSE,
        abundance.weighted = TRUE)
}
\arguments{
  \item{comm}{matrix or data.frame, community data matrix, rownames are sample names, colnames are OTU ids.}
  \item{pd.desc}{character, the name to describe bigmemory file of phylogenetic distance matrix, default is "pd.desc".}
  \item{pd.spname}{vector, the OTU ids (species names) in exactly the same order as the phylogenetic matrix rows or columns}
  \item{pd.wd}{the path of the folder saving the phylogenetic distance matrix.}
  \item{sp.limit}{numeric, the limitation of species number in each round of the looping inside this funciton, default is 10000, generally fine for most computer, if returns error, set lower value.}
  \item{spname.check}{logic, whether to check the OTU ids (species names) in community matrix and phylogenetic distance matrix are the same.}
  \item{abundance.weighted}{logic, whether weighted by species abundance, default is TRUE, means weighted.}
}
\details{
mean pairwise distance (MPD) in each community, which is the same index as 'mpd' in package 'picante', but calculated by matrix multiplication and using bigmemory to handle big dataset.
}
\value{
result is a numeric vector with sample names
}
\references{
Webb C, Ackerly D, McPeek M, and Donoghue M. (2002). Phylogenies and community ecology. Annual Review of Ecology and Systematics 33:475-505.

Kembel, S.W., Cowan, P.D., Helmus, M.R., Cornwell, W.K., Morlon, H., Ackerly, D.D. et al. (2010). Picante: R tools for integrating phylogenies and ecology. Bioinformatics, 26, 1463-1464.
}
\author{
Daliang Ning (ningdaliang@ou.edu)
}
\note{
Version 2: 2023.9.4.
Version 1: 2017.3.13.
}
\seealso{
\code{\link{mpdn}}, \code{\link{NRI.p}}, \code{\link{NRI.cm}}
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
save.wd=paste0(tempdir(),"/pdbig.mpd.big")
# you may change save.wd to the folder you want to save the pd.big output.
nworker=2 # parallel computing thread number
pd.big=pdist.big(tree = tree, wd=save.wd, nworker = nworker)
mpd.wt=mpd.big(comm=comm, pd.desc = pd.big$pd.file,
               pd.spname = pd.big$tip.label, pd.wd = pd.big$pd.wd,
               sp.limit = 10000, spname.check = FALSE,
               abundance.weighted = TRUE)
setwd(wd0)
}
}
\keyword{phylogenetic}
