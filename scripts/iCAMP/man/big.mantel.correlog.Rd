\name{big.mantel.correlog}
\alias{big.mantel.correlog}
\title{
Mantel Correlogram and stepwise Mantel test for big data
}
\description{
This function is to perform Mantel Correlogram and stepwise Mantel test for phylogenetic signal test or other similar analysis on big data.
}
\usage{
big.mantel.correlog(x.desc, x.name, x.wd, y.desc, y.name, y.wd,
                    nworker = 4, interval = 0.02, break.pts = NULL,
                    permutations = 999, strata = NULL, diag.in = FALSE,
                    temp.wd = getwd(), stepw = 10^7, padjust.method = "fdr",
                    cleanup = FALSE, parallelset = TRUE)
}
\arguments{
  \item{x.desc}{character, the name to describe bigmemory file of the distance matrix X, which defines values on the x-axis of correlogram. For phylogenetic signal test, this is the phylogenetic distance matrix file name, e.g., "pd.desc".}
  \item{x.name}{vector, the rownames of the distance matrix X. For phylogenetic signal test, these are species name or taxa IDs.}
  \item{x.wd}{the path of the folder saving the distance matrix X. For phylogenetic signal test, this is the folder saving phylogenetic distance matrix file.}
  \item{y.desc}{character, the name to describe bigmemory file of the distance matrix Y. For phylogenetic signal test, this is the niche preference difference matrix file name, e.g., "pH.desc".}
  \item{y.name}{vector, the rownames of the distance matrix Y. For phylogenetic signal test, these are species name or taxa IDs.}
  \item{y.wd}{the path of the folder saving the distance matrix Y. For phylogenetic signal test, this is the folder saving niche difference matrix file.}
  \item{nworker}{for parallel computing. Either a character vector of host names on which to run the worker copies of R, or a positive integer (in which case that number of copies is run on localhost). default is 4, means 4 threads will be run.}
  \item{interval}{the interval of distance X. If break.pts is NULL, this defines how to seperate the data into different classes.}
  \item{break.pts}{Vector containing the break points of the distance X distribution. Provide (n.class+1) breakpoints, that is, a list with a beginning and an ending point. Default: break.pts=NULL.}
  \item{permutations}{Number of permutations for the tests of significance. Default: nperm=999. For large data files, permutation tests are rather slow.}
  \item{strata}{An integer vector or factor specifying the strata for permutation. If supplied, observations are permuted only within the specified strata.}
  \item{diag.in}{logic, whether count in the diagonal elements when calculating r. default is FALSE.}
  \item{temp.wd}{path of a folder that temporarily stores the big matrix file for intermediate results.}
  \item{stepw}{step width. When processing a large dataset, need to divid the long vector to different steps to avoid exceed size or space limitation.}
  \item{padjust.method}{P value correction method, a character string. options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". see function p.adjust. default is 'fdr'.}
  \item{cleanup}{logic, if TRUE, it will try to remove the big matrix file for intermediate results.}
  \item{parallelset}{logic, if TRUE, it will use parallel computing in the first step to set up values.}
}
\details{
The Mantel Correlogram analysis is the same as the function mantel.correlog in package 'vegan', but applicable to large datasets.
The stepwise Mantel test is to perform Mantel test using all distance values within each class or using all distance values smaller than the end point of each class, to directly answer the question whether there is significant phylogenetic signal when phylogenetic distances are within a range or lower than a certain value. 
}
\value{
The output is a list,
  \item{mantel.res}{A table with the start, end, and middle points of each class (Break.point.start, Break.point.end, Class.midpoint), the number of distance values within each class (n.withinclass), the number of distance values smaller than the end point of each class (n.beforeend), various statistics of X and Y values within each class, correlogram correlation coefficient r value for each class (r.correlog), r value within each class (r.within.class), r value using distance values smaller than the end point of each class (r.before.endpt), corresponding p values (p.correlog, p.within.class, p.before.endpt) and adjusted p values (p.adj.correlog, p.adj.within.class, p.adj.before.endpt).}
  \item{n.class}{The number of distance classes.}
  \item{break.pts}{The break points provided by the user or computed by the program.}
  \item{padjust.method}{P value correction method.}
  \item{permutations}{permutation number.}
  \item{strata}{strata setting.}
  \item{diag.in}{whether include diagonal elements.}
  \item{temp.wd}{folder that stores intermediate results.}
}
\references{
# update later
}
\author{
Daliang Ning
}
\note{
version 1: 2023.3.4
}
\seealso{
\code{\link{dniche}}, \code{\link{ps.bin}}
}
\examples{
data("example.data")
comm=example.data$comm
env=example.data$env
tree=example.data$tree

# it is designed for large data, you need to use big.memory
# since big.memory need to specify a certain folder,
# it is set as 'not test'.
# but you may test the code on your computer after change the path for 'save.wd'.
\donttest{
  wd0=getwd()
  save.wd=paste0(tempdir(),"/dnichewd")
  # please change to the folder you want to save the big niche difference matrix.
  
  nworker=2 # parallel computing thread number
  nichedi=dniche(env = env, comm = comm,
                   method = "niche.value", nworker = nworker,
                   out.dist=FALSE,bigmemo=TRUE,nd.wd = save.wd)
  
  pd.big=pdist.big(tree = tree, wd=save.wd, nworker = nworker)
  
  mcrl1=big.mantel.correlog(x.desc=pd.big$pd.file,x.name=pd.big$tip.label,x.wd=pd.big$pd.wd,
                          y.desc=nichedi$nd[[1]],y.name=nichedi$names,y.wd=nichedi$nd.wd,
                          nworker=nworker,interval=0.02,break.pts=NULL,
                          permutations=999,strata=NULL,
                          diag.in=FALSE,temp.wd=save.wd,stepw=10^6,
                          padjust.method="fdr")
  
  setwd(wd0)
}
}
\keyword{phylogenetic signal}
