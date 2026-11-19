\name{metaset}
\alias{metaset}
\title{
Metacommunity settings
}
\description{
This function is to make the IDs, frequence, and abundnce in metacommunities consistent among different elements.
}
\usage{
metaset(comm = NULL, meta.group = NULL, meta.com = NULL, meta.frequency = NULL, meta.ab = NULL)
}
\arguments{
  \item{comm}{matrix or data.frame, community data, each row is a sample or site, each colname is a taxon (a species or OTU or ASV), thus rownames should be sample IDs, colnames should be taxa IDs.}
  \item{meta.group}{matrix or data.frame, a one-column (n x 1) matrix indicating which metacommunity each sample belongs to. Rownames are sample IDs. The first column is metacommunity names. Such that different samples can belong to different metacommunities. If input a n x m matrix, only the first column is used. NULL means all samples belong to the same metacommunity. Default is NULL, means all samples from the same metacommunity.}
  \item{meta.com}{a list object, each element is a matrix or data.frame to define abundance (or relative abundance) of taxa in a metacommunity (regional pool). The element names indicate metacommunity names, which should be consistent with the metacommunity names defined in meta.group. If there is only one metacommunity, meta.com can be a matrix or data.frame to define taxa abundance (or relative abundance) in the metacommunity. Default is NULL, means to calculate metacommunity structure from comm according to metacommunities defined in meta.group.}
  \item{meta.frequency}{a list object, each element is a matrix or data.frame to define abundance (or relative abundance) of taxa in a metacommunity (regional pool). The element names indicate metacommunity names, which should be consistent with the metacommunity names defined in meta.group. If there is only one metacommunity, meta.com can be a matrix or data.frame to define taxa abundance (or relative abundance) in the metacommunity. Default is NULL, means to calculate metacommunity structure from comm according to metacommunities defined in meta.group, or calculate from meta.com.}
  \item{meta.ab}{matrix or data.frame, each column represents a taxon, each row represents a metacommunity (regional pool), to define the aubndance (or relative abundance) of each taxon in each metacommunity. The rownames indicate metacommunity names, which should be the same as the metacommunity names in meta.group. Default setting is NULL, means to calculate meta.ab as average relative abundance of each taxon in comm across the samples within each metacommunity defined by meta.group, or calculate from meta.com.}
}
\details{
This is a handy tool for taxa.binphy.big.cm and other .cm (multiple metacommunities) functions.
}
\value{
Output is a list object. comm=comm,meta.group=meta.group,meta.com=meta.com,meta.frequency=meta.frequency,meta.ab=meta.ab,meta.spool=meta.spool,comms=comms
  \item{comm}{community data table after matching IDs.}
  \item{meta.group}{talbe showing metacommunity of each sample, after matchin IDs.}
  \item{meta.com}{a list object showing each metacommunity data.}
  \item{meta.frequence}{table showing observed frequency of each taxon in each metacommunity.}
  \item{meta.ab}{table showing abundances of each taxon in each metacommunity.}
  \item{meta.spool}{a list object showing detected taxon IDs in each metacommunity.}
  \item{comms}{a list object splitting the observed community data table to different metacommunities.}
}
\references{
# not yet
}
\author{
Daliang Ning
}
\note{
Version 1: 2024.12.22.
}
\examples{
library(iCAMP)
data("example.data")
comm=example.data$comm

# in this example, 10 samples from one metacommunity,
# the other 10 samples from another metacommunity.
meta.group=data.frame(meta.com=c(rep("meta1",10),rep("meta2",10)))
rownames(meta.group)=rownames(comm)

metas=metaset(comm=comm,meta.group=meta.group)
metas$meta.spool
}
\keyword{tool}
