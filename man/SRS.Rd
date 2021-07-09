\name{Scaling with ranked subsampling (SRS)}
\alias{SRS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scaling with ranked subsampling (SRS)
}
\description{
Scaling with ranked subsampling (SRS) for the normalization of ecological count data.
It is recommended to use \link[SRS]{SRS.shiny.app} for the determination of Cmin.
}
\usage{
SRS(data, Cmin, set_seed = TRUE, seed = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data frame (species count or OTU table) in which columns are samples and rows are the counts of species or OTUs. Only integers are accepted as data.
}
  \item{Cmin}{
The number of counts to which all samples will be normalized. Typically, the total OTU count of the sample with the lowest sequencing depth is chosen as Cmin. Samples with sequencing depth lower than the chosen Cmin will be discarded.
}
  \item{set_seed}{
Logical, if TRUE, a seed is set to enable reproducibility of SRS if OTUs with identical Cfrag as well as Cint are sampled randomly without replacement. See \link[base]{set.seed} for details. Default is TRUE.
}
  \item{seed}{
Integer, specifying the seed. See \link[base]{set.seed} for details. Default is 1.
}
}
\details{
It is recommended to use \link[SRS]{SRS.shiny.app} for the determination of Cmin.
SRS consists of two steps. In the first step, the total counts for all OTUs (operational taxonomic units) or species in each sample are divided by a scaling factor chosen in such a way that the sum of the scaled counts Cscaled equals Cmin. In the second step, the non-integer Cscaled values are converted into integers by an algorithm that we dub ranked subsampling. The Cscaled value for each OTU or species is split into the integer part Cint  (\eqn{Cint = floor(Cscaled)}) and the fractional part Cfrac (\eqn{Cfrac = Cscaled - Cint}). Since \eqn{\Sigma Cint \le Cmin} , additional \eqn{\Delta C = Cmin - \Sigma Cint} counts have to be added to the library to reach the total count of Cmin. This is achieved as follows. OTUs are ranked in the descending order of their Cfrac values. Beginning with the OTU of the highest rank, single count per OTU is added to the normalized library until the total number of added counts reaches \eqn{\Delta C} and the sum of all counts in the normalized library equals Cmin. When the lowest Cfrag involved in picking \eqn{\Delta C} counts is shared by several OTUs, the OTUs used for adding a single count to the library are selected in the order of their Cint values. This selection minimizes the effect of normalization on the relative frequencies of OTUs. OTUs with identical Cfrag as well as Cint are sampled randomly without replacement.
}
\value{
Data frame normalized to Cmin.
}
\references{
Beule L, Karlovsky P. 2020. Improved normalization of species count data in ecology by scaling with ranked subsampling (SRS): application to microbial communities. PeerJ 8:e9593

<https://doi.org/10.7717/peerj.9593>
}
\author{
Lukas Beule, Vitor Heidrich, Devon O'rourke, Petr Karlovsky
}
\examples{
##Samples should be arranged columnwise.
##Input data should not contain any categorial
##data such as taxonomic assignment or barcode sequences.
##An example of the input data can be found below:

example_input_data <- matrix(c(sample(1:20, 100, replace = TRUE),
sample(1:30, 100, replace = TRUE),sample(1:40, 100, replace = TRUE)), nrow = 100)
colnames(example_input_data) <- c("sample_1","sample_2","sample_3")
example_input_data <- as.data.frame(example_input_data)
example_input_data

##Selection of the desired number of counts
##(e.g., total OTU counts of the sample with the lowest sequencing depth):

Cmin <- min(colSums(example_input_data))
Cmin

##Running the SRS function
SRS_output <- SRS(example_input_data, Cmin)
SRS_output

##Samples that have a total number of counts < Cmin will be discarded:
SRS_output <- SRS(example_input_data, Cmin+1)
SRS_output
}
