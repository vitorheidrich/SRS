\name{Scaling with ranked subsampling curve (SRScurve)}
\alias{SRScurve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scaling with ranked subsampling curve (SRScurve)
}
\description{
For each column of the input data, draws a line plot of alpha diversity indices (see \link[SRS:SRScurve]{metric}) at different sample sizes (specified by \link[SRS:SRScurve]{step}) normalized by scaling with ranked subsampling (using \link[SRS]{SRS}). Minimum sample size (cutoff-level) can be evaluated by specifying \link[SRS:SRScurve]{sample}. The function further allows to visualize trade-offs between cutoff-level and alpha diversity and enables direct comparison of SRS and repeated rarefying.

See Beule & Karlovsky (2020) <doi:10.7717/peerj.9593> for details regarding SRS.
}
\usage{
SRScurve(data, metric = "richness", step = 50, sample = 0, max.sample.size = 0,
        rarefy.comparison = FALSE, rarefy.repeats = 10, 
        rarefy.comparison.legend = FALSE, xlab = "sample size", 
        ylab = "richness", label = FALSE, col, lty, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data frame (species count or OTU table) in which columns are samples and rows are the counts of species or OTUs. Only integers are accepted as data.
}
  \item{metric}{
Character, "richness" (using \link[vegan]{specnumber}) for species richness or "\link[vegan:diversity]{shannon}", "\link[vegan:diversity]{simpson}" or "\link[vegan:diversity]{invsimpson}" (using \link[vegan]{diversity}) for common diversity indices. Default is "richness".
}
  \item{step}{
Numeric, specifying the step used to vary the sample size. Default is 50.
}
  \item{sample}{
Numeric, specifying the cutoff-level to visualize trade-offs between cutoff-level and alpha diversity.
}
  \item{max.sample.size}{
Numeric, specifying the maximum sample size to which SRS curves are drawn. Default is 0 which does not limit the maximum sample size.
}
  \item{rarefy.comparison}{
Logical, if TRUE, median values of rarefy with n repeats (specified by \link[SRS:SRScurve]{rarefy.repeats}) will be drawn for comparison. Default is FALSE.
}
  \item{rarefy.repeats}{
Numeric, specifying the number of repeats used to obtain median values for rarefying. Default is 10.
}
  \item{rarefy.comparison.legend}{
Logical, if TRUE, a legend for the comparison between SRS and rarefy is plotted. Default is FALSE.
}
  \item{xlab, ylab, label, col, lty, ...}{
Graphical parameters.
}
}
\details{
See Beule & Karlovsky (2020) <doi:10.7717/peerj.9593> for details regarding scaling with ranked subsampling.
}
\value{
Returns a line plot visualizing the change in alpha diversity indices with changing sample size.
}
\references{
Beule L, Karlovsky P. 2020. Improved normalization of species count data in ecology by scaling with ranked subsampling (SRS): application to microbial communities. PeerJ 8:e9593

<https://doi.org/10.7717/peerj.9593>
}
\author{
Vitor Heidrich, Petr Karlovsky, Lukas Beule
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

##Default settings of SRScurve.
SRScurve(example_input_data, metric = "richness", step = 50,
         ylab = "richness", 
         col = c("#000000", "#E69F00", "#56B4E9"))

##Limit the compution of SRS curves to a sample size of 200.
SRScurve(example_input_data, metric = "richness", step = 50,
         max.sample.size = 200, ylab = "richness", 
         col = c("#000000", "#E69F00", "#56B4E9"))

##SRScurve with comparison of SRS (solid lines) and repeated rarefying (dashed lines).
##Different colors correspond to indiviual samples. Cuttoff-level set to 200.
SRScurve(example_input_data, metric = "richness", step = 50,
         sample = 200, max.sample.size = 200,
         rarefy.comparison = TRUE, rarefy.repeats = 10, rarefy.comparison.legend = TRUE,
         ylab = "richness", 
         col = c(rep(c("#000000", "#E69F00", "#56B4E9"),2)),
         lty = c(1,2))
}
