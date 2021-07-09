\name{Scaling with ranked subsampling (SRS) Shiny app}
\alias{SRS.shiny.app}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Shiny app for scaling with ranked subsampling (SRS) 
}
\description{
Shiny app for the determination of Cmin for scaling with ranked subsampling (SRS).
}
\usage{
SRS.shiny.app(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data frame (species count or OTU table) in which columns are samples and rows are the counts of species or OTUs. Only integers are accepted as data.
}
}
\details{
Shiny app that generates a visualization of retained samples, summary statistics, SRS curves, and an interactive table in response to varying minimum sample size (Cmin).
}
\value{
Launches Shiny app for SRS in the default web browser.
}
\references{
Beule L, Karlovsky P. 2020. Improved normalization of species count data in ecology by scaling with ranked subsampling (SRS): application to microbial communities. PeerJ 8:e9593

<https://doi.org/10.7717/peerj.9593>
}
\author{
Vitor Heidrich, Devon O'rourke, Petr Karlovsky, Lukas Beule
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

##Launching the SRS shiny app with example_input_data as input
if (interactive()) {SRS.shiny.app(example_input_data)}
}
