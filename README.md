# SRS
[R package](https://cran.r-project.org/package=SRS) for microbiome count data normalization by scaling with ranked subsampling (SRS)

![CRAN Badge](http://www.r-pkg.org/badges/version/SRS)  ![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/SRS)

Read more about this normalization method in the [SRS paper](https://doi.org/10.7717/peerj.9593) (Beule and Karlovsky, PeerJ 2020).

## Installing

Install the latest version from CRAN by running:
```R
install.packages("SRS")
```
## Using

The SRS R package features three R functions:
* `SRS()` - performs SRS normalization at a user-defined number of reads per sample (C<sub>min</sub>)
* `SRScurve()` -  draws alpha diversity rarefaction curves for SRS-normalized data (instead of rarefied data)
* `SRS.shiny.app()` - generates a visualization of retained samples, summary statistics, SRS curves, and an interactive table in response to varying C<sub>min</sub>

Refer to the SRS [reference manual](https://cran.r-project.org/web/packages/SRS/SRS.pdf) for usage details.

##### Citation
If you use this package in your research paper, please cite as:

Heidrich V, Karlovsky P, Beule L. 2021. ‘SRS’ R package and ‘q2-srs’ QIIME 2 plugin: Normalization of Microbiome Data Using Scaling with Ranked Subsampling (SRS). [*Appl. Sci.* 11(23), 11473](https://doi.org/10.3390/app112311473).

When referencing the SRS algorithm itself, please cite:

Beule L, Karlovsky P. 2020. Improved normalization of species count data in ecology by scaling with ranked subsampling (SRS): application to microbial communities. [*PeerJ* 8:e9593](https://doi.org/10.7717/peerj.9593).
