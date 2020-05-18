# MarSuPial
Estimating recombinant fraction, genetic distance, and recombination rate from Marker Selected Pools

Marker selected pools are pooled-sequencing of marker selected offsprings from a two generation recombinant backcross scheme. The MarSuPial package takes filtered and polarized allele-specific count data and infers the recombinant fraction, genetic distance and recombination rates from the allele-frequency decay across the chromosome. MarSuPial does not genotype, filter, or polarize counts. We recommend following GATK best practices for genotyping.
MarSuPial also simulates read count data for allele-frequency decay provided the chromosome-wide recombination rate.  

Installation

Download the package MarSuPial_1.0.tar.gz and install using command line:
    <p>$ R CMD INSTALL MarSuPial_1.0.tar.gz</p>

Then, in R, load the package:
    <p>> library(MarSuPial)</p>
