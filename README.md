# MarSuPial
**Estimating recombinant fraction, genetic distance, and recombination rate from Marker Selected Pools**

Marker selected pools are pooled-sequencing of marker selected offsprings from a two generation recombinant backcross scheme. The MarSuPial package takes filtered and polarized allele-specific count data and infers the recombinant fraction, genetic distance and recombination rates from the allele-frequency decay across the chromosome. MarSuPial does not genotype, filter, or polarize counts. We recommend following GATK best practices for genotyping.
MarSuPial also simulates read count data for allele-frequency decay provided the chromosome-wide recombination rate.  

## Installation

Download the package MarSuPial_1.0.tar.gz and install using command line:
    <p>$ R CMD INSTALL MarSuPial_1.0.tar.gz</p>

Then, in R, load the package:
    <p>> library(MarSuPial)</p>

OR Install in R by first extracting the tar ball and:
    <p>> install.packages(MarSuPial)</p>
    <p>> library(MarSuPial)</p>

## Usage
The MarSuPial package contains a collection of functions to analyzing, estimating, and simulating recombination in a two generation backcross.

### Analytical relationships:
Forcing a line-break\s\s
Next line in the list
r2d.locus(rate, l, size, pos, start, end)\s\s
Returns a list of the chromosome-wide windows and the recombination rate and genetic distance in cM from a locus (l). Rate (in cM/Mb) must be either a r function describing the relationship between recombination rate or chromosome position or a bed file (loaded as a dataframe) of equally sized chromosome windows. If pos (in Mb) is a single value, it is used as the window size of evenly spaced non-overlapping windows. If pos (in Mb) is a vector of numbers, they are used as positions to determine rate and distance. If a recombination rate map function is provided, the start and end positions of the recombination rate can also be inputted; positions < start and > end will have rates of 0 to emulate heterochromatic/centromeric/telomeric suppression.

d2D(d, method)
Converts genetic distance (d) into recombinant fraction (D), using either the "haldane" or "kosambi" mapping functions (method).

D2d(D, method)
Converts recombinant fraction (D) into genetic distance (d), using either the "haldane" or "kosambi" mapping functions (method).

D2AF(D, fitness, ploidy, mendel_rate)
Converts recombinant fraction (D) to allele frequency (AF) after taking taking into account the fitness differential (fitness) of the selected locus, the ploidy, and the expected mendelian ratio (mendel_rate).

AF2D(AF, fitness, ploidy, mendel_rate)
Converts allele frequency (AF) to recombinant fraction (D) after taking taking into account the fitness differential (fitness) of the selected locus, the ploidy, and the expected mendelian ratio (mendel_rate).

r2d.3pt(ratefunc, l1, l2, size, winsize = 0.01, start, end)
Same as r2d.locus, but for two loci under selection. Rate must be a function.

d2AF.3pt(list3pt, fitness1, fitness2, method)
Takes an r2d.3pt outlut (list), and generates the alleles frequency given the fitness differential of the two loci.
