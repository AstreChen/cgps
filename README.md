# cgps

###############
# preparation #
###############
install R package:
* Biobase
* limma
* edgeR
* gage
* GSVA
* EnrichmentBrowser

install python modules:
* numpy
* pandas
* sklearn

run under the code folder

Notes: 
1. only support human species for now
2. only support NCBI Entrez Gene ID for now

################
# steps        #
################
Rscript combined_methods.R expfile phefile datatype=[ma/rseq] output_directory
