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

Data format
Please refer to 'test' directory
the test data is microarray 

################
#     steps    #
################
Rscript combined_methods.R expfile phefile datatype=[ma/rseq] datadir
python predict.py datadir outdir

* expfile: file to save expression data
* phefile: file to save the experiment category of expfile
* data type : for expression data, "rseq" for RNA-Seq data, 'ma' for microarray data
* datadir: directory to save the results of individual methods
* outdir : directory to save the combined results of CGPS

