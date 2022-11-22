# DeconvolutionAtheroscleroticPlaques

You can find here all the code used for my project "Deconvolution of atherosclerotic plaques and its clinical association".

1st PART OF THE PROJECT, BENCHMARKING:

To reproduce the results, you will first need acess to the Athero-Express data (save it in input folder). Then, you should follow this workflow:
1. Run DataProcessing.R: processing of the bulk and single-cell RNA-seq datasets and creation of the files needed for the deconvolution algorithms (create folder output/[ourSC|ourBulk] to save the results).
2. Follow the analysis with pseudo-bulk artificially created data with pseudoBulk.R.
3. Deconvolution with MuSiC, Bisque and NNLS can be done following MusicBisqueNNLS_deconvolutions.R, but for Scaden you will need to work with Python on the command line (see https://scaden.readthedocs.io/en/latest/), and for CibersortX you need to use their web server (https://cibersortx.stanford.edu) (save all the results in results/[deconvolutionAlgorithm] folder).
4. Integrate all the deconvolution results and analyze them with IntegrationResults.R. 

2nd PART OF THE PROJECT, CLINICAL ASSOCIATION:
1. Obtain final proportions (read report).
2. Run clinical analysis with ClinicalAssociation.R.

To obtain the same graphs I used in my report, you can refer to reportFigures.R.

* in scadenModification you can find the script I modified to make the Scaden simulation step reproducible (changes commented with #gbb) with Scaden's original license.
