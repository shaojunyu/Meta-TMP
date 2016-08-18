
FSWeight and Other Predictions
Version 2.1
==============================



USAGE: 

predict.pl -i interactions \
         [-s annotscheme] [-a annotations] [-o outputfile] \
         [-p action] [-f folds] [-x informative] [-k] \
         [-m prediction method] [-w weighting scheme] \
         [-t protein list]


options:

General Parameters:

        -o      Output File (Default = output.txt)
			This specifies an output file for which results are printed.	

        -p      Action (Default = p)
                	p = prediction of novel annotations
				This will print predictions for all unannotated proteins to the specified output file.
				The format of the output is as follows:
				Protein Name\tGO Term\tScore
				E.g.:
				SGD_S000001924  0016481 0.121650140097154

			c = cross validation
				This will perform a N-Fold cross validation for predictions using Informative GO Terms.
				N is specified using the -f specifier. When N is set to 1 (Default), LOOCV is performed.

				Informative threshold is set using the -x specifier. (Default = 30)
				Predictions for each fold will be printed to the specified output file.
				Cross Validation results (Precision-Recall analysis and ROC analysis) are printed to the standard output.

				The format of the output is as follows:
				a) Precision-Recall Analysis:
					1) First Line:
					No. of Annotations\tPredicted Annotations\tNo. of proteins predicted
					E.g.:
					5299    60079   4889
					2) Subsequent Lines:
					Score Threshold\tPredictions\tRecall\tPrecision
					E.g.:
					0.743805825614686       112     0.0200037742970372      0.946428571428571
				
				b) ROC (Receiver Operating Characteristics) Analysis
					GO Namespace|GO Level|Description|Parent GO Terms\tGO Term\tROC
					E.g.:
					biological_process|6|cell cycle checkpoint|0000074      0000075  0.956277508936605
			
				For more information on Cross-Validation, GO Level, Informative GO Terms and validation measures, please refer to [1] & [2]

			w = weight interactions
				This will print the topological weight of network edges to the specified output file.
				Weighting Scheme is specified using the -w specifier. (Default = FSWEIGHT)
				The format of the output is as follows:
				Protein A\tProtein B\tType\tWeight
				E.g.:
				SGD_S000005552  SGD_S000005939  L1      0.107143836623493

Functional Annotations:

        -s      Annotation Scheme
			This is a flat text file describing GO Terms.
			The format of the file is as follows:
			GO ID\tGO Namespace|GO Level|Description|Parent GO Terms
			E.g.:
			0051177 biological_process|8|meiotic sister chromatid cohesion|0007062,0007127,0045132
			Processed GO Scheme files can be found at http://srs2.bic.nus.edu.sg/~kenny/integration/
			This input is necessary when p or c are selected for the -p specifier.
			Non-GO annotation schemes are also supported. See TIPS below.

        -a      Annotations
			This is a flat file describing GO Annotations.
			The format of the file is as follows:
			Protein Name|GO Term|
			E.g.:
			UniProt_O77783|0050508|
			A GO Annotation file can be found at http://srs2.bic.nus.edu.sg/~kenny/integration/
			This input is necessary when p or c are selected for the -p specifier.
			Non-GO annotation schemes are also supported. See TIPS below.

	-v	Vague List
			This is a comma-delimited list of vague annotation IDs which will be ignored during prediction or cross validation.
			Default list is 0000004,0005554,0008372,0003674,0005575,0008150 for GO Annotations
			If no vague terms are defined, set -v 0.

	-r	Minimum Meaningful Annotation Level
			This is the minimum annotation level (for a hierarchical annotation scheme) in which the sharing of function is considered meaningful.
			The default level is 4 for GO annotations.
			For non-GO annotation scheme, -r 1 usually suffices.

	TIP:	If an annotation scheme other than Gene ontology is used:

		a)	Non-Hierarchical annotation scheme
			i)	Supply annotation scheme file in the following format:
			ID\tName
			ii)	Set options -v 0 and -r 1

		b) 	Hierarchical annotation scheme
			i)	Supply annotation scheme file in the following format:
				GO ID\t|Level|Description|Parent Terms
			ii)	Propagate all annotations so that true path annotation is ensured.
			iii)	Set options -v based on vague terms in the scheme e.g. 799, 96, 98, 99 for MPS FunCat 2.0
			iv)	Set options -r 1

Protein Binary Pairs:

        -i      Interactions
		This is a flat file describing binary protein relationship (e.g. interactions)
		The format of the file is as follows:
		Protein A\tProtein B\tDatasource Name\tWeight
		E.g.:
		SGD_S000001474	SGD_S000005922	L1	0.0520662450872069
		Weight is optional and is only used for the PROB_COMBINE prediction method.

Cross Validation:

        -f      Folds (Default = 1)
		This specifies the number of folds for cross validation. (Only applies when -p c is set)
		Defaults to 1, which will perform a Leave-One-Out cross validation (LOOCV).

        -x      Threshold for Informative GO Terms (Default = 30)
		The minimum number of proteins annotated to a GO Term to define Informative GO Terms. See [1] & [2] for more details.
		
        -k      Ignore isolated proteins (no neighbours) (Default = Off)
		If set, isolated proteins (i.e. not involved in any interactions or relationship specified by -i) will be ignored when 
		computing informative GO Terms and during cross validation.

Function Prediction:

        -m      Prediction Method (MAJORITY_VOTE, CHI_SQUARE, WEIGHTED_AVG, PROB_COMBINE) (Default = PROB_COMBINE)
		Specify prediction method to use.
		MAJORITY_VOTE - Neighbour Counting Method [3]
		CHI_SQUARE - Chi-Square Method [4]
		WEIGHTED_AVG - Weighted Averaging [1] & [2]
		PROB_COMBINE - Data Fusion Methods [5]

        -w      Weighting Scheme (CD_DIST, GEOMETRIC, FSWEIGHT) (Default = FSWEIGHT)
		Specify topological weight for weighting edges.
		CD_DIST	- 1 - Czekanowski-Dice distance [6]
		GEOMETRIC - Geometric Distance [7]
		FSWEIGHT - FS-Weight [1] & [2]

        -o      Protein List File (Perform prediction/validation only on selected proteins)
		List of proteins to predict or validate. If this is not specified, all proteins will be used.



Credits:

This program is implemented by CHUA Hon Nian, and is supported in part by
an NGS scholarship, and URC grant "R-252-000-274-112: Graph-Based 
Protein Function Prediction". If you use this program, please cite Refs 1, 2, 5 below.




References:

[1]	H. N. Chua, W.K. Sung, L. Wong. (2006)
	Exploiting indirect neighbours and topological weight to predict protein function from protein-protein interactions. 
	Bioinformatics, 22:1623-1630.
	
[2]	H. N. Chua, W.K. Sung, L. Wong. (2007)
	Using Indirect Protein Interactions for the Prediction of Gene Ontology Functions. 
	BMC Bioinformatics, 8(Suppl 4):S8.

[3]	B. Schwikowski, P. Uetz, S. Fields. (2000)
	A network of interacting proteins in yeast. 
	Nat Biotechnol, 18:1257-1261.

[4]	H. Hishigaki, K. Nakai, T. Ono, A. Tanigami, T. Takagi. (2001)
	Assessment of prediction accuracy of protein function from protein-protein interaction data.
	Yeast, 18:525-531.

[5]	H. N. Chua, W.K. Sung, L. Wong. (2007)
	An efficient strategy for extensive integration of diverse biological data for protein function prediction.
	Submitted to Bioinformatics.

[6]	C. Brun, F. Chevenet, D. Martin, J. Wojcik, A. Guénoche, B. Jacq. (2003)
	Functional classification of proteins for the prediction of cellular function from a protein-protein interaction network.
	Genome Biol, 5:R6.

[7]	D.S. Goldberg and F.P. Roth. (2002)
	Assessing experimentally derived interactions in a small world.
	Proc Natl Acad Sci U S A 2003 Apr 15 100(8):4372-6.
