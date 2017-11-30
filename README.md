# Documentation v. 1.0

To run workflow, open terminal, cd to parent directory and run: bash runworkflow.sh

Be sure to have all dependencies installed and config.json configured to your needs.

## Dependencies

numpy
matplotlib
scipy
biopython
clustalW

## Config File

Config Dictionary Structure (config.json):

@params:
	desc 			        : (Str) 	- description of config file
	email 			        : (Str) 	- email used for Entrez
	db 			            : (Str) 	- Entrez db to query
	outputDirectory 	    : (Str) 	- output folder name
	geneDirectory 	    	: (Str) 	- gene folder name
	scoreDirectory 	    	: (Str) 	- score folder name 
	statsDirectory 		    : (Str) 	- stats folder name 
	groups 			        : (List) 	- query groups
	organisms 		        : (List) 	- query organisms
	targetOddsThreshold 	: (Float) 	- Odds pathway thresholds 
	diffExpressedDataPath 	: (Str) 	- Path to file containing data on differentially expressed gene probes
	universeDataPath	    : (Str) 	- Path to file containing data on all gene probes
	upathwaysDataPath	    : (Str) 	- Path to file containing data on KEGG pathways
    evidenceLevels          : (List)    - Strings representing levels of evidence to use in query, acceptable values include: known, reviewed, validated, providional, predicted, inferred

Default Settings of config.json:

{
    "desc": "config",
    "email": "anonymous@ohsu.edu",
    "db": "nuccore",
    "geneDirectory": "genes",
    "outputDirectory": "output",
    "scoreDirectory": "scores",
    "statsDirectory": "statistics",
    "groups": [
        "de",
        "nonde"
    ],
    "organisms": [
        "homo sapiens",
        "mus musculus",
        "gallus gallus"
    ],
    "targetOddsThreshold": 1.5,
    "diffExpressedDataPath": "data/H5N1_VN1203_DE_Probes.txt",
    "universeDataPath": "data/H5N1_VN1203_UNIVERSE_Probes.txt",
    "pathwaysDataPath": "data/KEGG_Pathway_Genes.txt"
}