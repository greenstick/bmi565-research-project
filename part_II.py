#! /usr/bin/env python3

"""
Developed With Python Version 2.7.8 by Benjamin Cordier

Refer to config.json dictionary for setup variables

Part II of workflow:
	- Parses targetpathways.json to dictionary
	- Formats Entrez queries to database for organisms and genes
	- Retrieves gene sequence IDs using formatted queries at specified evidence levels
	- Retrieves gene sequence files using retrieved sequence IDs
	- Performs multiple sequence alignment on genes using clustalW2
	- Computes Hamming distances and a normalized score for each gene
	- Performs statistical analysis using Mann-Whitney U test
	- Saves boxplots of distribution of scores by group (i.e. differentially expressed or not)

Files created per gene:
	- sequences file (fasta)
	- alignment file (aln)

Analysis Files Created:
	- hamming statistics files (json)
	- analysis file (json)
	- boxplots file (png)
"""

import sys 													as S
import csv 													as C
import re 													as Rgx
import numpy 												as Npy
import matplotlib.pyplot 									as Plot
import json 												as JSON
import scipy.stats 											as Stats
from multiprocessing        	import Pool, cpu_count
from functools              	import partial
from Bio 						import Entrez 				as E
from Bio.Align.Applications 	import ClustalwCommandline	as ClustalW
from Bio 						import AlignIO

import time
# 
# Utility Methods
# 

# Prompt user input from command line
def getUserInput(valid, prompt):
	"""
	Prompts user for input and validates
	@params:
		prompt 		- Required 	: verbose user prompt (Str)
		valid 		- Required 	: regex to validate against (Rgx)
	Returns: dicts (List)
	"""
	response = input(prompt)
	if Rgx.match(valid, response):
		return response
	else:
		print("Error: invalid input")
		getUserInput(valid, prompt)


# Parse JSON to list of dictionaries
def parseJSONToDicts(path):
	"""
	Parse JSON to list of dictionaries
	@params:
		path 		- Required 	: input path (Str)
	Returns: dicts (List)
	"""
	dicts = []
	with open(path, "r") as file:
		dicts = JSON.load(file)
	return dicts

# Write JSON output to file.
def writeJSON(data, path="output/output.json"):
	"""
	Write JSON output to file.
	@params:
		data 		- Required	: data to write (Str, Dict, List, Set)
		path 		- Optional 	: output path (Str)
	Returns: path (Str)
	"""
	with open(path, "w") as file:
		file.write(JSON.dumps(data, indent=4, separators=(',', ': '), sort_keys=False))
		file.close()
	return path

# General File Writing Function
def writeOutput(data, path="output/output.txt", stdout=False):
	"""
	General File Writing Function
	@params:
		data 		- Required	: data to write (Str, Dict, List, Set)
		path 		- Optional 	: output path (Str)
		stdout 		- Optional 	: write output to system (Boolean)
	Returns: path (Str)
	"""
	with open(path, "w") as file:
		if stdout == True:
			return S.stdout.write(file)
		else:
			for line in data:
				file.write(str(line) + "\n")
			file.close()
	return path

# Apply function to list in parallel
def applyParallelList (data = None, function = None, arguments = None, processes = None):
    """
    Parallelize function across list
    @params:
        data        : Required  - Pandas DataFrame  - List to apply function to
        function    : Required  - Reference         - Pre-defined function to call on data
        arguments   : Optional  - Tuple             - Tuple of ordered arguments for function
        processes   : Optional  - Integer           - Number of processors to use (defaults to cpu_count())
    @return:
        result      : Recombined data with the function performed on it
    """
    processes   = cpu_count() if processes is None else processes
    length      = len(data)
    grouped     = [data[i * (length // processes): (i + 1) * (length // processes)] for i in range(0, processes)]
    with Pool(processes) as pool:
        if arguments is None:
            resultList  = pool.map(partial(function), [group for group in grouped])
        else:
            resultList  = pool.map(partial(function, **arguments), [group for group in grouped])
    result      = sum([list(result) for result in resultList], [])
    return result

# 
# Database Query Methods
# 

# Format queries for Entrez
def formatQuery(pathways, geneKey="genes", selection="*", organisms=["homo sapiens"], type="DNA", seqFilter="Refseq"):
	"""
	Format queries for Entrez
	@params:
		parsed 				- Required 	: target pathways (List of Dicts)
		geneKey 			- Optional 	: key to genes in target pathway dict (Str)
		pathway 			- Optional 	: index of selected target pathway (Int)
		organisms			- Optional 	: organisms to query
		type 				- Optional 	: sequence type to query
		seqFilter 			- Optional	: additional filter for query
	Returns: queries (List of Dicts)
	"""
	queries = []
	print("\nStatus: formatting database queries")
	for organism in organisms:
		if selection == "*":
			for pathway in pathways:
				genes = {gene for gene in pathway[geneKey]}
				for gene in genes:
					queries.append({
						"gene"		: gene, 
						"desc"		: pathway["desc"], 
						"geneKey"	: geneKey, 
						"organisms"	: organisms, 
						"organism"	: organism, 
						"type"		: type, 
						"seqFilter"	: seqFilter, 
						"query"		: gene + "[Gene] AND " + organism + "[Organism] AND " + type + "[Filter] AND " + seqFilter + "[Filter]"
					})
		else:
			selection = int(selection)
			genes = {gene for gene in pathways[selection][geneKey]}
			for gene in genes:
				queries.append({
					"gene"		: gene, 
					"desc"		: pathways[selection]["desc"], 
					"geneKey" 	: geneKey, 
					"organisms" : organisms, 
					"organism" 	: organism, 
					"type" 		: type, 
					"seqFilter" : seqFilter, 
					"query" 	: gene + "[Gene] AND " + organism + "[Organism] AND " + type + "[Filter] AND " + seqFilter + "[Filter]"
				})
	# Sort Queries By Gene
	queries.sort(key = lambda x: x["gene"])
	print("Status: done\n")
	return queries

# Query database for sequence records
def getSequenceRecords(queries, evidenceLevels=["known"], db="nuccore", retmax=1):
	"""
	Query database for sequence records
	@params:
		queries 			- Required 	: database queries (List of Dicts)
		db 					- Optional 	: database to query
		retmax 				- Optional 	: max number of records to retrieve
	Returns: records (List of Dicts)
	"""
	tempRecords, validate, records = [], [], []
	# parallelArgs = {
	# 	"evidenceLevels" 	: evidenceLevels,
	# 	"db" 				: db,
	# 	"retmax" 			: retmax,
	# 	"i" 				: 0
	# }
	# tempRecords = applyParallelList(queries, function=__getSequenceRecord__, arguments=parallelArgs, processes=ncores)
	for query in queries:
		for evidenceLevel in evidenceLevels:
			fullQuery = query["query"] + " AND srcdb_refseq_" + evidenceLevel + "[Prop]"
			try:
				print("Status: querying",db,"database -",query["gene"],query["organism"],evidenceLevel)
				search = E.esearch(db=db, term=fullQuery, retmax=retmax)
				record = E.read(search)
				search.close()
				if len(record["IdList"]) > 0:
					tempRecords.append({
						"gene" 		: query["gene"], 
						"organism"	: query["organism"], 
						"record" 	: record
					})
					print("Status: query returned",record["Count"],evidenceLevel,"record(s)")
				else:
					print("Status: no",evidenceLevel,"records")
			except:
				print("Status: request failed on query: " + fullQuery)
		print("Status: done\n")
	print("Status: requests complete, cleaning up records . . .")
	# Removing gene records that are incomplete across organisms
	for organism in query["organisms"]:
		unique = {record["gene"] for record in tempRecords if record["organism"] == organism}
		validate.append(unique)
	shared = set.intersection(*validate)
	for record in tempRecords:
		for gene in shared:
			if record["gene"] == gene:
				records.append(record)
	if len(records) == 0:
		print("Warning: No gene records shared between selected organisms")
	print("Status: done\n")
	return records

def __getSequenceRecord__(queries, evidenceLevels = ["known"], db = "nuccore", retmax = 1, i = 0):
	records = []
	for query in queries:
		for evidenceLevel in evidenceLevels:
			fullQuery = query["query"] + " AND srcdb_refseq_" + evidenceLevel + "[Prop]"
			try:
				print("Status: querying",db,"database -",query["gene"],query["organism"],evidenceLevel)
				search = E.esearch(db=db, term=fullQuery, retmax=retmax)
				record = E.read(search)
				search.close()
				if len(record["IdList"]) > 0:
					records.append({
						"gene" 		: query["gene"], 
						"organism"	: query["organism"], 
						"record" 	: record
					})
					print("Status: query returned",record["Count"],evidenceLevel,"record(s)")
				else:
					print("Status: no",evidenceLevel,"records")
			except:
				print("Status: request failed on query: " + fullQuery)
		print("Status: done\n")
	return records

# Query database for sequences using record IDs
def getSequences(records, db="nuccore", path="outputs/genes/"):
	"""
	Query database for sequences using record IDs
	@params:
		records 			- Required 	: sequence records (List of Dicts)
		db 					- Optional 	: database to query (Str)
		path 				- Optional 	: output path for gene files (Str)
	Returns: sequences (List of Dicts), errors (List of Strings)
	"""
	total, i 			= 0, 1
	genes 				= set()
	errors, sequences 	= [], []
	print("Status: retrieving sequences from", db, "database")
	for record in records:
		genes.add(record["gene"])
		total += len(record["record"]["IdList"])
	for gene in genes:
		geneSequences = set()
		for record in records:
			if record["gene"] == gene:
				# The commented code will 
				# geneSequences = geneSequences + applyParallelList(record["record"]["IdList"], function = __getSequenceID__, arguments = (i), processes = ncores)
				for Id in record["record"]["IdList"]:
					try:
						print("Status: queue position: ",i,"of",total,"queries\nStatus: requesting sequence ID:",Id,". . .")
						response 	= E.efetch(db=db, id=Id, rettype="fasta", retmode="text")
						fasta 		= response.read()
						response.close()
						geneSequences.add(fasta)
						print("Status: request success")
					except:
						print("Status: request failed on ID: \"",Id,"\"")
						errors.append(Id)
					print("Status: request done\n")
					i += 1
		writeOutput(geneSequences, path=path + gene + ".fasta")
		sequences.append({
			"gene": gene, 
			"path": path + gene + ".fasta"
		})
	if len(errors) == 0:
		print("Status: request series completed successfully\n")
	else:
		print("Status: request series completed with", len(errors), "errors\n")
	return sequences, errors

def __getSequenceByID__ (ids, i = 0):
	"""
	Internal Function: Query database for sequences using record IDs
	To be applied in parallel to a list
	@params:
		id 					- Required  : list of ids (List)
		i 					- Required  : index
	Returns: geneSequences (List fasta files)
	"""
	for ID in ids:
		try:
			print("Status: queue position:", i, "of", total, "queries\nStatus: requesting sequence ID:", ID, ". . .")
			response 	= E.efetch(db=db, id=ID, rettype="fasta", retmode="text")
			fasta 		= response.read()
			response.close()
			geneSequences.add(fasta)
			print("Status: request success")
		except:
			print("Status: request failed on ID: \"",ID,"\"")
		print("Status: request done\n")
		# Increment i reference (side effect)
		i += 1
	return geneSequences

# 
# Alignment Methods
# 

# Perform multiple sequence alignment using ClustalW
def msaSequences(sequences):
	"""
	Perform multiple sequence alignment using ClustalW
	@params:
		sequences 			- Required 	: sequence files (List of Dicts)
	Returns: outputs (List of Dicts), errors (List of ClustalW Errors)
	"""
	outputs, errors = [], []
	i, l			= 1, len(sequences)
	print("Status: starting multiple sequence alignment\n")
	for sequence in sequences:
		command 		= ClustalW("clustalw2", infile=sequence["path"])
		output, error 	= command()
		outputs.append({
			"gene" 	: sequence["gene"], 
			"output": output
		})
		if len(error) > 0:
			errors.append(error)
		print("Status:", sequence["gene"], "alignment complete -", i, "of", l)
		i += 1
	if len(errors) > 0:
		print("\nStatus: alignments completed with", len(errors), "errors\n")
	else: 
		print("\nStatus: alignments completed without errors\n")
	return outputs, errors

# 
# Analysis Methods
# 

# Compute Hamming Distance and normalized Hamming Score
def getHammingScore(outputs, path="outputs/scores/"):
	"""
	Compute Hamming Distance and normalized Hamming Score
	@params:
		outputs 			- Required 	: alignment outputs (List of Dicts)
		path 				- Optional 	: file output path (Str)
	Returns: hammingStatistics (List of Dicts)
	"""
	hammingStatistics = []
	print("Status: calculating hamming distances\n")
	for alignments in outputs: 
		for aligned in alignments:
			print("Status: calculating hamming distance for", aligned["gene"])
			alignment = AlignIO.read(path + aligned["gene"] + ".aln", "clustal")
			columns = []
			i 		= 0
			count 	= len(alignment)
			while i < len(alignment[0].seq):
				columns.append(alignment[:, i])
				i += 1
			sequenceHamming = []
			length 			= float(len(columns))
			for column in columns:
				l 		= set(column)
				hamming = len(l) - 1
				sequenceHamming.append(hamming)
			hammingDistance = sum(sequenceHamming)
			hammingScore 	= (hammingDistance / length) / count
			hammingStatistics.append({
				"gene" 		: aligned["gene"],
				"hamming"	: hammingDistance,
				"score" 	: hammingScore
			})
			print("Status: done\n")
	print("Status: hamming calculations done\n")
	return hammingStatistics

# Compute p-value using Mann-Whitney U text
def computeStatistics(statistics, path="outputs/"):
	"""
	Compute p-value using Mann-Whitney U text
	@params:
		hammingStatistics 	- Required 	: hamming statistics by group (List of Dicts)
		path 				- Optional 	: output directory (Str)
	"""
	scores = []
	for group in statistics:
		groupScores = [gene["score"] for gene in group["genes"]]
		scores.append({
			"group": group["group"],
			"score": groupScores
		})
	u, pValue = Stats.mannwhitneyu(scores[0]["score"], scores[1]["score"])
	return u, pValue

# Generate box plots from hamming scores
def generateBoxPlots(analysis, path="outputs/"):
	"""
	Generate box plots from hamming scores
	@params:
		analysis 			- Required 	: analysis data (List of Dicts)
		path 				- Optional 	: output directory (Str)
	"""
	u, pValue, statistics 	= analysis["u"], analysis["pValue"], analysis["statistics"]
	boxplots, labels 		= [], []
	n 						= len(statistics)
	figure 					= Plot.figure()
	axes 					= figure.add_subplot(111)
	for group in statistics:
		groupScores = [gene["score"] for gene in group["genes"]]
		boxplots.append(Npy.array(groupScores))
		labels.append(group["group"])
	data = Npy.array(boxplots)
	axes.boxplot(data)
	axes.set_title("Box Plots of Normalized Hamming Scores by Gene Expression Group\nOne-Sided p-Value = " + str(pValue) + ", U Statistic = " + str(u), fontsize=12)
	axes.set_xticklabels(labels)
	Plot.xlabel("Gene Expression Groups")
	Plot.ylabel("Normalized Hamming Score (0.0 - 1.0)")
	Plot.savefig(path + "boxplots.png", format="PNG")

# 
# Main Routine
# 

if __name__ == "__main__":
	print("Status: part_II.py initialized from commandline\n")
	exitCode 	= 0
	E.email		= getUserInput(valid=r"[^@]+@[^@]+\.[^@]+", prompt="Input: enter email for Entrez queries: ")
	print("Status: configuring")
	try:
		# Parse Config & Set Global Variables
		config 				= parseJSONToDicts("config.json")
		db 					= config["db"]
		outputDirectory 	= config["outputDirectory"]
		geneDirectory 		= config["geneDirectory"]
		scoreDirectory 		= config["scoreDirectory"]
		statsDirectory 		= config["statsDirectory"]
		organisms 			= config["organisms"]
		evidenceLevels 		= config["evidenceLevels"]
		debug 				= config["debug"]
		pathways 			= parseJSONToDicts(outputDirectory + "/targetpathways.json")
		statistics 			= []
		option 				= 0
		numOptions 			= str(len(pathways) - 1)
		print("Status: done\n")
	except:
		print(("Error: unable to configure part II of workflow\n"))
		S.exit(exitCode)
	# User Target Pathway Selection
	print("Input: select pathway for analysis (0 - " + numOptions + "):")
	for pathway in pathways:
		print("	Option:", "[", option, "]", pathway["id"], "-", pathway["desc"])
		option += 1
	selection 	= getUserInput(valid=r"([0-9]{1,2}|\*)", prompt="Input: enter [ n ] to perform analysis on pathway or [ * ] to analyze all pathways: ")
	ncores 		= int(getUserInput(valid=r"[1-9]{1}[0-9]{0,1}", prompt="Input: enter [ n ] number of cores: "))
	print("Status: input received")
	# Core Workflow
	for group in config["groups"]:
		groupDirectory 		= group
		queries 			= formatQuery(pathways, geneKey=group, selection=selection, organisms=organisms, type="mRNA", seqFilter="RefSeq")
		sequenceRecords 	= getSequenceRecords(queries, evidenceLevels=evidenceLevels, db=db)
		sequences, errs 	= getSequences(sequenceRecords, db=db, path=outputDirectory + "/" + geneDirectory + "/" + groupDirectory + "/")
		alignments 			= msaSequences(sequences)
		hammingStatistics 	= getHammingScore(alignments, path=outputDirectory + "/" + geneDirectory + "/" + groupDirectory + "/")
		# Write Scores 
		writeJSON(hammingStatistics, path=outputDirectory + "/" + scoreDirectory + "/" + db + "-" + groupDirectory + "-hammingstatistics.json")
		print("Status: hamming statistics saved to", outputDirectory + "/" + scoreDirectory + "/" + db + "-" + groupDirectory + "-hammingstatistics.json")
		statistics.append({
			"group" : group,
			"genes": hammingStatistics
		})
	u, pValue = computeStatistics(statistics, outputDirectory + "/" + statsDirectory + "/")
	analysis = {
		"pValue"		: pValue,
		"u" 			: u, 
		"statistics" 	: statistics
	}
	writeJSON(analysis, path=outputDirectory + "/" + statsDirectory + "/analysis.json")
	print("Status: analysis saved to", outputDirectory + "/" + statsDirectory + "/analysis.json\n")
	generateBoxPlots(analysis, path=outputDirectory + "/" + statsDirectory + "/")
	print("Status: boxplots saved to", outputDirectory + "/" + statsDirectory + "/boxplots.png\n")
	print("Status: workflow completed successfully")
	exitCode = 1
	S.exit(exitCode)
else:
	pass