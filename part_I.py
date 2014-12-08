#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Developed With Python Version 2.7.8 by Benjamin Cordier

Refer to config.txt dictionary for setup variables

Part I of workflow:
	- Configures directory structure
	- Retrieves input data
	- Performs odds ratio calculations on genes 
	- Writes targetpathways.json output file

Files created:
	- target pathways dictionary file (txt)
"""

import sys 		as S
import csv 		as C
import os 		as OS
import json 	as JSON
import errno

# 
# Utility Methods
# 

# Parse List of Strings to Dicts
def parseJSONToDicts(path):
	"""
	Parse from strings to dictionaries
	@params:
		path 			- Required 	: input path (Str)
	Returns: dicts (List)
	"""
	dicts = []
	with open(path, "r") as file:
		dicts = JSON.load(file)
	return dicts

# JSON File Writing Function
def writeJSON(data, path="output/output.json"):
	"""
	Write JSON output to file.
	@params:
		data 			- Required	: data to write (Dict)
		path 			- Optional 	: output path (Str)
	Returns: path (Str)
	"""
	with open(path, "w") as file:
		file.write(JSON.dumps(data, indent=4, separators=(',', ': '), sort_keys=False))
		file.close()
	return path

# Make Directories From Path
def makeDirectories(path):
	"""
	Create directories from path
	@params:
		path 			- Required 	: directory path (Str)
	Returns: path (Str)
	"""
	try:
		OS.makedirs(path)
	except OSError as exception:
		if exception.errno != errno.EEXIST:
			raise
	return path

# General Data Setting Function
def setData(file, delimiter="\t"):
	"""
	General data setting function.
	@params:
		file:  			- Required 	: file path (str)
		delimiter: 		- Optional 	: file delimiter (Str)
	Returns: dict (Dict)
	"""
	with open(file) as data:
		dict = list(C.DictReader(data, delimiter=delimiter))
	return dict

# Set Kegg Pathway Members Data
def setPathwayData(file, delimiter="\t"):
	"""
	Kegg Pathway data setting function.
	@params:
		file:  			- Required 	: file path (str)
		delimiter: 		- Optional 	: file delimiter (Str)
	Returns: dictList (List of Dicts), columns (List)
	"""
	dictList = []
	with open(file, "r") as data:
		columns = data.readline().rstrip().split(delimiter)
		lines = [line.rstrip().split(delimiter) for line in data]
	data.close()
	for line in lines:
		dict = {}
		i = 0
		for column in columns:
			marker = len(columns) - 1
			if i == marker:
				dict[column] = line[marker:]
				dictList.append(dict)
				break
			else:
				dict[column] = line[i]
				i += 1
	return dictList, columns

# 
# Data Manipulation / Computation Methods
# 

# Get Unique Genes From File, Data Param Is List of Dict
def getUniqueGenes(data, key, returnList=False):
	"""
	Kegg Pathway data setting function.
	@params:
		data:  			- Required 	: file path (str)
		key: 			- Required 	: key to gene name (Str)
		returnList 		- Optional 	: (Bool)
	Returns: out (List)
	"""
	uniqueGenes = {dict[key] for dict in data}
	if returnList == True:
		uniqueGenes = list(uniqueGenes)
	return uniqueGenes

def generateMetaData(deGenes, universeGenes, pathwayData, pathwayIdKey, pathwayDescKey, pathwayMemKey):
	"""
	Generate pathway meta data.
	@params:
		degGenes:  		- Required 	: differentially expressed genes (List, Set)
		universeGenes: 	- Required 	: all probed genes (List, Set)
		pathwayData 	- Required 	: pathway data set (List of Dicts)
		pathwayIdKey 	- Required 	: key value to access pathway ID (Str)
		pathwayDescKey 	- Required 	: key value to access pathway Description (Str)
		pathwayMemKey 	- Required 	: key value to access pathway Members (Str)
	Returns: pathwayMetaData (List of Dicts)
	"""
	deGenes, universeGenes 	= set(deGenes), set(universeGenes)
	deUnionUniverse 		= set.intersection(deGenes, universeGenes)
	nonDeIntersectUniverse 	= set.difference(universeGenes, deUnionUniverse)
	pathwayMetaData 		= []
	for pathway in pathwayData:
		pathwayGenes 				= []
		path 						= {}
		pathwayGenes 				= set(pathway[pathwayMemKey])
		A 							= float(len(list(set.intersection(pathwayGenes, deUnionUniverse))))
		B 							= float(len(list(set.intersection(pathwayGenes, nonDeIntersectUniverse))))
		C 							= float(len(list(set.difference(deUnionUniverse, pathwayGenes))))
		D 							= float(len(list(set.difference(nonDeIntersectUniverse, pathwayGenes))))
		oddsRatio 					= (A * D) / (B * C)
		path["id"] 					= pathway[pathwayIdKey]
		path["desc"] 				= pathway[pathwayDescKey]
		path["pathwayGenes"] 		= pathway[pathwayMemKey]
		path["targetDeGenes"] 		= A
		path["targetNonDeGenes"] 	= B
		path["nonTargetDeGenes"] 	= C
		path["nonTargetNonDeGenes"] = D
		path["oddsRatio"] 			= oddsRatio
		path["DE"]					= list(set.intersection(pathwayGenes, deUnionUniverse))
		path["NonDE"]				= list(set.intersection(pathwayGenes, nonDeIntersectUniverse))
		pathwayMetaData.append(path)
	return pathwayMetaData

def oddsRatioThreshold(pMetaData, threshold=1.0):
	"""
	Returns pathways that meet or exceed odds ratio threshold
	@params:
		pMetaData:  	- Required 	: pathways meta data (List of Dicts)
		threshold:		- Optional 	: threshold for output inclusion (Float)
	Returns: passedThreshold (List of Dicts)
	"""
	passedThreshold = []
	for pathway in pMetaData: 
		if pathway["oddsRatio"] >= threshold:
			pathway["threshold"] = threshold
			passedThreshold.append(pathway)
	for pathway in passedThreshold:
		print "Status: target pathway -", pathway["id"], "-", pathway["desc"]
	return passedThreshold

# 
# Main Routine
# 

if __name__ == "__main__":
	print "\nStatus: part_I.py initialized from commandline\n"
	exitCode 						= 0
	try:
		# Parse Config & Set Global Variables
		print "Status: configuring"
		config 					= parseJSONToDicts("config.json")
		diffExpressedDataPath 	= config["diffExpressedDataPath"]
		universeDataPath 		= config["universeDataPath"]
		pathwaysDataPath 		= config["pathwaysDataPath"]
		outputDirectory  		= config["outputDirectory"]
		geneDirectory 			= config["geneDirectory"]
		scoreDirectory 			= config["scoreDirectory"]
		statsDirectory 			= config["statsDirectory"]
		groups 					= config["groups"]
		targetOddsThreshold 	= config["targetOddsThreshold"]
		print "Status: differentially expressed genes data location - ", diffExpressedDataPath
		print "Status: universe genes data location - ", universeDataPath
		print "Status: KEGG pathways data location - ", pathwaysDataPath
		print "Status: done\n"
	except:
		print ("Error: unable to configure part I of workflow; try validating config.json\n")
		S.exit(exitCode)
	print "Status: building directory structure"
	try: 
		# Scaffold Directory
		makeDirectories(outputDirectory)
		makeDirectories(outputDirectory + "/" + geneDirectory)
		makeDirectories(outputDirectory + "/" + scoreDirectory)
		makeDirectories(outputDirectory + "/" + statsDirectory)	
		for group in groups:
			makeDirectories(outputDirectory + "/" + geneDirectory + "/" + group)
		print "Status: done\n"
	except:
		print ("Error: unable scaffold directory structure\n")
		S.exit(exitCode)
	print "Status: beginning initial assessment"
	try:
		# Set Data
		deData 					= setData(diffExpressedDataPath)
		universeData 			= setData(universeDataPath)
		pathwayData, columns 	= setPathwayData(pathwaysDataPath)
	except:
		print ("Error: unable to retrieve data files, ensure directory paths in config.json are correct\n")
		S.exit(exitCode)
	try:
		# Perform Set Computations & Generate Pathway Meta Data
		deGenes 				= getUniqueGenes(deData, key="gene", returnList=True)
		universeGenes 			= getUniqueGenes(universeData, key="gene", returnList=True)
		pathwayMetaData 		= generateMetaData(deGenes=deGenes, universeGenes=universeGenes, pathwayData=pathwayData, pathwayIdKey="KEGG PATHWAY ID", pathwayDescKey="KEGG PATHWAY TITLE", pathwayMemKey="PATHWAY MEMBERS")
	except:
		print ("Error: unable to generate pathway meta analysis data\n")
		S.exit(exitCode)
	try: 
		# Set Target Pathways and Save
		targetPathways 			= oddsRatioThreshold(pathwayMetaData, targetOddsThreshold)
		writeJSON(targetPathways, outputDirectory + "/targetpathways.json")
		print "Status: target pathways saved to", outputDirectory + "/targetpathways.json"
		print "Status: done\n"
	except: 
		print ("Error: unable to save target pathways data\n")
		S.exit(exitCode)
	exitCode = 1
	S.exit(exitCode)
else:
	pass