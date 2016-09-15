#!/usr/bin/env python
"""
This script takes a fasta file (including the reference sequence) and produces a specified number of filtered, subsampled fasta files.  
The filter step is done by phydms_prepalign. 
The subsampling is down by randomizing the sequences and then taking the first instance of each year. 

Run script to get help message with info on command::

    python randomize_divpressures.py
    
SKH 08-18-2016
"""

import argparse
import sys
from Bio import SeqIO
import os
import random
import re

class ArgumentParserNoArgHelp(argparse.ArgumentParser):
	"""Like *argparse.ArgumentParser*, but prints help when no arguments."""
	def error(self, message):
		sys.stderr.write('error: %s\n\n' % message)
		self.print_help()
		sys.exit(2)


def ParseArguments():
	"""Argument parser for script."""
	parser = ArgumentParserNoArgHelp(
			description='Filter and subsample a fasta file',
			formatter_class=argparse.ArgumentDefaultsHelpFormatter,
			)
	parser.add_argument('alignment', help='File containing the alignment (including the reference sequence).')
	parser.add_argument('referenceheader', help='Name of reference sequence used for the DMS')
	parser.add_argument('--purge', help = "File with one substring per line of sequences to be purged")
	parser.add_argument('alignmentnumber', help = "The final number of filter, subsampled fasta files.")
        return parser

def oneSequencePerYear(fastaFile, seed):
	"""
	This function randomize the order of sequences in a fasta file and then selects one sequence per year. 
	"""
	random.seed(seed)
	sequences = SeqIO.parse(open(fastaFile),'fasta')
	sequences = [x for x in sequences]
	random.shuffle(sequences)	
	
	dates = []
	finalSeq = []
	noRegex = []
	for seq in sequences:
		regex = r"(/\d+_+\d+/)"
		match = re.search(regex, seq.id)
		if match:
			match = re.search(r"(\d\d\d\d)", (match.group(0)))
			year =  match.group(0)
		else:
			noRegex.append(seq.id)
		if int(year) not in range(1918,2017):
			noRegex.append(seq.id)
		elif int(year) not in dates:
			dates.append(int(year))
			finalSeq.append(seq)
	dates.sort(key=int)
	print "Couldn't parse ", noRegex
	
	outputName = os.path.basename(fastaFile)+"_%s_final.fasta"%(seed)
	print "The output file %s has %s sequences from %s to %s" %(outputName,str(len(finalSeq)),str(dates[0]), str(dates[-1]))
	print
	output_handle = open(outputName, "w")
	SeqIO.write(finalSeq, output_handle, "fasta")
	output_handle.close()
	
def main():
	args = vars(ParseArguments().parse_args())
	print args
	
	##phydms_prepalign
	prepalignOutputName = os.path.basename(args["alignment"])+"_prep.fasta"
	print prepalignOutputName
	if args["purge"]:
		command = "phydms_prepalignment %s %s %s --purgeseqs %s" %(args["alignment"], prepalignOutputName,args["referenceheader"],args["purge"])
	else:
		command = "phydms_prepalignment %s %s %s" %(args["alignment"], prepalignOutputName,args["referenceheader"])
	print command
	os.system(command)
	
	## subsample, 1 per year
	for seed in range(int(args["alignmentnumber"])):
		oneSequencePerYear(prepalignOutputName,seed)

	
	
	
if __name__ == '__main__':
	main()
