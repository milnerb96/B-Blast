#!/usr/bin/env python
#Blaze Milner
#DNA Sequence Project

from Bio.Seq import Seq #Seq allows for the creation of a sequence object to store the DNA string
from Bio.Blast import NCBIWWW #Allows access to NCBI wesbite
from Bio.Blast import NCBIXML #Allows access to BLAST database

from flask import Flask, render_template, request

import csv
import os
import json

class blastSequence:
	def __init__(self, alignment = " ", length = " ", eValue = " ", startValue = " ", endValue = " ", score = " ",
				 num_alignments = " ", identities = " ", positives = " ", strand = " ",
				 frame = " ", query = " ", match = " ", subject = " "):

		self.alignment = alignment
		self.length = length
		self.eValue = eValue
		self.startValue = startValue
		self.endValue = endValue
		self.score = score
		self.num_alignments = num_alignments
		self.identities = identities
		self.positives = positives
		self.strand = strand
		self.frame = frame
		self.query = query
		self.match = match
		self.subject = subject


def BlastParameterMenu(parameterBool):
	print("\nYou can modify your BLAST search in a variety of ways. Below is a list of different parameters.\nIf you would like to add a particular item, please choose the corresponding number.")

	while (True):
		print("\n~~~~~ Parameter List~~~~~")
		print("1. Hitlist Size")
		print("2. Use Megablast")
		print("3. Alignments")
		print("4. Exit")

		userChoice = int(raw_input("Enter your choice: "))

		if (userChoice == 1):
			print("\n~~~~~ HITLIST ~~~~~")
			print("1. Modify the Hitlist Size")
			print("2. Exit")

			userModification = int(raw_input("Enter your choice: "))

			if (userModification == 1):
				newHitlistSize = int(raw_input("Enter a size: "))
				parameterBool[0] = newHitlistSize
				continue

			if (userModification == 2):
				continue

			else:
				continue


		if (userChoice == 2):
			print("\n~~~~~ MEGABLAST ~~~~~")
			print("1. Use Megablast")
			print("2. Do not use Megablast")
			print("3. Exit")

			userModification = int(raw_input("Enter your choice: "))

			if (userModification == 1):
				parameterBool[1] = True
				continue

			if (userModification == 2):
				parameterBool[1] = False
				continue

			if (userModification == 3):
				continue

			else:
				continue

		if (userChoice == 3):
			print("\n~~~~~ ALIGNMENT ~~~~~")
			print("1. Change Alignment Return Size")
			print("2. Exit")

			userModification = int(raw_input("Enter you choice: "))

			if (userModification == 1):
				newAttributeSize = int(raw_input("Enter a new size: "))
				parameterBool[2] = newAttributeSize
				continue

			if (userModification == 2):
				continue

			else:
				continue

		else:
			break

	return parameterBool




def KeyBoardSequence(userSequence):

	mySeq = Seq(userSequence)

	return mySeq

def FastaFileSequence(userSeqeunce):

	record = open(userSeqeunce).read()

	return record

def CheckFileIsEmpty(fileName):

	isEmpty = os.stat(fileName).st_size == 0

	return isEmpty


def RunNCBIComparison(mySeq, parameterList):
	result_handle = NCBIWWW.qblast("blastn", "nt", mySeq, hitlist_size=parameterList[0], megablast=parameterList[1], alignments=parameterList[2])
	# Searches the NSCBI BLAST database with the nucleotide sequence in my_seq

	return result_handle

def PrintResults(resultHandle):
	blast_records = NCBIXML.parse(resultHandle)
	# The returned information from result_handle is a messy XML file.
	# This returns a cleaner version that is easier to read.

	blastList = []
	for blast_record in blast_records:
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				if hsp.expect is not None:
					# print("\n********************** Alignment **********************")
					# print("Sequence: {}".format(alignment.title))
					# print("Length: {}".format(alignment.length))
					# print("e value: {}\n".format(hsp.expect))
					# print("Start Value: {}".format(hsp.query_start)),
					# print("End Value: {}".format((hsp.query_start + hsp.align_length) - 1))

					blast = blastSequence(alignment.title, alignment.length, hsp.expect, hsp.query_start, hsp.query_start + hsp.align_length - 1,
										  hsp.score, hsp.num_alignments, hsp.identities, hsp.positives,
										  hsp.strand, hsp.frame, hsp.query, hsp.match, hsp.sbjct)
					blastList.append(blast)
				else:
					print("No matches!")

	return blastList

def WriteHeaderToFile(fileName, sequence):

		f = open(fileName, "w+")
		dnaWriter = csv.writer(f, delimiter=",")
		dnaWriter.writerow([sequence, "Alignment", "Length", "E Value", "Start Value", "End Vale", "Score", "Bits", "Number of Alignments", "Identities", "Positives", "Gaps", "Strand", "Frame", "Query", "Match", "Subject", "Subject Start"])

		f.close()


def WriteResultsToFile(fileName, blastList):

	with open(fileName, mode="a+") as dnaFile:
		dnaWriter = csv.writer(dnaFile, delimiter=",")
		#WriteHeaderToFile(fileName, blastList)
		for blast in blastList:
			dnaWriter.writerow([" ", blast.alignment, blast.length, blast.eValue, blast.startValue,blast.startValue + blast.length - 1, blast.score, blast.num_alignments,blast.identities, blast.positives, blast.strand, blast.frame, blast.query,blast.match, blast.subject])
	dnaFile.close()

def main():


	# with open('/Users/BlazeMilner 1/PycharmProjects/Capstone/sequence.json') as json_file:
	# 	data = json.load(json_file)
	# 	seq = data['sequence']
	# 	print(seq)


	userChoice = int(raw_input("Would you like to enter DNA Sequence by keyboard(1) or FASTA file(2)? "))

	if (userChoice == 1):
		#userSeq = raw_input("Enter the sequence: ")
		print("\nEnter the sequence. Please press once you are finished.")
		lines = []
		while True:
			line = raw_input()
			if line:
				lines.append(line)
			else:
				break

		text = "".join(lines)
		print(text)
		seq = KeyBoardSequence(text)

	if (userChoice == 2):
		userSeq = raw_input('\nEnter the path of the FASTA file.\nTo do this, locate your FASTA file, right click, hold option, and choose Copy "File Name" as Pathname.\n') or "/Users/BlazeMilner 1/Desktop/T2SOfl9 Reference ReZero (1).fasta"
		seq = FastaFileSequence(userSeq)

	# else:
	# 	print("Please restart program. Your choice is invalid.")

	userFile = raw_input('\nEnter the path of the csv file you would like to store the results in.\nTo do this, create a csv file, locate it in Finder, right click, hold option, and choose Copy "File Name" as Pathname.\n') or "/Users/BlazeMilner 1/Desktop/myResults.csv"

	parameterList = [999, False, 100]

	parameterList = BlastParameterMenu(parameterList)

	print("\nYour search is underway.")

	blastXMl = RunNCBIComparison(seq, parameterList)
	blastList = PrintResults(blastXMl)

	# for blast in blastList:
	# 	print blast.alignment

	WriteHeaderToFile(userFile, seq)
	WriteResultsToFile(userFile, blastList)


	# data = {}
	# data['blast_sequence'] = []
	#
	# for blast in blastList:
	# 	data['blast_sequence'].append({
	# 		'alignment': blast.alignment
	# 	})

	# with open('/Users/BlazeMilner 1/PycharmProjects/Capstone/results.json', 'w') as outfile:
	# 	json.dump(data, outfile)

	#print("\nYour BLAST search is complete!")

main()

print("Your search is complete. All results are located in the file you provided.")

# app = Flask(__name__)
#
# @app.route('/')
# def index():
#     return render_template('queryPage.html')
#
# @app.route('/foo.html', methods=['GET', 'POST'])
# def foo():
#     # execute whatever code you want when the button gets clicked here
# 	main()
# 	return "Your search is complete!"
#
# if __name__ == '__main__':
#     app.run()