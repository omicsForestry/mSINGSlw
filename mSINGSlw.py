#!/usr/bin/env python3

# version 0.1.1

import argparse
import pysam
import numpy
import sys
import os

parser = argparse.ArgumentParser(description='Lightweight implementation of MSIngs')
parser.add_argument('-c', '--controls', help = 'Input list of control bam files')
parser.add_argument('-p', '--prebuilt', help = 'Prebuilt control file, created by -b option. Either -c or -p must be specified.')
parser.add_argument('-b', '--build', help = 'Filename to save control information from files provided by -c option.')
parser.add_argument('-s', '--samples', help = 'Sample bam files to test in tab separated: file, sampleID format.')
parser.add_argument('-r', '--regions', help = 'Regions to test. Must be specified if -c option is used.')
parser.add_argument('-o', '--output', help = 'Output file to save if -s option is used to specify samples.')
parser.add_argument('-d', '--depth', type=int, default=30, help = 'Minimum depth required for processing repeat')

args = parser.parse_args()

# check that arguments do not contradict each other
def checkArgs(theseArgs):
	if not theseArgs.controls and not theseArgs.prebuilt:
		print('Either --controls or --prebuilt option must be stated', file=sys.stderr)
		quit()
	if theseArgs.build and not theseArgs.controls:
		print('--build option requires --controls option as well', file=sys.stderr)
		quit()
	if (theseArgs.controls and not theseArgs.regions) or (not theseArgs.controls and theseArgs.regions):
		print('--controls and --regions options must be specified together', file=sys.stderr)
		quit()
	if theseArgs.output and not theseArgs.samples:
		print('--output option requires --samples option as well', file=sys.stderr)
		quit()
	if theseArgs.samples and not theseArgs.output:
		print('--samples option requires --output option as well', file=sys.stderr)
		quit()
	if theseArgs.prebuilt and not theseArgs.samples:
		print('Loading a prebuilt control file, but no samples to process', file=sys.stderr)
		quit()
	if theseArgs.controls and not theseArgs.samples and not theseArgs.build:
		print('Control file specified, but no output region file to make or samples to process', file=sys.stderr)
		quit()


# function to report the length of a repeat given the surrounding sequence
def findRep(seq, l_flank, r_flank, rep):
	found=False
	if l_flank in seq:
		if r_flank in seq.split(l_flank)[1]:
			n=seq.split(l_flank)[1].split(r_flank)[0].count(rep)
			if l_flank+(n*rep)+r_flank in seq:
				found=True
	if found:
		return(n)
	else:
		return('NA')

# function to measure the number and distribution of a repeat in a bam file
def countRep(samfile, chr, pos, l_flank, r_flank, rep, genCount):
	d={}
	tot=0
	for read in samfile.fetch(chr, (pos-2), pos+(len(rep)*genCount)+2):
		seq=read.seq
		n=findRep(seq, l_flank, r_flank, rep)
		if n!='NA':
			if n not in d:
				d[n]=0
			d[n]+=1
			tot+=1
	maxTot=-1
	keyMax=0
	for i in d:
		if d[i]>maxTot:
			maxTot=d[i]
			keyMax=i
	l=[]
	for i in d:
		if d[i] > maxTot*0.05:
			l.append(d[i]/float(maxTot))
	if l!=[]:
		return([tot, float(maxTot)/tot, len(l), numpy.mean(l), numpy.std(l), l])
	else:
		return([0, 0, 0, 0, 0, l])

# function to collect all the repeats for one file
def countAllRepsBam(testBam, regionDict, depth):
	regRes={}
	samfile=pysam.AlignmentFile(testBam, 'rb')
	for i in regionDict:
		info=regionDict[i]
		res=countRep(samfile, info[0], int(info[1]), info[8], info[9], info[7], int(info[4]))
		if res[0] >= depth and res[1] >= 0.05:
			regRes[i]=res[2]
		else:
			regRes[i]='NA'
	samfile.close()
	return(regRes)


# collect all the repeats for multiple samples
def countAllRepsAllSamples(sampleDict, regionDict, depth):
	repsRes={}
	for sample in sampleDict:
		testBam=sampleDict[sample]
		repsRes[sample]=countAllRepsBam(testBam, regionDict, depth)
	return(repsRes)

# collect controls and remove failing samples and regions. Return distributions, ready for saving
def collectControls(conListFile, regionDict, depth):
	conDict={}
	infile=open(conListFile)
	for line in infile:
		bam=line.split()[0]
		conDict[bam]=bam
	infile.close()
	conReps=countAllRepsAllSamples(conDict, regionDict, depth)
	regRes={}
	for i in regionDict:
		regRes[i]=[]
	for i in conDict:
		for j in regRes:
			if conReps[i][j]!='NA':
				regRes[j].append(conReps[i][j])
	maxControls=0
	for i in regionDict:
		if len(regRes[i])>maxControls:
			maxControls=len(regRes[i])
	for i in regionDict:
		if len(regRes[i])<(maxControls/2):
			del regRes[i]
	regSumRes={}
	for i in regRes:
		regSumRes[i]=[numpy.mean(regRes[i]), numpy.std(regRes[i])]
	return(regSumRes)

# save collected control information to a file
def saveControls(conRes, regionDict, controlFile):
	outfile=open(controlFile, 'w')
	print('chromosome\tlocation\trepeat_unit_length\trepeat_unit_binary\trepeat_times\tleft_flank_binary\tright_flank_binary\trepeat_unit_bases\tleft_flank_bases\tright_flank_bases\tmean\tSD', file=outfile)
	for region in regionDict:
		if region in conRes:
			print('\t'.join(map(str,regionDict[region]+conRes[region])), file=outfile)
		else:
			print('\t'.join(regionDict[region])+'\tNA\tNA', file=outfile)
	outfile.close()


# test one sample
def testOneSample(conRes, repsDict, sample):
	tot=0
	unstable=0
	testSites=0
	for i in conRes:
		testSites+=1
		thisCutOff = conRes[i][0]+3*conRes[i][1]
		if i in repsDict[sample]:
			if repsDict[sample][i]!='NA':
				tot+=1
				if repsDict[sample][i] > thisCutOff:
					unstable+=1
	if tot > 0.1*testSites:
		res=100*(unstable/tot)
	else:
		res='NA'
	return(res)

# test list of samples from file, saving to new file
def testAllSamples(sampleFile, outputFile, regionDict, conRes, depth):
	infile=open(sampleFile)
	sampleDict={}
	for line in infile:
		field=line.split()
		sampleDict[field[1]]=field[0]
	infile.close()
	repsDict=countAllRepsAllSamples(sampleDict, regionDict, depth)
	outfile=open(outputFile,'w')
	print('sample\tscore', file=outfile)
	for sample in repsDict:
		thisRes=testOneSample(conRes, repsDict, sample)
		print(sample+'\t'+str(thisRes), file=outfile)
	outfile.close()

# collect region information plus controls if needed. File has header:
# chromosome\tlocation\trepeat_unit_length\trepeat_unit_binary\trepeat_times\tleft_flank_binary\tright_flank_binary\trepeat_unit_bases\tleft_flank_bases\tright_flank_bases

def collectRegions(regFile, getControls=False):
	infile=open(regFile)
	regionDict={}
	conRes={}
	line=infile.readline()
	for line in infile:
		field=line.split()
		id=field[0]+':'+field[1]
		if getControls:
			regionDict[id]=field[:-2]
			if field[-2]!='NA':
				conRes[id]=list(map(float,field[-2:]))
		else:
			regionDict[id]=field
	infile.close()
	return([regionDict, conRes])


checkArgs(args)

if args.controls:
	regionDict=collectRegions(args.regions)[0]
	controlRes=collectControls(args.controls, regionDict, args.depth)
elif args.prebuilt:
	regConRes=collectRegions(args.prebuilt, getControls=True)
	regionDict=regConRes[0]
	controlRes=regConRes[1]

if args.build:
	saveControls(controlRes, regionDict, args.build)

if args.output:
	testAllSamples(args.samples, args.output, regionDict, controlRes, args.depth)
