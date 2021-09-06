#!/usr/bin/env python
import warnings
import MYUTILS
import sys
import argparse
parser = argparse.ArgumentParser(description='Takes as input a file containing optimized PWMs (one file), and outputs multiple files each with a single matrix')
parser.add_argument('-i',dest='inFP',	metavar='<inFile>',help='Input file of PKdMs', required=True);
parser.add_argument('-in',dest='inFPNames',	metavar='<inFileNames>',help='Input file of names for the matrices', required=False);
parser.add_argument('-d',dest='hasHead' , help='inFileNames has a header', action='count');
parser.add_argument('-o',dest='outFPre', metavar='<outFPre>',help='Where to output results', required=True);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings [default=stderr]', required=False);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);
#args = lambda: None;  args.inFP = "thisScript.py -i ...".split() # debug settings

args = parser.parse_args();


inFile=MYUTILS.smartGZOpen(args.inFP,'r');


if (args.logFP is not None):
	logFile=MYUTILS.smartGZOpen(args.logFP,'w');
	sys.stderr=logFile;

#raise Exception("Reached bad state=%d for '%s.%d' '%s' at line '%s'" %(state,mid,ver,tfid,line));
if args.inFPNames is not None:
	inFileNames=MYUTILS.smartGZOpen(args.inFPNames,'r');
	motNames= [];
	for line in inFileNames:
		if line is None or line == "": continue
		if (args.hasHead>0):
			args.hasHead=0;
			continue;
		data=line.rstrip().split("\t");
		motNames.append(data[0]);
	inFileNames.close();
	if args.verbose>0: sys.stderr.write(", ".join(motNames)+"\n");
z=0
state=0
for line in inFile:
	if line is None or line == "": continue
	if line[0:1]!="#": continue;
	line = line[1:]; #remove the #
	if state==0:
		if args.inFPNames is not None:
			sys.stderr.write("z=%i; line=%s"%(z,line));
			curName = motNames[z];
		else:
			curName = line.rstrip().replace(" ","_");
		outFile = MYUTILS.smartGZOpen(args.outFPre+"/"+curName+".pkdm",'w');
	else:
		outFile.write(line);
	state = state+1;
	if state==5:
		state=0;
		outFile.close();
		z=z+1

inFile.close();
if (args.logFP is not None):
	logFile.close();
