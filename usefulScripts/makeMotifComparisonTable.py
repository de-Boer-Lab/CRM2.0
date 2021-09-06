#!/usr/bin/env python
import warnings
import MYUTILS
import sys
import argparse
parser = argparse.ArgumentParser(description='Takes as input a list of motif annotations (file), and a comma separated list of prefixes and another of suffixes, to PNG files to put in the table..')
parser.add_argument('-i',dest='inFP',	metavar='<inFile>',help='Input file of motif annotations', required=True);
parser.add_argument('-d',dest='header',	metavar='<inFile>',help='Comma-separated list of header names - otherwise taken from inFile frist line', required=False);
parser.add_argument('-ip',dest='inFPre',	metavar='<inFilePre>',help='Input PNG prefixes, comma separated', required=True);
parser.add_argument('-is',dest='inFSuf',	metavar='<inFileSuf>',help='Input PNG suffixes, comma separated', required=True);
parser.add_argument('-o',dest='outFP', metavar='<outFile>',help='Where to output results (html) [default=stdout]', required=False);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings [default=stderr]', required=False);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);
#args = lambda: None;  args.inFP = "thisScript.py -i ...".split() # debug settings

args = parser.parse_args();


inFile=MYUTILS.smartGZOpen(args.inFP,'r');


if (args.logFP is not None):
	logFile=MYUTILS.smartGZOpen(args.logFP,'w');
	sys.stderr=logFile;

if (args.outFP is None):
	outFile= sys.stdout;
else:
	if args.verbose>0: sys.stderr.write("Outputting to file "+args.outFP+"\n");
	outFile = MYUTILS.smartGZOpen(args.outFP,'w');

#raise Exception("Reached bad state=%d for '%s.%d' '%s' at line '%s'" %(state,mid,ver,tfid,line));
motData = [];
headerLine = [];
if (args.header is not None):
	headerLine = args.header.split(",");
for line in inFile:
	if line is None or line == "" or line[0]=="#": continue
	data=line.rstrip().split("\t");
	if (args.header is None):
		headerLine = data
		args.header=True;
		continue;
	motData.append(data);
inFile.close();

inPrefixes = args.inFPre.split(",");
inSuffixes = args.inFSuf.split(",");
outFile.write(" <!DOCTYPE html><html> <head> <script type=\"text/javascript\" src=\"../../sorttable.js\"></script> </head> <body> <table border=\"1\"  id=\"motifTable\" class=\"sortable\"> <thead> \n")
outFile.write("<tr><th>"+"</th><th>".join(headerLine)+"</th>\n")
for i in range(0,len(inPrefixes)):
	outFile.write("<th>%s*%s</th>"%(inPrefixes[i], inSuffixes[i]));
outFile.write("</thead> <tbody> \n")
outFile.write("</tr>\n");
for j in range(0,len(motData)):
	outFile.write("<tr><td>"+"</td><td>".join(motData[j])+"</td>");
	for i in range(0,len(inPrefixes)):
		outFile.write("<td><img src=\"%s%s%s\" height = 100/></td>"%(inPrefixes[i],motData[j][0],inSuffixes[i]));
	outFile.write("</tr>\n");
outFile.write("</tbody></table> </body></html>\n")
outFile.close();
if (args.logFP is not None):
	logFile.close();
