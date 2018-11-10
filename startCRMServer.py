#!/usr/bin/env python
import warnings
import sys
import argparse
parser = argparse.ArgumentParser(description='Starts a server that takes streams of sequences (list of lists) and returns model predictions.')
parser.add_argument('-M',dest='loadModel',	metavar='<restoreSavedModel>',help='Restore graph and parameters from file', required=True);
parser.add_argument('-p',dest='port',	metavar='<port>',help='What port to use [default=0]', default=0, required=False);
parser.add_argument('-s',dest='host',	metavar='<host>',help='What host to use [default=localhost]', default="localhost", required=False);
parser.add_argument('-t',dest='threads',	metavar='<threads>',help='Number of threads to make use of [default=1]',default = "1", required=False);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings [default=stderr]', required=False);
parser.add_argument('-b',dest='batch',	metavar='<batchSize>',help='Batch size to use to feed into the neural net [default = 1000] ', required=False, default = "1000");
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);
#sys.argv = "predictThermodynamicEnhancosomeModel.py	-res ../../PBound_learned_params_pTpA_ACP_OHC.txt.ckpt -i 20160613_RefinedTSSs_Scer_20160203_R64_goodR_TSS-before-ATG_-190_to_-80.OHC.txt.gz -o test.txt	-v -v -v	-t 1 -b 1024	-ntm -po -nm 244 -sl 110".split();
#sys.argv = "predictThermodynamicEnhancosomeModel.py -i 20161024_average_promoter_ELs_per_seq_3p1E7_Gly_ALL.shuffled_OHC_test.txt.gz -res EBound_progressive_learning_pTpA_Gly.ACPM.ckpt -o test.txt -v -v -v -t 1 -b 1024 -sl 110 -nm 245 -ml 25 -po".split();
#sys.argv = "predictThermodynamicEnhancosomeModel.py -i 20161024_average_promoter_ELs_per_seq_3p1E7_Gly_ALL.shuffled_OHC_test.txt.gz -res EBound_progressive_learning_pTpA_Gly.A.ckpt -o test.txt -v -v -v -t 1 -b 1024 -sl 110 -nm 245 -ml 25 -ntm -ntc".split();
#sys.argv = "predictThermodynamicEnhancosomeModel.py -i ../../../20160525_NextSeq_pTpA_and_Abf1TATA/analysis/tensorflow/20160503_average_promoter_ELs_per_seq_atLeast100Counts.OHC.txt.gz -res EBound_progressive_learning_pTpA_Gly.A.ckpt -o test.txt -v -v -v -t 1 -b 1024 -sl 110 -nm 245 -ntc -ntm -ic /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_minKds_polyA_and_FZF1_justKds.txt".split();
args = parser.parse_args();
import MYUTILS
import PWM;
import tensorflow as tf
import numpy as np;

args.batch = int(args.batch);
args.port = int(args.port);
args.threads = int(args.threads);


import SETUPOHCENHANCOSOMEMODEL;
#import CisRegModels;

import SocketServer
import json
import numpy as np
BASEORDER = list("ACGT")

def seqsToOHC(seqX):
	ohcX = np.zeros((len(seqX),4,len(seqX[0]),1))
	for i in range(0,4):#bases
		for j in range(0,len(seqX)):
			ohcX[j,i,[x==BASEORDER[i] for x in seqX[j][:]],:]=1;
	return (ohcX)

MAXBYTES=10000000

class PredictionServer(SocketServer.BaseRequestHandler):
	global args
	global myCRM;
	global MAXBYTES
	def getPredictions(self,batchX):
		global args
		global myCRM;
		ohcX = seqsToOHC(batchX);
		numSeqs=np.shape(ohcX)[0]
		if args.verbose>1: sys.stdout.write("Predicting on %i seqs\n"%numSeqs);
		predY = np.zeros(numSeqs);
		z=0
		while z < numSeqs:
			curLast = min(numSeqs,(z+args.batch));
			curPredY = myCRM.predictThese(ohcX[z:curLast,:,:,:])
			predY[z:curLast] = curPredY; 
			z=curLast;
		return (predY);
	def handle(self):
		curRecv = self.request.recv(MAXBYTES).strip()
		while True:
			try:
				self.data = json.loads(curRecv)
				break;
			except ValueError:
				if args.verbose>2: sys.stdout.write("combining packets\n")
				curRecv = curRecv + self.request.recv(MAXBYTES).strip()
		#print self.data
		pred = self.getPredictions(self.data)
		#print pred
		self.response = json.dumps(list(pred))
		self.request.sendall(self.response)
	#def __init__(self, args):
	#def __init__(self, request, client_address, server):
	#	super().__init__(request, client_address, server)
	

#myCRM = CisRegModels.SETUPOHCENHANCOSOMEMODEL.CRM(args);
myCRM = SETUPOHCENHANCOSOMEMODEL.CRM(args);
myCRM.loadModel()
server = SocketServer.TCPServer((args.host, args.port), PredictionServer)
print "Server starting on %s on port %i. Press Ctrl+C to interrupt..."%server.server_address
server.serve_forever()



##TO USE:
#s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
#s.connect((args.host, args.port))
#s.send(json.dumps(offspring))
#preds = json.loads(s.recv(1024))
