#!/usr/bin/env python
import warnings
import MYUTILS
import sys
import argparse
parser = argparse.ArgumentParser(description='Takes as input one or more lists of ports of servers, weights for training, and parameters for optimization and outputs a collection of sequences and outputs.')
parser.add_argument('-ids',dest='ids',	metavar='<ids>',help='Model IDs, comma separated', required=True);
parser.add_argument('-p',dest='ports',	metavar='<ports>',help='Model ports, comma separated, same order as -ids', required=True);
parser.add_argument('-hosts',dest='hosts',	metavar='<hosts>',help='Model hosts, comma separated, same order as -ids - can also be one entry for all [default="127.0.0.1"]', default = "127.0.0.1", required=False);
parser.add_argument('-w',dest='weights',	metavar='<weights>',help='Model weights, comma separated, same order as -ids', required=True);
parser.add_argument('-b',dest='batches',	metavar='<batches>',help='Batches of sequences', required=True);
parser.add_argument('-nb',dest='seqsPerBatch',	metavar='<seqsPerBatch>',help='Sequences taken per batch', required=True);
parser.add_argument('-pop',dest='popSize',	metavar='<popSize>',help='Population size [default=100]',default="100", required=False);
parser.add_argument('-gen',dest='gen',	metavar='<generations>',help='Number of generations [default=100]',default="100", required=False);
parser.add_argument('-mr',dest='mutationRate',	metavar='<mutationRate>',help='Mutation rate [default=0.02]',default="0.02", required=False);
parser.add_argument('-rf',dest='recombFreq',	metavar='<recombFreq>',help='Recombination rate [default=0.1]',default="0.1", required=False);
parser.add_argument('-ws',dest='weightedSum', action='count',help='Use a weighted sum for fitness rather than NSGA-II (multi-objective)?', required=False, default=0);
parser.add_argument('-o',dest='outFP', metavar='<outFile>',help='Where to output results [default=stdout]', required=False);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings [default=stderr]', required=False);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);
#args = lambda: None;	args.inFP = "thisScript.py -i ...".split() # debug settings

args = parser.parse_args();

import numpy as np;
from deap import creator, base, tools, algorithms
from tqdm import tqdm
import socket
import json
args.ports = [int(x) for x in args.ports.split(",")]
args.ids = args.ids.split(",")
args.hosts = args.hosts.split(",")
args.weights = [float(x) for x in args.weights.split(",")]

assert len(args.ids)==len(args.weights), "Must be same number of IDs as weights" 
assert len(args.ids)==len(args.ports), "Must be same number of IDs as ports" 
assert len(args.ids)==len(args.weights), "Must be same number of IDs as weights" 
assert len(args.ids)==len(args.hosts) or len(args.hosts)==1, "Must be same number of IDs as hosts, or only one host" 



if (args.logFP is not None):
	logFile=MYUTILS.smartGZOpen(args.logFP,'w');
	sys.stderr=logFile;

if (args.outFP is None):
	outFile= sys.stdout;
else:
	if args.verbose>0: sys.stderr.write("Outputting to file "+args.outFP+"\n");
	outFile = MYUTILS.smartGZOpen(args.outFP,'w');


argsGA    = {'sequence_length' : 80 , 'nucleotide_frequency' :[0.25,0.25,0.25,0.25] } 
LHS = list("TGCATTTTTTTCACATC")
RHS = list("GGTTACGGCTGTT");
MAXBYTES=1000000

args.batches = int(args.batches);
args.mutationRate = float(args.mutationRate);
args.recombFreq = float(args.recombFreq);
args.popSize = int(args.popSize);
args.seqsPerBatch = int(args.seqsPerBatch);
args.gen = int(args.gen);

randomizer=np.random
def random_sequence_generator(randomizer,args) :
	return randomizer.choice(list('ACGT') , p=args['nucleotide_frequency'] )


def mutation(individual, indpb):
	for i in xrange(len(individual)):
		if np.random.random() < indpb:
			if individual[i]=='A' :
				individual[i] = (randomizer.choice(list('CGT') , p=[ argsGA['nucleotide_frequency'][1]/(1- argsGA['nucleotide_frequency'][0]) , argsGA['nucleotide_frequency'][2]/(1- argsGA['nucleotide_frequency'][0]) , argsGA['nucleotide_frequency'][3]/(1- argsGA['nucleotide_frequency'][0]) ] ) )
			elif individual[i]=='C' :
				individual[i] = (randomizer.choice(list('AGT') , p=[ argsGA['nucleotide_frequency'][0]/(1- argsGA['nucleotide_frequency'][1]) , argsGA['nucleotide_frequency'][2]/(1- argsGA['nucleotide_frequency'][1]) , argsGA['nucleotide_frequency'][3]/(1- argsGA['nucleotide_frequency'][1]) ] ) )
			elif individual[i]=='G' :
				individual[i] = (randomizer.choice(list('CGT') , p=[ argsGA['nucleotide_frequency'][2]/(1- argsGA['nucleotide_frequency'][2]) , argsGA['nucleotide_frequency'][1]/(1- argsGA['nucleotide_frequency'][2]) , argsGA['nucleotide_frequency'][3]/(1- argsGA['nucleotide_frequency'][2]) ] ) )
			elif individual[i]=='T' : 
				individual[i] = (randomizer.choice(list('CGT') , p=[ argsGA['nucleotide_frequency'][0]/(1- argsGA['nucleotide_frequency'][3]) , argsGA['nucleotide_frequency'][1]/(1- argsGA['nucleotide_frequency'][3]) , argsGA['nucleotide_frequency'][2]/(1- argsGA['nucleotide_frequency'][3]) ] ) )
	return individual,

def getPrediction(sequences, socketInfo):
	global MAXBYTES
	s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	s.connect(socketInfo)
	s.send(json.dumps(sequences))
	predStr = s.recv(MAXBYTES);
	nerr=0
	while True:
		try:
			preds = json.loads(predStr)
			break;
		except ValueError as e:
			nerr=nerr+1
			if nerr>10:
				print predStr;
				raise e
			predStr = predStr + s.recv(MAXBYTES).strip()
	s.close();
	return preds;

def fitness(sequences):
	global socketInfo, RHS, LHS
	seqsFull= []
	for ind in range(len(sequences)) :
		seqsFull.append(LHS + sequences[ind] + RHS )
	allPreds = [[] for x in range(0,len(seqsFull))]
	for i in range(0,len(socketInfo)):
		curPreds = getPrediction(seqsFull, socketInfo[i])
		for j in range(0,len(seqsFull)):
			allPreds[j].append(curPreds[j]);
	for j in range(0,len(seqsFull)):
		allPreds[j] = tuple(allPreds[j])
	return(allPreds)

def weightedSum(sequences):
	preds = fitness(sequences);
	global args
	preds2 = [];
	for i in range(0, len(preds)):
		preds2.append((sum(np.array(preds[i])*np.array(args.weights)),))
	return preds2;

if args.weightedSum>0:
	creator.create("FitnessMax", base.Fitness, weights=((1.0,)))
else:
	creator.create("FitnessMax", base.Fitness, weights=(tuple(args.weights)))
#creator.create("FitnessMax", base.Fitness, weights=(0,))
creator.create("Individual", list , fitness=creator.FitnessMax)

toolbox = base.Toolbox()

toolbox.register("base", random_sequence_generator , randomizer , argsGA)
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.base, n=argsGA['sequence_length'])
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

if args.weightedSum>0:
	toolbox.register("evaluate", weightedSum)
else:
	toolbox.register("evaluate", fitness)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", mutation, indpb=args.mutationRate)
if len(args.weights)==1 or args.weightedSum>0:
	toolbox.register("select", tools.selTournament, tournsize=3)
else:
	toolbox.register("select", tools.selNSGA2)



outFile.write("seq\tbatch\t%s\n"%("\t".join(args.ids)));

if len(args.hosts)< len(args.ports):
	args.hosts = args.hosts*len(args.ports);

socketInfo = zip(args.hosts,args.ports)
args.weights = np.array(args.weights)
if args.verbose>1: sys.stderr.write("Running optimization\n")
for b in tqdm(range(args.batches)):
	population = toolbox.population(n=args.popSize)
	fits = toolbox.evaluate(population)
	for fit, ind in zip(fits, population): # below is where Fitness.setValues is actually run
		ind.fitness.setValues(fit)# somehow this sets ind.fitness to the unweighted values,
	if args.verbose>1: sys.stderr.write("batch=%i\n"%(b))
	for gen in range(args.gen):
		if args.verbose>1: sys.stderr.write("\tgen=%i\n"%(gen)) 
		#make offspring from survivors of last round
		offspring = algorithms.varAnd(population, toolbox, cxpb=args.recombFreq, mutpb=1.0)
		#print "".join(offspring[0])
		#calculate fitness of offspring
		fits = toolbox.evaluate(offspring)
		if args.verbose>1:
			if args.weightedSum>0:
				sys.stderr.write("\tBest fitness: %g\n"%(max([f for f in fits])));
			else:
				sys.stderr.write("\tBest fitness: %g\n"%(max([np.sum(np.array(f)*args.weights) for f in fits])));
		for fit, ind in zip(fits, offspring): # below is where Fitness.setValues is actually run
			ind.fitness.setValues(fit)# somehow this sets ind.fitness to the unweighted values,
		population = toolbox.select(population+offspring, k=len(population))
	if len(args.weights)==1 or args.weightedSum>0:
		top = tools.selBest(population, k=args.seqsPerBatch)
	else:
		top = tools.selNSGA2(population, k=args.seqsPerBatch)
	topFitnesses = fitness(top)
	top = ["".join(x) for x in top]
	for i in range(0, len(top)):
		outFile.write("%s\t%i\t%s\n"%(top[i], b, "\t".join([str(x) for x in topFitnesses[i]])))

	
#raise Exception("Reached bad state=%d for '%s.%d' '%s' at line '%s'" %(state,mid,ver,tfid,line));
outFile.close();
if (args.logFP is not None):
	logFile.close();
