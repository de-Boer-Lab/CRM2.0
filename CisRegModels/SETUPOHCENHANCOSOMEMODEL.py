import MYUTILS
import PWM;
import TFHELP;
import tensorflow as tf
import numpy as np;
from datetime import datetime
import os;
import sys

BASES = ['A','C','G','T']

class CRM:
	#global self.args, self.potentiation, self.constantPot, self.positionalActivityBias, self.sess, self.train_step, self.myLoss, self.mseTF, self.mergedSummaries, self.global_step, self.ohcX, self.realELY, self.summaryWriter,self.epBoundTensor, self.epBoundTensorRC, self.seqPotentialTensor, self.motifsTensor2, self.logConcs, self.activities, self.constant, self.saver, self.saveParams, self.positionalActivityBiasRC, self.activityDiffs, BASES
	
	def __init__(self, args2):
		self.args=args2;
	
	def loadModel(self):
		if self.args.loadModel is None:
			raise Exception("Must have -M set for loadModel to work");
		sys.stderr.write("Loading saved model: %s\n"%(self.args.loadModel));
		self.sess = tf.Session()
		tf.saved_model.loader.load(self.sess, ['main'], self.args.loadModel)
		varNames = [v.name for v in tf.global_variables()];
		self.ohcX = tf.get_default_graph().get_tensor_by_name("ohcX:0")
		try: 
			self.predELY = tf.get_default_graph().get_tensor_by_name("predELY:0")
		except KeyError:
			self.predELY = tf.get_default_graph().get_tensor_by_name("predELY")
		if not hasattr(self.args, 'outputKdMat'):
			self.args.outputKdMat=0;
		if  self.args.outputKdMat>0:
			self.motifScan = tf.get_default_graph().get_tensor_by_name("motifScan:0")
			self.motifScanRC = tf.get_default_graph().get_tensor_by_name("motifScanRC:0")
		if not hasattr(self.args, 'outputBinding'):
			self.args.outputBinding=0;
		if  self.args.outputBinding>0:
			self.epBoundTensor = tf.get_default_graph().get_tensor_by_name("epBoundTensor:0")
			if "epBoundTensorRC:0" in varNames:
				self.args.trainStrandedActivities=1
				self.epBoundTensorRC = tf.get_default_graph().get_tensor_by_name("epBoundTensorRC:0")
		if "seqPotentialTensor:0" in varNames:
			self.args.potentiation=1
			self.seqPotentialTensor = tf.get_default_graph().get_tensor_by_name("seqPotentialTensor:0")
		if self.args.dropoutTF is not None:
			if self.args.VARIABLE:
				try: 
					self.logConcs =  [var for var in tf.global_variables() if var.op.name=="Variable:0"][0] 
				except IndexError:
					self.logConcs =  [var for var in tf.global_variables() if var.op.name=="Variable"][0] 
			else:
				try: 
					self.logConcs =  [var for var in tf.global_variables() if var.op.name=="concentrations:0"][0]
				except IndexError:
					self.logConcs =  [var for var in tf.global_variables() if var.op.name=="concentrations"][0] 

	def predictThese(self, batchX): 
		return(self.predELY.eval(session=self.sess, feed_dict={self.ohcX: batchX}))
	def testModel(self):
		if (self.args.logFP is not None):
			logFile=MYUTILS.smartGZOpen(self.args.logFP,'w');
			sys.stderr=logFile;
		
		if (self.args.outFP is not None):
			if self.args.verbose>0: sys.stderr.write("Outputting to file "+self.args.outFP+"*\n");
			outFile=MYUTILS.smartGZOpen(self.args.outFP,'w');
		else:
			outFile = sys.stdout;
		if self.args.loadModel is None and self.args.restore is None:
			raise Exception("Must specify either -res or -M")
		if self.args.loadModel is not None:
			self.loadModel();
		else:
			self.makeGraph();
		if self.args.dropoutTF is not None:
			self.args.dropoutTF = int(self.args.dropoutTF);
			newConcs = self.logConcs.eval(session=self.sess);
			newConcs[0,0,0,self.args.dropoutTF]=-99999;
			self.logConcs.load(newConcs, self.sess)
		outFile.write("actual\tpredicted");
		if self.args.potentiation>0:
			outFile.write("\tpred_openness");
		
		if self.args.outputBinding > 0:
			if self.args.trainPositionalActivities > 0:
				raise Exception("Cannot output binding values while using positional self.activities");
			outFile.write("\t" + "\t".join(["Binding_%i"%x for x in range(0,self.args.numMotifs)]))
			if self.args.trainStrandedActivities>0:
				outFile.write("\t" + "\t".join(["RCBinding_%i"%x for x in range(0,self.args.numMotifs)]))
		if self.args.outputKdMat>0:
			sys.stderr.write("Outputting Kd header\n")
			for j in range(0,self.args.seqLen):
				outFile.write("\t" + "\t".join(["Kd_f%i_p%i"%(x,j) for x in range(0,self.args.numMotifs)]))
			for j in range(0,self.args.seqLen):
				outFile.write("\t" + "\t".join(["KdRC_f%i_p%i"%(x,j) for x in range(0,self.args.numMotifs)]))
			
		outFile.write("\n");
		
		b = 0
		batchX = np.zeros((self.args.batch,4,self.args.seqLen,1))
		batchY = np.zeros((self.args.batch))
		inFile=MYUTILS.smartGZOpen(self.args.inFP,'r');
		for line in inFile:
			if line is None or line == "" or line[0]=="#": continue
			curData = np.fromstring(line, dtype=float, sep="\t")
			batchY[b]=curData[0];
			batchX[b,:,:,0] = curData[1:].reshape((4,self.args.seqLen))
			b+=1
			if b==self.args.batch:
				#curPredY  = self.sess.run([predELY], feed_dict={self.ohcX: batchX})
				curPredY = self.predELY.eval(session=self.sess, feed_dict={self.ohcX: batchX})
				if self.args.outputBinding>0:
					bindingAmount = self.epBoundTensor.eval(session=self.sess, feed_dict={self.ohcX: batchX});
					if self.args.trainStrandedActivities>0:
						bindingAmountRC = self.epBoundTensorRC.eval(session=self.sess, feed_dict={self.ohcX: batchX});
				if self.args.potentiation>0:
					curPredOpenness = self.seqPotentialTensor.eval(session=self.sess, feed_dict={self.ohcX: batchX})
				if self.args.outputKdMat>0:
					## [None,1,seqLen,numMotifs]
					curMotifScan = self.motifScan.eval(session=self.sess, feed_dict={self.ohcX: batchX})
					curMotifScanRC = self.motifScanRC.eval(session=self.sess, feed_dict={self.ohcX: batchX})
				for i in range(0,batchY.shape[0]):
					outFile.write("%g\t%g"%(batchY[i],curPredY[i]));
					if self.args.potentiation>0:
						outFile.write("\t%g"%(curPredOpenness[i]));
					if self.args.outputBinding > 0:
						outFile.write("\t%s"%("\t".join(["%g"%ba for ba in bindingAmount[i,]])));
						if self.args.trainStrandedActivities>0:
							outFile.write("\t%s"%("\t".join(["%g"%ba for ba in bindingAmountRC[i,]])));
					if self.args.outputKdMat>0:
						for j in range(0,self.args.seqLen):
							outFile.write("\t" + "\t".join(["%g"%(curMotifScan[i,0,j,x]) for x in range(0,self.args.numMotifs)]))
						for j in range(0,self.args.seqLen):
							outFile.write("\t" + "\t".join(["%g"%(curMotifScanRC[i,0,j,x]) for x in range(0,self.args.numMotifs)]))
					outFile.write("\n");
				b=0;
		inFile.close();

		#test with remaining data, but remove anything past b
		if b > 0:
			batchX = batchX[0:b,:,:,:]
			batchY = batchY[0:b];
			curPredY = self.predELY.eval(session=self.sess, feed_dict={self.ohcX: batchX})
			if self.args.outputBinding>0:
				bindingAmount = self.epBoundTensor.eval(session=self.sess, feed_dict={self.ohcX: batchX});
				if self.args.trainStrandedActivities>0:
					bindingAmountRC = self.epBoundTensorRC.eval(session=self.sess, feed_dict={self.ohcX: batchX});
			if self.args.potentiation>0:
				curPredOpenness = self.seqPotentialTensor.eval(session=self.sess, feed_dict={self.ohcX: batchX})
			if self.args.outputKdMat>0:
				## [None,1,seqLen,numMotifs]
				curMotifScan = self.motifScan.eval(session=self.sess, feed_dict={self.ohcX: batchX})
				curMotifScanRC = self.motifScanRC.eval(session=self.sess, feed_dict={self.ohcX: batchX})
			for i in range(0,batchY.shape[0]):
				outFile.write("%g\t%g"%(batchY[i],curPredY[i]));
				if self.args.potentiation>0:
					outFile.write("\t%g"%(curPredOpenness[i]));
				if self.args.outputBinding > 0:
					outFile.write("\t%s"%("\t".join(["%g"%ba for ba in bindingAmount[i,]])));
					if self.args.trainStrandedActivities>0:
						outFile.write("\t%s"%("\t".join(["%g"%ba for ba in bindingAmountRC[i,]])));
				if self.args.outputKdMat>0:
					for j in range(0,self.args.seqLen):
						outFile.write("\t" + "\t".join(["%g"%(curMotifScan[i,0,j,x]) for x in range(0,self.args.numMotifs)]))
					for j in range(0,self.args.seqLen):
						outFile.write("\t" + "\t".join(["%g"%(curMotifScanRC[i,0,j,x]) for x in range(0,self.args.numMotifs)]))
				outFile.write("\n");
		self.sess.close()
		sys.stderr.write("Done!\n")
		if (self.args.logFP is not None):
			logFile.close();
		outFile.close()
		
	
	def makeModel(self):
		if (self.args.logFP is not None):
			logFile=MYUTILS.smartGZOpen(self.args.logFP,'w');
			sys.stderr=logFile;
		
		if self.args.dropoutTFTF is not None:
			self.args.dropoutTFTF = float(self.args.dropoutTFTF);
			sys.stderr.write("Using dropout on TF-TF interactions with probability (%f)\n"%(self.args.dropoutTFTF));
		
		if (self.args.outFPre is not None):
			if self.args.verbose>0: sys.stderr.write("Outputting to file "+self.args.outFPre+"*\n");
		self.makeGraph()
		
		if self.args.verbose>1: sys.stderr.write("saving session... ");
		save_path = self.saver.save(self.sess, self.args.outFPre+".ckpt")
		self.saveParams(self.sess)
		if self.args.verbose>1: sys.stderr.write("done. saved in %s\n"%save_path);
		
		#remove the file that labels the model as complete
		if os.path.isfile(self.args.outFPre+".done"):
			if self.args.verbose>1: sys.stderr.write("Deleting old \"done\" file %s\n"%(self.args.outFPre+".done"));
			os.remove(self.args.outFPre+".done")
		
		if self.args.tensorboard is not None:
			nextTBF = 1;
		
		globalStep=0;
		if self.args.verbose>1: sys.stderr.write("Running through data %i times\n"%self.args.runs)
		
		batchGetter = TFHELP.BatchGetterOneHot(self.args.inFP, self.args.batch, self.args.runs, self.args.seqLen)
		batchX, batchY, runsLeft = batchGetter.getNextBatch();
		runningMeanMSE = np.zeros(self.args.runningAverageWindow);#init to nans
		runningMeanMSE.fill(np.nan)
		runningMSE=0
		lastTime = datetime.now();
		curTimeNonTF = datetime.now()
		#main loop of training and occasionally saving sessions/summaries
		while runsLeft>0:
			#train with current data
			if self.args.verbose>3: sys.stderr.write("  Time for last non-TF code: %f\n"%(datetime.now()-curTimeNonTF).total_seconds())
			curTimeTF = datetime.now()
			if self.args.trace >0 and globalStep>100: #don't take the first one
				run_metadata = tf.RunMetadata()
				_, curLoss, curMSE, summaries, globalStep = self.sess.run([self.train_step, self.myLoss, self.mseTF, self.mergedSummaries, self.global_step], feed_dict={self.ohcX: batchX, self.realELY: batchY}, options=tf.RunOptions(trace_level=tf.RunOptions.FULL_TRACE), run_metadata=run_metadata)
				from tensorflow.python.client import timeline
				trace = timeline.Timeline(step_stats=run_metadata.step_stats)
				trace_file = open('%s.trace.json'%self.args.outFPre, 'w')
				trace_file.write(trace.generate_chrome_trace_format())
				quit();
			elif self.args.tensorboard is not None and ((self.args.tensorboardFrequency is None and nextTBF == (globalStep+1)) or (self.args.tensorboardFrequency is not None and (globalStep+1) % self.args.tensorboardFrequency==0 )):
				sys.stderr.write("  Getting tensorboard summaries\n")
				
				_, curLoss, curMSE, summaries, globalStep = self.sess.run([self.train_step, self.myLoss, self.mseTF, self.mergedSummaries, self.global_step], feed_dict={self.ohcX: batchX, self.realELY: batchY})
			else:
				_, curLoss, curMSE, globalStep = self.sess.run([self.train_step, self.myLoss, self.mseTF, self.global_step], feed_dict={self.ohcX: batchX, self.realELY: batchY})
			if self.args.verbose>3: sys.stderr.write("  Time for last TF code: %f\n"%(datetime.now()-curTimeTF).total_seconds())
			curTimeNonTF = datetime.now()
			runningMeanMSE[globalStep % self.args.runningAverageWindow]=curMSE;
			runningMSE = np.nanmean(runningMeanMSE);
			if self.args.verbose>3 or (self.args.verbose>2 and (globalStep % 100) ==0) or (self.args.verbose>1 and (globalStep % 10000) ==0): 
				curTime = datetime.now()
				deltaTime = (curTime-lastTime).total_seconds()
				lastTime=curTime;
				if self.args.L1 is not None:
					pl1Sum, pl1Num = self.sess.run([self.paramPenaltyL1Tensor, self.paramNumActivityTensor], feed_dict={})
					sys.stderr.write("	Batch = %i; examples = %i;\trunning MSE = %.3f;\tlast MSE = %.2f;\tlast time = %.2f\tloss = %.2f;\t paramsum = %.2f;\tparamnum = %i\n"%(globalStep, globalStep * self.args.batch,  runningMSE,curMSE, deltaTime,curLoss,pl1Sum, pl1Num));
				else:
					sys.stderr.write("	Batch = %i; examples = %i;\trunning MSE = %.3f;\tlast MSE = %.2f\tlast time = %.2f\n"%(globalStep, globalStep * self.args.batch, runningMSE,curMSE, deltaTime));
			#save state, if desired
			if self.args.tensorboard is not None and ((self.args.tensorboardFrequency is None and nextTBF == globalStep) or (self.args.tensorboardFrequency is not None and globalStep % self.args.tensorboardFrequency==0 )):
				nextTBF = globalStep*2;
				if self.args.verbose>1: sys.stderr.write("	saving tensorboard summary %s... "%globalStep);
				self.summaryWriter.add_summary(summaries, globalStep)
				if self.args.verbose>1: sys.stderr.write("done\n");
			if globalStep % self.args.saveEvery==0:
				if np.isnan(curMSE):
					raise Exception("ERROR: reached nan MSE - quitting without saving.");
				if self.args.verbose>1: sys.stderr.write("	saving session... ");
				save_path = self.saver.save(self.sess, self.args.outFPre+".ckpt")
				self.saveParams(self.sess)
				if self.args.verbose>1: sys.stderr.write("done. saved in %s\n"%save_path);
			curTimeB = datetime.now()
			batchX, batchY, runsLeft = batchGetter.getNextBatch();
			if self.args.verbose>3: sys.stderr.write("  Time To Get Batch: %f\n"%(datetime.now()-curTimeB).total_seconds())
			
		#train on the final batch
		if self.args.tensorboard is not None:
			_, curLoss, curMSE, summaries = self.sess.run([self.train_step, self.myLoss, self.mseTF, self.mergedSummaries], feed_dict={self.ohcX: batchX, self.realELY: batchY})
			if self.args.verbose>1: sys.stderr.write("Final saving tensorboard summary... ");
			self.summaryWriter.add_summary(summaries, globalStep)
			if self.args.verbose>1: sys.stderr.write("done\n");
		else:
			self.sess.run(self.train_step, feed_dict={self.ohcX: batchX, self.realELY: batchY})
		
		#final session saving
		if self.args.verbose>1: sys.stderr.write("Final saving session... ");
		save_path = self.saver.save(self.sess, self.args.outFPre+".ckpt")
		self.saveParams(self.sess)
		
		self.saveModel()
		# Add a second MetaGraphDef for inference.
		
		open(self.args.outFPre+".done", 'a').close()
		
		self.sess.close()
		if (self.args.logFP is not None):
			logFile.close();
	def describeModel(self):
		if self.args.loadModel is not None:
			self.loadModel();
		else:
			self.makeGraph();
		graph = self.sess.graph
		print([node.name for node in graph.as_graph_def().node])
	def saveModel(self):
		if self.args.verbose>1: sys.stderr.write("saving model... ");
		modelSaveDir = "%s.Model"%self.args.outFPre;
		if os.path.exists(modelSaveDir):
			import shutil
			if self.args.verbose>1: sys.stderr.write(" deleting old model directory... ");
			shutil.rmtree(modelSaveDir)
		builder = tf.saved_model.builder.SavedModelBuilder(modelSaveDir)
		builder.add_meta_graph_and_variables(self.sess, ["main"])#, signature_def_map=foo_signatures, assets_collection=foo_assets)
		builder.save()
		if self.args.verbose>1: sys.stderr.write("done. saved in %s\n"%modelSaveDir);
		
	def makeGraph(self):
		sys.stderr.write("Input arguments: " + " ".join(sys.argv)+"\n")
		
		if self.args.accIsAct>0 and self.args.potentiation==0:
			raise Exception("Cannot include accessibility influence on EL without including self.potentiation layer!")
		
		if not hasattr(self.args, 'activityClasses'):
			self.args.activityClasses=1;
		if not hasattr(self.args, 'VARIABLE'):
			self.args.VARIABLE=None;
		if not hasattr(self.args, 'tensorboard'):
			self.args.tensorboard=None;
		if not hasattr(self.args, 'learningRateED'):
			self.args.learningRateED=1.0;
		if not hasattr(self.args, 'cycleLearningRate'):
			self.args.cycleLearningRate=None;
		if not hasattr(self.args, 'trace'):
			self.args.trace=0;
		if not hasattr(self.args, 'useMomentumOptimizer'):
			self.args.useMomentumOptimizer=0;
			self.args.useNesterov=0;
		if not hasattr(self.args, 'useRMSProp'):
			self.args.useRMSProp=0;
			if not hasattr(self.args, 'rmsEpsilon'):
				self.args.rmsEpsilon=1e-10;
			if not hasattr(self.args, 'momentum'):
				self.args.momentum=0.0;
			if not hasattr(self.args, 'rmsDecay'):
				self.args.rmsDecay=0.9;
		if not hasattr(self.args, 'Abeta1'):
			self.args.Abeta1=0.9;
		if not hasattr(self.args, 'Abeta2'):
			self.args.Abeta2=0.999;
		if not hasattr(self.args, 'Aepsilon'):
			self.args.Aepsilon=1E-8;
		if not hasattr(self.args, 'clearAdamVars'):
			self.args.clearAdamVars=0;
		if not hasattr(self.args, 'L2Pos'):
			self.args.L2Pos=0.001;
		if not hasattr(self.args, 'L2'):
			self.args.L2=0.001;
		if not hasattr(self.args, 'L1int'):
			self.args.L1int=0.001;
		if not hasattr(self.args, 'L1'):
			self.args.L1=0.001;
		if not hasattr(self.args, 'learningRate'):
			self.args.learningRate=0.001;
		if not hasattr(self.args, 'meanSDFile'):
			self.args.meanSDFile=None;
		
		### set up tensorflow graph
		#useful: tf.Tensor.get_shape(predELY)
		#sequence layer
		self.ohcX = tf.placeholder(tf.float32, [None,4,self.args.seqLen,1], name="ohcX") # here, [None, 4,...] indicate that the first dimention is of unknown length
		#motif parameters:
		#allow for the user to provide a list of PWMs that we initialize the motifs tensor to
		if self.args.motifsFP is not None:
			if self.args.verbose > 0: sys.stderr.write("Reading in default set of TF Motifs from %s\n"%(self.args.motifsFP));
			#motifsFP is a one-per-line list of PWM paths
			motifsFile=MYUTILS.smartGZOpen(self.args.motifsFP,'r');
			defaultMotifs = np.random.normal(0,1,[4,self.args.motifLen,self.args.numMotifs])
			i=0
			for line in motifsFile:
				if line is None or line == "" or line[0]=="#": continue
				curPWM = PWM.loadPWM(line.rstrip());
				pwmLen = curPWM.len();
				pwmStart = 0;
				pwmEnd = pwmLen; # actually end index+1
				if pwmLen > self.args.motifLen: # can't fit the motif; 
					#trim the terminal bases until pwmEnd ==self.args.motifLen; greedy
					while (pwmEnd - pwmStart) > self.args.motifLen:
						startInfo = 0.0;
						endInfo = 0.0; 
						for b in BASES:
							startInfo += curPWM.mat[b][pwmStart]**2
							endInfo += curPWM.mat[b][pwmEnd-1]**2
						if startInfo > endInfo:
							pwmEnd-=1;
						else:
							pwmStart+=1;
				paddingFront = int((self.args.motifLen - (pwmEnd - pwmStart))/2); #truncated 
				for bi in range(0,len(BASES)):
					b=BASES[bi];
					defaultMotifs[bi,:,i] = 0.0 # set other entries to 0
					defaultMotifs[bi,paddingFront:(paddingFront+pwmEnd-pwmStart),i] = curPWM.mat[b][pwmStart:pwmEnd];
				i+=1;
			if self.args.noTrainMotifs>0:
				motifsTensor = tf.constant(defaultMotifs.astype(np.float32),name="motifs") 
			else:
				motifsTensor = tf.Variable(tf.convert_to_tensor(defaultMotifs.astype(np.float32)),name="motifs") 
		else:
			if self.args.noTrainMotifs>0:
				motifsTensor = tf.constant((np.random.standard_normal([4,self.args.motifLen,self.args.numMotifs])/10).astype(np.float32),name="motifs") 
			else:
				motifsTensor = tf.Variable(tf.random_normal([4,self.args.motifLen,self.args.numMotifs],stddev=0.1),name="motifs") 
		
		motifsRCTensor = tf.image.flip_left_right(tf.image.flip_up_down(motifsTensor)) #flips along first two dims; this works because order of bases is ACGT=rc(TGCA)
		
		self.motifsTensor2 = tf.reshape(motifsTensor, [4,self.args.motifLen,1,self.args.numMotifs]);#patch height, patch width, input channels, output channels
		motifsRCTensor2 = tf.reshape(motifsRCTensor, [4,self.args.motifLen,1,self.args.numMotifs]);
		#output: [4, motifLen, 1, numMotifs]
		
		############## INPUT ANY PROVIDED INITIAL PARAMS
		## ACTIVITIES
		initActiv = np.random.standard_normal(self.args.numMotifs).astype(np.float32);
		#input default concentrations if applicable
		if self.args.initActiv is not None:
			if self.args.verbose>0: sys.stderr.write("Initializing activities to %s\n"%self.args.initActiv)
			try: 
				vals = np.loadtxt(self.args.initActiv, dtype = np.str);
				initActiv[0:vals.shape[0]] = vals[:,1].astype("float32");
			except IOError:
				initActiv = np.zeros([self.args.numMotifs]).astype("float32")+float(self.args.initActiv)
		
		
		if self.args.noTrainActivities==0:
			if self.args.verbose>0: sys.stderr.write("Training activities\n")
			self.activities = tf.Variable(tf.convert_to_tensor(initActiv.reshape([self.args.numMotifs])),name="activities");
		else:
			if self.args.verbose>0: sys.stderr.write("Not training activities\n")
			self.activities = tf.constant(initActiv.reshape([self.args.numMotifs]),name="activities");
		
		if self.args.trainStrandedActivities>0:
			if self.args.verbose>0: sys.stderr.write("Training stranded activities\n")
			self.activityDiffs = tf.Variable(tf.zeros(self.args.numMotifs), name="activityDiffs")
			activitiesRC = tf.add(self.activities, self.activityDiffs);
		
		if self.args.trainPositionalActivities>0:
			if self.args.verbose>0: sys.stderr.write("Training positionally biased activities\n")
			self.positionalActivityBias = tf.Variable(tf.ones([self.args.seqLen,self.args.numMotifs]), name="positionalActivities") #[seqLen, numMotifs]
			positionalActivity = tf.multiply(self.positionalActivityBias, self.activities) #[seqLen, numMotifs]
			if self.args.trainStrandedActivities:
				self.positionalActivityBiasRC = tf.Variable(tf.ones([self.args.seqLen,self.args.numMotifs]), name="positionalActivitiesRC") #[seqLen, numMotifs]
				positionalActivityRC = tf.multiply(self.positionalActivityBiasRC, activitiesRC) #[seqLen, numMotifs]
		
		## CONCENTRATIONS
		initConcs = np.zeros(self.args.numMotifs).astype(np.float32);
		#input default concentrations if applicable
		if self.args.initConcs is not None:
			if self.args.verbose>0: sys.stderr.write("Initializing concentrations to %s\n"%self.args.initConcs)
			vals = np.loadtxt(self.args.initConcs, dtype = np.str);
			initConcs[0:vals.shape[0]] = vals[:,1].astype("float32");
		
		
		if self.args.noTrainConcs==0:
			if self.args.verbose>0: sys.stderr.write("Training concentrations\n")
			if not self.args.VARIABLE:
				self.logConcs = tf.Variable(tf.convert_to_tensor(initConcs.reshape((1,1,1,self.args.numMotifs))), name="concentrations"); 
			else:
				self.logConcs = tf.Variable(tf.convert_to_tensor(initConcs.reshape((1,1,1,self.args.numMotifs))), name="Variable"); 
		else:
			if self.args.verbose>0: sys.stderr.write("Not training concentrations\n")
			self.logConcs = tf.constant(initConcs.reshape((1,1,1,self.args.numMotifs)), name="concentrations"); ###TODO name="concentrations"
		
		## binding limits
		if self.args.bindingLimits:
			initBL = np.ones(self.args.numMotifs).astype(np.float32);
			#input default concentrations if applicable
			if self.args.initBindLim is not None:
				if self.args.verbose>0: sys.stderr.write("Initializing binding limits to %s\n"%self.args.initBindLim)
				vals = np.loadtxt(self.args.initBindLim, dtype = np.str);
				initBL[0:vals.shape[0]] = np.log(vals[:,1].astype("float32"));  #thus they are input as non-log
		
			if self.args.noTrainBL==0:
				if self.args.verbose>0: sys.stderr.write("Training binding limits\n")
				logBindingLimits = tf.Variable(tf.convert_to_tensor(initBL.reshape((self.args.numMotifs))), name="logBindingLimits");  #[numMotifs]
			else:
				if self.args.verbose>0: sys.stderr.write("Not training binding limits\n")
				logBindingLimits = tf.constant(initBL.reshape((self.args.numMotifs)), name="logBindingLimits");  #[numMotifs]
			self.bindingLimits=tf.exp(logBindingLimits, name="self.bindingLimits");
		
		#motif layer: conv layer 1
			#nodes: motifs * orientations * positions
			#params: motifs * motifLens * 4
		#strides all =1 - this is the step size
		# zero padding =SAME makes output dims =input; valid does not pad with 0s
		#motifScanTensor = tf.nn.conv2d(self.ohcX, self.motifsTensor2, strides = [1,1,1,1], padding='VALID', name="motifScan") #VALID so that the output dimensions are 1 * seqLen-motifLen+1, ...
		#motifScanRCTensor= tf.nn.conv2d(self.ohcX,motifsRCTensor2, strides = [1,1,1,1], padding='VALID', name="motifScanRC") #VALID so that the output dimensions are 1 * seqLen-motifLen+1, ...
		#### ##outputs [None,1,seqLen-motifLen+1,numMotifs]
		#these are log(Kds)
		motifScanTensor = tf.nn.conv2d(self.ohcX, self.motifsTensor2, strides = [1,4,1,1], padding='SAME', name="motifScan") 
		motifScanRCTensor= tf.nn.conv2d(self.ohcX,motifsRCTensor2, strides = [1,4,1,1], padding='SAME', name="motifScanRC") 
		##outputs [None,1,seqLen,numMotifs]
		
		logKdConcRatioTensor = tf.subtract(self.logConcs,motifScanTensor) # [None, 1, seqLen,numMotifs] 
		logKdConcRatioRCTensor = tf.subtract(self.logConcs,motifScanRCTensor) # [None, 1, seqLen,numMotifs] 
		pNotBoundTensor1 = tf.div(1.0,tf.add(1.0,tf.exp(logKdConcRatioTensor))); # size: [None,1,seqLen,numMotifs]
		pNotBoundRCTensor1 = tf.div(1.0,tf.add(1.0,tf.exp(logKdConcRatioRCTensor))); # size: [None,1,seqLen,numMotifs]
		
		if self.args.interactions>0:
			#the following code implements TF cooperativity/competition; see lab notebook "How best to model TF-TF interactions"
			#for cooperativity, a Kd scale is added in proportion to the bidning of the cofactor (log scaled)
			#for competition, a negative Kd scale is added that scales with the Kd of the competing factor (standard space)
			#negative coeficients indicate competition between the factors at the specified distance
			#positive coeficients indicate cooperativity
			#pBoundX = 1-(1/(1+exp(logXs-CxyComp*PBY*(logY))*(1+(exp(CxyCoop)-1)*PBY)))
			#these are the kd interaction filters
			#here, first F/R refers to current motif X [dim 3], second F/R refers to interacting motif Y [0];
			#for these parameters, Ccoop<0 implies competition, Ccoop>0 implies cooperation
			#self.kdBonusInteractionsFilterFF = tf.Variable(tf.random_normal([self.args.numMotifs, self.args.interactions,1,self.args.numMotifs],stddev=0.01),name="tftfInteractionsFF") #[numMotifs,interactions,1,numMotifs] = [1,X,interactions,1,Y]
			#self.kdBonusInteractionsFilterFR = tf.Variable(tf.random_normal([self.args.numMotifs, self.args.interactions,1,self.args.numMotifs],stddev=0.01),name="tftfInteractionsFR") 
			#self.kdBonusInteractionsFilterRF = tf.Variable(tf.random_normal([self.args.numMotifs, self.args.interactions,1,self.args.numMotifs],stddev=0.01),name="tftfInteractionsRF") 
			#self.kdBonusInteractionsFilterRR = tf.Variable(tf.random_normal([self.args.numMotifs, self.args.interactions,1,self.args.numMotifs],stddev=0.01),name="tftfInteractionsRR") 
			self.kdBonusInteractionsFilterFF = tf.Variable(tf.fill([self.args.numMotifs, self.args.interactions,1,self.args.numMotifs],-0.001),name="tftfInteractionsFF") #[numMotifs,interactions,1,numMotifs] = [1,X,interactions,1,Y]
			self.kdBonusInteractionsFilterFR = tf.Variable(tf.fill([self.args.numMotifs, self.args.interactions,1,self.args.numMotifs],-0.001),name="tftfInteractionsFR") 
			self.kdBonusInteractionsFilterRF = tf.Variable(tf.fill([self.args.numMotifs, self.args.interactions,1,self.args.numMotifs],-0.001),name="tftfInteractionsRF") 
			self.kdBonusInteractionsFilterRR = tf.Variable(tf.fill([self.args.numMotifs, self.args.interactions,1,self.args.numMotifs],-0.001),name="tftfInteractionsRR") 
			if self.args.dropoutTFTF is not None:
				kdBonusInteractionsFilterFF2 = tf.nn.dropout(self.kdBonusInteractionsFilterFF, 1.0-self.args.dropoutTFTF);
				kdBonusInteractionsFilterFR2 = tf.nn.dropout(self.kdBonusInteractionsFilterFR, 1.0-self.args.dropoutTFTF);
				kdBonusInteractionsFilterRF2 = tf.nn.dropout(self.kdBonusInteractionsFilterRF, 1.0-self.args.dropoutTFTF);
				kdBonusInteractionsFilterRR2 = tf.nn.dropout(self.kdBonusInteractionsFilterRR, 1.0-self.args.dropoutTFTF);
			else:
				kdBonusInteractionsFilterFF2 = self.kdBonusInteractionsFilterFF;
				kdBonusInteractionsFilterFR2 = self.kdBonusInteractionsFilterFR;
				kdBonusInteractionsFilterRF2 = self.kdBonusInteractionsFilterRF;
				kdBonusInteractionsFilterRR2 = self.kdBonusInteractionsFilterRR;
			#weights for cooperativity
			coopWeightsFilterFF = tf.subtract(tf.exp(tf.nn.relu(kdBonusInteractionsFilterFF2)),1, name="coopParamsFF")
			coopWeightsFilterFR = tf.subtract(tf.exp(tf.nn.relu(kdBonusInteractionsFilterFR2)),1, name="coopParamsFR")
			coopWeightsFilterRF = tf.subtract(tf.exp(tf.nn.relu(kdBonusInteractionsFilterRF2)),1, name="coopParamsRF")
			coopWeightsFilterRR = tf.subtract(tf.exp(tf.nn.relu(kdBonusInteractionsFilterRR2)),1, name="coopParamsRR")
			#weights for competition
			compWeightsFilterFF = tf.nn.relu(tf.negative(kdBonusInteractionsFilterFF2), name="compParamsFF") # all positive, so I will subtract these later
			compWeightsFilterFR = tf.nn.relu(tf.negative(kdBonusInteractionsFilterFR2), name="compParamsFR") # all positive, so I will subtract these later
			compWeightsFilterRF = tf.nn.relu(tf.negative(kdBonusInteractionsFilterRF2), name="compParamsRF") # all positive, so I will subtract these later
			compWeightsFilterRR = tf.nn.relu(tf.negative(kdBonusInteractionsFilterRR2), name="compParamsRR") # all positive, so I will subtract these later
			#pBoundX = 1-(1/(1+exp(logXs-compWF[X,I,1,Y]*PBY[,1,L,Y]*(logY[,1,L,Y]))*(1+(coopWF[X,I,1,Y]*PBY[,1,L,Y]))))
			#pBounds
			pBoundTensor1 = tf.subtract(1.0, pNotBoundTensor1)#[None,1,seqLen,numMotifs]
			pBoundRCTensor1 = tf.subtract(1.0, pNotBoundRCTensor1) #[None,1,seqLen,numMotifs] 
			#logY*PBY:
			logKDConcPBYTensor = tf.transpose(tf.multiply(pBoundTensor1, logKdConcRatioTensor), perm = [0,3,2,1])#[None,1,seqLen,numMotifs] -> [None, numMot, seqLen, 1]
			logKDConcPBYRCTensor = tf.transpose(tf.multiply(pBoundRCTensor1, logKdConcRatioRCTensor), perm = [0,3,2,1])#[None,1,seqLen,numMotifs] -> [None, numMot, seqLen, 1]
			pBoundTensor1 = tf.transpose(pBoundTensor1, perm = [0,3,2,1])     #[None, numMot, seqLen, 1]
			pBoundRCTensor1 = tf.transpose(pBoundRCTensor1, perm = [0,3,2,1]) #[None, numMot, seqLen, 1]
			#pBoundX = 1-(1/(1+exp(logXs-compWF[Y,I,1,X]*logYPBY[,Y,L,1]))*(1+(coopWF[Y,I,1,X]]*PBY[,Y,L,1]))))
			#the following contain the weighted cooperative and competitive effects for all potential binding sites and at all interaction distances between X and Y
			#apply convolution here to make these for i in seqLen refer to X positions
			#compEffects = compWF[Y,I,1,X]*logYPBY[,Y,L,1]
			compEffectsFFTensor = tf.nn.conv2d(logKDConcPBYTensor,   compWeightsFilterFF, strides = [1,self.args.numMotifs,1,1], padding='SAME', name="compEffectsFF")  #[None,1,seqLen,X]
			compEffectsRFTensor = tf.nn.conv2d(logKDConcPBYTensor,   compWeightsFilterRF, strides = [1,self.args.numMotifs,1,1], padding='SAME', name="compEffectsRF")  #[None,1,seqLen,X]
			compEffectsFRTensor = tf.nn.conv2d(logKDConcPBYRCTensor, compWeightsFilterFR, strides = [1,self.args.numMotifs,1,1], padding='SAME', name="compEffectsFR")  #[None,1,seqLen,X]
			compEffectsRRTensor = tf.nn.conv2d(logKDConcPBYRCTensor, compWeightsFilterRR, strides = [1,self.args.numMotifs,1,1], padding='SAME', name="compEffectsRR")  #[None,1,seqLen,X]
			#coopEffects = coopWF[Y,I,1,X]]*PBY[,Y,L,1]
			coopEffectsFFTensor = tf.nn.conv2d(pBoundTensor1,   coopWeightsFilterFF, strides = [1,self.args.numMotifs,1,1], padding='SAME', name="coopEffectsFF")  #[None,1,seqLen,X]
			coopEffectsRFTensor = tf.nn.conv2d(pBoundTensor1,   coopWeightsFilterRF, strides = [1,self.args.numMotifs,1,1], padding='SAME', name="coopEffectsRF")  #[None,1,seqLen,X]
			coopEffectsFRTensor = tf.nn.conv2d(pBoundRCTensor1, coopWeightsFilterFR, strides = [1,self.args.numMotifs,1,1], padding='SAME', name="coopEffectsFR")  #[None,1,seqLen,X]
			coopEffectsRRTensor = tf.nn.conv2d(pBoundRCTensor1, coopWeightsFilterRR, strides = [1,self.args.numMotifs,1,1], padding='SAME', name="coopEffectsRR")  #[None,1,seqLen,X]
			#pBoundX = 1-(1/(1+exp(logXs-compEffects)*(1+coopEffects)))
			#these now contain the cooperative and competititve effects for each TF in each position - use the logX to calculate PBound_X
			#add the bonus Kds from neighboring motifs to original motifs
			pNotBoundTensor = tf.div(1.0,tf.add(1.0,tf.multiply(tf.exp(tf.subtract(logKdConcRatioTensor, tf.add(compEffectsFFTensor,compEffectsFRTensor))),tf.add(1.0,tf.add(coopEffectsFFTensor,coopEffectsFRTensor)))))
			pNotBoundRCTensor = tf.div(1.0,tf.add(1.0,tf.multiply(tf.exp(tf.subtract(logKdConcRatioRCTensor, tf.add(compEffectsRFTensor,compEffectsRRTensor))),tf.add(1.0,tf.add(coopEffectsRFTensor,coopEffectsRRTensor)))))
		else:
			pNotBoundTensor = pNotBoundTensor1;
			pNotBoundRCTensor = pNotBoundRCTensor1
		
		
		if self.args.useEBound>0: 
			if self.args.verbose>0: sys.stderr.write("Using E-bound\n")
			self.epBoundTensor = tf.add(tf.reduce_sum(tf.subtract(1.0,pNotBoundRCTensor), reduction_indices=[1,2]),tf.reduce_sum(tf.subtract(1.0,pNotBoundTensor), reduction_indices=[1,2])) # size: [None, numMotifs] #expected amount of binding
		else:
			if self.args.verbose>0: sys.stderr.write("Using P-bound\n")
			self.epBoundTensor = tf.subtract(1.0,tf.multiply(tf.reduce_prod(pNotBoundRCTensor,reduction_indices=[1,2]),tf.reduce_prod(pNotBoundTensor, reduction_indices=[1,2]))) # size: [None, numMotifs] # p(bound)
		
		## POTENTIATION
		if self.args.potentiation>0:
			if self.args.verbose>0: sys.stderr.write("Using potentiation layer\n")
			initPotent = np.random.standard_normal(self.args.numMotifs).astype(np.float32);
			#input default concentrations if applicable
			if self.args.initPotent is not None:
				if self.args.verbose>0: sys.stderr.write("Initializing potentiations to %s\n"%self.args.initPotent)
				try: 
					vals = np.loadtxt(self.args.initPotent, dtype = np.str);
					initPotent[0:vals.shape[0]] = vals[:,1].astype("float32");
				except IOError:
					#initPotent = tf.fill([self.args.numMotifs],float(self.args.initPotent))
					initPotent = np.zeros([self.args.numMotifs]).astype("float32") + float(self.args.initPotent)
			if self.args.noTrainPotentiations==0:
				if self.args.verbose>0: sys.stderr.write("Training potentiations\n")
				self.potentiation = tf.Variable(tf.convert_to_tensor(initPotent.reshape([self.args.numMotifs])),name="potents");
			else:
				if self.args.verbose>0: sys.stderr.write("Not training potentiations\n")
				self.potentiation = tf.constant(initPotent.reshape([self.args.numMotifs]),name="potents");
			seqPotentialByTFTensor = tf.multiply(self.epBoundTensor, self.potentiation); #size: [None,numMotifs]
			self.constantPot = tf.Variable(tf.zeros(1),name="constantPot")
			self.seqPotentialTensor = tf.sigmoid(tf.add(tf.reduce_sum(seqPotentialByTFTensor,reduction_indices=[1]), self.constantPot), name="seqPotentialTensor") #[None, 1]
		else:
			if self.args.verbose>0: sys.stderr.write("Not using potentiation layer\n")
		
		if self.args.trainPositionalActivities>0: # account for positional activity with linear scaling of activity
			if self.args.trainStrandedActivities>0: # account for strand-specific activity biases
				pBoundPerPos = tf.subtract(1.0, pNotBoundTensor) # size: [None,1,seqLen,numMotifs]
				pBoundPerPosRC = tf.subtract(1.0, pNotBoundRCTensor) # size: [None,1,seqLen,numMotifs]
				if self.args.potentiation>0:
					pBoundPerPos = tf.transpose(tf.multiply(tf.transpose(pBoundPerPos, perm=(1,2,3,0)), self.seqPotentialTensor), perm = (3,0,1,2)) # size: None,1,seqLen,numMotifs]
					pBoundPerPosRC = tf.transpose(tf.multiply(tf.transpose(pBoundPerPosRC, perm=(1,2,3,0)), self.seqPotentialTensor), perm = (3,0,1,2)) # size: None,1,seqLen,numMotifs]
				#print(tf.Tensor.get_shape(pBoundPerPos))
				#print(tf.Tensor.get_shape(positionalActivity))
				if self.args.bindingLimits>0:
					expectedActivitySense = tf.reduce_sum(tf.multiply(pBoundPerPos, positionalActivity),reduction_indices=[1,2]) # size: [None,numMotifs]
					expectedActivityRC = tf.reduce_sum(tf.multiply(pBoundPerPosRC, positionalActivityRC),reduction_indices=[1,2]) # size: [None,numMotifs]
					expectedActivityPerTF = tf.add( #min of positive self.activities and max of negative self.activities accounting for binding limits.
						tf.nn.relu(                        tf.minimum(tf.reshape(tf.multiply(self.bindingLimits,self.activities),(1,self.args.numMotifs)),tf.add(expectedActivitySense, expectedActivityRC))), #positive
						tf.negative(tf.nn.relu(tf.negative(tf.maximum(tf.reshape(tf.multiply(self.bindingLimits,self.activities),(1,self.args.numMotifs)),tf.add(expectedActivitySense, expectedActivityRC))))) #negative
						) #[None,numMotifs]
					
				else:
					#expectedActivitySense = tf.multiply(tf.reshape(pBoundPerPos, (-1, self.args.seqLen*self.args.numMotifs)), tf.reshape(positionalActivity, (self.args.seqLen*self.args.numMotifs,1))) # size: [None,numMotifs]
					#expectedActivityRC = tf.multiply(tf.reshape(pBoundPerPosRC, (-1, self.args.seqLen*self.args.numMotifs)), tf.reshape(positionalActivityRC, (self.args.seqLen*self.args.numMotifs,1))) # size: [None,1]
					expectedActivitySense = tf.reduce_sum(tf.multiply(pBoundPerPos, positionalActivity), reduction_indices=[1,2]) # size: [None,numMotifs]
					#print(tf.Tensor.get_shape(expectedActivitySense))
					expectedActivityRC = tf.reduce_sum(tf.multiply(pBoundPerPosRC, positionalActivityRC), reduction_indices=[1,2])  # size: [None,numMotifs]
					expectedActivityPerTF = tf.add(expectedActivitySense, expectedActivityRC);
			else:
				if self.args.useEBound>0:
					pBoundPerPos = tf.add(tf.subtract(1.0,pNotBoundTensor), tf.subtract(1.0,pNotBoundRCTensor)) # size: [None,1,seqLen,numMotifs]
				else:
					pBoundPerPos = tf.subtract(1.0,tf.multiply(pNotBoundTensor, pNotBoundRCTensor)) # size: [None,1,seqLen,numMotifs]
				#print(tf.Tensor.get_shape(pBoundPerPos))
				if self.args.potentiation>0:
					pBoundPerPos = tf.transpose(tf.multiply(tf.transpose(pBoundPerPos, perm=(1,2,3,0)), self.seqPotentialTensor), perm = (3,0,1,2)) # size: [None,1,seqLen,numMotifs]
				if self.args.bindingLimits>0:
					pBoundPerPos = tf.minimum(pBoundPerPos, self.bindingLimits); #[None,1,seqLen,numMotifs]
				#print(tf.Tensor.get_shape(pBoundPerPos))
				#print(tf.Tensor.get_shape(positionalActivity))
				if self.args.bindingLimits>0:
					expectedActivitySense = tf.reduce_sum(tf.multiply(pBoundPerPos, positionalActivity),reduction_indices=[1,2]) # size: [None,numMotifs]
					expectedActivityPerTF = tf.add( #min of positive self.activities and max of negative self.activities accounting for binding limits.
						tf.nn.relu(                        tf.minimum(tf.reshape(tf.multiply(self.bindingLimits,self.activities),(1,self.args.numMotifs)),expectedActivitySense)), #positive
						tf.negative(tf.nn.relu(tf.negative(tf.maximum(tf.reshape(tf.multiply(self.bindingLimits,self.activities),(1,self.args.numMotifs)),expectedActivitySense)))) #negative
						) #[None,numMotifs]
				else:
					expectedActivityPerTF = tf.multiply(tf.reshape(pBoundPerPos, (-1, self.args.seqLen*self.args.numMotifs)), tf.reshape(positionalActivity, (self.args.seqLen*self.args.numMotifs,1))) # size: [None,numMotifs]
		else: #no positional self.activities
			if self.args.trainStrandedActivities>0: # account for strand-specific activity biases
				if self.args.useEBound>0: 
					self.epBoundTensor   = tf.reduce_sum(tf.subtract(1.0,  pNotBoundTensor), reduction_indices=[1,2], name="epBoundTensor") # size: [None, numMotifs] #expected amount of binding
					self.epBoundTensorRC = tf.reduce_sum(tf.subtract(1.0,pNotBoundRCTensor), reduction_indices=[1,2], name="epBoundTensorRC") # size: [None, numMotifs] #expected amount of binding
				else:
					self.epBoundTensor = tf.subtract(1.0,tf.reduce_prod(pNotBoundTensor,reduction_indices=[1,2]), name="epBoundTensor") # numMotifs: [None, numMotifs] # p(bound)
					self.epBoundTensorRC = tf.subtract(1.0,tf.reduce_prod(pNotBoundRCTensor,reduction_indices=[1,2]), name="epBoundTensorRC") # size: [None, numMotifs] # p(bound)
				if self.args.potentiation>0:
					self.epBoundTensor = tf.transpose(tf.multiply(tf.transpose(self.epBoundTensor), self.seqPotentialTensor), name="epBoundTensor"); # [None, numMotifs]
					self.epBoundTensorRC = tf.transpose(tf.multiply(tf.transpose(self.epBoundTensorRC), self.seqPotentialTensor), name="epBoundTensorRC"); # [None, numMotifs]
					#print(tf.Tensor.get_shape(self.epBoundTensor))
					#print(tf.Tensor.get_shape(self.epBoundTensorRC))
				if self.args.bindingLimits>0:#note that when adding strand-specific self.activities when there are binding limits, the output will change without changing the params because you can't limit both simultaneously now.
					expectedActivityPerTF = tf.add( #min of positive self.activities and max of negative self.activities accounting for binding limits.
						tf.nn.relu(                        tf.minimum(tf.reshape(tf.multiply(self.bindingLimits,self.activities),(1,self.args.numMotifs)),tf.add(tf.multiply(self.epBoundTensor, self.activities), tf.multiply(self.epBoundTensorRC, activitiesRC)))), #positive
						tf.negative(tf.nn.relu(tf.negative(tf.maximum(tf.reshape(tf.multiply(self.bindingLimits,self.activities),(1,self.args.numMotifs)),tf.add(tf.multiply(self.epBoundTensor, self.activities), tf.multiply(self.epBoundTensorRC, activitiesRC)))))) #negative
						)#[None,numMotifs]
					#sys.stderr.write(" ".join([str(x) for x in tf.Tensor.get_shape(expectedActivity)])+"\n")
				else:
					expectedActivityPerTF = tf.add(tf.multiply(self.epBoundTensor, tf.reshape(self.activities,(self.args.numMotifs,1))), tf.multiply(self.epBoundTensorRC, tf.reshape(activitiesRC,(self.args.numMotifs,1)))) #[None,numMotifs]
			else: #no positional or strand effects
				if self.args.potentiation>0:
					self.epBoundTensor = tf.transpose(tf.multiply(tf.transpose(self.epBoundTensor),self.seqPotentialTensor), name="epBoundTensor"); # [None,numMotifs]
				if self.args.bindingLimits>0:
					self.epBoundTensor = tf.minimum(self.epBoundTensor, self.bindingLimits, name="epBoundTensor"); #[None, numMotifs]
				#expectedActivityPerTF = tf.multiply(self.epBoundTensor, tf.reshape(self.activities,(self.args.numMotifs,1))); #size: [None,numMotifs]
				expectedActivityPerTF = tf.multiply(self.epBoundTensor, self.activities); #size: [None,numMotifs]
		
		
		if self.args.activityClasses>1:
			if self.args.verbose>0: sys.stderr.write("Using activity classes\n")
			initActivClassWeights = np.random.normal(0,0.01,[self.args.numMotifs,self.args.activityClasses]).astype("float32")+0.01
			#initActivClassWeights = np.zeros([self.args.numMotifs,self.args.activityClasses]).astype("float32")+1.0
			self.activityClassWeights = tf.Variable(tf.convert_to_tensor(initActivClassWeights.reshape([self.args.numMotifs,self.args.activityClasses])),name="activityClassWeights");
			self.activityClassActivationBiases =  tf.Variable(tf.convert_to_tensor(np.random.normal(0,1,[self.args.activityClasses]).astype("float32").reshape([self.args.activityClasses])),name="activityActivationBiases");
			self.activityClassExpressionBiases = tf.Variable(tf.convert_to_tensor(np.random.normal(0,1,[self.args.activityClasses,1]).astype("float32").reshape([self.args.activityClasses,1])),name="activityExpressionBiases");
			# first layer of activity types
			activityClassValues = tf.matmul(expectedActivityPerTF, self.activityClassWeights) #[None, self.args.activityClasses] 
			activityClassActivations = tf.sigmoid(tf.add(activityClassValues, self.activityClassActivationBiases)) #[None, self.args.activityClasses]
			expectedActivity = tf.matmul(activityClassActivations, self.activityClassExpressionBiases) #[None, 1]
		else:
			expectedActivity = tf.reduce_sum(expectedActivityPerTF, reduction_indices=[1])
		
		self.constant = tf.Variable(tf.zeros(1),name="constant")
		
		if self.args.accIsAct>0:
			self.accActivity = tf.Variable(tf.zeros(1),name="accessActiv")
			accActivELTensor = tf.multiply(self.seqPotentialTensor, self.accActivity) # [None,1]
			self.predELY= tf.add(tf.add(tf.reshape(expectedActivity, [-1]),tf.reshape(accActivELTensor, [-1])), self.constant, name="predELY") #size: [None]
		else:
			self.predELY= tf.add(tf.reshape(expectedActivity, [-1]),self.constant, name="predELY") #size: [None]
		self.realELY = tf.placeholder(tf.float32, [None]);
		
		EPSILON=0.0001
		
		self.mseTF = tf.reduce_mean(tf.square(self.realELY - self.predELY))
		if self.args.meanSDFile is not None:
			vals = np.loadtxt(self.args.meanSDFile, dtype = np.str);
			llMeans = vals[:,0].astype("float32");
			llSDs = vals[:,1].astype("float32");
			maxllMean = llMeans.max()
			llSDLen = llSDs.shape[0];
			if self.args.verbose>0: sys.stderr.write("Minimizing the negative log liklihood with an SD vector of length %i and a maximum value of %f\n"%(llSDLen,maxllMean))
			llSDTensor = tf.constant(llSDs, name="llSDs")
			predELLLIndeces = tf.minimum(llSDLen-1, tf.maximum(0,tf.to_int32(tf.round(tf.multiply(self.predELY,(llSDLen/maxllMean))))), name="predELLLInd") 
			predELYSDs = tf.nn.embedding_lookup(llSDTensor, predELLLIndeces, name="predELYSDs") #[None]
			predZ = tf.div(tf.subtract(self.predELY, self.realELY), predELYSDs, name="predZ")
			negLogLik = tf.reduce_sum(tf.square(predZ), name="negLogLik")
			self.myLoss = negLogLik;
		else:
			self.myLoss = self.mseTF;
		if self.args.L2 is not None and self.args.noTrainMotifs==0:
			if self.args.verbose>0: sys.stderr.write("Using L2 regularization of PWMs with lambda=%s\n"%(self.args.L2))
			self.args.L2 = float(self.args.L2);
			paramPenaltyL2Tensor = tf.nn.l2_loss(motifsTensor) #doesn't make sense to l2 loss the concentrations since this brings them to 0; however, bringing the PWM entries to 0 does make sense.
			self.myLoss = tf.add(self.myLoss, tf.multiply(paramPenaltyL2Tensor,self.args.L2));
		
		if self.args.L1 is not None:
			if self.args.verbose>0: sys.stderr.write("Using L1 regularization of activities with lambda=%s\n"%(self.args.L1))
			self.args.L1 = float(self.args.L1);
			self.paramPenaltyL1Tensor = tf.reduce_sum(tf.abs(self.activities))
			self.paramNumActivityTensor = tf.reduce_sum(tf.cast(tf.greater(tf.abs(self.activities),EPSILON),tf.int32))
			if self.args.activityClasses>1:
				self.paramPenaltyL1Tensor = self.paramPenaltyL1Tensor + tf.reduce_sum(tf.abs(self.activityClassWeights))
				self.paramNumActivityTensor = self.paramNumActivityTensor +tf.reduce_sum(tf.cast(tf.greater(tf.abs(self.activityClassWeights),EPSILON),tf.int32))
			if self.args.potentiation>0:
				self.paramPenaltyL1Tensor = self.paramPenaltyL1Tensor + tf.reduce_sum(tf.abs(self.potentiation))
				self.paramNumActivityTensor = self.paramNumActivityTensor +tf.reduce_sum(tf.cast(tf.greater(tf.abs(self.potentiation),EPSILON),tf.int32))
			if self.args.trainStrandedActivities>0:
				self.paramPenaltyL1Tensor = self.paramPenaltyL1Tensor + tf.reduce_sum(tf.abs(self.activityDiffs))
				self.paramNumActivityTensor = self.paramNumActivityTensor +tf.reduce_sum(tf.cast(tf.greater(tf.abs(self.activityDiffs),EPSILON),tf.int32))
			if self.args.accIsAct>0:
				self.paramPenaltyL1Tensor = self.paramPenaltyL1Tensor + tf.abs(self.accActivity)
				self.paramNumActivityTensor = self.paramNumActivityTensor +tf.reduce_sum(tf.cast(tf.greater(tf.abs(self.accActivity),EPSILON),tf.int32))
			if self.args.trainPositionalActivities>0:
				self.paramPenaltyL1Tensor = self.paramPenaltyL1Tensor + tf.reduce_sum(tf.abs(tf.slice(self.positionalActivityBias,[0,0],[1,self.args.numMotifs]))) #penalize only the first column of positional self.activities
				#So only one position per TF is L1; the differences between the others are L2
				self.paramNumActivityTensor = self.paramNumActivityTensor +tf.reduce_sum(tf.cast(tf.greater(tf.abs(tf.slice(self.positionalActivityBias,[0,0],[1,self.args.numMotifs])),EPSILON),tf.int32))
				if self.args.trainStrandedActivities>0:
					self.paramPenaltyL1Tensor = self.paramPenaltyL1Tensor + tf.reduce_sum(tf.abs(tf.slice(self.positionalActivityBiasRC,[0,0],[1,self.args.numMotifs]))) 
					self.paramNumActivityTensor = self.paramNumActivityTensor +tf.reduce_sum(tf.cast(tf.greater(tf.abs(tf.slice(self.positionalActivityBiasRC,[0,0],[1,self.args.numMotifs])),EPSILON),tf.int32))
			self.myLoss = tf.add(self.myLoss, tf.multiply(self.paramPenaltyL1Tensor,self.args.L1));
		
		if self.args.L2Pos is not None and self.args.trainPositionalActivities>0:
			self.args.L2Pos = float(self.args.L2Pos);
			paramPenaltyL2PosTensor = tf.reduce_sum(tf.abs(tf.subtract(tf.slice(self.positionalActivityBias,[0,0],[self.args.seqLen-1,self.args.numMotifs]),tf.slice(self.positionalActivityBias,[1,0],[self.args.seqLen-1,self.args.numMotifs]))))
			if self.args.trainStrandedActivities>0:
				paramPenaltyL2PosTensor = paramPenaltyL2PosTensor + tf.nn.l2_loss(tf.subtract(tf.slice(self.positionalActivityBiasRC,[0,0],[self.args.seqLen-1,self.args.numMotifs]),tf.slice(self.positionalActivityBiasRC,[1,0],[self.args.seqLen-1,self.args.numMotifs])))
			self.myLoss = tf.add(self.myLoss, tf.multiply(paramPenaltyL2PosTensor, self.args.L2Pos));
			
		
		if self.args.interactions>0:
			if self.args.verbose>0: sys.stderr.write("Using L1 regularization of interaction coeficients with lambda=%s\n"%(self.args.L1int))
		 	paramPenaltyL1intTensor = tf.reduce_sum(tf.abs(self.kdBonusInteractionsFilterFF)) + tf.reduce_sum(tf.abs(self.kdBonusInteractionsFilterFR)) +tf.reduce_sum(tf.abs(self.kdBonusInteractionsFilterRF)) +tf.reduce_sum(tf.abs(self.kdBonusInteractionsFilterRR));
		 	paramNumIntTensor = tf.reduce_sum(tf.cast(tf.greater(tf.abs(self.kdBonusInteractionsFilterFF),EPSILON),tf.int32))+tf.reduce_sum(tf.cast(tf.greater(tf.abs(self.kdBonusInteractionsFilterFR),EPSILON),tf.int32))+tf.reduce_sum(tf.cast(tf.greater(tf.abs(self.kdBonusInteractionsFilterRF),EPSILON),tf.int32))+tf.reduce_sum(tf.cast(tf.greater(tf.abs(self.kdBonusInteractionsFilterRR),EPSILON),tf.int32))
			self.myLoss = tf.add(self.myLoss, tf.multiply(paramPenaltyL1intTensor,self.args.L1int));
		
		self.global_step = tf.Variable(0, trainable=False, name="global_step")
		learningRate = self.args.learningRate
		if self.args.cycleLearningRate is not None:
			if self.args.verbose>0: sys.stderr.write("Cycling learning rate low,high,period %s\n"%self.args.cycleLearningRate)
			clrParams = self.args.cycleLearningRate.split(",")
			learningRate = TFHELP.cycle_learning_rate(float(clrParams[0]), float(clrParams[1]), self.global_step, int(clrParams[2]), name="learningRate")
		else:
			learningRate =tf.train.exponential_decay(learningRate, self.global_step, 175, self.args.learningRateED, staircase=True, name="learningRate")
		if self.args.useMomentumOptimizer > 0 or self.args.useRMSProp>0:
			if self.args.useMomentumOptimizer > 0 and self.args.useRMSProp>0:
				raise Exception("Cannot use both MomentumOptimizer and RMSProp! pick one!");
			if isinstance(self.args.momentum, basestring):
				momentum = tf.constant(float(self.args.momentum),name="momentum");
			else: 
				momentum=self.args.momentum
			if self.args.useMomentumOptimizer > 0:
				if self.args.verbose>0: sys.stderr.write("Using MomentumOptimizer instead of Adam\n")
				opt = tf.train.MomentumOptimizer(self.args.learningRate, momentum, use_nesterov = self.args.useNesterov>0);
		
			elif self.args.useRMSProp>0:
				if self.args.verbose>0: sys.stderr.write("Using RMSProp instead of Adam\n")
				opt = tf.train.RMSPropOptimizer(self.args.learningRate, decay = self.args.rmsDecay, momentum=momentum, epsilon=self.args.rmsEpsilon);
		else:
			if self.args.verbose>0: sys.stderr.write("Using AdamOptimizer\n")
			opt = tf.train.AdamOptimizer(self.args.learningRate, epsilon=self.args.Aepsilon, beta1=self.args.Abeta1, beta2=self.args.Abeta2);
		
		self.train_step = opt.minimize(self.myLoss, global_step=self.global_step);
		
		#raise Exception("Reached bad state=%d for '%s.%d' '%s' at line '%s'" %(state,mid,ver,tfid,line));
		if self.args.trace > 2:
			self.sess = tf.Session(config=tf.ConfigProto(intra_op_parallelism_threads=self.args.threads-1, log_device_placement=True));
		else:
			self.sess = tf.Session(config=tf.ConfigProto(intra_op_parallelism_threads=self.args.threads-1));
		
		
		if self.args.tensorboard is not None:
			#learning rate, loss
			tf.summary.scalar('learning rate', learningRate)
			tf.summary.scalar('MSE', self.mseTF)
			tf.summary.scalar('loss', self.myLoss)
			if self.args.interactions>0:
				tf.summary.scalar('TFTF L1 interaction penalty', paramPenaltyL1intTensor)
				tf.summary.scalar('number of non-zero TFTF interactions', paramNumIntTensor)
				tf.summary.histogram('TFTF interactions FF', self.kdBonusInteractionsFilterFF)
				tf.summary.histogram('TFTF interactions FR', self.kdBonusInteractionsFilterFR)
				tf.summary.histogram('TFTF interactions RF', self.kdBonusInteractionsFilterRF)
				tf.summary.histogram('TFTF interactions RR', self.kdBonusInteractionsFilterRR)
			if self.args.trainPositionalActivities>0:
				tf.summary.histogram('TF positional activity bias', self.positionalActivityBias)
				if self.args.trainStrandedActivities>0:
					tf.summary.histogram('TF positional self.activities (RC)', self.positionalActivityBiasRC)
				if self.args.L2Pos is not None:
					tf.summary.scalar('positional L2 penalty', paramPenaltyL2PosTensor)
			if self.args.L1 is not None:
				tf.summary.scalar('general L1 penalty', self.paramPenaltyL1Tensor)
				tf.summary.scalar('number of non-zero L1 parameters', self.paramNumActivityTensor)
			# learned parameters
			tf.summary.histogram('TF self.activities', self.activities)
			if self.args.potentiation>0:
				tf.summary.histogram('TF potentiations', self.potentiation)
			if self.args.trainStrandedActivities>0:
				tf.summary.histogram('TF activity strand difference', self.activityDiffs)
			self.mergedSummaries = tf.summary.merge_all()
			self.summaryWriter = tf.summary.FileWriter(self.args.tensorboard, self.sess.graph)
		
		init_op = tf.global_variables_initializer()
		self.saver = tf.train.Saver()
		self.sess.run(init_op)
		if self.args.restore is not None:
			if self.args.verbose>1:
				sys.stderr.write("Loading initial parameter settings: %s\n"%(self.args.restore))
			reader = tf.train.NewCheckpointReader(self.args.restore);
			restoringThese = {v.name[0:-2]: v for v in tf.global_variables()};
			#print(list(restoringThese.keys()));
			if "concentrations" in reader.get_variable_to_shape_map().keys()  and  self.args.VARIABLE:
				raise Exception("BACKWARD COMPATIBILITY ERROR: loaded graph has variable named concentrations - re-run withOUT -VARIABLE option\n");
			if "Variable" in reader.get_variable_to_shape_map().keys() and not self.args.VARIABLE:
				raise Exception("BACKWARD COMPATIBILITY ERROR: re-run with -VARIABLE option\n");
				#restoringThese["Variable"] = restoringThese["concentrations"]
				#del restoringThese["concentrations"]
			restoringThese = {k: restoringThese[k] for k in list(restoringThese.keys()) if k in reader.get_variable_to_shape_map().keys()}
			import re;
			restoringThese = {k: restoringThese[k] for k in list(restoringThese.keys()) if re.search("global_step", k) is None}
			if self.args.clearAdamVars>0:
				restoringThese = {k: restoringThese[k] for k in list(restoringThese.keys()) if re.search("/Adam", k) is None}
				restoringThese = {k: restoringThese[k] for k in list(restoringThese.keys()) if re.search("beta[12]_power", k) is None}
			if self.args.verbose>1:
				sys.stderr.write("Loading these variables: %s\n"%(", ".join([k for k in list(restoringThese.keys())])))
			restoringThese = [restoringThese[k] for k in list(restoringThese.keys())]
			restorer = tf.train.Saver(restoringThese);
			restorer.restore(self.sess, self.args.restore)
		
		#make sure there are no unnamed variables
		import re;
		#if [k.name for k in tf.trainable_variables() if re.search("Variable", k.name) is not None].len() > 0:
		if len([k.name for k in tf.trainable_variables() if re.search("Variable", k.name) is not None]) > 0 and not self.args.VARIABLE:
			print([k.name for k in tf.trainable_variables()])
			raise(Exception("Error: one or more variables with default names: %s"%", ".join([k.name for k in tf.trainable_variables() if re.search("Variable", k.name) is not None])));
		
	def saveParams(self, sess):
		if self.args.activityClasses>1:
			if (self.args.outFPre is None):
				outFile= sys.stdout;
			else:
				outFile = MYUTILS.smartGZOpen(self.args.outFPre+".activityClasses.gz",'w');
			classWeights = self.activityClassWeights.eval(session=self.sess); #[4,self.args.motifLen,1,self.args.numMotifs]);[numMotifs, activityClasses]
			sigmoidBiases = self.activityClassActivationBiases.eval(session=self.sess); #[activityClasses]
			sigmoidWeights = self.activityClassExpressionBiases.eval(session=self.sess); #[activityClasses, 1]
			outFile.write("i");
			for i in range(0, self.args.activityClasses):
				outFile.write("\tAC%i"%i);
			outFile.write("\n");
			outFile.write("-2")
			for i in range(0, self.args.activityClasses):
				outFile.write("\t%g"%sigmoidBiases[i]);
			outFile.write("\n");
			outFile.write("-1")
			for i in range(0, self.args.activityClasses):
				outFile.write("\t%g"%sigmoidWeights[i,0]);
			outFile.write("\n");
			for i in range(0, self.args.numMotifs):
				outFile.write("%i"%i);
				for j in range(0, self.args.activityClasses):
					outFile.write("\t%g"%classWeights[i,j]);
				outFile.write("\n");
			outFile.close();
			
		if self.args.noTrainMotifs==0:
			if (self.args.outFPre is None):
				outFile= sys.stdout;
			else:
				outFile = MYUTILS.smartGZOpen(self.args.outFPre+".pkdms",'w');
			#print out weight matrices
			pkdms = self.motifsTensor2.eval(session=self.sess); #[4,self.args.motifLen,1,self.args.numMotifs]);
			for i in range(0, self.args.numMotifs):
				outFile.write("#Motif %i\n"%(i));
				for j in range(0,4):
					outFile.write("#%s"%BASES[j]);
					for k in range(0,self.args.motifLen):
						outFile.write("\t%g"%(pkdms[j,k,0,i]));
					outFile.write("\n");
			outFile.write("\n");
			outFile.close();
		if self.args.interactions>0:
			if (self.args.outFPre is None):
				outFile= sys.stdout;
			else:
				outFile = MYUTILS.smartGZOpen(self.args.outFPre+".interactions.gz",'w');
			outFile.write("TF1\tTF2\tstrands\tdist\tcoef\n")
			strands = ["++","+-","-+","--"];
			interactionTensors = [self.kdBonusInteractionsFilterFF, self.kdBonusInteractionsFilterFR, self.kdBonusInteractionsFilterRF, self.kdBonusInteractionsFilterRR];
			for s in range(0,4):
				strand =strands[s]
				tfInteractions = interactionTensors[s].eval(session=self.sess); #[numMot, distance, 1, numMot]
				for i in range(0, self.args.numMotifs): # primary motif
					for j in range(0, self.args.numMotifs): # interacting motif
						for d in range(0, self.args.interactions):
							outFile.write("%i\t%i\t%s\t%i\t%g\n"%(i,j,strand,d-(self.args.interactions//2),tfInteractions[j,d ,0,i]));
			outFile.close();
		if self.args.trainPositionalActivities>0:
			if (self.args.outFPre is None):
				outFile= sys.stdout;
			else:
				outFile = MYUTILS.smartGZOpen(self.args.outFPre+".positional.gz",'w');
			outFile.write("TF\tposition\tpositionalActivityBias");
			positionalActivityBiasVals = self.positionalActivityBias.eval(session=self.sess).reshape((self.args.seqLen,self.args.numMotifs));
			if self.args.trainStrandedActivities>0:
				outFile.write("\tpositionalActivityBiasRC");
				positionalActivityBiasRCVals = self.positionalActivityBiasRC.eval(session=self.sess).reshape((self.args.seqLen,self.args.numMotifs));
			outFile.write("\n");
			for j in range(0,self.args.numMotifs):
				for i in range(0,self.args.seqLen):
					outFile.write("%i\t%i\t%g"%(j,i,positionalActivityBiasVals[i,j]));
					if self.args.trainStrandedActivities>0:
						outFile.write("\t%g"%positionalActivityBiasRCVals[i,j]);
					outFile.write("\n");
			outFile.close();
				
		if (self.args.outFPre is None):
			outFile= sys.stdout;
		else:
			outFile = MYUTILS.smartGZOpen(self.args.outFPre+".params",'w');
		#print params
		concs = self.logConcs.eval(session=self.sess).reshape((self.args.numMotifs));
		activityVals = self.activities.eval(session=self.sess).reshape((self.args.numMotifs));
		outFile.write("i\tlogConc\tactivity");
		if self.args.potentiation>0:
			outFile.write("\tpotentiations");
			potentiationVals = self.potentiation.eval(session=self.sess).reshape((self.args.numMotifs));
		if self.args.bindingLimits>0:
			outFile.write("\tbindingLimits");
			bindingLimitVals = self.bindingLimits.eval(session=self.sess).reshape((self.args.numMotifs));
		if self.args.trainStrandedActivities>0:
			outFile.write("\tactivityDiffs");
			activityDiffVals = self.activityDiffs.eval(session=self.sess).reshape((self.args.numMotifs));
		outFile.write("\n");
		if self.args.accIsAct>0:
			outFile.write("-1\t%g\t%g"%(self.accActivity.eval(session=self.sess),self.constant.eval(session=self.sess)));#intercepts
		else:
			outFile.write("-1\tNA\t%g"%self.constant.eval(session=self.sess));#intercepts
		if self.args.potentiation>0:
			outFile.write("\t%g"%(self.constantPot.eval(session=self.sess)));
		if self.args.bindingLimits>0:
			outFile.write("\tNA");
		if self.args.trainStrandedActivities>0:
			outFile.write("\tNA");
		outFile.write("\n");
		for i in range(0,self.args.numMotifs):
			outFile.write("%i\t%g\t%g"%(i, concs[i], activityVals[i]));#intercepts
			if self.args.potentiation>0:
				outFile.write("\t%g"%(potentiationVals[i]));
			if self.args.bindingLimits>0:
				outFile.write("\t%g"%(bindingLimitVals[i]));
			if self.args.trainStrandedActivities>0:
				outFile.write("\t%g"%(activityDiffVals[i]));
			outFile.write("\n");
		outFile.close();
	
