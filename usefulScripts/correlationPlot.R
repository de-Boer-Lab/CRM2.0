#!/usr/bin/env Rscript
allArgs = commandArgs(TRUE);

if (length(allArgs)<1){
  print("Usage: correlationPlot.R <inFile> <outFPre> <field1> <field2> \n");
  quit("no");
}

inFile= allArgs[1];
outFPre= allArgs[2];
f1= allArgs[3];
f2= allArgs[4];

#library(ggplot2, lib.loc=.libPaths()[-1]) # RHEL6
library(ggplot2)
library(ggrastr)
#library(reshape)

dataTable = read.csv(inFile, header=TRUE, sep="\t", row.names=NULL);

r = cor(dataTable[[f1]], dataTable[[f2]], use="pairwise.complete.obs")

p1 = ggplot(dataTable, aes_string(x=f1, y=f2)) + geom_point_rast(alpha=0.01)+ggtitle(sprintf("r2 = %g",r^2))+theme_classic();
p2 = ggplot(dataTable, aes_string(x=f1, y=f2)) + geom_point_rast(alpha=0.1)+ggtitle(sprintf("r2 = %g",r^2))+theme_classic();
p3 = ggplot(dataTable, aes_string(x=f1, y=f2)) + geom_point_rast(alpha=1)+ggtitle(sprintf("r2 = %g",r^2))+theme_classic();
if (nrow(dataTable)<=1E6){
	p1 = p1 + geom_smooth(); 
	p2 = p2 + geom_smooth(); 
	p3 = p3 + geom_smooth(); 
}

ggsave(sprintf('%s.pdf',outFPre),width=5, height=5, plot=p1)
#ggsave(sprintf('%s.png',outFPre),width=5, height=5, plot=p1)
ggsave(sprintf('%s_a0.1.pdf',outFPre),width=5, height=5, plot=p2)
#ggsave(sprintf('%s_a0.1.png',outFPre),width=5, height=5, plot=p2)
ggsave(sprintf('%s_a1.pdf',outFPre),width=5, height=5, plot=p3)
#ggsave(sprintf('%s_a1.png',outFPre),width=5, height=5, plot=p3)

write.table(data.frame(id=outFPre, r=r), file=sprintf('%s.r',outFPre), col.names=F, row.names=F, quote=F, sep="\t")
