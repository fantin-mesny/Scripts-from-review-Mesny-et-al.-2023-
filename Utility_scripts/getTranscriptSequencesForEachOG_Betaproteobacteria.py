from Bio import SeqIO
import os
scriptDir=os.path.dirname(os.path.realpath(__file__))

with open(scriptDir+'/../Genomic_data/Betaproteobacteria/cds/all.ffn','r') as inp:
	allt=SeqIO.to_dict(SeqIO.parse(inp, 'fasta'))

Dir=scriptDir+'/../Genomic_data/Betaproteobacteria/OrthoFinder/Results/Orthogroup_Sequences/'
for f in os.listdir(Dir):
	tseqs=[]
	with open(Dir+f, 'r') as inp:
		seqs=SeqIO.to_dict(SeqIO.parse(inp,'fasta'))
	for seq in seqs:
		tseqs.append(allt[seq])
	with open(scriptDir+'/../Genomic_data/Betaproteobacteria/OrthoFinder/Results/Orthogroup_Sequences_Transcripts/'+f.replace('fa','ffn'),'w+') as outp:
		SeqIO.write(tseqs,outp,'fasta')
