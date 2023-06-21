from Bio import SeqIO
import os

with open('Genomic_data/Betaproteobacteria/cds/all.ffn','r') as inp:
	allt=SeqIO.to_dict(SeqIO.parse(inp, 'fasta'))

Dir='Genomic_data/Betaproteobacteria/OrthoFinder/Results/Orthogroup_Sequences/'
for f in os.listdir(Dir):
	tseqs=[]
	with open(Dir+f, 'r') as inp:
		seqs=SeqIO.to_dict(SeqIO.parse(inp,'fasta'))
	for seq in seqs:
		tseqs.append(allt[seq])
	with open('Genomic_data/Betaproteobacteria/OrthoFinder/Results/Orthogroup_Sequences_Transcripts/'+f.replace('fa','ffn'),'w+') as outp:
		SeqIO.write(tseqs,outp,'fasta')