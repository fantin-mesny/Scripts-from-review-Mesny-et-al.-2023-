from Bio import SeqIO
import pandas as pd
import os
scriptDir=os.path.dirname(os.path.realpath(__file__))

sorda=['Chafu1','Chagl1','Chame1','Cylol1','Daces1','Dacma1','Fusco1','Fuseq1','Fusoxy1','Fusoxys1','Fusre1','Fusso1','Fustr1','Fustri1','Fusven1','Ilyeu1','Mictri1','Neora1','Plecuc1','Plecucu1','Sorhu1','Stael1','Truan1','Verdah1']

## pId to tId
gffFiles={f.split('_')[0]:f for f in os.listdir(scriptDir+'/../Genomic_data/Sordariomycetes/gff') if f.split('_')[0] in sorda}
t2p={s:{} for s in sorda}
for f in gffFiles:
    print(f, gffFiles[f])
    gff=pd.read_csv(scriptDir+'/../Genomic_data/Sordariomycetes/gff/'+gffFiles[f],sep='\t',header=None,comment="#")
    gff=gff[gff[2]=='mRNA'][[8]].reset_index(drop=True)
    for ind in gff.index:
#        p2t[f][gff.loc[ind,8].split('proteinId=')[1].split(';')[0]]=gff.loc[ind,8].split('transcriptId=')[1]
        t2p[f][gff.loc[ind,8].split('transcriptId=')[1]]=gff.loc[ind,8].split('proteinId=')[1].split(';')[0]

with open(scriptDir+'/../Genomic_data/Sordariomycetes/cds/all.ffn','r') as inp:
	allt=SeqIO.to_dict(SeqIO.parse(inp, 'fasta'))
allt_renamed={t.replace('|'+t.split('|')[2]+'|','|'+t2p[t.split('|')[1]][t.split('|')[2]]+'|'):allt[t] for t in allt}


Dir=scriptDir+'/../Genomic_data/Sordariomycetes/protein/OrthoFinder/Results/Orthogroup_Sequences/'
for f in os.listdir(Dir):
	tseqs=[]
	with open(Dir+f, 'r') as inp:
		seqs=SeqIO.to_dict(SeqIO.parse(inp,'fasta'))
	for seq in seqs:
		tseqs.append(allt_renamed[seq])
		tseqs[-1].id=seq
		tseqs[-1].name=seq
		tseqs[-1].description=seq
		tseqs[-1].seq=tseqs[-1].seq.replace('*','')
	with open(scriptDir+'/../Genomic_data/Sordariomycetes/protein/OrthoFinder/Results/Orthogroup_Sequences_Transcripts/'+f.replace('fa','ffn'),'w+') as outp:
		SeqIO.write(tseqs,outp,'fasta')