import pandas as pd
import os
import re
scriptDir=os.path.dirname(os.path.realpath(__file__))

allHyphyLogs=scriptDir+'/../Fubar_outputs/allHyphyLogs_Betaproteobacteria.txt'
ogTable=scriptDir+'/../Genomic_data/Betaproteobacteria/OrthoFinder/Results/Orthogroups/Orthogroups.tsv'
eggnog=scriptDir+'/../Genomic_data/Betaproteobacteria/eggnogAnnotations_genes.tsv'
dbCan=scriptDir+'/../Functional_annotations/Betaproteobacteria/cazymesAnnotation_genes.tsv'
signalP=scriptDir+'/../Functional_annotations/Betaproteobacteria/signalP_genes.tsv'
positiveSelectionThreshold=1 #superior to, meaning >=2

def parseAllHyphyLogs(allHyphyLogs):
    df=pd.read_csv(allHyphyLogs,sep=' ',header=None)
    df=df[[0,3]]
    df[0]=df[0].str.split('.hyphy').str[0]
    df[3]=df[3].str.replace('no','0')
    df[3]=df[3].astype(int)
    df=df.rename(columns={0:'OG',3:'Number of sites under positive selection'})
    return df.set_index('OG')

def parseAnnotations(eggnog,signalP,dbCan):
    allEgg=pd.read_csv(eggnog,sep='\t', header=None,comment='#')
    allEgg=allEgg[[0,6,20]].set_index(0).rename(columns={6:'GOs',20:'COGs'})
    
    signalP=pd.read_csv(signalP,sep='\t',comment="#",header=None).set_index(0)
    signalP['Effectors']=signalP[1]!='OTHER'
    signalP=signalP[['Effectors']]

    allCaz=pd.read_csv(dbCan,sep='\t').set_index('Gene ID')
    allCaz['CAZannot']=allCaz['HMMER']+' '+allCaz['eCAMI']+' '+allCaz['DIAMOND']
    allCaz=allCaz[['CAZannot','#ofTools']].rename(columns={'#ofTools':'CAZymes'})
    
    annots=signalP.merge(allEgg,how='left',left_index=True,right_index=True)
    annots=annots.merge(allCaz,how='left',left_index=True,right_index=True)
    annots['GOs']=annots['GOs'].fillna('')
    annots['COGs']=annots['COGs'].fillna('')
    annots['GOs']=annots['GOs'].str.replace(',',';')
    annots['CAZymes']=annots['CAZymes']>0
    return annots



df=parseAllHyphyLogs(allHyphyLogs)
annots=parseAnnotations(eggnog,signalP,dbCan)


og=pd.read_csv(ogTable,sep='\t').set_index('Orthogroup')
og=og.fillna('')




OGtoGO={}
OGtoCaz={}
OGtoEff={}
OGtoCOG={}
OGtoCazFam={}
OGtoProtNumber={}
for ind in df.index:
    OGtoGO[ind]=[]
    OGtoCOG[ind]=''
    OGtoCaz[ind]=0
    OGtoEff[ind]=0
    OGtoProtNumber[ind]=0
    OGtoCazFam[ind]=[]
    for col in og.columns:
        prots=og.loc[ind, col].split(', ')
        for prot in prots:
            if prot!='':
                OGtoProtNumber[ind]+=1
                OGtoGO[ind]+=annots.loc[prot,'GOs'].split(';')
                OGtoCaz[ind]+=annots.loc[prot,'CAZymes'].astype(int)
                OGtoEff[ind]+=annots.loc[prot,'Effectors'].astype(int)
                OGtoCOG[ind]+=annots.loc[prot,'COGs']
                if str(annots.loc[prot,'CAZannot'])!='nan':
                    OGtoCazFam[ind].append(';'.join(list(set([fam for fam in re.sub(r'\([^)]*\)', '', annots.loc[prot,'CAZannot']).replace('-','').replace(' ','+').split('+') if fam not in ['+','-',' ']]))))


OGtoGO_filtered={o:';'.join([prot for prot in set(OGtoGO[o]) if prot!='']) for o in OGtoGO}

OGtoCazFam2={OG:list(set([a for a in ';'.join(OGtoCazFam[OG]).split(';') if a!=''])) for OG in OGtoCazFam}
OGtoCazFam_f={OG:';'.join(list(set(OGtoCazFam2[OG]))) for OG in OGtoCazFam2}
OGtoCOG={o:';'.join(set(list(OGtoCOG[o]))) for o in OGtoCOG}

df['Number of proteins in OG']=df.index.map(OGtoProtNumber)
df['OG under positive selection']=(df['Number of sites under positive selection']>positiveSelectionThreshold).astype(str).str.replace('True','Yes').str.replace('False','No')
df['Secreted']=df.index.map(OGtoEff)
df['Secreted proteins in OG']=(df['Secreted']>0).astype(str).str.replace('True','Yes').str.replace('False','No')

df['CAZyme']=df.index.map(OGtoCaz)
df['CAZymes in OG']=(df['CAZyme']>0).astype(str).str.replace('True','Yes').str.replace('False','No')
df['CAZyme annotation']=df.index.map(OGtoCazFam_f)
df['COG annotation']=df.index.map(OGtoCOG)
df['GO annotation']=df.index.map(OGtoGO_filtered)

df=df[['Number of proteins in OG','OG under positive selection','Number of sites under positive selection','Secreted proteins in OG','CAZyme annotation','COG annotation','GO annotation']]
df.to_csv('table_betaproteobacteria.csv')

