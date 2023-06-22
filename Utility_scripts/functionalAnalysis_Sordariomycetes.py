import pandas as pd
import os
import re
scriptDir=os.path.dirname(os.path.realpath(__file__))

allHyphyLogs=scriptDir+'/../Fubar_outputs/allHyphyLogs_Sordariomycetes.txt'
ogTable=scriptDir+'/../Genomic_data/Sordariomycetes/OrthoFinder/Results/Orthogroups/Orthogroups.tsv'
dbCan=scriptDir+'/../Functional_annotations/Sordariomycetes/cazymesAnnotation_genes.tsv'
signalP=scriptDir+'/../Functional_annotations/Sordariomycetes/signalP_genes.tsv'

goAnn='/netscratch/dep_psl/grp_hacquard/Fantin/REVIEW/Sordario/GOAnnotation.csv'
positiveSelectionThreshold=1 #superior to

def parseAllHyphyLogs(allHyphyLogs):
    df=pd.read_csv(allHyphyLogs,sep=' ',header=None)
    df=df[[0,3]]
    df[0]=df[0].str.split('.hyphy').str[0]
    df[3]=df[3].str.replace('no','0')
    df[3]=df[3].astype(int)
    df=df.rename(columns={0:'OG',3:'Number of sites under positive selection'})
    return df.set_index('OG')

def parseAnnotations(goAnn,signalP,dbCan):
    goAnn=pd.read_csv(goAnn).rename(columns={'goAcc':'GOs','#proteinId':'shorterId'})
    
    signalP=pd.read_csv(signalP,sep='\t',skiprows=1).set_index('# ID')
    signalP['hasSignal']=signalP['Prediction']!='OTHER'
    signalP=signalP[['hasSignal']]
    
    allCaz=pd.read_csv(dbCan,sep='\t').set_index('Gene ID')
    allCaz['CAZannot']=allCaz['HMMER']+' '+allCaz['eCAMI']+' '+allCaz['DIAMOND']
    allCaz=allCaz[['CAZannot','#ofTools']].rename(columns={'#ofTools':'CAZymes'})
    
    annots=signalP
    annots=annots.merge(allCaz,how='left',left_index=True,right_index=True)
    annots['shorterId']=annots.index.str.split('|').str[0:-1]
    annots['shorterId']=annots['shorterId'].apply('|'.join)
    annots=annots.reset_index().merge(goAnn,on='shorterId',how='left').set_index('# ID').drop(columns='shorterId')
    annots['GOs']=annots['GOs'].fillna('')
    annots['CAZymes']=annots['CAZymes']>0
    return annots



df=parseAllHyphyLogs(allHyphyLogs)

annots=parseAnnotations(goAnn,signalP,dbCan)

og=pd.read_csv(ogTable,sep='\t').set_index('Orthogroup')
og=og.fillna('')



OGtoGO={}
OGtoCaz={}
OGtoEff={}
OGtoCazFam={}
OGtoProtNumber={}
for ind in df.index:
    OGtoGO[ind]=[]
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
                OGtoEff[ind]+=annots.loc[prot,'hasSignal'].astype(int)


                if str(annots.loc[prot,'CAZannot'])!='nan':
                    OGtoCazFam[ind].append(';'.join(list(set([fam for fam in re.sub(r'\([^)]*\)', '', annots.loc[prot,'CAZannot']).replace('-','').replace(' ','+').split('+') if fam not in ['+','-',' ']]))))
OGtoGO_filtered={o:';'.join([prot for prot in set(OGtoGO[o]) if prot!='']) for o in OGtoGO}
OGtoCazFam2={OG:list(set([a for a in ';'.join(OGtoCazFam[OG]).split(';') if a!=''])) for OG in OGtoCazFam}
OGtoCazFam_f={OG:';'.join(list(set(OGtoCazFam2[OG]))) for OG in OGtoCazFam2}

df['Number of proteins in OG']=df.index.map(OGtoProtNumber)
df['OG under positive selection']=(df['Number of sites under positive selection']>positiveSelectionThreshold).astype(str).str.replace('True','Yes').str.replace('False','No')
df['Secreted']=df.index.map(OGtoEff)
df['Secreted proteins in OG']=(df['Secreted']>0).astype(str).str.replace('True','Yes').str.replace('False','No')

df['CAZyme']=df.index.map(OGtoCaz)
df['CAZymes in OG']=(df['CAZyme']>0).astype(str).str.replace('True','Yes').str.replace('False','No')
df['CAZyme annotation']=df.index.map(OGtoCazFam_f)

df['GO Annotation']=df.index.map(OGtoGO_filtered)

df=df[['Number of proteins in OG','OG under positive selection','Number of sites under positive selection','Secreted proteins in OG','CAZyme annotation','GO Annotation']]
df.to_csv('table_sordariomycetes.csv')

