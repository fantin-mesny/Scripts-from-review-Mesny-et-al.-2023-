#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

## 0. Unzip all compressed files in the repository
gunzip $SCRIPT_DIR/Genomic_data/*/*/*.gz
gunzip $SCRIPT_DIR/Functional_annotations/*/*.gz
tar -xf $SCRIPT_DIR/Orthology_prediction/*.tar.gz

## 1. Orthology prediction
# Main output files from these orthology predictions are provided in folder "Orthology_prediction"
orthofinder -f $SCRIPT_DIR/Genomic_data/Sordariomycetes/protein -S blast
orthofinder -f $SCRIPT_DIR/Genomic_data/Betaproteobacteria/protein -S blast
# For simplicity and future scripts, let's rename the Orthofinder output folder 'Results_{DATE}' into 'Results'
mv $SCRIPT_DIR/Genomic_data/Betaproteobacteria/protein/OrthoFinder/$(ls $SCRIPT_DIR/Genomic_data/Betaproteobacteria/protein/OrthoFinder/) $SCRIPT_DIR/Genomic_data/Betaproteobacteria/protein/OrthoFinder/Results
mv $SCRIPT_DIR/Genomic_data/Sordariomycetes/protein/OrthoFinder/$(ls $SCRIPT_DIR/Genomic_data/Sordariomycetes/protein/OrthoFinder/) $SCRIPT_DIR/Genomic_data/Sordariomycetes/protein/OrthoFinder/Results


## 2. Align sequences within orthogroups
# This step uses FAMSA to align sequences clustered into orthogroups by OrthoFinder.
# It uses the OrthoFinder output folder "Orthogroup_Sequences", which was not uploaded on GitHub. However, the protein sequences and the file 'Orthogroups.tsv' can be used to retrieve individual sequences in each orthogroup.
for class in "Sordariomycetes" "Betaproteobacteria"
do
    cd $SCRIPT_DIR/Genomic_data/$class/protein/OrthoFinder/Results/Orthogroup_Sequences #Output data from OrthoFinder, not uploaded on GitHub
    ls * | while read file; do famsa-1.6.1-linux $file ../Orthogroup_Sequences_Aligned/$file.aln ;done
done


## 3. Organize CDS sequences by orthogroups
# This step is necessary to then be able to obtain codon alignments from protein alignments in step 4.
cat $SCRIPT_DIR/Genomic_data/Betaproteobacteria/cds/*.ffn > $SCRIPT_DIR/Genomic_data/Betaproteobacteria/cds/all.ffn
python $SCRIPT_DIR/Utility_scripts/getTranscriptSequencesForEachOG_Betaproteobacteria.py
rm $SCRIPT_DIR/Genomic_data/Betaproteobacteria/cds/all.ffn

cat $SCRIPT_DIR/Genomic_data/Sordariomycetes/cds/*.ffn > $SCRIPT_DIR/Genomic_data/Sordariomycetes/cds/all.ffn
python $SCRIPT_DIR/Utility_scripts/getTranscriptSequencesForEachOG_Sordariomycetes.py
rm $SCRIPT_DIR/Genomic_data/Sordariomycetes/cds/all.ffn


## 3. Obtain codon alignments using PAL2NAL
# In this step, we use PAL2NAL to obtain codon alignments from protein alignments and CDS sequences
for class in "Sordariomycetes" "Betaproteobacteria"
do
    cd $SCRIPT_DIR/Genomic_data/$class/protein/OrthoFinder/Results/Orthogroup_Sequences_Aligned
    for f in *.aln
    do
	   og=$(python -c "print('$f'.split('.')[0])");\
	   pal2nal.pl $f $SCRIPT_DIR/Genomic_data/$class/protein/OrthoFinder/Results/Orthogroup_Sequences_Transcripts/$og.ffn -output fasta > $f.pal2nal.fasta;\
	   python $SCRIPT_DIR/Utility_scripts/renameInPal2Nal_$class.py -i $f.pal2nal.fasta;\
    done;
done


## 4. FUBAR
# Run the FUBAR algorithm on all orthogroups
# To do so, we use the Gene Trees generated for each orthogroup by OrthoFinder 
for class in "Sordariomycetes" "Betaproteobacteria"
do
    cd $SCRIPT_DIR/Genomic_data/$class/protein/OrthoFinder/Results/Orthogroup_Sequences_Aligned
    mkdir hyphy
    for f in *.aln
    do
	   og=$(python -c "print('$f'.split('.')[0])");\
       hyphy fubar $f.2.pal2nal.fasta $SCRIPT_DIR/Genomic_data/$class/protein/OrthoFinder/Results/Gene_Trees/$og\_tree.txt > hyphy/$og.hyphy.log;
    done;
    grep '## FUBAR inferred ' hyphy/*.hyphy.log > $SCRIPT_DIR/Fubar_outputs/allHyphyLogs_$class.txt
done


## 5. Functional analysis
# This part of the pipeline is used to generate the data tables provided as Supplementary Material of the review
# These tables were also used to generate Figure 1
python $SCRIPT_DIR/Utility_scripts/functionalAnalysis_Betaproteobacteria.py
python $SCRIPT_DIR/Utility_scripts/functionalAnalysis_Sordariomycetes.py