#Align genomic sequences using stretcher

#Install EMBOSS

conda activate mamba
conda install mamba
mamba install emboss

#Run stretcher
#This software is used to align long genomic sequences. The online version has size restrictions for the sequences

stretcher -asequence sequence_1.fasta -bsequence sequence_2.fasta -outfile alignment.fasta -aformat fasta

#Generate the dotplot

#Import the alignment into Geneious software to generate the dotplot