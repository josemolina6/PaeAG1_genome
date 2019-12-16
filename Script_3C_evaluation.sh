# Pipeline for assembly assessment by contiguity, completeness and correctness (3C) criterion 
echo " "
echo "PIPELINE FOR ASSEMBLY ASSESSMENT BY CONTIGUITY, COMPLETENESS AND CORRECTNESS (3C) CRITERION"
echo "----------------------------------------------"
echo "Implemented by Jose Arturo Molina Mora"
echo "University of Costa Rica, San José, Costa Rica"
echo "----------------------------------------------"
echo " "

echo "++++++++++++++++++++++++++++++++++++++++++++++"
echo " 3C criterion analysis for assemblies"
echo "++++++++++++++++++++++++++++++++++++++++++++++"
echo "Required files: Assembly.fasta, Reference.fasta, shortreads1.fastq, shortreads2.fastq, longreads.fastq, known_sequences.fasta (Sanger sequencing for example)"

echo "******* Contiguity analysis *********"
echo " "

#"USAGE: Create a folder with all the assemblies to evaluate. Call the script indicating the number of assemblies (fasta files), the reference genome and then the names of the assemblies files. 

mkdir contiguity

assembly-stats -t *.fasta > ./Contiguity/Genenal_metrics.txt

quast *.fasta -o ./contiguity/Quast_results
quast *.fasta -o ./contiguity/Quast_results_with_reference -r Reference.fasta
# Note: QUAST will also be run for the reference. 


echo " "
echo "******* Completeness analysis *********"
echo " "

mkdir completeness
cd completeness

# Short and long reads must be also included, and here we call them shortreads1.fastq, shortreads2.fastq and longreads1.fastq. In addition, annotation files (GFF, GBK and other are also required).

# Single copy ortholog gene sets
echo "Access gVolante plataform (https://gvolante.riken.jp) to run Single copy ortholog gene sets"
# Single copy ortholog gene sets:  were searched (expected: 100%) in the assemblies using the BUSCO tool (Simão, Waterhouse, Ioannidis, Kriventseva, & Zdobnov, 2015) within the gVolante plataform (https://gvolante.riken.jp)

# ---------------------------------------------------------------------------------------

# Evaluation of circularization
echo "Running Circlator"
#USAGE: For each Assembly.fasta, run the command. Due short reads files needs to be provided in single fasta file, we first will convert them using fq2fa.  

fq2fa --merge --filter shortreads1.fastq shortreads2.fastq short_reads_together.fa

circlator all Assembly.fasta short_reads_together.fa Circlator_output

# ---------------------------------------------------------------------------------------

echo "Mapping reads and evaluation with QUALIMAP"
# USAGE: For each assembly, short and long reads files need to be mapped back to the assembly. Then, the bam file will be analyzed by qualimap.

bwa index Assembly.fasta

#Mapping short reads
bwa mem -t 4 Assembly.fasta shortreads1.fastq shortread2.fastq | samtools sort > assembly_map_short-reads.bam

#Evaluation of mapping (QUALIMAP)
samtools index assembly_map_short-reads.bam
qualimap bamqc -bam assembly_map_short.bam -outfile Qualimap_map_short.pdf

#Mapping long reads
bwa mem -t 4 Assembly.fasta longReads.fastq | samtools sort > assembly_map_long-reads.bam

#Evaluation of mapping (QUALIMAP)
samtools index assembly_map_long-reads.bam
qualimap bamqc -bam assembly_map_long-reads.bam -outfile Qualimap_map_long.pdf

# ---------------------------------------------------------------------------------------

echo "Running Integron_Finder"
integron_finder Assembly.fasta --func_annot --local_max --path_func_annot bank_hmm

cd ..

echo " "
echo "******* Correctness analysis *********"
echo " "

mkdir correctness
cd correctness 

#See results of Contiguity of QUALIMAP report.
quast *.fasta -o Quast_results_with_reference -r Reference.fasta

# ---------------------------------------------------------------------------------------

#BLAST of known Sanger sequences: For each Assembly.fasta
makeblastdb –in knonw_sequences.fasta –dbtype nucl –parse_seqids
blastn –db knonw_sequences.fasta –query Assembly.fasta –out blast_results.out

cd ..

# END
