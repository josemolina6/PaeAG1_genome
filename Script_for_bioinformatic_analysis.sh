# HIGH QUALITY DE NOVO ASSEMBLY OF A MULTIDRUG RESISTANCE ST-111 P. AERUGINOSA GENOME: BECNHMARK OF HYBRID AND NON-HYBRID ASSEMBLERS

echo "SCRIPT FOR BIOINFORMATICS ANALYSIS"
echo "----------------------------------------------"
echo "Implemented by Jose Arturo Molina Mora"
echo "University of Costa Rica"
echo "----------------------------------------------"


#SCRIPT 1: SHORT READS ASSEMBLIES

echo "Starting genome assembly (short reads) with" $1 $2 "data"

$1  Read1.fastq
$2  Read2.fastq

# BLOCK 1: QC by fastQC and trimming by trimmomatic:

fastqc -o QC_raw-data/ $1 $2

trimmomatic PE $1 $2 -baseout ./02_Trimming_QC/AG1_trim.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:5 SLIDINGWINDOW:4:20 MINLEN:101 #due high coverage, full quality complete reads were kept 

fastqc -o 02_Trimming_QC ./02_Trimming_QC/*.fastq

# MultiQC: all files of fastQC in same directory and run:

multiqc .

# processed clean data:

#$1  ShortRead1.fastq
#$2  ShortRead2.fastq

#BLOCK 2: Assemblies using 6 algorithms

mkdir All_assemblies

echo "Running Velvet (Kmergenie to evaluate multiple kmer values)"

#Kmergenie searching between 11 and 101 kmers, each 10 values. 

ls -1 *.fastq > Short-reads_files
kmergenie Short-reads_files -k 101 -l 11 -s 10
# kmer=51 was selected as best value

# Assembly
velveth Velvet/ 51 -shortPaired -separate -fastq $1 $2
velvetg Velvet/ -ins_length 101 -exp_cov auto -cov_cutoff auto -amos_file yes


echo "Running SPAdes"

spades.py -o 02_SPAdes -1 $1 -2 $2 --careful

echo "Running IDBA"

fq2fa --merge --filter $1 $2 Reads_together.fa

idba_ud -r Reads_together.fa -o 03_IDBA_UD

echo "Running Megahit"

megahit -1 $1 -2 $2 -o 04_Megahit

echo "Running SKESA"


mkdir 05_SKESA
skesa --fastq $1,$2 > ./05_SKESA/Skesa_contigs.fa

echo "Running UNICYCLER"

unicycler -1 $1 -2 $2 -o 06_Unicycler


#Calculating metrics: assembly-stats *.fasta

assembly-stats -t ./All_assemblies/*.fasta > ./All_assemblies/00_Metrics_all_assemblies.txt

# ---------------------------------------------------------------------------------------

#SCRIPT 2: LONG READS ASSEMBLIES

echo "Starting genome assembly (long reads) with" $1 "data"

# BLOCK 1: QC by Poretools and filtering by Porechop and Filtlong:

cd Nanopore_data

poretools fastq Long_data/RawlongAG1.fast5
poretools yield_plot --plot-type reads Long_data/
poretools squiggle Long_data/RawlongAG1.fast5
poretools stats Long_data/
poretools qualpos Long_data/

porechop -i RawlongAG1.fastq.gz -o preLongReads.fastq.gz

filtlong --min_length 1000 --mean_q_weight 10 --keep_percent 90 --target_bases 500000000 preLongReads.fastq.gz | gzip > LongReads.fastq.gz


$1  LongReads.fastq.gz

#BLOCK 2: Assemblies using 3 algorithms

mkdir All_assemblies

echo "Running Canu"

canu -p canu -d 01_Canu genomeSize=7.2m -nanopore-raw $1

echo "Running Flye"

flye --nano-raw $1 --genome-size 7.2m --out-dir 02_flye

echo "Running Unicycler"

unicycler -l $1 -o 03_unicycler

# ---------------------------------------------------------------------------------------

#SCRIPT 3: HYBRID ASSEMBLIES

echo "Starting genome assembly (hybrid) with" $1 $2 $3 "data"

# Filtered data by quality:
$1  ShortRead1.fastq
$2  ShortRead2.fastq
$3  LongReads.fastq.gz

echo "Running Unicycler"

unicycler -1 $1 -2 $2 -l $3 -o hyb_unicycler

echo "Running IDBA-hyb"

fq2fa --merge --filter $1 $2  Short_reads_together.fa
idba_hybrid -r Short_reads_together.fa -o hyb_idba -l $3


echo "Running SPAdes-hyb"
spades.py -o hyb_spades -1 $1 -2 $2 --nanopore $3


#NOTE: Next analysis were run for each algorithm, but here we call it "Assembly.fasta" to refer as generic file for each fasta file from assemblers.  

# ---------------------------------------------------------------------------------------

#SCRIPT 4: Filtering out contigs (less than 1000 bp) and statistics

echo "Starting contig filtering by size"

# For each assembly: 
faidx Assembly.fasta -a 1000,50000000 > ../Assembly.fasta

#statistics:
assembly-stats *.fasta

# ---------------------------------------------------------------------------------------

#SCRIPT 5: Scaffolding 

echo "Scaffolding"

mkdir Scaffolds

# For each assembly:
java -jar medusa.jar -f reference_genomes/ -i ./Assembly.fasta -v -random 5 -o ./03_Scaffolds/

# ---------------------------------------------------------------------------------------

#SCRIPT 6: Evaluation of assemblies 

#Each assembly at scaffold level was evaluated, here we call it "Assembly_scaf.fasta" as a generic name

# Evaluation of circularization
echo "Run Circlator"

circlator all Assembly_scaf.fasta Short_reads_together.fa Circlator_output

echo "Running QUAST"

# Comparison of metrics for each assembly using QUAST

quast ./Scaffolds/*.fasta -o Quast_scaffolds
quast ./Scaffolds/*.fasta -o Quast_with_reference -r ./reference_genomes/Reference.fasta
# Note: QUAST was run again after final assembly (post polishing)

# ---------------------------------------------------------------------------------------

#SCRIPT 7: Annotation using Prokka

echo "Annotation using Prokka"

# For each scaffold assembly:
prokka --compliant --proteins Pseudomonas_aeruginosa_PAO1_107.gbk Assembly_scaf.fasta --outdir Assembly_Annotation
# NOte: Prokka was also run after polishing of seleted genome

# ---------------------------------------------------------------------------------------

#SCRIPT 8: Polishing

# This step was run for the "winner assembly": Unicycler_hybrid_scaf.fasta

echo "Starting genome polishing"

#First: BWA for aligment reads to genome:
bwa index Unicycler_hybrid_scaf.fasta
#Align the short reads:
bwa mem -t 4 Unicycler_hybrid_scaf.fasta ShortRead1.fastq ShortRead2.fastq | samtools sort > pilon_aln.bam

#Index the files:
samtools index pilon_aln.bam
samtools faidx Unicycler_hybrid_scaf.fasta

pilon --genome Unicycler_hybrid_scaf.fasta --frags pilon_aln.bam --output PaeAG1_genome_assembly_final.fasta --fix all --mindepth 0.1 --changes --verbose  -Xmx4096m 

# ---------------------------------------------------------------------------------------

#SCRIPT 9: Integrons searching

# Run Integron_finder

echo "Starting Integron identification"

integron_finder ../PaeAG1_genome_assembly_final.fasta --func_annot --local_max --path_func_annot bank_hmm


# ---------------------------------------------------------------------------------------

#SCRIPT 10: Mapping short reads to final assembly

echo "Mapping short reads and evaluation with QUALIMAP"

#Mapping
bwa index PaeAG1_genome_assembly_final.fasta
bwa mem -t 4 PaeAG1_genome_assembly_final.fasta ShortRead1.fastq ShortRead2.fastq | samtools sort > assembly_map_short-reads.bam

#Evaluation of mapping (QUALIMAP)
samtools index assembly_map_short-reads.bam
qualimap bamqc -bam assembly_map_short.bam -outfile Qualimap_map_short.pdf

# ---------------------------------------------------------------------------------------

#SCRIPT 11: Mapping long reads to final assembly

echo "Mapping long reads and evaluation with QUALIMAP"

#Mapping
bwa index PaeAG1_genome_assembly_final.fasta
bwa mem -t 4 PaeAG1_genome_assembly_final.fasta LongReads.fastq.gz | samtools sort > assembly_map_long-reads.bam

#Evaluation of mapping (QUALIMAP)
samtools index assembly_map_long-reads.bam
qualimap bamqc -bam assembly_map_long-reads.bam -outfile Qualimap_map_long.pdf


# ---------------------------------------------------------------------------------------

#SCRIPT 12 EXTRA: Mapping of RNASeq reads to final assembly (data not fot this study but mapping of RNAseq to evaluate assembly completeness) 

echo "Mapping RNASeq reads and evaluation with QUALIMAP"

#Mapping using HISAT2: 12 samples labeled as A, B, C and D (3 replicates each)

for i in A1 A2 A3 B1 B2 B3 C1 C2 C3 D1 D2 D3
do
hisat2 -x PaeAG1_genome_assembly_final.fasta  -1 ${i}_1P.fq.gz -2 ${i}_2P.fq.gz| samtools sort > ${i}_sort.bam
samtools index ${i}_sort.bam
#QUALIMAP
qualimap bamqc -bam ${i}_sort.bam -outdir qualimap_bamqc/${i}_bamqc -outformat PDF:HTML --java-mem-size=10G
echo "${i} processed"
done

# All analysis in a single report with MultiQC
multiqc ./qualimap_bamqc/. -o qualimap_bamqc

# ---------------------------------------------------------------------------------------

# END
