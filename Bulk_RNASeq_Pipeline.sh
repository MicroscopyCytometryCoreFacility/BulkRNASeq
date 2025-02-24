#!/bin/bash
# Author:		Michael de Kok
# Created:		Feb 20, 2023
# Last Updated: 	Feb 24, 2025
# Title:		Modular and Fully Automatic Linux Bash Pipieline for RNA-Bulk Seq Analysis
# Modules:		SRA FastqDump, Host Genome Downloading & STAR Indexing, STAR Aligning & Samtools Sorting/Indexing, and FeatureCounts
# Software Versions:	SRA Toolkit 3.1.1, STAR 2.7.10a, samtools 1.13, and featureCounts v2.0.3 
# Prerequisites:	Make sure all of the softwares except SRA Toolkit are installed and on the PATH

echo "Starting BulkRNA Sequencing Pipeline by Michael 'MikedeKokkie' de Kok, for more information please reach out through our GitHub: https://github.com/MicroscopyCytometryCoreFacility/BulkRNASeq"

# -------------- GENERAL PARAMETERS --------------

RUN_SRA=true				 			# Can be "true" or "false"
OBTAIN_AND_INDEX_HOST_GENOMES=true	 			# Can be "true" or "false"
RUN_STAR=true				 			# Can be "true" or "false"
RUN_COUNTS=true 			 			# Can be "true" or "false"
PAIRED_ENDED=true			 			# Can be "true" or "false"
HOST_GENOME="Human"			 			# Can be "Human" or "Mouse"
NUM_THREADS=8				 			# Can be any integer 
SRA_FASTQS_TO_GET="SRR00000000 SRR00000001"			# Can be any length, only needed if RUN_SRA=true
# Where the SRR codes correspond to the SRA Accession Codes of the Dataset Runs that you want to download the .fastq sequences of.

# -------------- FOLDER SPECIFICATIONS ---------------

SRA_FOLDER="/mnt/Data2/mike/Algorithms/sratoolkit.3.1.1-ubuntu64" 	# Only needed if RUN_SRA=true

MAIN__FOLDER="/mnt/Data2/mike/RNA_Datasets/DatasetName"			# Main Dataset Folder, always needed
FASTQ_FOLDER="${MAIN__FOLDER}/Fastq"
BAMS__FOLDER="${MAIN__FOLDER}/Bams"
COUNT_FOLDER="${MAIN__FOLDER}/Results"

REF_FOLDER="/mnt/Data2/mike/Reference_Genomes_113"			# Main Reference Genome Folder, always needed
if [ $HOST_GENOME = "Human" ]
	then
		RAW_HOST="${REF_FOLDER}/Human/"
		INDXHOST="${REF_FOLDER}/Human/STAR-Index"
		GTF_FILE="${REF_FOLDER}/Human/Homo_sapiens.GRCh38.113.gtf"
fi
if [ $HOST_GENOME = "Mouse" ]
	then
		RAW_HOST="${REF_FOLDER}/Mouse/"
		INDXHOST="${REF_FOLDER}/Mouse/STAR-Index"
		GTF_FILE="${REF_FOLDER}/Mouse/Mus_musculus.GRCm39.113.gtf"
fi

if [ ! -d "$MAIN__FOLDER" ]; then
    mkdir "$MAIN__FOLDER"
    echo "Directory $MAIN__FOLDER created."
fi
if [ ! -d "$FASTQ_FOLDER" ]; then
    mkdir "$FASTQ_FOLDER"
    echo "Directory $FASTQ_FOLDER created."
fi
if [ ! -d "$BAMS__FOLDER" ]; then
    mkdir "$BAMS__FOLDER"
    echo "Directory $BAMS__FOLDER created."
fi
if [ ! -d "$COUNT_FOLDER" ]; then
    mkdir "$COUNT_FOLDER"
    echo "Directory $COUNT_FOLDER created."
fi
if [ ! -d "$RAW_HOST" ]; then
    mkdir "$RAW_HOST"
    echo "Directory $RAW_HOST created."
fi
if [ ! -d "$INDXHOST" ]; then
    mkdir "$INDXHOST"
    echo "Directory $INDXHOST created."
fi

# --------------- SRA SCRIPT ----------------

if [ $RUN_SRA = true ]
	then
		echo "Obtaining Fasta files from SRA Database..." 
		cd $FASTQ_FOLDER
		for i in $SRA_FASTQS_TO_GET
		do
		  echo $i
		  #${SRA_FOLDER}/bin/fastq-dump -X 5 -Z $i
		  ${SRA_FOLDER}/bin/fastq-dump $i --gzip 
		  #${SRA_FOLDER}/bin/fastq-dump -I --split-files $i
		done
fi

# --------------- STAR GENOME MODULE ----------------

if [ $OBTAIN_AND_INDEX_HOST_GENOMES = true ]
	then
		if [ $HOST_GENOME = "Human" ]
			then
				cd $RAW_HOST
				echo "Obtaining Raw Genome..." 
				curl -o Homo_sapiens.GRCh38.113.gtf.gz ftp://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
				gunzip Homo_sapiens.GRCh38.113.gtf.gz

				curl -o ensemble113genome.fa.gz ftp://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
				gunzip ensemble113genome.fa.gz
		fi
		if [ $HOST_GENOME = "Mouse" ]
			then
				cd $RAW_HOST
				echo "Obtaining Raw Genome..." 
				curl -o Mus_musculus.GRCm39.113.gtf.gz ftp://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz
				gunzip Mus_musculus.GRCm39.113.gtf.gz

				curl -o ensemble113genome.fa.gz ftp://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz
				gunzip ensemble113genome.fa.gz
		fi

		echo "Indexing Genome..."
		STAR --runThreadN $NUM_THREADS \
		--runMode genomeGenerate \
		--genomeDir $INDXHOST \
		--genomeFastaFiles ${RAW_HOST}/ensemble113genome.fa \
		--sjdbGTFfile $GTF_FILE \
		--sjdbOverhang 99
fi

# --------------- STAR ALIGNMENT MODULE ----------------

if [ $RUN_STAR = true ]
	then
		echo "Loading Indexed Genome..."
		STAR --genomeLoad LoadAndExit --genomeDir $INDXHOST

		echo "Starting Main Aligning Loop..."
		for fq in `find -L ${FASTQ_FOLDER} -iname "*_R1_001.fastq.gz"`; do
		echo "File name: $fq" 
		bn=`basename $fq _R1_001.fastq.gz`
		echo "Sample name: $bn" 
		if [ ! -e ${bn}.sam ]
		        then
				if [ $PAIRED_ENDED = true ]
					then
						echo "Aligning... "
						STAR --runThreadN $NUM_THREADS \
							--genomeLoad LoadAndKeep \
							--genomeDir $INDXHOST \
							--readFilesIn <(gunzip -c ${FASTQ_FOLDER}/${bn}_R1_001.fastq.gz) <(gunzip -c ${FASTQ_FOLDER}/${bn}_R2_001.fastq.gz) \
							--outFileNamePrefix ${BAMS__FOLDER}/${bn}. \
							--outFilterMismatchNoverLmax 0.04 \
							#> ${COUNT_FOLDER}/STAR.screen-output.txt
						echo "Sorting..."
						samtools view -bS ${BAMS__FOLDER}/${bn}.Aligned.out.sam | samtools sort -o ${BAMS__FOLDER}/${bn}.sorted.bam # unsure if.bam 
						echo "Indexing..."		
						samtools index ${BAMS__FOLDER}/${bn}.sorted.bam
				fi
				if [ $PAIRED_ENDED = false ]
					then
						echo "Aligning... "
						STAR --runThreadN $NUM_THREADS \
							--genomeLoad LoadAndKeep \
							--genomeDir $INDXHOST \
							--readFilesIn <(gunzip -c ${FASTQ_FOLDER}/${bn}_R1_001.fastq.gz) \
							--outFileNamePrefix ${BAMS__FOLDER}/${bn}. \
							--outFilterMismatchNoverLmax 0.04 \
							#> ${COUNT_FOLDER}/STAR.screen-output.txt
						echo "Sorting..."
						samtools view -bS ${BAMS__FOLDER}/${bn}.Aligned.out.sam | samtools sort -o ${BAMS__FOLDER}/${bn}.sorted.bam # unsure if.bam 
						echo "Indexing..."		
						samtools index ${BAMS__FOLDER}/${bn}.sorted.bam
				fi
		fi
		
		done
		
		echo "Unloading Indexed Genome..."
		STAR --genomeLoad Remove --genomeDir $INDXHOST
fi

# --------------- FEATURECOUNTS MODULE ----------------

if [ $RUN_COUNTS = true ]
	then
		if [ $PAIRED_ENDED = true ]
			then
				echo "Mapping Alignment into Counts..."
				featureCounts -T $NUM_THREADS \
					-p \
					-a $GTF_FILE \
					-o ${COUNT_FOLDER}/counttable.txt \
					${BAMS__FOLDER}/*.sorted.bam \
					2> /${COUNT_FOLDER}/featurecounts.screen-output.log
		fi
		if [ $PAIRED_ENDED = false ]
			then
				echo "Mapping Alignment into Counts..."
				featureCounts -T $NUM_THREADS \
					-a $GTF_FILE \
					-o ${COUNT_FOLDER}/counttable.txt \
					${BAMS__FOLDER}/*.sorted.bam \
					2> /${COUNT_FOLDER}/featurecounts.screen-output.log
		fi
fi		

# --------------- END OF PIPELINE ----------------

echo "Pipeline Finished."
read -n1 -r -p "Press any key to continue..." key
