#!/bin/bash

##### Data cleaning ATAC-seq data ######

########################################
## Initial quality assessment
########################################
## Assess read quality using FastQC. 
## Save results in `results_fastqc/fastqc_raw/`

mkdir -p results/results_fastqc/fastqc_raw/

for filename in ./data/*fastq.gz ;
do
echo "Starting FastQC analysis of" $filename

fastqc $filename \
       -o results/results_fastqc/fastqc_raw/ \
       -t 15
done

## Download
aws s3 sync ~/project/results/ s3://kadm-results/
#aws s3 sync s3://kadm-results ./results

########################################
## Adapter removal (`fastq_trim/`)
########################################
## Remove adapters (14 bp)
## Remove reads with > 1 ambiguous base
## Trim ends until reach base with quality 30+
## Remove reads < 15 bp

mkdir -p results/fastq_trim

paste <(ls data/*R1_001.fastq.gz) <(ls data/*R2_001.fastq.gz) |

while read file1 file2;
do
  name1=$(paste -d '\0' \
            <(echo 'results/fastq_trim/') \
            <(awk -F'[_]S' '{print $1}' <(basename $file1)))
  
  AdapterRemoval --file1 $file1 --file2 $file2 \
    --basename $name1 --gzip \
    --trim5p 14 --maxns 1 --minlength 15 \
    --trimqualities --minquality 30 \
    --threads 15
done

########################################
## Re-assess quality
########################################
## Assess trimmed read quality using FastQC.
mkdir -p results/results_fastqc/fastqc_trim/

for filename in results/fastq_trim/*pair[12].truncated.gz ;
do
echo "Starting FastQC analysis of" $filename

fastqc $filename \
       -o results/results_fastqc/fastqc_trim/ \
       -t 15
done

## Download
aws s3 sync ~/project/results/ s3://kadm-results/
#aws s3 sync s3://kadm-results ./results

########################################
## Alignment (`bam/`)
########################################
## Get ref data files
mkdir -p ref
sudo chmod 777 -R ref

s3fs kadm-ref ref -o passwd_file=~/.passwd-s3fs \
    -o default_acl=public-read -o uid=1000 -o gid=1000 -o umask=0007
     
## Align with STAR
mkdir -p results/bam/

paste <(ls results/fastq_trim/*pair1.truncated.gz) \
      <(ls results/fastq_trim/*pair2.truncated.gz) |

while read file1 file2;
do
    echo "Aligning" $(basename  -- "$file1");
    
    name=$(paste -d '\0' \
            <(echo 'results/bam/') \
            <(awk -F'[.]pair' '{print $1}' <(basename $file1)) \
            <(echo '_'))
    
    STAR --genomeDir ~/project/ref/STARindex \
         --readFilesIn $file1 $file2 \
         --readFilesCommand zcat \
         --outFileNamePrefix $name \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN 15 \
         --runRNGseed 8756
done

########################################
## Assess aligned data
########################################
## median CV of gene model coverage

mkdir -p results/results_cleaning/

for bam_file in results/bam/*sortedByCoord.out.bam ;
do
    java -XX:ParallelGCThreads=15 \
        -jar apps/anaconda3/share/picard-2.22.3-0/picard.jar \
        CollectRnaSeqMetrics \
        REF_FLAT=ref/PICARDref/refFlat.ensembl.txt \
        INPUT=$bam_file  \
        OUTPUT=results/results_cleaning/temp.tsv \
        ASSUME_SORTED=true \
        STRAND_SPECIFICITY=NONE \
        MINIMUM_LENGTH=500
    
    #Append results
    echo $bam_file >> results/results_cleaning/bam.metrics.tsv
    cat results/results_cleaning/temp.tsv >> results/results_cleaning/bam.metrics.tsv
    #Remove this iteration
    rm results/results_cleaning/temp.tsv
done
        
## mapped_reads_w_dups aka alignment percentage

for bam_file in results/bam/*sortedByCoord.out.bam ;
do
    echo "Processing" $bam_file
    echo $bam_file >> results/results_cleaning/summary.alignment.tsv
    
    samtools flagstat -@ 15 $bam_file \
    >> results/results_cleaning/summary.alignment.tsv
done

## Insert sizes
## Calculate insert lengths for all sequences in all sorted BAM files.

for bam_file in results/bam/*sortedByCoord.out.bam ;
do
  echo "Calculating inserts for" $bam_file
  
  sampID=$(awk -F'[_]' '{print $1}' <(basename $bam_file))

  samtools stats $bam_file --threads 14 \
    | grep ^IS \
    | cut -f 2- \
    > temp.insert.txt
  
  awk -v sampID=$sampID '{print sampID"\t"$0}' temp.insert.txt \
    >> results/results_cleaning/insert.size.tsv
  
  rm temp.insert.txt

done

########################################
## Peak calls
########################################
# -j ATACseq mode
# -r remove PCR dups
# -e remove specific chromosomes
# -d expand cut sites to X bp. Default is 100
# Default removes unpaired alignments and adjusts for Tn5 insertion
# Multi-mapping alignments are counted as fractions 
#   if secondary alignment has equal MAPQ to primary

#Re-sort by query name

for bam_file in results/bam/*sortedByCoord.out.bam ;
do
    echo "Sorting" $bam_file
    
    name=$(paste -d '\0' \
            <(echo 'results/bam/') \
            <(awk -F'[_]' '{print $1}' <(basename $bam_file)) \
            <(echo '_Aligned.sortedBySamp.out.bam'))
  
    samtools sort -n -o $name -@ 15 $bam_file
done

#PEAKS

mkdir -p results/peaks/

for bam_file in results/bam/*sortedBySamp.out.bam ;
do
    echo "Determining Genrich peaks in" $bam_file
    
    name1=$(paste -d '\0' \
            <(echo 'results/peaks/') \
            <(awk -F'[_]' '{print $1}' <(basename $bam_file)) \
            <(echo '.narrowPeak'))
    name2=$(paste -d '\0' \
            <(echo 'results/peaks/') \
            <(awk -F'[_]' '{print $1}' <(basename $bam_file)) \
            <(echo '.bed'))
    #Call peaks
    Genrich -t $bam_file \
            -o $name1 \
            -j -r -e chrM,chrY,chrY_KI270740v1_random \
            -d 100 -b $name2
done

########################################
## Save to S3 storage
########################################

aws s3 sync ~/project/results/ s3://kadm-results

################# END ##################