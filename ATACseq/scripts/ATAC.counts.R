"Count reads in paired-end BAM files in ATAC-seq peaks

#################

Kim Dill-Mcfarland
University of Washington, kadm@uw.edu
Copyright (C) 2019 Kim Dill-Mcfarland
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Input parameters:
REQUIRED
  peaks = data frame of peaks in annotation style with 
          peak, Chr, Start, End, Strand
  bam = character list of file paths to bam files
  minF = integer giving the minimum fragment length for paired-end reads
  maxF = integer giving the maximum fragment length for paired-end reads
  filename = file path for output to be written to
  
OPTIONAL
  cores = set number of cores to use for parallel computing. Default is 1
  
Example
  peaks <- data.frame(
    peak = c('peak1', 'peak2'),
    Chr = c('chr1', 'ch2'),
    Start = c(100, 1000),
    End = c(200, 2000),
    Strand = c('*', '*')
  )
  
  bam <- dir('data/data_trim_align/', 
                pattern = '*.bam$', full.names = TRUE) %>% 
                gsub('//', '/', .)

  count.atac(peaks=peaks, bam=bam, minF=50, maxF=100, filename='results.csv', cores=5)
"

#################

count.atac <- function(peaks, bam, minF, maxF, filename, cores=1) 
{
  #####Load packages#####
  require(tidyverse)
  require(Rsubread)
  require(parallel)
  
  #Set seed
  set.seed(589)
  
  #####Format data#####
  ### Peaks
  if(!is.data.frame(peaks)){
    stop("ERROR: Peaks must be a data frame")}
  if(any(colnames(peaks) != c("peak","Chr","Start","End","Strand"))){
    stop("ERROR: Peaks must contain peak, Chr, Start, End, Strand")}
  ### Bam files
  if(!is.character(bam)){
    stop("ERROR: BAM file list must be a character")}
  
  #####Counting#####
  print("Counting reads in peaks")
  
  Fcounts <- featureCounts(bam, annot.ext = peaks,
                           nthreads=cores, 
                           #Fragment size
                           minFragLength=minF, maxFragLength = maxF,
                           #Quality cutoffs
                           isPairedEnd = TRUE,            ## Must be concordant pairs
                           requireBothEndsMapped = TRUE,
                           minMQS = 20,                   ## MAPQ > 20
                           countMultiMappingReads = TRUE, ## Include multi-mapping reads
                           fraction = TRUE,               ## as fractional counts
                           allowMultiOverlap = TRUE       ## Reads can be assigned to >1 peak
  )
  
  print("Counting complete. Saving results")
  
  #Extract count table
  Fcounts.df <- as.data.frame(Fcounts$counts) %>% 
    rownames_to_column("peak")
  #Format column names
  samp <- colnames(Fcounts.df)[2:length(colnames(Fcounts.df))] %>% 
    gsub(., pattern = "_Aligned.sortedByCoord.out.bam", 
         replacement="")
  
  #Fix short RSID
  samp2 <- sapply(samp, FUN = function(x) ifelse(nchar(x) == 7,
                                   gsub("RS102", "RS1020", x), x))
  colnames(Fcounts.df) <- c("peak", samp2)
  
  #Save to disk
  write_csv(Fcounts.df, filename)
  
  print("FIN")
}
