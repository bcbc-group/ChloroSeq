#!/bin/sh
###############################################################
# pipeline for running RNAseq DE on multiple single end files #
#                                                             #
# usage:                                                      #
#                                                             #
#      RNAseq_pe_pipeline.sh $indir $CPU $bin                 #
#                                                             #
###############################################################
#get commandline options
cd $1        #input dir
CPU=$2       #no CPUS
bin=$3       #full path to executables
ref=$4       #full path and file name of indexed reference fasta
gtf=$5       #full path and file name of gtf file
chloro=$6    #chloroplast name, ex: NC_022850.1

##Attention!!! Go to bottom of script and adjust chloroseq inputs

#filter rRNA
ls *.fastq |parallel -j 1 $bin/sortmerna-2.1b/sortmerna --ref $bin/sortmerna-2.1b/rRNA_databases/silva-euk-18s-id95.fasta,$bin/sortmerna-2.1b/rRNA_databases/silva-euk-18s-id95-db:$bin/sortmerna-2.1b/rRNA_databases/silva-euk-28s-id98.fasta,$bin/sortmerna-2.1b/rRNA_databases/silva-euk-28s-id98-db --reads {} --aligned {.}_rRNA --other {.}_other --num_alignments 1 --fastx -a $CPU --log

#clean reads
for file in `dir -d *other.fastq`; do
   out1=`echo "$file" |sed 's/.fastq/_cln.fastq/'`

   java -jar $bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE $file $out1 ILLUMINACLIP:$bin/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:50

done

#get read count
nohup wc *_cln.fastq | awk '{print $1/4 "\t" $4}' > filtered_read_count.txt 2>&1 &

#map reads to genome
for file in `dir -d *_cln.fastq` ; do

    #create output file name                                                                                          
    samfile=`echo "$file" | sed 's/.fastq/.sam/'`
    
    hisat2 --max-intronlen 120000 --dta --rna-strandness RF -p $CPU -x $ref -U $file -S $samfile

done

#convert files and get metrics
ls *.sam |parallel --gnu -j $CPU samtools view -Sb -o {.}.bam {}
ls *.bam |parallel --gnu -j $CPU samtools sort {} {.}.sort.bam
ls *.sort.bam |parallel --gnu -j $CPU samtools flagstat {} ">" {.}.flagstat

#calc mapping efficiency
ls *sort.bam |parallel -j $CPU samtools view -F 260 {} ">" {.}_mapped

#get fpkm and counts
mkdir ballgown

for file in `dir -d *.sort.bam` ; do
    
    outfile=`echo "$file" | sed 's/.bam/.gtf/'`  
    outdir=`echo "$file" |sed 's/.bam//'`
    
    stringtie  -e -B -p $CPU -G $gtf -o ballgown/$outdir/$outfile $file

done

#create count table for import to R with DESeqDataSetFromMatrix (DEseq) and DGEList (EdgeR)  
python3 $bin/prepDE.py -i ballgown -g gene_count_matrix.csv -t transcript_count_matrix.csv 

#extract chloroplast bams                                                                                          
ls *sort.bam |parallel -j $CPU bamtools split -in {} -reference

#calc chloro mapping efficiency
ls *$chloro.bam |parallel -j 20 samtools view -F 260 {} ">" {.}_mapped

#filter reads that are spliced >= 2500 bp                                                                          
ls *$chloro.bam | parallel -j $CPU samtools view -o {.}.sam {}
ls *$chloro.sam | parallel -j $CPU $bin/cigar_filter.pl {} ">" {.}.filt.sam
ls *filt.sam |parallel -j 10 samtools view -SbT /home/srs57/Leila4888_2017_07_01/fifth/reference/NC_022850.1.fasta -o {.}.bam {}
ls *filt.bam |parallel -j $CPU samtools index {}

#run chloroseq (output can be plotted with accompanying R scripts                                                  
ls *filt.bam |parallel -j $CPU $bin/ChloroSeq_for_bin/chloroseq_v2.pl -a all -b {} -e /home/srs57/Leila4888_2017_07_01/fifth/reference/NC_022850.1_exon.gff3 -i /home/srs57/Leila4888_2017_07_01/fifth/reference/NC_022850.1_intron.gff3 -g 135516 -n $chloro -f /home/srs57/Leila4888_2017_07_01/ref_genome/Sitalica_312_v2_mscaff14_NC.fasta -s /home/srs57/Leila4888_2017_07_01/fifth/reference/Si_plastid_splice_sites_corrected.txt -v /home/srs57/Leila4888_2017_07_01/fif
th/reference/all_samples_edited.gff3 -k

#make bedgraphs for plotting
ls *coverage_* |parallel -j 20 $CPU \'{print \$1 \"\\t\" \$2 - 1 \"\\t\" \$3 \"\\t\" \$4}\' {} ">" {.}.bedgraph
