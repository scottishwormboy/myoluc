#Map Thomas data
#consists of 6 individuals at fairly low coverage to get high quality SNPs
#plus 5 populations

#ada10, chos 7

workdir=$HOME/Datafiles/Bats/myoluc
mkdir $HOME/Datafiles/Bats/labbook

#indvs
readdir1=/agf/illum/LIMS12904_3194059ea08c3479/Trimmed

#pops
readdir2=/agf/illum/LIMS12896_3487f003e8df1fd6/Trimmed

mkdir -p $workdir
cd $workdir
mkdir $workdir/ref
mkdir -p $workdir/indvs/bams
mkdir -p $workdir/pops/bams
mkdir -p $workdir/indvs/bamstats
mkdir -p $workdir/pops/bamstats
mkdir -p $workdir/indvs/coverage
mkdir -p $workdir/pops/coverage
mkdir -p $workdir/indvs/logs
mkdir -p $workdir/pops/logs
mkdir -p $workdir/indvs/vcfs
mkdir -p $workdir/pops/vcfs


#versions of software
bowtie2 --version > $workdir/version.txt
samtools --version >> $workdir/version.txt
java -version  2>> $workdir/version.txt
ls -l ~/bin/picard >> $workdir/version.txt
java -jar ~/bin/GATK3_4/GenomeAnalysisTK.jar -version >> $workdir/version.txt



wget ftp://ftp.broadinstitute.org/pub/assemblies/mammals/microbat/myoluc2/assembly.bases.gz 
wget ftp://ftp.broadinstitute.org/pub/assemblies/mammals/microbat/myoluc2/assembly.agp
wget https://github.com/WormBase/wormbase-pipeline/blob/master/for_cvs/agp2fasta.pl
gunzip $workdir/assembly.bases.gz 
mv assembly.bases contig.fa
perl agp2fasta.pl assembly.agp contig.fa 
mv supercontigs.fa $workdir/ref/myoluc2.0.fa

ref=$workdir/ref/myoluc2.0.fa

bowtie2-build $workdir/ref/myoluc2.0.fa $workdir/ref/myoluc2.0.fa

java -jar ~/bin/picard/picard.jar CreateSequenceDictionary R= $workdir/ref/myoluc2.0.fa O= $workdir/ref/myoluc2.0.dict

find $readdir1 -name '*R1*' > $workdir/indvsR1.txt
find $readdir1 -name '*R2*' > $workdir/indvsR2.txt
awk 'BEGIN {FS="/"} {print $6}' $workdir/indvsR1.txt > $workdir/indvs_samples.txt
paste $workdir/indvs_samples.txt $workdir/indvsR1.txt $workdir/indvsR2.txt > $workdir/indvs_sheet.txt

mkdir /scratch/stevep11bam

#put first sample through manually, change this back
#cat $workdir/indvs_sheet.txt | while read samp R1 R2
awk 'NR>1 {print $0}' $workdir/indvs_sheet.txt | while read samp R1 R2
do
  echo $samp
  bowtie2 -x $ref --sensitive-local -p 64 -1 $R1 -2 $R2 | samtools view -b - -o /scratch/stevep11bam/1.${samp}.bam
  samtools sort /scratch/stevep11bam/1.${samp}.bam /scratch/stevep11bam/2.${samp}
  samtools flagstat /scratch/stevep11bam/2.${samp}.bam > $workdir/indvs/logs/${samp}.flagstat &
  java -jar ~/bin/picard/picard.jar AddOrReplaceReadGroups I=/scratch/stevep11bam/2.${samp}.bam O=/scratch/stevep11bam/3.${samp}.bam RGID=${samp} RGLB=${samp} RGPL=Illumina RGSM=${samp} RGPU=${samp}
  java -Xmx5000M -jar ~/bin/picard/picard.jar MarkDuplicates INPUT=/scratch/stevep11bam/3.${samp}.bam OUTPUT=$workdir/indvs/bams/${samp}.rmdup.bam METRICS_FILE=$workdir/indvs/logs/${samp}.rmdupMetrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
  samtools index $workdir/indvs/bams/${samp}.rmdup.bam
  #~/cgrpipe/Modules/BT2_Map -a $R1 -b $R2 -o /scratch/stevep11bam/${samp}.bam -r $ref -t 64 -s $samp -e $workdir/indvs/logs/bt2.${samp}.log
    
done

#call SNPs in individuals
java -jar ~/bin/GATK3_4/GenomeAnalysisTK.jar -R $ref -T HaplotypeCaller \
-I ${workdir}/indvs/bams/Sample_1-602849.rmdup.bam \
-I ${workdir}/indvs/bams/Sample_2-602892.rmdup.bam \
-I ${workdir}/indvs/bams/Sample_3-2025UPMIPRE.rmdup.bam \
-I ${workdir}/indvs/bams/Sample_4-2091UPMIPRE.rmdup.bam \
-I ${workdir}/indvs/bams/Sample_6-26MPRENY.rmdup.bam \
-I ${workdir}/indvs/bams/Sample_8-29FPRENY.rmdup.bam \
-nct 60 -o ${workdir}/indvs/vcfs/combined.indv.v1.vcf

pigz -p 32 ${workdir}/indvs/vcfs/combined.indv.v1.vcf

mkdir /scratch/stevep11bam

#ada08
#cat $workdir/pop_sheet.txt | while read samp R1 R2
#do
#  echo $samp
#  bowtie2 -x $ref --sensitive-local -p 100 -1 $R1 -2 $R2 | samtools view -b - -o /scratch/stevep11bam/1.${samp}.bam
#  samtools sort /scratch/stevep11bam/1.${samp}.bam /scratch/stevep11bam/2.${samp}
#  samtools flagstat /scratch/stevep11bam/2.${samp}.bam > $workdir/pops/logs/${samp}.flagstat &
#  java -jar ~/bin/picard/picard.jar AddOrReplaceReadGroups I=/scratch/stevep11bam/2.${samp}.bam O=/scratch/stevep11bam/3.${samp}.bam RGID=${samp} RGLB=${samp} RGPL=Illumina RGSM=${samp} RGPU=${samp}
#  java -Xmx5000M -jar ~/bin/picard/picard.jar MarkDuplicates INPUT=/scratch/stevep11bam/3.${samp}.bam OUTPUT=$workdir/indvs/bams/${samp}.rmdup.bam METRICS_FILE=$workdir/pops/logs/${samp}.rmdupMetrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
#  samtools index $workdir/pops/bams/${samp}.rmdup.bam
  #~/cgrpipe/Modules/BT2_Map -a $R1 -b $R2 -o /scratch/stevep11bam/${samp}.bam -r $ref -t 64 -s $samp -e $workdir/indvs/logs/bt2.${samp}.log
    
#done


#retry, hangs at the end of mapping
find $readdir2 -name '*R1*' > $workdir/popsR1.txt
sed 's/R1/R2/' $workdir/popsR1.txt > $workdir/popsR2.txt
awk 'BEGIN {FS="/"} {print $6}' $workdir/popsR1.txt > $workdir/pops_samples.txt
paste $workdir/pops_samples.txt $workdir/popsR1.txt $workdir/popsR2.txt > $workdir/pops_sheet2.txt

#cat $workdir/pops_sheet2.txt | while read samp R1 R2
#awk 'NR>1 {print $0}' $workdir/pops_sheet2.txt | while read samp R1 R2
#failed to read some files for unknown reason, check carefully at end (sample_2, lane2 not run, sample 1 not run
#awk '$1 ~ /Sample_1/ {print $0}' $workdir/pops_sheet2.txt | while read samp R1 R2
awk '$1 ~ /Sample_[345]/ {print $0}' $workdir/pops_sheet2.txt | while read samp R1 R2
#awk 'NR>=22 && NR<=25 {print $0}' $workdir/pops_sheet2.txt | while read samp R1 R2
do  
  #randbam=`echo $RANDOM`
  #echo ${samp}.${randbam}
  lane=`echo $R1 | sed 's/.*_L00/L00/'|sed 's/_R1.*//'`
  bowtie2 -x $ref --sensitive-local -p 60 -1 $R1 -2 $R2 | samtools view -b - -o /scratch/stevep11bam/1.${samp}.${lane}.bam
  samtools sort /scratch/stevep11bam/1.${samp}.${lane}.bam /scratch/stevep11bam/2.${samp}.${lane}
  samtools flagstat /scratch/stevep11bam/2.${samp}.${lane}.bam > $workdir/pops/logs/${samp}.${lane}.flagstat &
  java -jar ~/bin/picard/picard.jar AddOrReplaceReadGroups I=/scratch/stevep11bam/2.${samp}.${lane}.bam O=/scratch/stevep11bam/3.${samp}.${lane}.bam RGID=${samp} RGLB=${samp} RGPL=Illumina RGSM=${samp} RGPU=${lane}
  java -Xmx5000M -jar ~/bin/picard/picard.jar MarkDuplicates INPUT=/scratch/stevep11bam/3.${samp}.${lane}.bam OUTPUT=$workdir/pops/bams/${samp}.${lane}.rmdup.bam METRICS_FILE=$workdir/pops/logs/${samp}.${lane}.rmdupMetrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
  samtools index $workdir/pops/bams/${samp}.${lane}.rmdup.bam
  #~/cgrpipe/Modules/BT2_Map -a $R1 -b $R2 -o /scratch/stevep11bam/${samp}.bam -r $ref -t 64 -s $samp -e $workdir/indvs/logs/bt2.${samp}.log
  rm /scratch/stevep11bam/*
done


#example of mpileup
#generate a test file
zcat ${workdir}/indvs/vcfs/combined.indv.v1.vcf.gz|grep -v "#" |awk 'BEGIN {OFS="\t"} $6>200 {print $1,$2}'|head -n 1000 > ${workdir}/indvs/vcfs/tst.sites

mkdir $workdir/pops/mpileup
find $workdir/pops/bams -name "*.bam" > $workdir/pops/mpileup/pop_bamlist.txt
samtools mpileup -b $workdir/pops/mpileup/pop_bamlist.txt \
-f $ref -l ${workdir}/indvs/vcfs/tst.sites -o $workdir/pops/mpileup/tst.pileup
#for a test it will hang until it gets to end of file, spools lines from memory, may produce 990 with last 10 hours later

samtools mpileup -f $ref -l ${workdir}/indvs/vcfs/tst.sites \
-o $workdir/pops/mpileup/tst.pileup /pub34/stevep11/Datafiles/Bats/myoluc/pops/bams/Sample_3-PA_POST.L005.rmdup.bam

#output file will be very large, this pipes to gzip
samtools mpileup -b $workdir/pops/mpileup/pop_bamlist.txt \
-f $ref -l ${workdir}/indvs/vcfs/tst.sites |gzip > $workdir/pops/mpileup/tst.pileup.gz

#merge lanes to a sample first. Too slow
grep 'Sample_1' $workdir/pops/mpileup/pop_bamlist.txt > $workdir/pops/mpileup/Sample_1list.txt
#samtools merge -b $workdir/pops/mpileup/Sample_1list.txt /scratch/stevep11bam/Sample1.bam

#just mpileup sample 1 and process with awk. Assume 5 files per sample
samtools mpileup -b $workdir/pops/mpileup/Sample_1list.txt \
-f $ref -l ${workdir}/indvs/vcfs/tst.sites | \
awk 'BEGIN {FS="\t";OFS="\t"} {print $1,$2,$3,$4+$7+$10+$13+$16,$5$8$11$14$17,$6$9$12$15$18}'|gzip > sample1.pileup.gz &


#cat your.trimmed.vcf | awk 'BEGIN {OFS="\t"} {print $1,$2}' > positions.txt

#Thomas files still have too many SNPs
zcat ~tmlill/Datafiles/Bats/myoluc/indvs/vcfs/combined.indv.v1.vcf.gz |grep -v "#" |\
awk 'BEGIN {OFS="\t";FS="\t"} $6>100 && length($4)==1 && length($5)==1 {sub(/.*AF=/,"",$8);sub(/;.*/,"",$8);if($8<1){print $0}}' >$workdir/indvs/vcfs/trim100.lt1.txt 

#create set of sites from individual data
zcat ~tmlill/Datafiles/Bats/myoluc/indvs/vcfs/combined.indv.v1.vcf.gz | grep -v "#"|\
awk 'BEGIN {OFS="\t";FS="\t"} $6>100 && length($4)==1 && length($5)==1 {
  AF=$8;DP=$8;sub(/.*AF=/,"",AF);sub(/;.*/,"",AF);sub(/.*DP=/,"",DP);sub(/;.*/,"",DP);
  AF=AF+0;DP=DP+0;if(AF<1 && DP>29 && DP<101) {print $1,$2}}' > $workdir/indvs/vcfs/trim100.lt1.txt
#score > 100; AF<1; 30>=DP>=100; simple SNPs

ref=~tmlill/Datafiles/Bats/myoluc/ref/myoluc2.0.fa
samtools mpileup -b ~tmlill/Datafiles/Bats/myoluc/pops/mpileup/Sample_1list.txt \
-f $ref -l $workdir/indvs/vcfs/trim100.lt1.txt | \
awk 'BEGIN {FS="\t";OFS="\t"} {print $1,$2,$3,$4+$7+$10+$13+$16,$5$8$11$14$17,$6$9$12$15$18}'|gzip > $workdir/pops/mpileup/sample1.pileup.gz &

samtools mpileup -b ~tmlill/Datafiles/Bats/myoluc/pops/mpileup/Sample_2list.txt \
-f $ref -l $workdir/indvs/vcfs/trim100.lt1.txt | \
awk 'BEGIN {FS="\t";OFS="\t"} {print $1,$2,$3,$4+$7+$10+$13+$16,$5$8$11$14$17,$6$9$12$15$18}'|gzip > $workdir/pops/mpileup/sample2.pileup.gz &

samtools mpileup -b ~tmlill/Datafiles/Bats/myoluc/pops/mpileup/Sample_3list.txt \
-f $ref -l $workdir/indvs/vcfs/trim100.lt1.txt | \
awk 'BEGIN {FS="\t";OFS="\t"} {print $1,$2,$3,$4+$7+$10+$13+$16,$5$8$11$14$17,$6$9$12$15$18}'|gzip > $workdir/pops/mpileup/sample3.pileup.gz &

samtools mpileup -b ~tmlill/Datafiles/Bats/myoluc/pops/mpileup/Sample_4list.txt \
-f $ref -l $workdir/indvs/vcfs/trim100.lt1.txt | \
awk 'BEGIN {FS="\t";OFS="\t"} {print $1,$2,$3,$4+$7+$10+$13+$16,$5$8$11$14$17,$6$9$12$15$18}'|gzip > $workdir/pops/mpileup/sample4.pileup.gz &

samtools mpileup -b ~tmlill/Datafiles/Bats/myoluc/pops/mpileup/Sample_5list.txt \
-f $ref -l $workdir/indvs/vcfs/trim100.lt1.txt | \
awk 'BEGIN {FS="\t";OFS="\t"} {print $1,$2,$3,$4+$7+$10+$13+$16,$5$8$11$14$17,$6$9$12$15$18}'|gzip > $workdir/pops/mpileup/sample5.pileup.gz &

#get counts for each
seq 1 5 | while read samp
do
    zcat $workdir/pops/mpileup/sample${samp}.pileup.gz| awk 'BEGIN {OFS="\t";FS="\t"} {
    c5=$5;
    ref=gsub(/[.,]/,"",c5);
    alt=gsub(/[atcgATCG]/,"",c5);
    print $1,$2,$3,$4,$5,ref,alt}' | gzip > $workdir/pops/counts/sample${samp}.counts.gz
done

#Thomas has run 6th population
zcat ~tmlill/Datafiles/Bats/myoluc/pops/mpileup/sample6.pileup.gz| awk 'BEGIN {OFS="\t";FS="\t"} {
    c5=$5;
    ref=gsub(/[.,]/,"",c5);
    alt=gsub(/[atcgATCG]/,"",c5);
    print $1,$2,$3,$4,$5,ref,alt}' | pigz -p 4 > $workdir/pops/counts/sample6.counts.gz