gtf='/mnt/sda/jwz/genome/mouse/mm39/gencode.vM33.annotation.gtf'
mm39='/mnt/sda/jwz/genome/mouse/mm39/GRCm39'

##QC
for i in `cat sample.txt`
do
fastqc rawdata/${i}_1.fq.gz -o rawdata/fastqc
done
multiqc rawdata/fastqc -o rawdata/fastqc

##Trim and Mapping
for i in `cat sample.txt`
do 
trim_galore --length 35 --stringency 3 --paired rawdata/${i}_1.fq.gz rawdata/${i}_2.fq.gz -o ./cleandata/ 
hisat2 -p 4 -x ${mm39} -1 cleandata/${i}_1_val_1.fq.gz -2 cleandata/${i}_2_val_2.fq.gz -S ./hisat2/${i}.sam 2> ./hisat2/${i}.log
samtools view -bS hisat2/${i}.sam | samtools sort - > hisat2/${i}.sorted.bam 
samtools index hisat2/${i}.sorted.bam
bamCoverage -b hisat2/${i}.sorted.bam -o bw/${i}.bw --normalizeUsing RPKM
stringtie hisat2/${i}.sorted.bam -p 4 -G ${gtf} -o ./stringtie/${i}/${i}.gtf -A ./stringtie/${i}/${i}.abund -e -B
done

###Count
featureCounts -T 10 -p -a ${gtf} -o read.count -t exon -g gene_id CON-3.sorted.bam CON-4.sorted.bam CON-5.sorted.bam ISO-4.sorted.bam ISO-5.sorted.bam ISO-6.sorted.bam



