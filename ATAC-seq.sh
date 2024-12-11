gtf='/mnt/sda/jwz/genome/mouse/mm39/gencode.vM33.annotation.gtf'
mm39='/mnt/sda/jwz/genome/mouse/mm39/GRCm39'

for i in `cat sample.txt`
do 
cutadapt -m 5 -e 0.2 -o cleandata/${i}_R1.fq.gz -p cleandata/${i}_R2.fq.gz ./rawdata/${i}_combined_R1.fastq.gz ./rawdata/${i}_combined_R2.fastq.gz
bowtie2 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --dovetail -I 10 -X 1000 -x ${mm39} -1 cleandata/${i}_R1.fq.gz -2 cleandata/${i}_R2.fq.gz -S ./bowtie2/${i}.sam 2> ./bowtie2/${i}.log
samtools view -bS ./bowtie2/${i}.sam | samtools sort - >./bowtie2/${i}.sorted.bam
samtools index ./bowtie2/${i}.sorted.bam
java -jar /mnt/sda/jwz/software/picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT I=./bowtie2/${i}.sorted.bam O=./bowtie2/${i}.markdup.bam M=./bowtie2/${i}.markdup_metrics.txt
samtools view -h -f 2 -q 30 ./bowtie2/${i}.markdup.bam -o ./bowtie2/${i}.last.bam
samtools index ./bowtie2/${i}.last.bam
samtools flagstat ./bowtie2/${i}.last.bam >./bowtie2/${i}.last.flastat
macs3 callpeak -f BAMPE -t ./bowtie2/${i}.last.bam -g mm -n ${i} --outdir ./macs3 -B -q 0.01
bamCoverage -b ./bowtie2/${i}.last.bam -o ./bw/${i}.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2725521370
done

multiBamSummary bins --bamfiles CON-3.last.bam CON-4.last.bam CON-5.last.bam ISO-4.last.bam ISO-5.last.bam ISO-6.last.bam --labels CON-3 CON-4 CON-5 ISO-4 ISO-5 ISO-6 -out readCounts.npz --outRawCounts readCounts.tab
plotPCA -in readCounts.npz -o ATAC-PCA.pdf
plotCorrelation -in readCounts.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of ATAC-seq" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o ATAC-spearman.pdf --outFileCorMatrix SpearmanCorr_readCounts.tab



##merge replicates
samtools merge ./bowtie2/CON.bam ./bowtie2/CON-3.last.bam ./bowtie2/CON-4.last.bam ./bowtie2/CON-5.last.bam
samtools index ./bowtie2/CON.bam
macs3 callpeak -f BAMPE -t ./bowtie2/CON.bam -g mm -n CON --outdir ./macs3 -B -q 0.01
samtools merge ./bowtie2/ISO.bam ./bowtie2/ISO-3.last.bam ./bowtie2/ISO-4.last.bam ./bowtie2/ISO-5.last.bam
samtools index ./bowtie2/ISO.bam
macs3 callpeak -f BAMPE -t ./bowtie2/ISO.bam -g mm -n ISO --outdir ./macs3 -B -q 0.01
bamCoverage -b ./bowtie2/CON.bam -o ./bw/CON.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2725521370
bamCoverage -b ./bowtie2/ISO.bam -o ./bw/ISO.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2725521370


##TSS profile
computeMatrix reference-point --referencePoint TSS -b 2000 -a 2000 -R ~/genome/mouse/mm39.bed -S CON.bw ISO.bw --skipZeros -o TSS.gz --outFileSortedRegions TSSsortedRegions.bed
plotHeatmap -m TSS.gz -out TSS.Heatmap.pdf --plotFileFormat pdf
plotProfile -m TSS.gz -out TSS.Profile.pdf --plotFileFormat pdf
computeMatrix scale-regions -b 2000 -a 2000 -R ~/genome/mouse/mm39.bed -S CON.bw ISO.bw --skipZeros -o body.gz --outFileSortedRegions genebodysortedRegions.bed
plotHeatmap -m body.gz -out body.Heatmap.pdf --plotFileFormat pdf
plotProfile -m body.gz -out body.Profile.pdf --plotFileFormat pdf

