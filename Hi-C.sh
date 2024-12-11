##Hi-C and HiCUT analysis

##消化参考基因组
digest_genome.py -r ^GATC C^TNAG -o MboI-DdeI_resfrag_GRCm39.bed ~/genome/mouse/GRCm39.fa

#Hicpro
for i in `cat sample.txt`
do 
nohup HiC-Pro -i fastq/${i} -o ${i} -c config-hicpro.txt &
done

##FitHiChIP calling loops for HiCUT data
for i in `cat sample.txt`; do nohup PeakInferHiChIP.sh -H ./${i}/hic_results/data/${i} -D ../fithichip/${i} -R mm -L 150 & done
bash FitHiChIP_HiCPro.sh -C CON.config
bash FitHiChIP_HiCPro.sh -C ISO.config

##HiC A/B compartment and TAD
#Format conversion
for i in `cat sample.txt`
do
hicpro2juicebox.sh -i ${i}/data/${i}/${i}.allValidPairs -g /mnt/sda/jwz/software/HiC-Pro/HiC-Pro_3.1.0/annotation/GRCm39.size -j /mnt/sda/jwz/software/juicer_tools.1.7.6_jcuda.0.8.jar -r /mnt/sda/jwz/software/HiC-Pro/HiC-Pro_3.1.0/annotation/MboI-DdeI_resfrag_GRCm39.bed -o ./
hicConvertFormat -m ${i}.allValidPairs.hic --inputFormat hic --outputFormat cool -o ${i}.cool --resolutions 25000
hicConvertFormat -m ${i}_25000.cool --inputFormat cool --outputFormat h5 -o ${i}-10k.h5 --resolutions 25000
done

##homer分析compartment
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' CON.allValidPairs >homer/CON.allValidPairs.homer
makeTagDirectory ./CON -format HiCsummary CON.allValidPairs.homer
runHiCpca.pl CON-25kb ./CON -cpu 20 -res 25000  -genome mm39 -pc 1
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' ISO.allValidPairs >homer/ISO.allValidPairs.homer
makeTagDirectory ./ISO -format HiCsummary ISO.allValidPairs.homer
runHiCpca.pl ISO-25kb ./ISO -cpu 20 -res 25000  -genome mm39 -pc 1

##HiCExplorer分析TAD
hicFindTADs -m CON-25kb.h5 --outPrefix CON-25kb --thresholdComparisons 0.05 --delta 0.01 --correctForMultipleTesting fdr -p 64
hicFindTADs -m ISO-25kb.h5 --outPrefix ISO-25kb --thresholdComparisons 0.05 --delta 0.01 --correctForMultipleTesting fdr -p 64

