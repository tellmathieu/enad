import os, sys, glob, pandas

############# --- INPUT VARIABLES --- #############
### --- Change file paths here to your data --- ###
##### --- Full filepath starting with "/"" --- ####
### --- No slashes at the end of a filepath --- ###

# from the experiment
fastqDir="/home/tmathieu/eNAD/fastq"

# from ncbi : https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/863/925/GCF_002863925.1_EquCab3.0
reference="/home/tmathieu/eNAD/resources/GCF_002863925.1_EquCab3.0_genomic.fna"
annotation="/home/tmathieu/eNAD/resources/GCF_002863925.1_EquCab3.0_genomic.gff"
refAssemblyReport="/home/tmathieu/eNAD/resources/GCF_002863925.1_EquCab3.0_assembly_report.txt"
threads = 20
pointsOfInterest="/home/tmathieu/eNAD/resources/points_of_interest_historical.bed"
tracksPOI = "/home/tmathieu/eNAD/resources/tracksPOI.txt"

#order of comparisons - based on coancestry coefficient and control or treat category
order_comps = ['horse_05_608', 'horse_04_56', 'horse_06_1212', 'horse_03_238', 'horse_01_200', 'horse_02_205']
PHASES = ['hap1','hap2']

###################### --- PIPELINE--- #######################
### --- Don't change unless you know what you're doing --- ###

### Using this tutorial https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9254108/
### And more

mainDir = os.getcwd()
outputGenomeDir = "outputGenomes"
controlprimary = os.path.join(mainDir, order_comps[0], order_comps[0] + '.scaffold_chr_only.fa')
controlphase1 = os.path.join(mainDir, order_comps[0] + '.hap1', order_comps[0] + '.hap1_scaffold_chr_only.fa')
controlphase2 = os.path.join(mainDir, order_comps[0] + '.hap2', order_comps[0] + '.hap2_scaffold_chr_only.fa')

SAMPLES, = glob_wildcards(os.path.join(fastqDir, '{sample}.fastq'))
PHASES = ['hap1','hap2']

hap = ['phasedSVs.compto5-1','phasedSVs.compto5-2']

f = open(os.path.join(mainDir, 'resources', 'chromosome.name.list.txt'))
chrList = f.readlines()
for i in range(0,len(chrList)):
  chrList[i] = chrList[i].strip()
f.close()

first_comp = order_comps[0:-1]
second_comp = order_comps[1:]

#print(sniffles_comp)
print(first_comp)
print(second_comp)
print(order_comps[0])
print(SAMPLES)
print(chrList)

rule all:
  input:
    os.path.join(mainDir, 'snifflesPrimary', 'candidates_filtered.bed'),
    expand(os.path.join(mainDir, 'plotsrGraphsWithSVs', 'syri_{chr}_combined_output_plot.png'), chr = chrList),
    expand(os.path.join(mainDir, '{sample}', '{sample}.scaffold_chr_only.fa'), sample = order_comps),
    expand(os.path.join(mainDir, 'sniffles{hap}', 'plotDone.txt'), hap = hap),
    expand(os.path.join(mainDir, 'sniffles', '{hap}_filtered.bed'), hap = hap)
    #os.path.join(mainDir, 'coverageStats.txt')


rule alignToFirstControl:
  resources:
    mem_mb=50000
  input:
    control1 = controlphase1,
    control2 = controlphase2,
    fastq = os.path.join(mainDir, 'fastq', '{sample}.fastq')
  params:
    snifflesDir = os.path.join(mainDir, 'sniffles')
  output:
    hap1 = os.path.join(mainDir, 'sniffles', '{sample}.compto5-1.sam'),
    hap2 = os.path.join(mainDir, 'sniffles', '{sample}.compto5-2.sam')
  shell: '''
    
  '''

rule convertSamToBam:
  resources:
    mem_mb=50000
  input:
     os.path.join(mainDir, 'sniffles', '{sample}.sam')
  output:
    os.path.join(mainDir, 'sniffles', '{sample}.bam')
  shell: '''
    samtools view -@ 20 -S -b {input} > {output}
    rm {input}
  '''

rule sortBam:
  resources:
    mem_mb=50000
  input:
    os.path.join(mainDir, 'sniffles', '{sample}.bam')
  output:
    os.path.join(mainDir, 'sniffles', '{sample}_sorted.bam')
  shell: '''
    samtools sort -@ 20 -o {output} {input}
    rm {input}
  '''

rule indexBam:
  resources:
    mem_mb=50000
  input:
    os.path.join(mainDir, 'sniffles', '{sample}_sorted.bam')
  output:
    os.path.join(mainDir, 'sniffles', '{sample}_sorted.bam.bai')
  shell: '''
    samtools index -@ 20 -b {input}
  '''

rule findSVshap1:
  resources:
    mem_mb=50000
  input:
    compareToControlBam = os.path.join(mainDir, 'sniffles', '{sample}.compto5-1_sorted.bam'),
    bai = os.path.join(mainDir, 'sniffles', '{sample}.compto5-1_sorted.bam.bai'),
    refChr = controlphase1
  output:
    snf = os.path.join(mainDir, 'sniffles', '{sample}.compto5-1.snf'),
    vcf = os.path.join(mainDir, 'sniffles', '{sample}.compto5-1.vcf')
  shell: '''
    sniffles \
      --input {input.compareToControlBam} \
      --snf {output.snf} \
      --vcf {output.vcf} \
      --reference {input.refChr} \
      --threads 16
  '''

rule findSVshap2:
  resources:
    mem_mb=50000
  input:
    compareToControlBam = os.path.join(mainDir, 'sniffles', '{sample}.compto5-2_sorted.bam'),
    bai = os.path.join(mainDir, 'sniffles', '{sample}.compto5-2_sorted.bam.bai'),
    refChr = controlphase2
  output:
    snf = os.path.join(mainDir, 'sniffles', '{sample}.compto5-2.snf'),
    vcf = os.path.join(mainDir, 'sniffles', '{sample}.compto5-2.vcf')
  shell: '''
    sniffles \
      --input {input.compareToControlBam} \
      --snf {output.snf} \
      --vcf {output.vcf} \
      --reference {input.refChr} \
      --threads 16
  '''

rule compareSVshap1:
  resources:
    mem_mb=50000
  input:
    snf = expand(os.path.join(mainDir, 'sniffles', '{sample}.compto5-1.snf'),sample = order_comps),
    refChr = controlphase1
  output:
    multiVCF = os.path.join(mainDir, 'sniffles', 'phasedSVs.compto5-1.vcf')
  shell: '''
    sniffles \
      --input {input.snf} \
      --vcf {output.multiVCF} \
      --reference {input.refChr} \
      --threads 16
  '''

rule compareSVshap2:
  resources:
    mem_mb=50000
  input:
    snf = expand(os.path.join(mainDir, 'sniffles', '{sample}.compto5-2.snf'),sample = order_comps),
    refChr = controlphase2
  output:
    multiVCF = os.path.join(mainDir, 'sniffles', 'phasedSVs.compto5-2.vcf')
  shell: '''
    sniffles \
      --input {input.snf} \
      --vcf {output.multiVCF} \
      --reference {input.refChr} \
      --threads 16
  '''

rule cleanAndFilterVCF:
  input:
    multiVCF = os.path.join(mainDir, 'sniffles', '{hap}.vcf')
  output:
    filteredVCF = os.path.join(mainDir, 'sniffles', '{hap}_filtered.vcf')
  shell: '''
    sed "s/GT:GQ:DR:DV:PS:ID/GT:GQ:DR:DV:ID/g" {input.multiVCF} > temp_cleaned.vcf
    bgzip temp_cleaned.vcf -k
    bcftools index temp_cleaned.vcf.gz
    bcftools view -i '(GT[0]=="0/0" && GT[1]=="0/0") && ((GT[2]=="0/1" || GT[2]=="1/0" || GT[2]=="1/1") && (GT[3]=="0/1" || GT[3]=="1/0" || GT[3]=="1/1") && (GT[4]=="0/1" || GT[4]=="1/0" || GT[4]=="1/1"))' temp_cleaned.vcf.gz -Ov -o {output.filteredVCF}
    rm temp_cleaned*
  '''

rule convert2Bed:
  input:
    filteredVCF = os.path.join(mainDir, 'sniffles', '{hap}_filtered.vcf')
  output:
    filteredBED = os.path.join(mainDir, 'sniffles', '{hap}_filtered.bed')
  shell: '''
    gvcf2bed -I {input} -O {output}
  '''

rule plotSVs:
  input:
    multiVCF = os.path.join(mainDir, 'sniffles', '{hap}.vcf')
  params:
    plotDir = os.path.join(mainDir, 'sniffles{hap}')
  output:
    os.path.join(mainDir, 'sniffles{hap}', 'plotDone.txt')
  shell: '''
    mkdir -p {params.plotDir}

    python3 -m sniffles2_plot -i {input.multiVCF} -o {params.plotDir}

    echo "Plots made" > {output}
  '''
