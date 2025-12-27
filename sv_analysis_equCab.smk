import os, sys, glob, pandas

############# --- INPUT VARIABLES --- #############
### --- Change file paths here to your data --- ###
##### --- Full filepath starting with "/"" --- ####
### --- No slashes at the end of a filepath --- ###

# from ncbi : https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/863/925/GCF_002863925.1_EquCab3.0
reference="/home/tmathieu/eNAD/resources/GCF_002863925.1_EquCab3.0_genomic.fna"
annotation="/home/tmathieu/eNAD/resources/GCF_002863925.1_EquCab3.0_genomic.gff"
refAssemblyReport="/home/tmathieu/eNAD/resources/GCF_002863925.1_EquCab3.0_assembly_report.txt"
threads = 20
ref="/home/tmathieu/eNAD/resources/reference_genome_chr_only.fa"
gff="/home/tmathieu/eNAD/resources/EquCab3.0_annotation_chr_only.gff"
refName="EquCab3.0"

###################### --- PIPELINE--- #######################
### --- Don't change unless you know what you're doing --- ###


mainDir = os.getcwd()

SAMPLES, = glob_wildcards(os.path.join(mainDir, 'fastq', '{sample}.fastq'))

SAMPLE = ['horse_05_608', 'horse_04_56', 'horse_06_1212', 'horse_03_238', 'horse_01_200', 'horse_02_205']

rule all:
  input:
    os.path.join(mainDir, 'ref-results', refName + '_candidates_expanded_unfiltered.bed'),
    os.path.join(mainDir, 'qc', 'multiqc_done.txt')

rule alignToFirstControl:
  resources:
    mem_mb=30000
  input:
    ref = ref,
    fastq = os.path.join(mainDir, 'fastq', '{sample}.fastq')
  params:
    refDir = os.path.join(mainDir, 'refAlign')
  output:
    sam = os.path.join(mainDir, 'refAlign', '{sample}.sam')
  shell: '''
    mkdir -p {params.refDir}

    minimap2 \
      -a \
      -x map-hifi \
      --eqx \
      -t 20 \
      {input.ref} \
      {input.fastq} \
      > {output.sam}
  '''

rule fastQC:
  input:
    fastq = os.path.join(mainDir, 'fastq', '{sample}.fastq')
  params:
    qcDir = os.path.join(mainDir, 'qc')
  output:
    zipqc = os.path.join(mainDir, 'qc', '{sample}_fastqc.zip'),
    html = os.path.join(mainDir, 'qc', '{sample}_fastqc.html')
  shell: '''
    mkdir -p {params.qcDir}

    fastqc -o {params.qcDir} -f fastq {input.fastq}
  '''

rule multiQC:
  input:
    zipqc = expand(os.path.join(mainDir, 'qc', '{sample}_fastqc.zip'), sample = SAMPLE),
    html = expand(os.path.join(mainDir, 'qc', '{sample}_fastqc.html'), sample = SAMPLE)
  params:
    qcDir = os.path.join(mainDir, 'qc')
  output:
    report = os.path.join(mainDir, 'qc', 'multiqc_done.txt')
  shell: '''
    multiqc {params.qcDir}/

    echo multiQC Done! > {output.report}
  '''

rule convertSamToBam:
  resources:
    mem_mb=30000
  input:
     os.path.join(mainDir, 'refAlign', '{sample}.sam')
  output:
    os.path.join(mainDir, 'refAlign', '{sample}.bam')
  shell: '''
    samtools view -@ 20 -S -b {input} > {output}
    rm {input}
  '''

rule sortBam:
  resources:
    mem_mb=30000
  input:
    os.path.join(mainDir, 'refAlign', '{sample}.bam')
  output:
    os.path.join(mainDir, 'refAlign', '{sample}_sorted.bam')
  shell: '''
    samtools sort -@ 20 -o {output} {input}
    rm {input}
  '''

rule indexBam:
  resources:
    mem_mb=30000
  input:
    os.path.join(mainDir, 'refAlign', '{sample}_sorted.bam')
  output:
    os.path.join(mainDir, 'refAlign', '{sample}_sorted.bam.bai')
  shell: '''
    samtools index -@ 20 -b {input}
  '''

rule findSVs:
  resources:
    mem_mb=30000
  input:
    compareToControlBam = os.path.join(mainDir, 'refAlign', '{sample}_sorted.bam'),
    bai = os.path.join(mainDir, 'refAlign', '{sample}_sorted.bam.bai'),
    refChr = ref
  output:
    snf = os.path.join(mainDir, 'refAlign', '{sample}.snf'),
    vcf = os.path.join(mainDir, 'refAlign', '{sample}.vcf')
  shell: '''
    sniffles \
      --input {input.compareToControlBam} \
      --snf {output.snf} \
      --vcf {output.vcf} \
      --reference {input.refChr} \
      --threads 16
  '''

rule compareSVs:
  resources:
    mem_mb=50000
  input:
    snf = expand(os.path.join(mainDir, 'refAlign', '{sample}.snf'),sample = SAMPLE),
    refChr = ref
  output:
    multiVCF = os.path.join(mainDir, 'refAlign', 'reference_candidates.vcf')
  shell: '''
    sniffles \
      --input {input.snf} \
      --vcf {output.multiVCF} \
      --reference {input.refChr} \
      --threads 16
  '''

rule cleanAndFilterVCF:
  input:
    multiVCF = os.path.join(mainDir, 'refAlign', 'reference_candidates.vcf')
  params:
    name = refName
  output:
    filteredVCF = os.path.join(mainDir, 'refAlign', 'reference_candidates_filtered.vcf')
  shell: '''
    sed "s/GT:GQ:DR:DV:PS:ID/GT:GQ:DR:DV:ID/g" {input.multiVCF} > temp_cleaned_{params.name}.vcf
    bgzip temp_cleaned_{params.name}.vcf -k
    bcftools index temp_cleaned_{params.name}.vcf.gz
     bcftools view -i '(GT[0]=="0/0" && GT[1]=="0/0" && GT[2]=="0/0") && ((GT[3]=="0/1" || GT[3]=="1/0" || GT[3]=="1/1") && (GT[4]=="0/1" || GT[4]=="1/0" || GT[4]=="1/1") && (GT[5]=="0/1" || GT[5]=="1/0" || GT[5]=="1/1"))' temp_cleaned_{params.name}.vcf.gz -Ov -o {output.filteredVCF}
    rm temp_cleaned_{params.name}*

    # order of the horses horse_03_238  horse_04_56 horse_01_200  horse_06_1212 horse_05_608  horse_02_205
  '''

rule convert2Bed:
  input:
    filteredVCF = os.path.join(mainDir, 'refAlign', 'reference_candidates_filtered.vcf')
  output:
    filteredBED = os.path.join(mainDir, 'refAlign', 'reference_candidates_filtered.bed')
  shell: '''
    gvcf2bed -I {input} -O {output}
  '''

rule plotSVs:
  input:
    multiVCF = os.path.join(mainDir, 'refAlign', 'reference_candidates.vcf')
  params:
    plotDir = os.path.join(mainDir, 'refAlign')
  output:
    os.path.join(mainDir, 'refAlign', 'plotDone.txt')
  shell: '''
    mkdir -p {params.plotDir}

    python3 -m sniffles2_plot -i {input.multiVCF} -o {params.plotDir}

    echo "Plots made" > {output}
  '''

rule intersectGFFandVCF:
  input:
    filteredVCF = os.path.join(mainDir, 'refAlign', 'reference_candidates_filtered.vcf'),
    gff = gff
  params:
    os.path.join(mainDir, 'ref-results')
  output:
    intersect = os.path.join(mainDir, 'ref-results', 'reference_candidates_vcf_gff_intersect.txt')
  shell: '''
    mkdir -p {params}

    bedtools intersect -wo -a {input.gff} -b {input.filteredVCF} > {output.intersect}
  '''

rule runGOanalysis:
  input:
    intersect = os.path.join(mainDir, 'ref-results', 'reference_candidates_vcf_gff_intersect.txt')
  params:
    main = os.path.join(mainDir),
    prog = os.path.join(mainDir, 'scripts', 'get_go_terms_auto.R'),
    name = refName
  output:
    os.path.join(mainDir, 'ref-results', refName + '_candidates_expanded_unfiltered.bed')
  shell: '''
    Rscript {params.prog} {params.main} {input.intersect} {params.name}
  '''







