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

SAMPLES, reads,  = glob_wildcards(os.path.join(mainDir, 'shortReads', '{sample}_{reads}_001.fastq.gz'))

SAMPLES = sorted(set(SAMPLES))

#SAMPLE = ['horse_05_608', 'horse_04_56', 'horse_06_1212', 'horse_03_238', 'horse_01_200', 'horse_02_205']

rule all:
  input:
    #os.path.join(mainDir, 'ref-results', refName + '_candidates_expanded_unfiltered_short_read.bed'),
    os.path.join(mainDir, 'qc', 'multiqc_done_short_read.txt'),
    os.path.join(mainDir, 'ref-results', 'reference_candidates_vcf_gff_intersect_short_read.txt'),
    os.path.join(mainDir, 'refAlign', 'reference_candidates_short_read.vcf.gz.tbi')

rule alignToFirstControl:
  resources:
    mem_mb=30000
  input:
    ref = ref,
    fastq1 = os.path.join(mainDir, 'shortReads', '{sample}_R1_001.fastq.gz'),
    fastq2 = os.path.join(mainDir, 'shortReads', '{sample}_R2_001.fastq.gz')
  params:
    refDir = os.path.join(mainDir, 'refAlign')
  output:
    sam = os.path.join(mainDir, 'refAlign', '{sample}.sam')
  shell: '''
    mkdir -p {params.refDir}

    minimap2 \
      -ax sr \
      -t 20 \
      {input.ref} \
      {input.fastq1} {input.fastq2} \
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
    zipqc = expand(os.path.join(mainDir, 'qc', '{sample}_fastqc.zip'), sample = SAMPLES),
    html = expand(os.path.join(mainDir, 'qc', '{sample}_fastqc.html'), sample = SAMPLES)
  params:
    qcDir = os.path.join(mainDir, 'qc')
  output:
    report = os.path.join(mainDir, 'qc', 'multiqc_done_short_read.txt')
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

rule indexBamMaybeAddReadGroups:
  resources:
    mem_mb=30000
  threads: threads
  input:
    bam = os.path.join(mainDir, 'refAlign', '{sample}_sorted.bam')
  output:
    rgbam = os.path.join(mainDir, 'refAlign', '{sample}_rg.bam'),
    rgbai = os.path.join(mainDir, 'refAlign', '{sample}_rg.bam.bai')
  run:
    import subprocess

    # Check for @RG in the header
    has_rg = subprocess.run(
      f"samtools view -H {input.bam} | grep -q '^@RG'",
      shell=True
    ).returncode == 0

    if has_rg:
      # if the check for read groups yielded a positive result
      # I just copy the file and do the indexing step
      shell("""
        cp {input.bam} {output.rgbam}
        samtools index -@ {threads} {output.rgbam}
      """)
    else:
      # Extract sample name, first part before first underscore
      sample_id = wildcards.sample.split("_")[0]

      # Make RG tags
      rgid = wildcards.sample
      rglb = f"lib_{sample_id}"
      rgpu = wildcards.sample

      shell("""
        gatk AddOrReplaceReadGroups \
          -I {input.bam} \
          -O {output.rgbam} \
          --RGID {rgid} \
          --RGLB {rglb} \
          --RGPL ILLUMINA \
          --RGPU {rgpu} \
          --RGSM {sample_id}

        samtools index -@ {threads} {output.rgbam}
      """)


# used this tutorial for direction in this section 
# https://omicstutorials.com/snp-calling-a-step-by-step-guide/#Step_3_Mark_Duplicates_Optional
rule markDups:
  conda: 'envs/markDups.yaml'
  resources:
    mem_mb=30000
  input:
    bam = os.path.join(mainDir, 'refAlign', '{sample}_rg.bam')
  output:
    dedup = os.path.join(mainDir, 'refAlign', '{sample}_rg_dedup.bam'),
    metrics = os.path.join(mainDir, 'refAlign', '{sample}_metrics_dedup.txt')
  shell: '''
    gatk MarkDuplicates -I {input.bam} -O {output.dedup} -M {output.metrics}
  '''

rule reIndex:
  resources:
    mem_mb=30000
  input:
    bam = os.path.join(mainDir, 'refAlign', '{sample}_rg_dedup.bam')
  output:
    bai = os.path.join(mainDir, 'refAlign', '{sample}_rg_dedup.bam.bai')
  shell: '''
    samtools index -@ {threads} {input.bam}
  '''

rule makeDict:
  input:
    refChr = ref
  output:
    dict1 = os.path.join(mainDir,'resources','reference_genome_chr_only.dict')
  shell: '''
    gatk CreateSequenceDictionary -R {input.refChr}
  '''

rule findSNPs:
  resources:
    mem_mb=30000
  input:
    dict1 = os.path.join(mainDir,'resources','reference_genome_chr_only.dict'),
    bam = os.path.join(mainDir, 'refAlign', '{sample}_rg_dedup.bam'),
    bai = os.path.join(mainDir, 'refAlign', '{sample}_rg_dedup.bam.bai'),
    refChr = ref
  output:
    vcf = os.path.join(mainDir, 'refAlign', '{sample}_raw.vcf')
  shell: '''
    gatk HaplotypeCaller \
      -R {input.refChr} \
      -I {input.bam} \
      -O {output.vcf}
  '''

rule filterVCFQuality:
  resources:
    mem_mb=50000
  input:
    vcf = os.path.join(mainDir, 'refAlign', '{sample}_raw.vcf'),
    refChr = ref
  output:
    filterVcf = os.path.join(mainDir, 'refAlign', '{sample}_filtered.vcf')
  shell: '''
    gatk VariantFiltration \
      -R {input.refChr} \
      -V {input.vcf} \
      --filter-expression "QUAL < 30.0 || DP < 10" \
      --filter-name "LowQual" \
      -O {output.filterVcf}

  '''

rule compressVCFs:
  input:
    filterVcf = os.path.join(mainDir, 'refAlign', '{sample}_filtered.vcf')
  output:
    gzVcf = os.path.join(mainDir, 'refAlign', '{sample}_filtered.vcf.gz')
  shell: '''
    bgzip {input.filterVcf}
  '''

rule indexVCFs:
  input:
    gzVcf = os.path.join(mainDir, 'refAlign', '{sample}_filtered.vcf.gz')
  output:
    index = os.path.join(mainDir, 'refAlign', '{sample}_filtered.vcf.gz.tbi')
  shell: '''
    tabix -p vcf {input.gzVcf}
  '''

rule combineVCFs:
  resources:
    mem_mb=50000
  input:
    expand(os.path.join(mainDir, 'refAlign', '{sample}_filtered.vcf.gz'), sample = SAMPLES)
  output:
    multiVCF = os.path.join(mainDir, 'refAlign', 'reference_candidates_short_read.vcf')
  shell: '''
    echo {input} | awk '{{for (i=1;i<=NF;i++) print $i}}' | sort -u > short_vcf_list.txt

    bcftools merge \
      -o {output.multiVCF} \
      $(cat short_vcf_list.txt)
  '''

rule compressVCFmerged:
  input:
    multiVCF = os.path.join(mainDir, 'refAlign', 'reference_candidates_short_read.vcf')
  output:
    multiVCFgz = os.path.join(mainDir, 'refAlign', 'reference_candidates_short_read.vcf.gz')
  shell: '''
    bgzip {input.multiVCF}
  '''

rule indexVCFmerged:
  input:
    multiVCFgz = os.path.join(mainDir, 'refAlign', 'reference_candidates_short_read.vcf.gz')
  output:
    index = os.path.join(mainDir, 'refAlign', 'reference_candidates_short_read.vcf.gz.tbi')
  shell: '''
    tabix -p vcf {input.multiVCFgz}
  '''

