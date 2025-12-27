import os, sys, glob, pandas

############# --- INPUT VARIABLES --- #############
### --- Change file paths here to your data --- ###
##### --- Full filepath starting with "/"" --- ####
### --- No slashes at the end of a filepath --- ###


#order of comparisons - based on coancestry coefficient and control or treat category
order_comps = ['horse_05_608', 'horse_04_56', 'horse_06_1212', 'horse_03_238', 'horse_01_200', 'horse_02_205']
PHASES = ['hap1','hap2']

###################### --- PIPELINE--- #######################
### --- Don't change unless you know what you're doing --- ###

### Using this tutorial https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9254108/
### And more

mainDir = os.getcwd()

PHASES = ['1','2']

first_comp = order_comps[0:-1]

rule all:
  input:
    expand(os.path.join(mainDir, 'results', '{sample}.hap{hap}.candidates_expanded_unfiltered.bam.bai'), sample = order_comps, hap = PHASES)

rule samtoolsView:
  resources:
    mem_mb=50000
  input:
    bam = os.path.join(mainDir, 'sniffles', '{sample}.compto5-{hap}_sorted.bam'),
    bed = os.path.join(mainDir, 'results', 'hap{hap}_candidates_expanded_unfiltered.bed')
  output:
    candidates = os.path.join(mainDir, 'results', '{sample}.hap{hap}.candidates_expanded_unfiltered.bam')
  shell: '''
    samtools view -b -L \
      {input.bed} \
      {input.bam} \
      > {output.candidates}
  '''

rule samtoolsIndex:
  resources:
    mem_mb=50000
  input:
    os.path.join(mainDir, 'results', '{sample}.hap{hap}.candidates_expanded_unfiltered.bam')
  output:
    os.path.join(mainDir, 'results', '{sample}.hap{hap}.candidates_expanded_unfiltered.bam.bai')
  shell: '''
    samtools index {input}
  '''
