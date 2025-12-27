import os, sys, glob, pandas

############# --- INPUT VARIABLES --- #############
### --- Change file paths here to your data --- ###
##### --- Full filepath starting with "/"" --- ####
### --- No slashes at the end of a filepath --- ###


#order of comparisons - based on coancestry coefficient and control or treat category
order_comps = ['horse_05_608', 'horse_04_56', 'horse_06_1212', 'horse_03_238', 'horse_01_200', 'horse_02_205']

###################### --- PIPELINE--- #######################
### --- Don't change unless you know what you're doing --- ###

### Using this tutorial https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9254108/
### And more

mainDir = os.getcwd()
first_comp = order_comps[0:-1]

rule all:
  input:
    expand(os.path.join(mainDir, 'ref-results', '{sample}.candidates_expanded_unfiltered.bam.bai'), sample = order_comps)

rule samtoolsView:
  resources:
    mem_mb=50000
  input:
    bam = os.path.join(mainDir, 'refAlign', '{sample}.sorted.bam'),
    bed = os.path.join(mainDir, 'ref-results', 'EquCab3.0_candidates_expanded_unfiltered.bed')
  output:
    candidates = os.path.join(mainDir, 'ref-results', '{sample}.candidates_expanded_unfiltered.bam')
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
    os.path.join(mainDir, 'ref-results', '{sample}.candidates_expanded_unfiltered.bam')
  output:
    os.path.join(mainDir, 'ref-results', '{sample}.candidates_expanded_unfiltered.bam.bai')
  shell: '''
    samtools index {input}
  '''
