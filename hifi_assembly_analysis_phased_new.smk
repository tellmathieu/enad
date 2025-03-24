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
controlphase = os.path.join(mainDir, order_comps[0] + '.hap1', order_comps[0] + '.hap1_scaffold_chr_only.fa')

SAMPLES, = glob_wildcards(os.path.join(fastqDir, '{sample}.fastq'))
PHASES = ['hap1','hap2']

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
    expand(os.path.join(mainDir, 'snifflesPrimary', '{comp_sample}_statsDone.txt'), comp_sample = order_comps),
    expand(os.path.join(mainDir, 'refAlign', '{sample}_statsDone.txt'), sample = SAMPLES)
    #os.path.join(mainDir, 'coverageStats.txt')


#this is de novo assembly, without a reference, we use the unzipped fastq files here, but this does not use the quality scores
rule initial_assembly:
  resources:
    mem_mb=200000
  input:
    fastq = os.path.join(fastqDir, '{sample}.fastq')
  params:
    threads = threads,
    outputDir = os.path.join(mainDir, outputGenomeDir)
  output:
    deNovogfa = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.asm.{phases}.bp.p_ctg.gfa')
  shell: '''
    mkdir -p {params.outputDir}/{wildcards.sample}
    
    #de novo assembly
    hifiasm -o {params.outputDir}/{wildcards.sample}/{wildcards.sample}.asm -t {params.threads} {input.fastq}
  '''

rule de_novo_gfa_to_fa:
  input:
    deNovogfa = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.asm.bp.{phases}.p_ctg.gfa')
  output:
    deNovoFa = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.asm.bp.{phases}.p_ctg.fa')
  shell: '''
    #convert the output gfa file to a traditional fasta file
    awk '/^S/{{print ">"$2;print $3}}' {input.deNovogfa} > {output.deNovoFa}
  '''

rule de_novo_gfa_to_fa_primary:
  input:
    deNovogfa = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.asm.bp.p_ctg.gfa')
  output:
    deNovoFa = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.asm.bp.p_ctg.fa')
  shell: '''
    #convert the output gfa file to a traditional fasta file
    awk '/^S/{{print ">"$2;print $3}}' {input.deNovogfa} > {output.deNovoFa}
  '''

#getting qualty stats and saving to a text file
rule de_novo_faassembly_stats:
  input:
    deNovoFa = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.asm.bp.{phases}.p_ctg.fa')
  output:
    deNovoAssemblyStatsTxt = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.asm.bp.{phases}.p_ctg.fa.txt')
  shell: '''
    assembly-stats {input.deNovoFa} > {output.deNovoAssemblyStatsTxt}
  '''

rule create_coverage_table:
  input:
    deNovoFa = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.asm.bp.{phases}.p_ctg.fa')
  output:
    length = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.{phases}.length.csv')
  shell: '''
    #create csv file with headers
    echo "line,length,type,coverage" > {output.length}
    
    #subcommand
    LEN1=`bioawk -c fastx '{{sum+=length($seq)}}END{{print sum}}' {input.deNovoFa}`
    
    #adding length lines to csv we created
    cat {input.deNovoFa} | bioawk -c fastx -v line="HIFI" '{{print line","length($seq)","length ($seq)}}' | sort -k3rV -t "," | awk -F "," -v len="$LEN1" -v type="contig" 'OFS=","{{print $1,$2,type,(sum+0)/len; sum+=$3 }}' >> {output.length}
  '''

rule create_coverage_plot:
  input:
    length = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.{phases}.length.csv')
  params:
    rscript = os.path.join(mainDir, 'scripts', 'coveragePlot.R')
  output:
    plot = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.{phases}.coverage.pdf')
  shell: '''
    Rscript {params.rscript} \
      {input.length} \
      {output.plot}
  '''

rule BUSCOanalysis:
  resources:
    mem_mb=70000
  input:
    deNovoFa = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.asm.bp.{phases}.p_ctg.fa')
  output:
    outputSummary = os.path.join(mainDir, '{sample}.{phases}', 'short_summary.specific.laurasiatheria_odb10.{sample}.{phases}.txt'),
    full_table = os.path.join(mainDir, '{sample}.{phases}','run_laurasiatheria_odb10', 'full_table.tsv')
  shell: '''
    busco -c 20 --metaeuk \
      --metaeuk_parameters="--disk-space-limit=1000M,--remove-tmp-files=1" \
      --metaeuk_rerun_parameters="--disk-space-limit=1000M,--remove-tmp-files=1" \
      -i {input.deNovoFa} \
      -o {wildcards.sample}.{wildcards.phases} \
      -m genome \
      -f \
      -l laurasiatheria_odb10
  '''

rule parseBUSCO:
  input:
    outputSummary = os.path.join(mainDir, '{sample}.{phases}', 'short_summary.specific.laurasiatheria_odb10.{sample}.{phases}.txt'),
    full_table = os.path.join(mainDir, '{sample}.{phases}','run_laurasiatheria_odb10', 'full_table.tsv')
  params:
    buscoDir = os.path.join(mainDir, '{sample}.{phases}')
  output:
    buscoCSV = os.path.join(mainDir, '{sample}.{phases}', 'busco.csv')
  shell: '''
    echo "Strain,Complete_single_copy,Complete_duplicated,Fragmented,Missing" > {output.buscoCSV}

    cat {input.outputSummary} | grep "(S)" | awk -v strain="{wildcards.sample}.{wildcards.phases}" '{{print strain","$1}}' > {params.buscoDir}/complete_single.txt

    cat {input.outputSummary} | grep "(D)" | awk '{{print $1}}' > {params.buscoDir}/complete_duplicated.txt

    cat {input.outputSummary} | grep "(F)" | awk '{{print $1}}' > {params.buscoDir}/fragmented.txt

    cat {input.outputSummary} | grep "(M)" | awk '{{print $1}}' > {params.buscoDir}/missing.txt

    paste -d "," {params.buscoDir}/complete_single.txt {params.buscoDir}/complete_duplicated.txt {params.buscoDir}/fragmented.txt {params.buscoDir}/missing.txt >> {output.buscoCSV}

    rm {params.buscoDir}/complete_single.txt {params.buscoDir}/complete_duplicated.txt {params.buscoDir}/fragmented.txt {params.buscoDir}/missing.txt
  '''

rule combineBUSCO:
  input:
    buscoCSV = expand(os.path.join(mainDir, '{sample}.{phases}', 'busco.csv'), sample = SAMPLES, phases = PHASES)
  output:
    combinedBuscoCSV = os.path.join(mainDir, 'combinedBusco.csv')
  shell: '''
    echo "Strain,Complete_single_copy,Complete_duplicated,Fragmented,Missing" > {output.combinedBuscoCSV}
    cat {input.buscoCSV} | sort | uniq -u >> {output.combinedBuscoCSV}
  '''

rule buscoRscript:
  input:
    combinedBuscoCSV = os.path.join(mainDir, 'combinedBusco.csv')
  params:
    buscoRscript = os.path.join(mainDir, 'scripts', 'busco.R')
  output:
    buscoPDF = os.path.join(mainDir, 'busco.pdf')
  shell: '''
    Rscript {params.buscoRscript} \
      {input.combinedBuscoCSV} \
      {output.buscoPDF}
  '''

# separates misassemblies - added this step for horse_01_200's chr20 which was unusually short
rule correct:
  resources:
    mem_mb=50000  
  input:
    buscoPDF = os.path.join(mainDir, 'busco.pdf'),
    deNovoFa = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.asm.bp.{phases}.p_ctg.fa'),
    ref = reference
  output:
    ragtagFasta = os.path.join(mainDir, '{sample}.{phases}', 'ragtag.correct.fasta')
  shell: '''
    # this ensures that all the code for this step is run again
    rm -f {wildcards.sample}.{wildcards.phases}/ragtag.correct.* 

    ragtag.py correct \
      -u \
      -o {wildcards.sample}.{wildcards.phases} \
      {input.ref} \
      {input.deNovoFa}
  '''

rule correctPrimary:
  resources:
    mem_mb=50000  
  input:
    buscoPDF = os.path.join(mainDir, 'busco.pdf'),
    deNovoFa = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.asm.bp.p_ctg.fa'),
    ref = reference
  output:
    ragtagFasta = os.path.join(mainDir, '{sample}', '{sample}.ragtag.correct.fasta')
  shell: '''
    # this ensures that all the code for this step is run again
    rm -f {wildcards.sample}/ragtag.correct.* 

    ragtag.py correct \
      -u \
      -o {wildcards.sample} \
      {input.ref} \
      {input.deNovoFa}

    mv {wildcards.sample}/ragtag.correct.fasta {output.ragtagFasta}
  '''

# This uses the reference genome to make scaffolds out of the contigs - assembling contigs into bigger contigs (the scaffolds)
rule scaffold:
  resources:
    mem_mb=50000  
  input:
    ragtagFasta = os.path.join(mainDir, '{sample}.{phases}', 'ragtag.correct.fasta'),
    ref = reference
  output:
    ragtagFasta = os.path.join(mainDir, '{sample}.{phases}', 'ragtag.scaffold.fasta')
  shell: '''
    # this ensures that all the code for this step is run again
    rm -f {wildcards.sample}.{wildcards.phases}/ragtag.scaffold.*

    ragtag.py scaffold \
      -u \
      -o {wildcards.sample}.{wildcards.phases} \
      {input.ref} \
      {input.ragtagFasta}
  '''

rule scaffoldPrimary:
  resources:
    mem_mb=50000  
  input:
    ragtagFasta = os.path.join(mainDir, '{sample}', '{sample}.ragtag.correct.fasta'),
    ref = reference
  output:
    ragtagFasta = os.path.join(mainDir, '{sample}', '{sample}.ragtag.scaffold.fasta')
  shell: '''
    # this ensures that all the code for this step is run again
    rm -f {wildcards.sample}/ragtag.scaffold.*

    ragtag.py scaffold \
      -u \
      -o {wildcards.sample} \
      {input.ref} \
      {input.ragtagFasta}

    mv {wildcards.sample}/ragtag.scaffold.fasta {output.ragtagFasta}
  '''

rule removeRagTag:
  input:
    ragtagFasta = os.path.join(mainDir, '{sample}.{phases}', 'ragtag.scaffold.fasta')
  output:
    ragtagRenamedFasta = os.path.join(mainDir, '{sample}.{phases}', '{sample}.{phases}_rename.scaffold.fasta')
  shell: '''
    sed 's/_RagTag//' {input.ragtagFasta} > {output.ragtagRenamedFasta}
  '''

rule removeRagTagPrimary:
  input:
    ragtagFasta = os.path.join(mainDir, '{sample}', '{sample}.ragtag.scaffold.fasta')
  output:
    ragtagRenamedFasta = os.path.join(mainDir, '{sample}', '{sample}.rename.scaffold.fasta')
  shell: '''
    sed 's/_RagTag//' {input.ragtagFasta} > {output.ragtagRenamedFasta}
  '''

rule chrNames:
  input:
    refAssemblyReport = refAssemblyReport
  output:
    chrNameList = os.path.join(mainDir, 'resources', 'chromosome.name.list.txt'),
    accessionList = os.path.join(mainDir, 'resources', 'chr.accesssion.list.txt')
  shell: '''
    grep "^chr" {input.refAssemblyReport} | awk '{{ print $1 }}' > {output.chrNameList}
    grep "^chr" {input.refAssemblyReport} | awk '{{ print $1, $7 }}' > {output.accessionList}
  '''

rule cleanRef:
  input:
    ref = reference
  output:
    refCleaned = os.path.join(mainDir, 'resources', 'cleaned_ref_genome.fa')
  shell: '''
    sed -e 's/\\(^[>].*$\\)/#\\1#/' /home/tmathieu/eNAD/resources/GCF_002863925.1_EquCab3.0_genomic.fna \
    | tr -d "\n" \
    | awk 'BEGIN {{RS="#"}} {{print}}' \
    > {output.refCleaned}
  '''

rule sameChrNames:
  input:
    ragtagRenamedFasta = os.path.join(mainDir, '{sample}.{phases}', '{sample}.{phases}_rename.scaffold.fasta'),
    chrNameList = os.path.join(mainDir, 'resources', 'chr.accesssion.list.txt'),
    refCleaned = os.path.join(mainDir, 'resources', 'cleaned_ref_genome.fa')
  output:
    ragtagChr = os.path.join(mainDir, '{sample}.{phases}', '{sample}.{phases}_scaffold_chr.fa')
  shell: '''
    {wildcards.sample}lines=$(wc -l {input.chrNameList} | awk '{{ print $1 }}')
    {wildcards.sample}currentline=1
    while [ "${wildcards.sample}currentline" -le "${wildcards.sample}lines" ]
    do
      {wildcards.sample}name=$(sed "${{{wildcards.sample}currentline}}q;d" {input.chrNameList} | awk '{{ print $2 }}')
      grep -A1 "${{{wildcards.sample}name}}" {input.ragtagRenamedFasta} >> {output.ragtagChr}
      {wildcards.sample}currentline=$((${wildcards.sample}currentline+1))
    done
  '''

rule sameChrNamesPrimary:
  input:
    ragtagRenamedFasta = os.path.join(mainDir, '{sample}', '{sample}.rename.scaffold.fasta'),
    chrNameList = os.path.join(mainDir, 'resources', 'chr.accesssion.list.txt'),
    refCleaned = os.path.join(mainDir, 'resources', 'cleaned_ref_genome.fa')
  output:
    ragtagChr = os.path.join(mainDir, '{sample}', '{sample}.scaffold_chr.fa')
  shell: '''
    {wildcards.sample}lines=$(wc -l {input.chrNameList} | awk '{{ print $1 }}')
    {wildcards.sample}currentline=1
    while [ "${wildcards.sample}currentline" -le "${wildcards.sample}lines" ]
    do
      {wildcards.sample}name=$(sed "${{{wildcards.sample}currentline}}q;d" {input.chrNameList} | awk '{{ print $2 }}')
      grep -A1 "${{{wildcards.sample}name}}" {input.ragtagRenamedFasta} >> {output.ragtagChr}
      {wildcards.sample}currentline=$((${wildcards.sample}currentline+1))
    done
  '''

rule refSameChrNames:
  input:
    chrNameList = os.path.join(mainDir, 'resources', 'chr.accesssion.list.txt'),
    refCleaned = os.path.join(mainDir, 'resources', 'cleaned_ref_genome.fa')
  output:
    refChr = os.path.join(mainDir, 'resources', 'reference_genome_chr.fa')
  shell: '''
    reflines=$(wc -l {input.chrNameList} | awk '{{ print $1 }}')
    refcurrentline=1
    while [ "$refcurrentline" -le "$reflines" ]
    do
      refname=$(sed "${{refcurrentline}}q;d" {input.chrNameList} | awk '{{ print $2 }}')
      grep -A1 "${{refname}}" {input.refCleaned} >> {output.refChr}
      refcurrentline=$(($refcurrentline+1))
    done
  '''

rule refChangeNCtoChr:
  input:
    refChr = os.path.join(mainDir, 'resources', 'reference_genome_chr.fa'),
    chrNameList = os.path.join(mainDir, 'resources', 'chr.accesssion.list.txt')
  output:
    refChr2 = os.path.join(mainDir, 'resources', 'reference_genome_chr_only.fa')
  shell: '''
    cp {input.refChr} {output.refChr2}
    reflines=$(wc -l {input.chrNameList} | awk '{{ print $1 }}')
    refcurrentline=1
    while [ "$refcurrentline" -le "$reflines" ]
    do
      refchr=$(sed "${{refcurrentline}}q;d" {input.chrNameList} | awk '{{ print $1 }}')
      refnc=$(sed "${{refcurrentline}}q;d" {input.chrNameList} | awk '{{ print $2 }}')
      sed -i "s/${{refnc}}/${{refchr}}/g" {output.refChr2}
      refcurrentline=$(($refcurrentline+1))
    done
  '''

rule changeNCtoChr:
  input:
    ragtagChr = os.path.join(mainDir, '{sample}.{phases}', '{sample}.{phases}_scaffold_chr.fa'),
    chrNameList = os.path.join(mainDir, 'resources', 'chr.accesssion.list.txt'),
    refChr = os.path.join(mainDir, 'resources', 'reference_genome_chr_only.fa')
  output:
    ragtagChr2 = os.path.join(mainDir, '{sample}.{phases}', '{sample}.{phases}_scaffold_chr_only.fa'),
    refChr2 = os.path.join(mainDir, '{sample}.{phases}', '{sample}.{phases}_ref_chr_only.fa')
  shell: '''
    cp {input.ragtagChr} {output.ragtagChr2}
    cp {input.refChr} {output.refChr2}
    {wildcards.sample}lines=$(wc -l {input.chrNameList} | awk '{{ print $1 }}')
    {wildcards.sample}currentline=1
    while [ "${wildcards.sample}currentline" -le "${wildcards.sample}lines" ]
    do
      {wildcards.sample}chr=$(sed "${{{wildcards.sample}currentline}}q;d" {input.chrNameList} | awk '{{ print $1 }}')
      {wildcards.sample}nc=$(sed "${{{wildcards.sample}currentline}}q;d" {input.chrNameList} | awk '{{ print $2 }}')
      sed -i "s/${{{wildcards.sample}nc}}/${{{wildcards.sample}chr}}/g" {output.ragtagChr2}
      {wildcards.sample}currentline=$((${wildcards.sample}currentline+1))
    done
  '''

rule changeNCtoChrPrimary:
  input:
    ragtagChr = os.path.join(mainDir, '{sample}', '{sample}.scaffold_chr.fa'),
    chrNameList = os.path.join(mainDir, 'resources', 'chr.accesssion.list.txt'),
    refChr = os.path.join(mainDir, 'resources', 'reference_genome_chr_only.fa')
  output:
    ragtagChr2 = os.path.join(mainDir, '{sample}', '{sample}.scaffold_chr_only.fa'),
    refChr2 = os.path.join(mainDir, '{sample}', '{sample}.ref_chr_only.fa')
  shell: '''
    cp {input.ragtagChr} {output.ragtagChr2}
    cp {input.refChr} {output.refChr2}
    {wildcards.sample}lines=$(wc -l {input.chrNameList} | awk '{{ print $1 }}')
    {wildcards.sample}currentline=1
    while [ "${wildcards.sample}currentline" -le "${wildcards.sample}lines" ]
    do
      {wildcards.sample}chr=$(sed "${{{wildcards.sample}currentline}}q;d" {input.chrNameList} | awk '{{ print $1 }}')
      {wildcards.sample}nc=$(sed "${{{wildcards.sample}currentline}}q;d" {input.chrNameList} | awk '{{ print $2 }}')
      sed -i "s/${{{wildcards.sample}nc}}/${{{wildcards.sample}chr}}/g" {output.ragtagChr2}
      {wildcards.sample}currentline=$((${wildcards.sample}currentline+1))
    done
  '''

rule alignFastQToReference:
  resources:
    mem_mb=50000 
  threads: 20
  input:
    ref = os.path.join(mainDir, 'resources', 'reference_genome_chr_only.fa'),
    fastq = os.path.join(mainDir, 'fastq', '{sample}.fastq')
  params:
    alignDir = os.path.join(mainDir, 'refAlign')
  output:
    sam = os.path.join(mainDir, 'refAlign', '{sample}.sam')
  shell: '''
    mkdir -p {params.alignDir}

    minimap2 \
      -a \
      -x map-hifi \
      --eqx \
      -t 20 \
      {input.ref} \
      {input.fastq} \
      > {output.sam}
  '''

rule convertSamToBamSortIndex:
  resources:
    mem_mb=50000 
  threads: 20
  input:
     os.path.join(mainDir, 'refAlign', '{sample}.sam')
  output:
    bam = os.path.join(mainDir, 'refAlign', '{sample}.sorted.bam'),
    bai = os.path.join(mainDir, 'refAlign', '{sample}.sorted.bam.bai')
  shell: '''
    samtools view -@ 20 -S -b {input} > temp{wildcards.sample}.bam
    samtools sort -@ 20 -o temp{wildcards.sample}.bam > {output.bam}
    samtools index -@ 20 -b {output.bam}
    rm {input}
    rm temp{wildcards.sample}.bam
  '''

rule coverageStatsInitial:
  input:
    os.path.join(mainDir, 'refAlign', '{sample}.sorted.bam')
  params:
    os.path.join(mainDir, 'coverageStatsRef.txt')
  output: 
    stats = os.path.join(mainDir, 'refAlign', '{sample}_statsDone.txt')
  shell: '''
    cov=$(samtools depth -a {input} | awk '{{sum+=$3}} END {{print sum/NR}}') 
    bamstats{wildcards.comp_sample}=$(bamtools stats -in {input})
    echo {wildcards.comp_sample}\tcoverage>\t${{cov}} >> {params}
    echo $bamstats{wildcards.comp_sample} >> {params}
    echo {wildcards.comp_sample}\tcoverage>\t${{cov}} >> {output.stats}
    echo $bamstats{wildcards.comp_sample} >> {output.stats}
  '''

rule alignToChainOfGenomes:
  resources:
    mem_mb=50000 
  input:
    comp1 = os.path.join(mainDir, '{comp1}', '{comp1}.scaffold_chr_only.fa'),
    comp2 = os.path.join(mainDir, '{comp2}', '{comp2}.scaffold_chr_only.fa')
  params:
    chainDir = os.path.join(mainDir, 'chainAlign')
  output:
    compareSam = os.path.join(mainDir, 'chainAlign', '{comp1}-{comp2}comparison.sam')
  shell: '''
    mkdir -p {params.chainDir}

    minimap2 \
      -a \
      -x map-hifi \
      --eqx \
      {input.comp1} \
      {input.comp2} \
      > {output.compareSam}
  '''

rule syntenyAnalysisForCombined:
  threads: 15
  input:
    comp1 = os.path.join(mainDir, '{comp1}', '{comp1}.scaffold_chr_only.fa'),
    comp2 = os.path.join(mainDir, '{comp2}', '{comp2}.scaffold_chr_only.fa'),
    compareSam = os.path.join(mainDir, 'chainAlign', '{comp1}-{comp2}comparison.sam')
  params:
    outputDir = os.path.join(mainDir, 'chainAlign', '{comp1}-{comp2}syri')
  output:
    outFile = os.path.join(mainDir, 'chainAlign', '{comp1}-{comp2}syri', '{comp1}-{comp2}syri.out')
  shell: '''
    mkdir -p {params.outputDir}

    syri \
      -c {input.compareSam} \
      -r {input.comp1} \
      -q {input.comp2} \
      -k -F S \
      --prefix {wildcards.comp1}-{wildcards.comp2}\
      --dir {params.outputDir} \
      --nc {threads} \
      --no-chrmatch
  '''

rule makeGenomeFileForPlotSRCombined:
  input:
    control = controlprimary
  params:
    control_name = order_comps[0],
    order_comps = second_comp,
    mainDir = mainDir
  output:
    genomeTxtFile = os.path.join(mainDir, 'syri_genomes.txt')
  shell: '''
    echo "#file  name" > {output.genomeTxtFile}
    echo "{input.control}\t{params.control_name}\tlw:6" >> {output.genomeTxtFile}
    COMPS=({params.order_comps})
    for sample in ${{COMPS[@]}}
      do
        echo "{params.mainDir}/$sample/${{sample}}.scaffold_chr_only.fa\t$sample\tlw:6" >> {output.genomeTxtFile}
      done
  '''

rule plotStructuralRearrangementCombined:
  resources:
    mem_mb=50000
  input:
    outFileChain = expand(os.path.join(mainDir, 'chainAlign', '{comp1}-{comp2}syri', '{comp1}-{comp2}syri.out'), zip, comp1 = first_comp, comp2 = second_comp),
    genomeTxtFile = os.path.join(mainDir, 'syri_genomes.txt'),
    svBed = os.path.join(mainDir, 'snifflesPrimary', 'candidates_filtered.bed')
  params:
    plotsrGraphs = os.path.join(mainDir, 'plotsrGraphsWithSVs')
  output:
    plot = os.path.join(mainDir, 'plotsrGraphsWithSVs', 'syri_{chr}_combined_output_plot.png')
  shell: '''
    mkdir -p {params.plotsrGraphs}

    COMPS=({input.outFileChain})
    SR=()
    for sample in ${{COMPS[@]}}
      do
        SR+=("--sr $sample ")
      done
    echo ${{SR[@]}}
    echo "plotsr \
      ${{SR[@]}} \
      --genomes {input.genomeTxtFile} \
      --markers {input.svBed} \
      -H 8 \
      -W 12 \
      --chr {wildcards.chr} \
      -o {output.plot}" \
      | bash
  '''

rule alignToFirstControl:
  resources:
    mem_mb=50000
  input:
    control = controlphase,
    ragtagChr2 = os.path.join(mainDir, '{comp_sample}.{phases}', '{comp_sample}.{phases}_scaffold_chr_only.fa')
  params:
    snifflesDir = os.path.join(mainDir, 'sniffles')
  output:
    compareToControl = os.path.join(mainDir, 'sniffles', '{comp_sample}.{phases}_compto5.sam')
  shell: '''
    mkdir -p {params.snifflesDir}

    minimap2 \
      -a \
      -x map-hifi \
      --eqx \
      {input.control} \
      {input.ragtagChr2} \
      > {output.compareToControl}
  '''

rule alignToFirstControlPrimary:
  resources:
    mem_mb=50000
  input:
    control = controlprimary,
    ragtagChr2 = os.path.join(fastqDir, '{sample}.fastq')
  params:
    snifflesDir = os.path.join(mainDir, 'snifflesPrimary')
  output:
    compareToControl = os.path.join(mainDir, 'snifflesPrimary', '{sample}.primary_compto5.sam')
  shell: '''
    mkdir -p {params.snifflesDir}

    minimap2 \
      -a \
      -x map-hifi \
      --eqx \
      -t 16 \
      {input.control} \
      {input.ragtagChr2} \
      > {output.compareToControl}
  '''

rule convertSamToBam:
  input:
     os.path.join(mainDir, 'sniffles', '{comp_sample}.{phases}_compto5.sam')
  output:
    os.path.join(mainDir, 'sniffles', '{comp_sample}.{phases}_compto5.bam')
  shell: '''
    samtools view -S -b {input} > {output}
    rm {input}
  '''

rule sortBam:
  input:
    os.path.join(mainDir, 'sniffles', '{comp_sample}.{phases}_compto5.bam')
  output:
    os.path.join(mainDir, 'sniffles', '{comp_sample}.{phases}_compto5_sorted.bam')
  shell: '''
    samtools sort -o {output} {input}
    rm {input}
  '''

rule indexBam:
  input:
    os.path.join(mainDir, 'sniffles', '{comp_sample}.{phases}_compto5_sorted.bam')
  output:
    os.path.join(mainDir, 'sniffles', '{comp_sample}.{phases}_compto5_sorted.bam.bai')
  shell: '''
    samtools index -b {input}
  '''

rule findSVs:
  resources:
    mem_mb=50000
  input:
    compareToControlBam = os.path.join(mainDir, 'sniffles', '{comp_sample}.{phases}_compto5_sorted.bam'),
    bai = os.path.join(mainDir, 'sniffles', '{comp_sample}.{phases}_compto5_sorted.bam.bai'),
    refChr = controlphase
  output:
    snf = os.path.join(mainDir, 'sniffles', '{comp_sample}.{phases}.snf'),
    vcf = os.path.join(mainDir, 'sniffles', '{comp_sample}.{phases}.vcf')
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
    snfhap2 = os.path.join(mainDir, 'sniffles', order_comps[0] + '.hap2.snf'),
    snf = expand(os.path.join(mainDir, 'sniffles', '{comp_sample}.{phases}.snf'), comp_sample = second_comp, phases = PHASES),
    refChr = controlphase
  output:
    multiVCF = os.path.join(mainDir, 'sniffles', 'phasedSVs1_compto5.vcf')
  shell: '''
    sniffles \
      --input {input.snf} {input.snfhap2} \
      --vcf {output.multiVCF} \
      --reference {input.refChr} \
      --threads 16
  '''

rule convertSamToBamPrimary:
  input:
     os.path.join(mainDir, 'snifflesPrimary', '{comp_sample}.primary_compto5.sam')
  output:
    os.path.join(mainDir, 'snifflesPrimary', '{comp_sample}.primary_compto5.bam')
  shell: '''
    samtools view -S -b {input} > {output}
    rm {input}
  '''

rule sortBamPrimary:
  input:
    os.path.join(mainDir, 'snifflesPrimary', '{comp_sample}.primary_compto5.bam')
  output:
    os.path.join(mainDir, 'snifflesPrimary', '{comp_sample}.primary_compto5_sorted.bam')
  shell: '''
    samtools sort -o {output} {input}
    rm {input}
  '''

rule indexBamPrimary:
  input:
    os.path.join(mainDir, 'snifflesPrimary', '{comp_sample}.primary_compto5_sorted.bam')
  output:
    os.path.join(mainDir, 'snifflesPrimary', '{comp_sample}.primary_compto5_sorted.bam.bai')
  shell: '''
    samtools index -b {input}
  '''

rule coverageStats:
  input:
    os.path.join(mainDir, 'snifflesPrimary', '{comp_sample}.primary_compto5_sorted.bam')
  params:
    os.path.join(mainDir, 'coverageStatsComp5.txt')
  output: 
    stats = os.path.join(mainDir, 'snifflesPrimary', '{comp_sample}_statsDone.txt')
  shell: '''
    cov=$(samtools depth -a {input} | awk '{{sum+=$3}} END {{print sum/NR}}') 
    bamstats{wildcards.comp_sample}=$(bamtools stats -in {input})
    echo {wildcards.comp_sample}\tcoverage>\t${{cov}} >> {params}
    echo $bamstats{wildcards.comp_sample} >> {params}
    echo {wildcards.comp_sample}\tcoverage>\t${{cov}} >> {output.stats}
    echo $bamstats{wildcards.comp_sample} >> {output.stats}
  '''

rule findSVsPrimary:
  resources:
    mem_mb=50000
  input:
    compareToControlBam = os.path.join(mainDir, 'snifflesPrimary', '{comp_sample}.primary_compto5_sorted.bam'),
    bai = os.path.join(mainDir, 'snifflesPrimary', '{comp_sample}.primary_compto5_sorted.bam.bai'),
    refChr = controlprimary
  output:
    snf = os.path.join(mainDir, 'snifflesPrimary', '{comp_sample}.primary.snf'),
    vcf = os.path.join(mainDir, 'snifflesPrimary', '{comp_sample}.primary.vcf')
  shell: '''
    sniffles \
      --input {input.compareToControlBam} \
      --snf {output.snf} \
      --vcf {output.vcf} \
      #--reference {input.refChr} \
      --threads 16
  '''

rule compareSVsPrimary:
  resources:
    mem_mb=50000
  input:
    snf = expand(os.path.join(mainDir, 'snifflesPrimary', '{comp_sample}.primary.snf'), comp_sample = second_comp),
    refChr = controlprimary
  output:
    multiVCF = os.path.join(mainDir, 'snifflesPrimary', 'primarySVs1_compto5.vcf')
  shell: '''
    sniffles \
      --input {input.snf} \
      --vcf {output.multiVCF} \
      --reference {input.refChr} \
      --phase \
      --threads 16
  '''

rule cleanAndFilterVCF:
  input:
    multiVCF = os.path.join(mainDir, 'snifflesPrimary', 'primarySVs1_compto5.vcf')
  output:
    filteredVCF = os.path.join(mainDir, 'snifflesPrimary', 'primarySVs_filtered.vcf')
  shell: '''
    sed "s/GT:GQ:DR:DV:PS:ID/GT:GQ:DR:DV:ID/g" {input.multiVCF} > temp_cleaned.vcf
    bgzip temp_cleaned.vcf -k
    bcftools index temp_cleaned.vcf.gz
    bcftools view -i '(GT[0]=="0/0" && GT[1]=="0/0") && ((GT[2]=="0/1" || GT[2]=="1/0" || GT[2]=="1/1") && (GT[3]=="0/1" || GT[3]=="1/0" || GT[3]=="1/1") && (GT[4]=="0/1" || GT[4]=="1/0" || GT[4]=="1/1"))' temp_cleaned.vcf.gz -Ov -o {output.filteredVCF}
    rm temp_cleaned*
  '''

rule convert2Bed:
  input:
    filteredVCF = os.path.join(mainDir, 'snifflesPrimary', 'primarySVs_filtered.vcf')
  output:
    filteredBED = os.path.join(mainDir, 'snifflesPrimary', 'primarySVs_filtered.bed')
  shell: '''
    gvcf2bed -I {input} -O {output}
  '''

rule makeCandidateBed:
  input:
    os.path.join(mainDir, 'snifflesPrimary', 'primarySVs_filtered.bed')
  output:
    candidates = os.path.join(mainDir, 'snifflesPrimary', 'candidates_filtered.bed')
  shell: '''
    sed -e "s/$/\thorse_03_238\tmt:|;mc:green;tp:0.0;tt:candidate/" {input} > {output}
  '''

rule plotSVs:
  input:
    multiVCF = os.path.join(mainDir, 'sniffles', 'phasedSVs1.vcf')
  params:
    plotDir = os.path.join(mainDir, 'snifflesPlot')
  output:
    os.path.join(mainDir, 'snifflesPlot', 'plotDone.txt')
  shell: '''
    python3 -m sniffles2_plot -i {input.multiVCF} -o {params.plotDir}

    echo "Plots made" > {output}
  '''

rule plotSVsPrimary:
  input:
    multiVCF = os.path.join(mainDir, 'snifflesPrimary', 'primarySVs1.vcf')
  params:
    plotDir = os.path.join(mainDir, 'snifflesPlotPrimary')
  output:
    os.path.join(mainDir, 'snifflesPlotPrimary', 'primaryPlotDone.txt')
  shell: '''
    python3 -m sniffles2_plot -i {input.multiVCF} -o {params.plotDir}

    echo "Plots made" > {output}
  '''

