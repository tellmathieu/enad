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

###################### --- PIPELINE--- #######################
### --- Don't change unless you know what you're doing --- ###

### Using this tutorial https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9254108/

mainDir = os.getcwd()
outputGenomeDir = "outputGenomes"

SAMPLES, = glob_wildcards(os.path.join(fastqDir, '{sample}.fastq'))

f = open(os.path.join(mainDir, 'resources', 'chromosome.name.list.txt'))
chrList = f.readlines()
for i in range(0,len(chrList)):
  chrList[i] = chrList[i].strip()
f.close()

first_comp = order_comps[0:-1]
second_comp = order_comps[1:]

sniffles_comp = order_comps[1:]

print(sniffles_comp)

'''print(first_comp)
print(second_comp)
print(order_comps[0])
print(SAMPLES)'''

rule all:
  input:
    expand(os.path.join(mainDir, '{sample}', 'syri', '{sample}_output_plot.png'), sample = SAMPLES),
    expand(os.path.join(mainDir, 'plotsrGraphs', 'syri_{chr}_combined_output_plot.png'), chr = chrList),
    expand(os.path.join(mainDir, 'plotsrGraphsPOI', 'syri_{chr}_combined_output_plot.png'), chr = chrList),
    os.path.join(mainDir, 'sniffles', 'nonphaseSVs.vcf')

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
    deNovogfa = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.asm.bp.p_ctg.gfa')
  shell: '''
    mkdir -p {params.outputDir}/{wildcards.sample}
    
    #de novo assembly
    hifiasm -o {params.outputDir}/{wildcards.sample}/{wildcards.sample}.asm -t {params.threads} {input.fastq}
  '''

rule de_novo_gfa_to_fa:
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
    deNovoFa = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.asm.bp.p_ctg.fa')
  output:
    deNovoAssemblyStatsTxt = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.asm.bp.p_ctg.fa.txt')
  shell: '''
    assembly-stats {input.deNovoFa} > {output.deNovoAssemblyStatsTxt}
  '''

rule create_coverage_table:
  input:
    deNovoFa = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.asm.bp.p_ctg.fa')
  output:
    length = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.length.csv')
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
    length = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.length.csv')
  params:
    rscript = os.path.join(mainDir, 'scripts', 'coveragePlot.R')
  output:
    plot = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.coverage.pdf')
  shell: '''
    Rscript {params.rscript} \
      {input.length} \
      {output.plot}
  '''

rule BUSCOanalysis:
  resources:
    mem_mb=70000
  input:
    deNovoFa = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.asm.bp.p_ctg.fa')
  output:
    outputSummary = os.path.join(mainDir, '{sample}', 'short_summary.specific.laurasiatheria_odb10.{sample}.txt'),
    full_table = os.path.join(mainDir, '{sample}','run_laurasiatheria_odb10', 'full_table.tsv')
  shell: '''
    busco -c 20 --metaeuk \
      --metaeuk_parameters="--disk-space-limit=1000M,--remove-tmp-files=1" \
      --metaeuk_rerun_parameters="--disk-space-limit=1000M,--remove-tmp-files=1" \
      -i {input.deNovoFa} \
      -o {wildcards.sample} \
      -m genome \
      -f \
      -l laurasiatheria_odb10
  '''

rule parseBUSCO:
  input:
    outputSummary = os.path.join(mainDir, '{sample}', 'short_summary.specific.laurasiatheria_odb10.{sample}.txt'),
    full_table = os.path.join(mainDir, '{sample}','run_laurasiatheria_odb10', 'full_table.tsv')
  params:
    buscoDir = os.path.join(mainDir, '{sample}')
  output:
    buscoCSV = os.path.join(mainDir, '{sample}', 'busco.csv')
  shell: '''
    echo "Strain,Complete_single_copy,Complete_duplicated,Fragmented,Missing" > {output.buscoCSV}

    cat {input.outputSummary} | grep "(S)" | awk -v strain="{wildcards.sample}" '{{print strain","$1}}' > {params.buscoDir}/complete_single.txt

    cat {input.outputSummary} | grep "(D)" | awk '{{print $1}}' > {params.buscoDir}/complete_duplicated.txt

    cat {input.outputSummary} | grep "(F)" | awk '{{print $1}}' > {params.buscoDir}/fragmented.txt

    cat {input.outputSummary} | grep "(M)" | awk '{{print $1}}' > {params.buscoDir}/missing.txt

    paste -d "," {params.buscoDir}/complete_single.txt {params.buscoDir}/complete_duplicated.txt {params.buscoDir}/fragmented.txt {params.buscoDir}/missing.txt >> {output.buscoCSV}

    rm {params.buscoDir}/complete_single.txt {params.buscoDir}/complete_duplicated.txt {params.buscoDir}/fragmented.txt {params.buscoDir}/missing.txt
  '''

rule combineBUSCO:
  input:
    buscoCSV = expand(os.path.join(mainDir, '{sample}', 'busco.csv'), sample = SAMPLES)
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
    deNovoFa = os.path.join(mainDir, outputGenomeDir, '{sample}', '{sample}.asm.bp.p_ctg.fa'),
    ref = reference
  output:
    ragtagFasta = os.path.join(mainDir, '{sample}', 'ragtag.correct.fasta')
  shell: '''
    # this ensures that all the code for this step is run again
    rm -f {wildcards.sample}/ragtag.correct.* 

    ragtag.py correct \
      -u \
      -o {wildcards.sample} \
      {input.ref} \
      {input.deNovoFa}
  '''

# This uses the reference genome to make scaffolds out of the contigs - assembling contigs into bigger contigs (the scaffolds)
rule scaffold:
  resources:
    mem_mb=50000  
  input:
    ragtagFasta = os.path.join(mainDir, '{sample}', 'ragtag.correct.fasta'),
    ref = reference
  output:
    ragtagFasta = os.path.join(mainDir, '{sample}', 'ragtag.scaffold.fasta')
  shell: '''
    # this ensures that all the code for this step is run again
    rm -f {wildcards.sample}/ragtag.scaffold.*

    ragtag.py scaffold \
      -u \
      -o {wildcards.sample} \
      {input.ref} \
      {input.ragtagFasta}
  '''

rule removeRagTag:
  input:
    ragtagFasta = os.path.join(mainDir, '{sample}', 'ragtag.scaffold.fasta')
  output:
    ragtagRenamedFasta = os.path.join(mainDir, '{sample}', '{sample}_rename.scaffold.fasta')
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
    ragtagRenamedFasta = os.path.join(mainDir, '{sample}', '{sample}_rename.scaffold.fasta'),
    chrNameList = os.path.join(mainDir, 'resources', 'chr.accesssion.list.txt'),
    refCleaned = os.path.join(mainDir, 'resources', 'cleaned_ref_genome.fa')
  output:
    ragtagChr = os.path.join(mainDir, '{sample}', '{sample}_scaffold_chr.fa')
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
    ragtagChr = os.path.join(mainDir, '{sample}', '{sample}_scaffold_chr.fa'),
    chrNameList = os.path.join(mainDir, 'resources', 'chr.accesssion.list.txt'),
    refChr = os.path.join(mainDir, 'resources', 'reference_genome_chr_only.fa')
  output:
    ragtagChr2 = os.path.join(mainDir, '{sample}', '{sample}_scaffold_chr_only.fa'),
    refChr2 = os.path.join(mainDir, '{sample}', '{sample}_ref_chr_only.fa')
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

rule alignToReference:
  resources:
    mem_mb=50000 
  input:
    ragtagChr2 = os.path.join(mainDir, '{sample}', '{sample}_scaffold_chr_only.fa'),
    refChr2 = os.path.join(mainDir, '{sample}', '{sample}_ref_chr_only.fa')
  output:
    compareSam = os.path.join(mainDir, '{sample}', '{sample}_comparison.sam')
  shell: '''
    minimap2 \
      -a \
      -x asm5 \
      --eqx \
      {input.refChr2} \
      {input.ragtagChr2} \
      > {output.compareSam}
  '''

rule syntenyAnalysis:
  input:
    ragtagChr2 = os.path.join(mainDir, '{sample}', '{sample}_scaffold_chr_only.fa'),
    refChr2 = os.path.join(mainDir, '{sample}', '{sample}_ref_chr_only.fa'),
    compareSam = os.path.join(mainDir, '{sample}', '{sample}_comparison.sam')
  params:
    outputDir = os.path.join(mainDir, '{sample}', 'syri')
  output:
    outFile = os.path.join(mainDir, '{sample}', 'syri', 'syri.out')
  shell: '''
    mkdir -p {params.outputDir}

    syri \
      -c {input.compareSam} \
      -r {input.refChr2} \
      -q {input.ragtagChr2} \
      -k -F S \
      --dir {params.outputDir}
  '''

rule makeGenomeFileForPlotSR:
  input:
    outFile = os.path.join(mainDir, '{sample}', 'syri', 'syri.out'),
    ragtagChr2 = os.path.join(mainDir, '{sample}', '{sample}_scaffold_chr_only.fa'),
    refChr2 = os.path.join(mainDir, '{sample}', '{sample}_ref_chr_only.fa')
  output:
    genomeTxtFile = os.path.join(mainDir, '{sample}', 'syri', 'genomes.txt'),
    renamedOutFile = os.path.join(mainDir, '{sample}', 'syri', 'reference_{sample}syri.out'),
  shell: '''
    echo "#file  name" > {output.genomeTxtFile}
    echo "{input.refChr2}\treference\tlw:1.5" >> {output.genomeTxtFile}
    echo "{input.ragtagChr2}\t{wildcards.sample}\tlw:1.5" >> {output.genomeTxtFile}

    cp {input.outFile} {output.renamedOutFile}
  '''

rule plotStructuralRearrangement:
  input:
    outFile = os.path.join(mainDir, '{sample}', 'syri', 'syri.out'),
    genomeTxtFile = os.path.join(mainDir, '{sample}', 'syri', 'genomes.txt')
  output:
    plot = os.path.join(mainDir, '{sample}', 'syri', '{sample}_output_plot.png')
  shell: '''
    plotsr \
      --sr {input.outFile} \
      --genomes {input.genomeTxtFile} \
      -H 8 \
      -o {output.plot}
  '''

rule alignToChainOfGenomes:
  resources:
    mem_mb=50000 
  input:
    comp1 = os.path.join(mainDir, '{comp1}', '{comp1}_scaffold_chr_only.fa'),
    comp2 = os.path.join(mainDir, '{comp2}', '{comp2}_scaffold_chr_only.fa')
  params:
    chainDir = os.path.join(mainDir, 'chainAlign')
  output:
    compareSam = os.path.join(mainDir, 'chainAlign', '{comp1}-{comp2}comparison.sam')
  shell: '''
    mkdir -p {params.chainDir}

    minimap2 \
      -a \
      -x asm5 \
      --eqx \
      {input.comp1} \
      {input.comp2} \
      > {output.compareSam}
  '''

rule syntenyAnalysisForCombined:
  threads: 15
  input:
    comp1 = os.path.join(mainDir, '{comp1}', '{comp1}_scaffold_chr_only.fa'),
    comp2 = os.path.join(mainDir, '{comp2}', '{comp2}_scaffold_chr_only.fa'),
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

rule getChrLen:
  input:
    refChr2 = os.path.join(mainDir, 'resources', 'reference_genome_chr_only.fa')
  output:
    chrLen = os.path.join(mainDir, 'resources', 'chr.len')
  shell: '''
    seqkit split -i {input.refChr2}

    seqkit stats {input.refChr2}.split/*.fa -b -T -S \
      | awk -F'\\t' 'BEGIN {{OFS = FS}} {{ split($1, chr, ".") }} {{ print chr[2], $5 }}' \
      | awk -F'\\t' 'BEGIN {{OFS = FS}} {{ split($1, chr, "_") }} {{ print chr[2], $2 }}' \
      | tail -n +2 > {output.chrLen}
  '''

rule getChrLen1:
  input:
    refChr2 = os.path.join(mainDir, order_comps[0], order_comps[0] + '_scaffold_chr_only.fa')
  output:
    chrLen = os.path.join(mainDir, 'resources', 'chr_' + order_comps[0] + '.len')
  shell: '''
    seqkit split -i {input.refChr2}

    seqkit stats {input.refChr2}.split/*.fa -b -T -S \
      | awk -F'\\t' 'BEGIN {{OFS = FS}} {{ split($1, chr, ".") }} {{ print chr[2], $5 }}' \
      | awk -F'\\t' 'BEGIN {{OFS = FS}} {{ split($1, chr, "_") }} {{ print chr[2], $2 }}' \
      | tail -n +2 > {output.chrLen}
  '''

rule compareChrLen:
  input:
    chrLen = os.path.join(mainDir, 'resources', 'chr.len'),
    chrLen1 = os.path.join(mainDir, 'resources', 'chr_' + order_comps[0] + '.len')
  output:
    chrLenFinal = os.path.join(mainDir, 'resources', 'chrFinal.len')
  shell: '''
    echo "" > {output.chrLenFinal}

    lines=$(wc -l {input.chrLen} | awk '{{ print $1 }}')
    currentline=1
    while [ "$currentline" -le "$lines" ]
    do
      chrname=$(sed -n ${{currentline}}p {input.chrLen} | awk '{{ print $1 }}')
      reflen=$(sed -n ${{currentline}}p {input.chrLen} | awk '{{ print $2 }}')
      onelen=$(sed -n ${{currentline}}p {input.chrLen1} | awk '{{ print $2 }}')
      if [ "$reflen" -le "$onelen" ]; then
        echo "$chrname\t$onelen" >> {output.chrLenFinal}
        currentline=$(($currentline+1))
      else
        echo "$chrname\t$reflen" >> {output.chrLenFinal}
        currentline=$(($currentline+1))
      fi
    done
  '''

rule makeGenomeFileForPlotSRCombined:
  input:
    refChr2 = os.path.join(mainDir, 'resources', 'reference_genome_chr_only.fa'),
    chrLen = os.path.join(mainDir, 'resources', 'chrFinal.len')
  params:
    order_comps = order_comps,
    mainDir = mainDir
  output:
    genomeTxtFile = os.path.join(mainDir, 'syri_genomes.txt')
  shell: '''
    echo "#file  name" > {output.genomeTxtFile}
    echo "{input.refChr2}\treference\tlw:6" >> {output.genomeTxtFile}
    COMPS=({params.order_comps})
    for sample in ${{COMPS[@]}}
      do
        echo "{params.mainDir}/$sample/${{sample}}_scaffold_chr_only.fa\t$sample\tlw:6" >> {output.genomeTxtFile}
      done
  '''

rule plotStructuralRearrangementCombined:
  input:
    first_comp1_File = os.path.join(mainDir, order_comps[0], 'syri', 'reference_' + order_comps[0] + 'syri.out'),
    outFileChain = expand(os.path.join(mainDir, 'chainAlign', '{comp1}-{comp2}syri', '{comp1}-{comp2}syri.out'), zip, comp1 = first_comp, comp2 = second_comp),
    genomeTxtFile = os.path.join(mainDir, 'syri_genomes.txt')
  params:
    plotsrGraphs = os.path.join(mainDir, 'plotsrGraphs')
  output:
    plot = os.path.join(mainDir, 'plotsrGraphs', 'syri_{chr}_combined_output_plot.png')
  shell: '''
    mkdir -p {params.plotsrGraphs}

    COMPS=({input.outFileChain})
    SR=()
    SR+=("--sr {input.first_comp1_File}")
    for sample in ${{COMPS[@]}}
      do
        SR+=("--sr $sample ")
      done
    echo ${{SR[@]}}
    echo "plotsr \
      ${{SR[@]}} \
      --genomes {input.genomeTxtFile} \
      -H 8 \
      -W 12 \
      --chr {wildcards.chr} \
      -o {output.plot}" \
      | bash
  '''

rule plotStructuralRearrangementCombinedWithPOI:
  resources:
    mem_mb=50000 
  input:
    first_comp1_File = os.path.join(mainDir, order_comps[0], 'syri', 'reference_' + order_comps[0] + 'syri.out'),
    outFileChain = expand(os.path.join(mainDir, 'chainAlign', '{comp1}-{comp2}syri', '{comp1}-{comp2}syri.out'), zip, comp1 = first_comp, comp2 = second_comp),
    genomeTxtFile = os.path.join(mainDir, 'syri_genomes.txt'),
    pointsOfInterest = pointsOfInterest,
    tracksPOI = tracksPOI,
    chrLen = os.path.join(mainDir, 'resources', 'chrFinal.len')
  params:
    plotsrGraphsPOI = os.path.join(mainDir, 'plotsrGraphsPOI')
  output:
    plot = os.path.join(mainDir, 'plotsrGraphsPOI', 'syri_{chr}_combined_output_plot.png')
  shell: '''
    mkdir -p {params.plotsrGraphsPOI}

    COMPS=({input.outFileChain})
    SR{wildcards.chr}=()
    SR{wildcards.chr}+=("--sr {input.first_comp1_File}")
    for sample in ${{COMPS[@]}}
      do
        SR{wildcards.chr}+=("--sr $sample ")
      done
    echo ${{SR{wildcards.chr}[@]}}
    echo "plotsr \
      ${{SR{wildcards.chr}[@]}} \
      --genomes {input.genomeTxtFile} \
      --markers {input.pointsOfInterest} \
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
    control = os.path.join(mainDir, order_comps[0], order_comps[0] + '_scaffold_chr_only.fa'),
    ragtagChr2 = os.path.join(mainDir, '{comp_sample}', '{comp_sample}_scaffold_chr_only.fa')
  params:
    snifflesDir = os.path.join(mainDir, 'sniffles')
  output:
    compareToControl = os.path.join(mainDir, 'sniffles', order_comps[0] + '_{comp_sample}_comp.sam')
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

rule convertSamToBam:
  input:
     os.path.join(mainDir, 'sniffles', order_comps[0] + '_{comp_sample}_comp.sam')
  output:
    os.path.join(mainDir, 'sniffles', order_comps[0] + '_{comp_sample}_comp.bam')
  shell: '''
    samtools view -S -b {input} > {output}
  '''

rule sortBam:
  input:
    os.path.join(mainDir, 'sniffles', order_comps[0] + '_{comp_sample}_comp.bam')
  output:
    os.path.join(mainDir, 'sniffles', order_comps[0] + '_{comp_sample}_comp_sorted.bam')
  shell: '''
    samtools sort -o {output} {input}
  '''

rule indexBam:
  input:
    os.path.join(mainDir, 'sniffles', order_comps[0] + '_{comp_sample}_comp_sorted.bam')
  output:
    os.path.join(mainDir, 'sniffles', order_comps[0] + '_{comp_sample}_comp_sorted.bam.bai')
  shell: '''
    samtools index -b {input}
  '''

rule findSVs:
  resources:
    mem_mb=50000
  input:
    compareToControlBam = os.path.join(mainDir, 'sniffles', order_comps[0] + '_{comp_sample}_comp_sorted.bam'),
    bai = os.path.join(mainDir, 'sniffles', order_comps[0] + '_{comp_sample}_comp_sorted.bam.bai')
  output:
    snf = os.path.join(mainDir, 'sniffles', order_comps[0] + '_{comp_sample}.snf')
  shell: '''
    sniffles \
      --input {input.compareToControlBam} \
      --snf {output.snf}
  '''

rule compareSVs:
  resources:
    mem_mb=50000
  input:
    snf = expand(os.path.join(mainDir, 'sniffles', order_comps[0] + '_{comp_sample}.snf'), comp_sample = sniffles_comp)
  output:
    multiVCF = os.path.join(mainDir, 'sniffles', 'nonphaseSVs.vcf')
  shell: '''
    sniffles \
      --input {input.snf} \
      --vcf {output.multiVCF}
  '''