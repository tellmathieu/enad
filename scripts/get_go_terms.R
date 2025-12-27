# This script

#args<-commandArgs(TRUE)

#wd <- args[0]
#intersectfile <- args[1]

#install.packages("gprofiler2")

#setwd(wd)
#intersectfile <- 'test_intersect.txt'
library(dplyr)
library(tidyverse)
library(openxlsx)
library(gprofiler2)
#install.packages("openxlsx")



#gostres <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"),
#                organism = "hsapiens")

#cns_terms <- c('spine','spinal','axon','dendritic','dendrite','glia','nervous','nerve','neuron','cholesterol','sterol','Vitamin E','bile','metabolism')

#hap1
hap1 <- read.csv("hap1_gProfiler_ecaballus_3-26-2025_9-54-08 AM__intersections.csv", header=TRUE)
hap1_filtered <- data.frame()
for (term in cns_terms) {
  hap1_filter <- grepl(term, hap1$term_name, fixed = TRUE)
  hap1_term <- hap1 %>% filter(hap1_filter == TRUE)
  hap1_filtered <- rbind(hap1_filtered, hap1_term)
}
hap1_filtered <- hap1_filtered %>% arrange(term_name)
hap1_filtered <- unique(hap1_filtered)
hap1list <- read.xlsx('phasedSVs.compto5-1_hap1_goi.xlsx',
                      sheet='phasedSVs.compto5-1_hap1_vcf_gf',
                      colNames = FALSE)
hap1_terms <- unique(unlist(strsplit(hap1_filtered$intersections,',')))

hap1_svs_filtered <- data.frame()
for (termq in hap1_terms) {
    filter1 <- grepl(termq, hap1list$X9, fixed = TRUE)
    term <- hap1list %>% filter(filter1 == TRUE)
    term$X26 <- termq
    hap1_svs_filtered <- rbind(hap1_svs_filtered, term)
  }


term_names <- data.frame()
for (term in hap1_terms) {
  hap1_termx <- hap1_filtered
  hap1_termx$forterm <- term
  BPlist <- hap1_termx %>% filter(hap1_termx$source == "GO:BP")
  BP <- grepl(term, BPlist$intersections, fixed = TRUE)
  BP_terms <- BPlist %>% select(forterm,term_name) %>% filter(BP == TRUE)
  BP_terms <- BP_terms %>% group_by(forterm) %>%
    summarise(BP=paste(term_name, collapse=', '))
  grouped_terms <- BP_terms
  #CC 
  CClist <- hap1_termx %>% filter(hap1_termx$source == "GO:CC")
  CC <- grepl(term, CClist$intersections, fixed = TRUE)
  CC_terms <- CClist %>% select(forterm,term_name) %>% filter(CC == TRUE)
  CC_terms <- CC_terms %>% group_by(forterm) %>%
    summarise(CC=paste(term_name, collapse=', '))
  grouped_terms <- merge(grouped_terms, CC_terms, by.x="forterm", by.y="forterm",all = TRUE)
  #HP 
  HPlist <- hap1_termx %>% filter(hap1_termx$source == "HP")
  HP <- grepl(term, HPlist$intersections, fixed = TRUE)
  HP_terms <- HPlist %>% select(forterm,term_name) %>% filter(HP == TRUE)
  HP_terms <- HP_terms %>% group_by(forterm) %>%
    summarise(HP=paste(term_name, collapse=', '))
  grouped_terms <- merge(grouped_terms, HP_terms, by.x="forterm", by.y="forterm",all = TRUE)
  term_names <- rbind(term_names,grouped_terms)
}

hap1_svs_filtered <- merge(hap1_svs_filtered,term_names, by.x='X26',by.y='forterm',all.y=TRUE)
colnames(hap1_svs_filtered) <- c('gene','seqname','source','feature','start','end','score','strand','frame','attribute','CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','horse_05_608','horse_04_56','horse_06_1212','horse_03_238','horse_01_200','horse_02_205','overlap','BP','CC','HP')
hap1_svs_filtered_genes <- hap1_svs_filtered %>% filter(feature == 'gene')
hap1_vcf <- hap1_svs_filtered_genes[c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','horse_05_608','horse_04_56','horse_06_1212','horse_03_238','horse_01_200','horse_02_205')]
write_delim(hap1_vcf,'hap1_filter.vcf',delim='\t',col_names = FALSE)
write_csv(hap1_svs_filtered,"hap1_svs.csv")
write_csv(hap1_svs_filtered_genes,"hap1_svs_genes.csv")

#hap1 unfiltered
hap1_terms_unfiltered <- unique(unlist(strsplit(hap1$intersections,',')))
hap1_terms_unfiltered <- unique(hap1_terms_unfiltered)

hap1_svs_unfiltered <- data.frame()
for (termq in hap1_terms_unfiltered) {
  filter1 <- grepl(termq, hap1list$X9, fixed = TRUE)
  term <- hap1list %>% filter(filter1 == TRUE)
  term$X26 <- termq
  hap1_svs_unfiltered <- rbind(hap1_svs_unfiltered, term)
}


term_names_unfiltered <- data.frame()
for (term in hap1_terms_unfiltered) {
  hap1x <- hap1
  hap1x$forterm <- term
  BPlist <- hap1x %>% filter(hap1x$source == "GO:BP")
  BP <- grepl(term, BPlist$intersections, fixed = TRUE)
  BP_terms <- BPlist %>% select(forterm,term_name) %>% filter(BP == TRUE)
  BP_terms <- BP_terms %>% group_by(forterm) %>%
    summarise(BP=paste(term_name, collapse=', '))
  grouped_terms <- BP_terms
  #CC 
  CClist <- hap1x %>% filter(hap1x$source == "GO:CC")
  CC <- grepl(term, CClist$intersections, fixed = TRUE)
  CC_terms <- CClist %>% select(forterm,term_name) %>% filter(CC == TRUE)
  CC_terms <- CC_terms %>% group_by(forterm) %>%
    summarise(CC=paste(term_name, collapse=', '))
  grouped_terms <- merge(grouped_terms, CC_terms, by.x="forterm", by.y="forterm",all = TRUE)
  #HP 
  HPlist <- hap1x %>% filter(hap1x$source == "HP")
  HP <- grepl(term, HPlist$intersections, fixed = TRUE)
  HP_terms <- HPlist %>% select(forterm,term_name) %>% filter(HP == TRUE)
  HP_terms <- HP_terms %>% group_by(forterm) %>%
    summarise(HP=paste(term_name, collapse=', '))
  grouped_terms <- merge(grouped_terms, HP_terms, by.x="forterm", by.y="forterm",all = TRUE)
  term_names_unfiltered <- rbind(term_names_unfiltered,grouped_terms)
}

hap1_svs_unfiltered <- merge(hap1_svs_unfiltered,term_names_unfiltered, by.x='X26',by.y='forterm',all.y=TRUE)
colnames(hap1_svs_unfiltered) <- c('gene','seqname','source','feature','start','end','score','strand','frame','attribute','CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','horse_05_608','horse_04_56','horse_06_1212','horse_03_238','horse_01_200','horse_02_205','overlap','BP','CC','HP')
hap1_svs_unfiltered_genes <- hap1_svs_unfiltered %>% filter(feature == 'gene')
hap1_vcf_unfiltered <- hap1_svs_unfiltered_genes[c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','horse_05_608','horse_04_56','horse_06_1212','horse_03_238','horse_01_200','horse_02_205')]
hap1_svs_unfiltered_exons <- hap1_svs_unfiltered %>% filter(feature == 'exon')
write_delim(hap1_vcf_unfiltered,'hap1_unfiltered.vcf',delim='\t',col_names = FALSE)
write_csv(hap1_svs_unfiltered,"hap1_svs_unfiltered.csv")
write_csv(hap1_svs_unfiltered_genes,"hap1_svs_genes_unfiltered.csv")
write_csv(hap1_svs_unfiltered_exons,"hap1_svs_exons_unfiltered.csv")




#hap2
hap2 <- read.csv("hap2_gProfiler_ecaballus_3-26-2025_9-57-50 AM__intersections.csv", header=TRUE)
hap2_filtered <- data.frame()
for (term in cns_terms) {
  hap2_filter <- grepl(term, hap2$term_name, fixed = TRUE)
  hap2_term <- hap2 %>% filter(hap2_filter == TRUE)
  hap2_filtered <- rbind(hap2_filtered, hap2_term)
}
hap2_filtered <- hap2_filtered %>% arrange(term_name)
hap2_filtered <- unique(hap2_filtered)
hap2list <- read.xlsx('phasedSVs.compto5-2_hap2_goi.xlsx',
                      sheet='phasedSVs.compto5-2_hap2_vcf_gf',
                      colNames = FALSE)
hap2_terms <- unique(unlist(strsplit(hap2_filtered$intersections,',')))
hap2_svs_filtered <- data.frame()
for (termq in hap2_terms) {
  filter1 <- grepl(termq, hap2list$X9, fixed = TRUE)
  term <- hap2list %>% filter(filter1 == TRUE)
  term$X26 <- termq
  hap2_svs_filtered <- rbind(hap2_svs_filtered, term)
}


term_names <- data.frame()
for (term in hap2_terms) {
  hap2_termx <-  hap2_filtered
  hap2_termx$forterm <- term
  BPlist <- hap2_termx %>% filter(hap2_termx$source == "GO:BP")
  BP <- grepl(term, BPlist$intersections, fixed = TRUE)
  BP_terms <- BPlist %>% select(forterm,term_name) %>% filter(BP == TRUE)
  BP_terms <- BP_terms %>% group_by(forterm) %>%
    summarise(BP=paste(term_name, collapse=', '))
  grouped_terms <- BP_terms
  #CC 
  CClist <- hap2_termx %>% filter(hap2_termx$source == "GO:CC")
  CC <- grepl(term, CClist$intersections, fixed = TRUE)
  CC_terms <- CClist %>% select(forterm,term_name) %>% filter(CC == TRUE)
  CC_terms <- CC_terms %>% group_by(forterm) %>%
    summarise(CC=paste(term_name, collapse=', '))
  grouped_terms <- merge(grouped_terms, CC_terms, by.x="forterm", by.y="forterm",all = TRUE)
  #HP 
  HPlist <- hap2_termx %>% filter(hap2_termx$source == "HP")
  HP <- grepl(term, HPlist$intersections, fixed = TRUE)
  HP_terms <- HPlist %>% select(forterm,term_name) %>% filter(HP == TRUE)
  HP_terms <- HP_terms %>% group_by(forterm) %>%
    summarise(HP=paste(term_name, collapse=', '))
  grouped_terms <- merge(grouped_terms, HP_terms, by.x="forterm", by.y="forterm",all = TRUE)
  term_names <- rbind(term_names,grouped_terms)
}

hap2_svs_filtered <- merge(hap2_svs_filtered,term_names, by.x='X26',by.y='forterm',all.y=TRUE)
colnames(hap2_svs_filtered) <- c('gene','seqname','source','feature','start','end','score','strand','frame','attribute','CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','horse_05_608','horse_04_56','horse_06_1212','horse_03_238','horse_01_200','horse_02_205','overlap','BP','CC','HP')
hap2_svs_filtered_genes <- hap2_svs_filtered %>% filter(feature == 'gene')
hap2_vcf <- hap2_svs_filtered_genes[c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','horse_05_608','horse_04_56','horse_06_1212','horse_03_238','horse_01_200','horse_02_205')]
write_delim(hap2_vcf,'hap2_filter.vcf',delim='\t',col_names = FALSE)
write_csv(hap2_svs_filtered,"hap2_svs.csv")
write_csv(hap2_svs_filtered_genes,"hap2_svs_genes.csv")


#hap2 unfiltered
hap2_terms_unfiltered <- unique(unlist(strsplit(hap2$intersections,',')))
hap2_terms_unfiltered <- unique(hap2_terms_unfiltered)

hap2_svs_unfiltered <- data.frame()
for (termq in hap2_terms_unfiltered) {
  filter1 <- grepl(termq, hap2list$X9, fixed = TRUE)
  term <- hap2list %>% filter(filter1 == TRUE)
  term$X26 <- termq
  hap2_svs_unfiltered <- rbind(hap2_svs_unfiltered, term)
}


term_names_unfiltered <- data.frame()
for (term in hap2_terms_unfiltered) {
  hap2x <- hap2
  hap2x$forterm <- term
  BPlist <- hap2x %>% filter(hap2x$source == "GO:BP")
  BP <- grepl(term, BPlist$intersections, fixed = TRUE)
  BP_terms <- BPlist %>% select(forterm,term_name) %>% filter(BP == TRUE)
  BP_terms <- BP_terms %>% group_by(forterm) %>%
    summarise(BP=paste(term_name, collapse=', '))
  grouped_terms <- BP_terms
  #CC 
  CClist <- hap2x %>% filter(hap2x$source == "GO:CC")
  CC <- grepl(term, CClist$intersections, fixed = TRUE)
  CC_terms <- CClist %>% select(forterm,term_name) %>% filter(CC == TRUE)
  CC_terms <- CC_terms %>% group_by(forterm) %>%
    summarise(CC=paste(term_name, collapse=', '))
  grouped_terms <- merge(grouped_terms, CC_terms, by.x="forterm", by.y="forterm",all = TRUE)
  #HP 
  HPlist <- hap2x %>% filter(hap2x$source == "HP")
  HP <- grepl(term, HPlist$intersections, fixed = TRUE)
  HP_terms <- HPlist %>% select(forterm,term_name) %>% filter(HP == TRUE)
  HP_terms <- HP_terms %>% group_by(forterm) %>%
    summarise(HP=paste(term_name, collapse=', '))
  grouped_terms <- merge(grouped_terms, HP_terms, by.x="forterm", by.y="forterm",all = TRUE)
  term_names_unfiltered <- rbind(term_names_unfiltered,grouped_terms)
}

hap2_svs_unfiltered <- merge(hap2_svs_unfiltered,term_names_unfiltered, by.x='X26',by.y='forterm',all.y=TRUE)
colnames(hap2_svs_unfiltered) <- c('gene','seqname','source','feature','start','end','score','strand','frame','attribute','CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','horse_05_608','horse_04_56','horse_06_1212','horse_03_238','horse_01_200','horse_02_205','overlap','BP','CC','HP')
hap2_svs_unfiltered_genes <- hap2_svs_unfiltered %>% filter(feature == 'gene')
hap2_vcf_unfiltered <- hap2_svs_unfiltered_genes[c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','horse_05_608','horse_04_56','horse_06_1212','horse_03_238','horse_01_200','horse_02_205')]
hap2_svs_unfiltered_exons <- hap2_svs_unfiltered %>% filter(feature == 'exon')
write_delim(hap2_vcf_unfiltered,'hap2_unfiltered.vcf',delim='\t',col_names = FALSE)
write_csv(hap2_svs_unfiltered,"hap2_svs_unfiltered.csv")
write_csv(hap2_svs_unfiltered_genes,"hap2_svs_genes_unfiltered.csv")
write_csv(hap2_svs_unfiltered_exons,"hap2_svs_exons_unfiltered.csv")


#combined
hap1_svs_filtered_genes$hap <- 'hap1'
hap2_svs_filtered_genes$hap <- 'hap2'
combined_svs <- rbind(hap1_svs_filtered_genes,hap2_svs_filtered_genes)
combined_svs <- combined_svs %>% arrange(gene)
duplicated_svs <- combined_svs[combined_svs$gene %in% combined_svs$gene[duplicated(combined_svs$gene)],]
colnames(duplicated_svs) <- c('gene','seqname','source','feature','start','end','score','strand','frame','attribute','CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','horse_05_608','horse_04_56','horse_06_1212','horse_03_238','horse_01_200','horse_02_205','overlap','BP','CC','HP','hap')
write_csv(duplicated_svs,"duplicated_svs.csv")

#getting bed files
hap1_bed <- data.frame()
for (i in 1:(nrow(hap1_svs_filtered_genes))) {
  print(i)
  chr <- hap1_svs_filtered_genes$CHROM[i]
  start <- hap1_svs_filtered_genes$POS[i] - 100000
  info_group <- unlist(strsplit(hap1_svs_filtered_genes$INFO[i],';'))
  end <- unlist(strsplit(info_group[grep('END=',info_group,fixed=TRUE)],'='))[2]
  end <- as.numeric(end) + 100000
  name <- paste(hap1_svs_filtered_genes$ID[i],hap1_svs_filtered_genes$gene[i],sep=';')
  bed <- data.frame(chr=chr,
    start=start,
    end=end,
    name=name)
  hap1_bed <- rbind(hap1_bed,bed)
}
write_delim(hap1_bed,'hap1_candidates_expanded.bed',delim='\t',col_names = FALSE)

hap2_bed <- data.frame()
for (i in 1:(nrow(hap2_svs_filtered_genes))) {
  chr <- hap2_svs_filtered_genes$CHROM[i]
  start <- hap2_svs_filtered_genes$POS[i] - 100000
  info_group <- unlist(strsplit(hap2_svs_filtered_genes$INFO[i],';'))
  end <- unlist(strsplit(info_group[grep('END=',info_group,fixed=TRUE)],'='))[2]
  end <- as.numeric(end) + 100000
  name <- paste(hap2_svs_filtered_genes$ID[i],hap2_svs_filtered_genes$gene[i],sep=';')
  bed <- data.frame(chr=chr,
                    start=start,
                    end=end,
                    name=name)
  hap2_bed <- rbind(hap2_bed,bed)
}
write_delim(hap2_bed,'hap2_candidates_expanded.bed',delim='\t',col_names = FALSE)


#getting bed files
hap1_bed_unfiltered <- data.frame()
for (i in 1:(nrow(hap1_svs_unfiltered_genes))) {
  chr <- hap1_svs_unfiltered_genes$CHROM[i]
  start <- hap1_svs_unfiltered_genes$POS[i] - 100000
  info_group <- unlist(strsplit(hap1_svs_unfiltered_genes$INFO[i],';'))
  end <- unlist(strsplit(info_group[grep('END=',info_group,fixed=TRUE)],'='))[2]
  end <- ifelse(is.null(end), start, as.numeric(end)) + 100000
  name <- paste(hap1_svs_unfiltered_genes$ID[i],hap1_svs_unfiltered_genes$gene[i],sep=';')
  bed <- data.frame(chr=chr,
                    start=start,
                    end=end,
                    name=name)
  hap1_bed_unfiltered <- rbind(hap1_bed_unfiltered,bed)
}
write_delim(hap1_bed_unfiltered,'hap1_candidates_expanded_unfiltered.bed',delim='\t',col_names = FALSE)

hap2_bed_unfiltered <- data.frame()
for (i in 1:(nrow(hap2_svs_unfiltered_genes))) {
  chr <- hap2_svs_unfiltered_genes$CHROM[i]
  start <- hap2_svs_unfiltered_genes$POS[i] - 100000
  info_group <- unlist(strsplit(hap2_svs_unfiltered_genes$INFO[i],';'))
  end <- unlist(strsplit(info_group[grep('END=',info_group,fixed=TRUE)],'='))[2]
  end <- ifelse(is.null(end), start, as.numeric(end)) + 100000
  name <- paste(hap2_svs_unfiltered_genes$ID[i],hap2_svs_unfiltered_genes$gene[i],sep=';')
  bed <- data.frame(chr=chr,
                    start=start,
                    end=end,
                    name=name)
  hap2_bed_unfiltered <- rbind(hap2_bed_unfiltered,bed)
}
write_delim(hap2_bed_unfiltered,'hap2_candidates_expanded_unfiltered.bed',delim='\t',col_names = FALSE)

