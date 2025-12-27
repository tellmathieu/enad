args <- commandArgs(TRUE)
wd <- args[1]
intersectfile <- args[2]
outputname1 <- args[3]

library(conflicted)
library(tidyverse)

conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")
library(dplyr)
library(gprofiler2)

#for testing purposes
#wd <- "/Users/tell/Desktop/eNAD"
#intersectfile <- "/Users/tell/Desktop/eNAD/test_intersect.txt"
#outputname <- "intersect1"

outputname <- gsub(" ", "", outputname1)
setwd(wd)

#load intersect file loaded from process
intersection <- read.delim(intersectfile,header=FALSE,sep="\t")

#isolate gene name to asearch using gprofiler
intersectiongene <- intersection %>% filter(V3=="gene")
intersectiongene[c('ID', 'DBXref','unk1','unk2','pregene','other')] <- str_split_fixed(intersectiongene$V9, ';', 6)
intersectiongene[c('pre', 'gene')] <- str_split_fixed(intersectiongene$pregene, '=', 2)
genes <- intersectiongene$gene
names(genes) <- genes



#running gost from gprofiler to get intersections and GO terms
gostres <- gost(query = genes,
                organism = "ecaballus",
                significant = FALSE,
                evcodes = TRUE)

#loading data into format that's easier for me to use -- data frame
intersect <- data.frame(source = gostres[["result"]][["source"]], 
                  term_name = gostres[["result"]][["term_name"]],
                  term_id = gostres[["result"]][["term_id"]],
                  p_value = gostres[["result"]][["p_value"]],
                  term_size = gostres[["result"]][["term_size"]],
                  query_size = gostres[["result"]][["query_size"]],
                  intersection_size = gostres[["result"]][["intersection_size"]],
                  effective_domain_size = gostres[["result"]][["effective_domain_size"]],
                  intersections = gostres[["result"]][["intersection"]]
                  )

#terms we're searching for
cns_terms <- c('spine','spinal','axon','dendritic','dendrite','glia','nervous','nerve','neuron','cholesterol','sterol','Vitamin E','bile','metabolism')


intersect_filtered <- data.frame()
for (term in cns_terms) {
  intersect_filter <- grepl(term, intersect$term_name, fixed = TRUE)
  intersect_term <- intersect %>% filter(intersect_filter == TRUE)
  intersect_filtered <- rbind(intersect_filtered, intersect_term)
}
intersect_filtered <- intersect_filtered %>% arrange(term_name)
intersect_filtered <- unique(intersect_filtered)
intersectlist <- intersection
intersect_terms <- unique(unlist(strsplit(intersect_filtered$intersections,',')))

intersect_svs_filtered <- data.frame()
for (termq in intersect_terms) {
  filter1 <- grepl(termq, intersectlist$V9, fixed = TRUE)
  term <- intersectlist %>% filter(filter1 == TRUE)
  term$X26 <- termq
  intersect_svs_filtered <- rbind(intersect_svs_filtered, term)
}


term_names <- data.frame()
for (term in intersect_terms) {
  intersect_termx <- intersect_filtered
  intersect_termx$forterm <- term
  BPlist <- intersect_termx %>% filter(intersect_termx$source == "GO:BP")
  BP <- grepl(term, BPlist$intersections, fixed = TRUE)
  BP_terms <- BPlist %>% select(forterm,term_name) %>% filter(BP == TRUE)
  BP_terms <- BP_terms %>% group_by(forterm) %>%
    summarise(BP=paste(term_name, collapse=', '))
  grouped_terms <- BP_terms
  #CC 
  CClist <- intersect_termx %>% filter(intersect_termx$source == "GO:CC")
  CC <- grepl(term, CClist$intersections, fixed = TRUE)
  CC_terms <- CClist %>% select(forterm,term_name) %>% filter(CC == TRUE)
  CC_terms <- CC_terms %>% group_by(forterm) %>%
    summarise(CC=paste(term_name, collapse=', '))
  grouped_terms <- merge(grouped_terms, CC_terms, by.x="forterm", by.y="forterm",all = TRUE)
  #HP 
  HPlist <- intersect_termx %>% filter(intersect_termx$source == "HP")
  HP <- grepl(term, HPlist$intersections, fixed = TRUE)
  HP_terms <- HPlist %>% select(forterm,term_name) %>% filter(HP == TRUE)
  HP_terms <- HP_terms %>% group_by(forterm) %>%
    summarise(HP=paste(term_name, collapse=', '))
  grouped_terms <- merge(grouped_terms, HP_terms, by.x="forterm", by.y="forterm",all = TRUE)
  term_names <- rbind(term_names,grouped_terms)
}

intersect_svs_filtered <- merge(intersect_svs_filtered,term_names, by.x='X26',by.y='forterm',all.y=TRUE)
colnames(intersect_svs_filtered) <- c('gene','seqname','source','feature','start','end','score','strand','frame','attribute','CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','horse_05_608','horse_04_56','horse_06_1212','horse_03_238','horse_01_200','horse_02_205','overlap','BP','CC','HP')
intersect_svs_filtered_genes <- intersect_svs_filtered %>% filter(feature == 'gene')
intersect_vcf <- intersect_svs_filtered_genes[c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','horse_05_608','horse_04_56','horse_06_1212','horse_03_238','horse_01_200','horse_02_205')]
write_delim(intersect_vcf,paste('ref-results/',outputname, '_filter.vcf',sep = ""),delim='\t',col_names = FALSE)
write_csv(intersect_svs_filtered,paste('ref-results/',outputname, "_svs.csv",sep = ""))
write_csv(intersect_svs_filtered_genes,paste('ref-results/',outputname, "_svs_genes.csv",sep = ""))

#intersect unfiltered
intersect_terms_unfiltered <- unique(unlist(strsplit(intersect$intersections,',')))
intersect_terms_unfiltered <- unique(intersect_terms_unfiltered)

intersect_svs_unfiltered <- data.frame()
for (termq in intersect_terms_unfiltered) {
  filter1 <- grepl(termq, intersectlist$V9, fixed = TRUE)
  term <- intersectlist %>% filter(filter1 == TRUE)
  term$X26 <- termq
  intersect_svs_unfiltered <- rbind(intersect_svs_unfiltered, term)
}


term_names_unfiltered <- data.frame()
for (term in intersect_terms_unfiltered) {
  intersectx <- intersect
  intersectx$forterm <- term
  BPlist <- intersectx %>% filter(intersectx$source == "GO:BP")
  BP <- grepl(term, BPlist$intersections, fixed = TRUE)
  BP_terms <- BPlist %>% select(forterm,term_name) %>% filter(BP == TRUE)
  BP_terms <- BP_terms %>% group_by(forterm) %>%
    summarise(BP=paste(term_name, collapse=', '))
  grouped_terms <- BP_terms
  #CC 
  CClist <- intersectx %>% filter(intersectx$source == "GO:CC")
  CC <- grepl(term, CClist$intersections, fixed = TRUE)
  CC_terms <- CClist %>% select(forterm,term_name) %>% filter(CC == TRUE)
  CC_terms <- CC_terms %>% group_by(forterm) %>%
    summarise(CC=paste(term_name, collapse=', '))
  grouped_terms <- merge(grouped_terms, CC_terms, by.x="forterm", by.y="forterm",all = TRUE)
  #HP 
  HPlist <- intersectx %>% filter(intersectx$source == "HP")
  HP <- grepl(term, HPlist$intersections, fixed = TRUE)
  HP_terms <- HPlist %>% select(forterm,term_name) %>% filter(HP == TRUE)
  HP_terms <- HP_terms %>% group_by(forterm) %>%
    summarise(HP=paste(term_name, collapse=', '))
  grouped_terms <- merge(grouped_terms, HP_terms, by.x="forterm", by.y="forterm",all = TRUE)
  term_names_unfiltered <- rbind(term_names_unfiltered,grouped_terms)
}

intersect_svs_unfiltered <- merge(intersect_svs_unfiltered,term_names_unfiltered, by.x='X26',by.y='forterm',all.y=TRUE)
colnames(intersect_svs_unfiltered) <- c('gene','seqname','source','feature','start','end','score','strand','frame','attribute','CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','horse_05_608','horse_04_56','horse_06_1212','horse_03_238','horse_01_200','horse_02_205','overlap','BP','CC','HP')
intersect_svs_unfiltered_genes <- intersect_svs_unfiltered %>% filter(feature == 'gene')
intersect_vcf_unfiltered <- intersect_svs_unfiltered_genes[c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','horse_05_608','horse_04_56','horse_06_1212','horse_03_238','horse_01_200','horse_02_205')]
intersect_svs_unfiltered_exons <- intersect_svs_unfiltered %>% filter(feature == 'exon')
write_delim(intersect_vcf_unfiltered,paste('ref-results/',outputname, '_unfiltered.vcf',sep = ""),delim='\t',col_names = FALSE)
write_csv(intersect_svs_unfiltered,paste('ref-results/',outputname, "_svs_unfiltered.csv",sep = ""))
write_csv(intersect_svs_unfiltered_genes,paste('ref-results/',outputname, "_unfiltered.csv",sep = ""))
write_csv(intersect_svs_unfiltered_exons,paste('ref-results/',outputname, "_svs_exons_unfiltered.csv",sep = ""))

#getting bed file
intersect_bed <- data.frame()
for (i in 1:(nrow(intersect_svs_filtered_genes))) {
  chr <- intersect_svs_filtered_genes$CHROM[i]
  start <- intersect_svs_filtered_genes$POS[i] - 100000
  info_group <- unlist(strsplit(intersect_svs_filtered_genes$INFO[i],';'))
  end <- unlist(strsplit(info_group[grep('END=',info_group,fixed=TRUE)],'='))[2]
  end <- ifelse(is.null(end), start, as.numeric(end)) + 100000
  name <- paste(intersect_svs_filtered_genes$ID[i],intersect_svs_filtered_genes$gene[i],sep=';')
  bed <- data.frame(chr=chr,
                    start=start,
                    end=end,
                    name=name)
  intersect_bed <- rbind(intersect_bed,bed)
}
write_delim(intersect_bed,paste('ref-results/',outputname,'_candidates_expanded.bed',sep = ""),delim='\t',col_names = FALSE)

#getting unfiltered bed file
intersect_bed_unfiltered <- data.frame()
for (i in 1:(nrow(intersect_svs_unfiltered_genes))) {
  chr <- intersect_svs_unfiltered_genes$CHROM[i]
  start <- intersect_svs_unfiltered_genes$POS[i] - 100000
  info_group <- unlist(strsplit(intersect_svs_unfiltered_genes$INFO[i],';'))
  end <- unlist(strsplit(info_group[grep('END=',info_group,fixed=TRUE)],'='))[2]
  end <- ifelse(is.null(end), start, as.numeric(end)) + 100000
  name <- paste(intersect_svs_unfiltered_genes$ID[i],intersect_svs_unfiltered_genes$gene[i],sep=';')
  bed <- data.frame(chr=chr,
                    start=start,
                    end=end,
                    name=name)
  intersect_bed_unfiltered <- rbind(intersect_bed_unfiltered,bed)
}
write_delim(intersect_bed_unfiltered,paste('ref-results/',outputname, '_candidates_expanded_unfiltered.bed',sep = ""),delim='\t',col_names = FALSE)


