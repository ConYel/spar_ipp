## load libraries -----
library(plyranges)
library(tidyverse)
## load pirnadb -----
piRNA_dbs_files <- list.files("/home/0/piRNA_DBs", full.names = TRUE)
# piRBase
pirbase <- piRNA_dbs_files[2] %>% 
  read_bed() %>% 
  as_tibble() %>% 
  rename(pirbase = "name") %>% 
  select(-score) %>% 
  mutate(pirbase_coor = str_c(.$pirbase, .$seqnames, .$start, .$end, .$strand, sep = "_")) %>% 
  as_granges()
# piRNADB
pirnadb <- piRNA_dbs_files[5] %>% 
  read_tsv(comment = "#", col_names = c("seqnames", "X2", "x3", "start", "end", "x6", "strand", "x7", "piRNAdb" )) %>% 
  select(-X2, -x3, -x6, -x7) %>% 
  mutate(seqnames = if_else(seqnames == "chrMT", "chrM", as.character(seqnames))) %>% 
  mutate(pirnadb_coor = str_c(.$piRNAdb, .$seqnames, .$start, .$end, .$strand, sep = "_")) %>% 
  as_granges() %>% 
  arrange(start)
# piRNADB cluster
pirnadb_cl <- piRNA_dbs_files[4] %>% 
  read_tsv() %>% 
  select(-Build_Code) %>% 
  mutate(Strand = if_else(Strand == "biDirectional","*", Strand),
                                                Chromosome = str_replace(.$Chromosome,"^", "chr")) %>% 
  as_granges(seqnames = Chromosome,
             start = Start,
             end = End, 
             strand = Strand
             ) %>% 
  arrange(start) %>% 
  keepStandardChromosomes(pruning.mode="coarse")
# dashr DB cluster
dashr_db <- piRNA_dbs_files[1] %>% 
  read_tsv(col_names = c("seqnames", 
                         "start", 
                         "end", 
                         "dashr_srna", 
                         "dashr_type",
                         "strand")) %>% 
  mutate(dashr_srna_coor = str_c(.$dashr_srna, .$seqnames, .$start, .$end, .$strand, sep = "_")) %>% 
  as_granges %>% 
  arrange(start)
# cluster db 
pirna_cl_db <- piRNA_dbs_files[3] %>% 
  read_gff() %>% 
  as_tibble %>% 
  mutate(score = 1:length(.$score),
         seqnames = str_replace(.$seqnames,"^", "chr")) %>% 
  unite(cl_db, type:score) %>% 
  select(-group, -phase, -source) %>% 
  as_granges() %>% 
  arrange(start) %>% 
  keepStandardChromosomes(pruning.mode="coarse")
# check all databases to find overlaps
pirna_cl_db %>% join_overlap_intersect_directed(pirnadb_cl)
#overlaps_3db <- dashr_db %>% 
  join_overlap_intersect_directed(pirnadb, minoverlap = 15) %>% 
  join_overlap_intersect_directed(pirbase, minoverlap = 15)

# create a union of all dbs ------
dashr_db_piRNA <- dashr_db %>% filter(dashr_type == "piRNA") %>% select(-dashr_type)

dashr_db_piRNA_red <- dashr_db_piRNA %>% reduce_ranges_directed()

pirbase_red <- pirbase %>% reduce_ranges_directed()

pirnadb_red <- pirnadb %>% reduce_ranges_directed()

# concat them and reduce
pirna_DB_union <- c(dashr_db_piRNA_red, pirbase_red, pirnadb_red) %>% 
  reduce_ranges_directed()
#rm(pirbase_red, pirnadb_red, dashr_db_piRNA_red)

# make the file for the 3 dbds
#p_DB_U_dashr <- pirna_DB_union %>% 
#  join_overlap_left_directed(dashr_db_piRNA) %>% 
#  filter(!is.na(dashr_type)) %>% select(-dashr_type)

p_DB_U_piRBase <- pirna_DB_union %>% 
  join_overlap_left_directed(pirbase) %>% 
  filter(!is.na(pirbase)) %>% select(-score)

p_DB_U_pirnaDB <- pirna_DB_union %>% 
  join_overlap_left_directed(pirnadb) %>% 
  filter(!is.na(piRNAdb)) 

# create the complete file
pirna_DB_union_1 <- pirna_DB_union %>% 
  as_tibble() %>% 
  left_join(as_tibble(p_DB_U_dashr))

pirna_DB_union_1 <- pirna_DB_union_1 %>% 
  left_join(as_tibble(p_DB_U_pirnaDB))

#pirna_DB_union_1 <- pirna_DB_union_1 %>% not working probably overflow
  left_join(as_tibble(p_DB_U_piRBase)) 

# filter for each db sequences are included in genes ----
# keep the primary chromosomes
#pirbase
pirbase_red %>% 
    keepStandardChromosomes(pruning.mode="coarse") %>% 
    as_tibble() %>% 
  group_by(seqnames) %>% 
  summarise_at(vars(width) ,list(min = min, Q1=~quantile(., probs = 0.25),
                 median=median, Q3=~quantile(., probs = 0.75),
                 max=max)) %>% 
  arrange(as.character(seqnames)) %>% tail(20)
# evaluating summary statistics we find that most of the regions
# that are around ~34 base pairs for that reason we will remove 
# all regions bigger than 39

pirbase_short <- pirbase_red %>% 
  keepStandardChromosomes(pruning.mode="coarse") %>% 
  as_tibble() %>% 
  filter(width < 40)
#dashr
dashr_db_piRNA_red %>% 
  keepStandardChromosomes(pruning.mode="coarse") %>% 
  as_tibble() %>% 
  group_by(seqnames) %>% 
  summarise_at(vars(width) ,list(min = min, Q1=~quantile(., probs = 0.25),
                                 median=median, Q3=~quantile(., probs = 0.75),
                                 max=max)) %>% 
  arrange(as.character(seqnames)) %>% tail(20)

dashr_db_piRNA_short <- dashr_db_piRNA_red %>% 
  keepStandardChromosomes(pruning.mode="coarse") %>% 
  as_tibble() %>% 
  filter(width < 40)
#pirnadb
pirnadb_red %>% 
  keepStandardChromosomes(pruning.mode="coarse") %>% 
  as_tibble() %>% 
  group_by(seqnames) %>% 
  summarise_at(vars(width) ,list(min = min, Q1=~quantile(., probs = 0.25),
                                 median=median, Q3=~quantile(., probs = 0.75),
                                 max=max)) %>% 
  arrange(as.character(seqnames)) %>% tail(20)

pirnadb_short <- pirnadb_red %>% 
  keepStandardChromosomes(pruning.mode="coarse") %>% 
  as_tibble() %>% 
  filter(width < 40)
# concat them and reduce
pirna_DB_union <- c(as_granges(dashr_db_piRNA_short), as_granges(pirbase_short), as_granges(pirnadb_short)) %>% 
  reduce_ranges_directed() %>%  
  keepStandardChromosomes(pruning.mode="coarse") %>% 
  as_tibble() %>%  
  mutate(rnaID = 1:length(.$strand) %>% str_c("RNA", .)) %>% 
  as_granges()

pirna_DB_union %>% as_tibble() %>% 
  group_by(seqnames) %>% 
  summarise_at(vars(width) ,list(min = min, Q1=~quantile(., probs = 0.25),
                                 median=median, Q3=~quantile(., probs = 0.75),
                                 max=max)) %>% 
  arrange(as.character(seqnames)) %>% tail(20)
# make the file for the 3 dbds
p_DB_U <- pirna_DB_union %>% 
  join_overlap_left_directed(dashr_db_piRNA) %>% 
  join_overlap_left_directed(pirnadb) %>% 
  join_overlap_left_directed(pirbase) %>% 
  join_overlap_left_directed(pirnadb_cl) %>%
  join_overlap_left_directed(pirna_cl_db)
# find gene / transcript regions -----
library(TxDb.Hsapiens.UCSC.hg38.knownGene);library(org.Hs.eg.db)
txd <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(bumphunter)
genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg38.knownGene,annotation="org.Hs.eg.db") %>% keepStandardChromosomes(pruning.mode="coarse")
#annotated_peaks <- matchGenes(p_DB_U, genes, type = "any", promoterDist = 2500, skipExons = FALSE, verbose = TRUE) %>% as_tibble()
#p_DB_U_chr1 <- pirna_DB_union %>% filter(seqnames == "chr1")
p_DB_U_chr <- map(seqlevels(pirna_DB_union), ~ pirna_DB_union  %>% 
                     filter(seqnames == .x))

annotated_peaks_chr11 <- matchGenes(p_DB_U_chr11, genes, type = "any", promoterDist = 2500, skipExons = FALSE, verbose = TRUE) %>% as_tibble()
chrY <- 
test1 <- list(p_DB_U_chr[[25]][1:100],  p_DB_U_chr[[23]])

library(BiocParallel)
SerialParam()
MulticoreParam()
mt_param <- MulticoreParam(workers = 6)
fun <- function(v){
  message("working")
  matchGenes(v, genes, type = "any", promoterDist = 2500, skipExons = FALSE, verbose = TRUE) %>% as_tibble()
}

test_genes_1 <- bplapply(test1, fun, BPPARAM = mt_param)  
test_genes_1 <- test_genes_1 %>% bind_rows()

test1 <- bind_ranges(test1)

test_genes_1 <- test_genes_1 %>% 
  bind_rows() %>% 
  bind_cols(as_tibble(test1)) %>% 
  dplyr::select(name:subjectHits,rnaID)
# find how many threads you are going to use -----
length(p_DB_U_chr[[23]]) %/% 15

(l_chr <- plyr::round_any(length(p_DB_U_chr[[23]]), 10, ceiling)) 
(l_chr %/% 14)

(chunks_chr <- seq(1,l_chr, by = l_chr %/% 14))

p_DB_U_chr_1 <- p_DB_U_chr[[23]][chunks_chr[1]:(chunks_chr[2]-1)]
p_DB_U_chr_2 <- p_DB_U_chr[[23]][chunks_chr[2]:(chunks_chr[3]-1)]
p_DB_U_chr_3 <- p_DB_U_chr[[23]][chunks_chr[3]:(chunks_chr[4]-1)]
p_DB_U_chr_4 <- p_DB_U_chr[[23]][chunks_chr[4]:(chunks_chr[5]-1)]
p_DB_U_chr_5 <- p_DB_U_chr[[23]][chunks_chr[5]:(chunks_chr[6]-1)]
p_DB_U_chr_6 <- p_DB_U_chr[[23]][chunks_chr[6]:(chunks_chr[7]-1)]
p_DB_U_chr_7 <- p_DB_U_chr[[23]][chunks_chr[7]:(chunks_chr[8]-1)]
p_DB_U_chr_8 <- p_DB_U_chr[[23]][chunks_chr[8]:(chunks_chr[9]-1)]
p_DB_U_chr_9 <- p_DB_U_chr[[23]][chunks_chr[9]:(chunks_chr[10]-1)]
p_DB_U_chr_10 <- p_DB_U_chr[[23]][chunks_chr[10]:(chunks_chr[11]-1)]
p_DB_U_chr_11 <- p_DB_U_chr[[23]][chunks_chr[11]:(chunks_chr[12]-1)]
p_DB_U_chr_12 <- p_DB_U_chr[[23]][chunks_chr[12]:(chunks_chr[13]-1)]
p_DB_U_chr_13 <- p_DB_U_chr[[23]][chunks_chr[13]:(chunks_chr[14]-1)]
p_DB_U_chr_14 <- p_DB_U_chr[[23]][chunks_chr[14]:length(p_DB_U_chr[[23]])]

#p_DB_U_chr_1_13 <- p_DB_U_chr[[1]][300001:325000]
#p_DB_U_chr_1_14 <- p_DB_U_chr[[1]][325001:338735]

test15_chrM <- list(p_DB_U_chr_1, p_DB_U_chr_2, 
                    p_DB_U_chr_3, p_DB_U_chr_4,
                    p_DB_U_chr_5, p_DB_U_chr_6,
                    p_DB_U_chr_7, p_DB_U_chr_8,
                    p_DB_U_chr_9, p_DB_U_chr_10,
                    p_DB_U_chr_11, p_DB_U_chr_12,
                    p_DB_U_chr_13, p_DB_U_chr_14)

rm(p_DB_U_chr_1, p_DB_U_chr_2, 
   p_DB_U_chr_3, p_DB_U_chr_4,
   p_DB_U_chr_5, p_DB_U_chr_6,
   p_DB_U_chr_7, p_DB_U_chr_8,
   p_DB_U_chr_9, p_DB_U_chr_10,
   p_DB_U_chr_11, p_DB_U_chr_12,
   p_DB_U_chr_13, p_DB_U_chr_14)
# chr1-----
mt_param <- MulticoreParam(workers = 14)
test_genes_3_chr1 <- bplapply(test3_chr1, fun, BPPARAM = mt_param)
test3_chr1 <- bind_ranges(test3_chr1)
test_genes_3_chr1 <- test_genes_3_chr1 %>%
  bind_rows() %>% 
  bind_cols(as_tibble(test3_chr1)) %>% 
  dplyr::select(name:subjectHits,rnaID)
chr_21_Y_1 <- bind_rows(test_genes_1,test_genes_2,test_genes_3_chr1)

p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_Y_21_1.txt")
# chr22 ----
mt_param <- MulticoreParam(workers = 12)
test_genes_4_chr22 <- bplapply(test4_chr22, fun, BPPARAM = mt_param)

test4_chr22 <- bind_ranges(test4_chr22)

test_genes_4_chr22 <- test_genes_4_chr22 %>%
  bind_rows() %>% 
  bind_cols(as_tibble(test4_chr22)) %>% 
  dplyr::select(name:subjectHits,rnaID)

chr_21_Y_1_22 <- bind_rows(test_genes_1,test_genes_2,test_genes_3_chr1,test_genes_4_chr22)

p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22.txt")

# chr18 ----
mt_param <- MulticoreParam(workers = 12)

test_genes_5_chr18 <- bplapply(test5_chr18, fun, BPPARAM = mt_param)

test5_chr18 <- bind_ranges(test5_chr18)

test_genes_5_chr18 <- test_genes_5_chr18 %>%
  bind_rows() %>% 
  bind_cols(as_tibble(test5_chr18)) %>% 
  dplyr::select(name:subjectHits,rnaID)

chr_21_Y_1_22_18 <- bind_rows(test_genes_1,test_genes_2,test_genes_3_chr1,test_genes_4_chr22, test_genes_5_chr18)

p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18.txt")
# chr20 ----
mt_param <- MulticoreParam(workers = 12)

test_genes_6_chr20 <- bplapply(test6_chr20, fun, BPPARAM = mt_param)

test6_chr20 <- bind_ranges(test6_chr20)

test_genes_6_chr20 <- test_genes_6_chr20 %>%
  bind_rows() %>% 
  bind_cols(as_tibble(test6_chr20)) %>% 
  dplyr::select(name:subjectHits,rnaID)

chr_21_Y_1_22_18_20 <- bind_rows(test_genes_1,
                              test_genes_2,
                              test_genes_3_chr1,
                              test_genes_4_chr22, 
                              test_genes_5_chr18,
                              test_genes_6_chr20)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20.txt")
# chr13 ----
mt_param <- MulticoreParam(workers = 14)

test_genes_7_chr13 <- bplapply(test7_chr13, fun, BPPARAM = mt_param)

test7_chr13 <- bind_ranges(test7_chr13)

test_genes_7_chr13 <- test_genes_7_chr13 %>%
  bind_rows() %>% 
  bind_cols(as_tibble(test7_chr13)) %>% 
  dplyr::select(name:subjectHits,rnaID)

chr_21_Y_1_22_18_20_13 <- bind_rows(test_genes_1,
                                 test_genes_2,
                                 test_genes_3_chr1,
                                 test_genes_4_chr22, 
                                 test_genes_5_chr18,
                                 test_genes_6_chr20,
                                 test_genes_7_chr13)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13.txt")
# chr14 ----
mt_param <- MulticoreParam(workers = 14)

test_genes_8_chr14 <- bplapply(test8_chr14, fun, BPPARAM = mt_param)

test8_chr14 <- bind_ranges(test8_chr14)

test_genes_8_chr14 <- test_genes_8_chr14 %>%
  bind_rows() %>% 
  bind_cols(as_tibble(test8_chr14)) %>% 
  dplyr::select(name:subjectHits,rnaID)

chr_21_Y_1_22_18_20_13_14 <- bind_rows(test_genes_1,
                                    test_genes_2,
                                    test_genes_3_chr1,
                                    test_genes_4_chr22, 
                                    test_genes_5_chr18,
                                    test_genes_6_chr20,
                                    test_genes_7_chr13,
                                    test_genes_8_chr14)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13_14.txt")
# chr19 ----
mt_param <- MulticoreParam(workers = 14)

test_genes_9_chr19 <- bplapply(test9_chr19, fun, BPPARAM = mt_param)

test9_chr19 <- bind_ranges(test9_chr19)

test_genes_9_chr19 <- test_genes_9_chr19 %>%
  bind_rows() %>% 
  bind_cols(as_tibble(test9_chr19)) %>% 
  dplyr::select(name:subjectHits,rnaID)

chr_21_Y_1_22_18_20_13_14_19 <- bind_rows(test_genes_1,
                                       test_genes_2,
                                       test_genes_3_chr1,
                                       test_genes_4_chr22, 
                                       test_genes_5_chr18,
                                       test_genes_6_chr20,
                                       test_genes_7_chr13,
                                       test_genes_8_chr14,
                                       test_genes_9_chr19)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14_19, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13_14_19.txt")
# chr16 ----
mt_param <- MulticoreParam(workers = 14)

test_genes_10_chr16 <- bplapply(test10_chr16, fun, BPPARAM = mt_param)

test10_chr16 <- bind_ranges(test10_chr16)

test_genes_10_chr16 <- test_genes_10_chr16 %>%
  bind_rows() %>% 
  bind_cols(as_tibble(test10_chr16)) %>% 
  dplyr::select(name:subjectHits,rnaID)

chr_21_Y_1_22_18_20_13_14_19_16 <- bind_rows(test_genes_1,
                                          test_genes_2,
                                          test_genes_3_chr1,
                                          test_genes_4_chr22, 
                                          test_genes_5_chr18,
                                          test_genes_6_chr20,
                                          test_genes_7_chr13,
                                          test_genes_8_chr14,
                                          test_genes_9_chr19,
                                          test_genes_10_chr16)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14_19_16, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13_14_19_16.txt")
# chr15 ----
mt_param <- MulticoreParam(workers = 14)

test_genes_11_chr15 <- bplapply(test11_chr15, fun, BPPARAM = mt_param)

test11_chr15 <- bind_ranges(test11_chr15)

test_genes_11_chr15 <- test_genes_11_chr15 %>%
  bind_rows() %>% 
  bind_cols(as_tibble(test11_chr15)) %>% 
  dplyr::select(name:subjectHits,rnaID)

chr_21_Y_1_22_18_20_13_14_19_16_15 <- bind_rows(test_genes_1,
                                             test_genes_2,
                                             test_genes_3_chr1,
                                             test_genes_4_chr22, 
                                             test_genes_5_chr18,
                                             test_genes_6_chr20,
                                             test_genes_7_chr13,
                                             test_genes_8_chr14,
                                             test_genes_9_chr19,
                                             test_genes_10_chr16,
                                             test_genes_11_chr15)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14_19_16_15, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13_14_19_16_15.txt")
# annotate for each ----
annotate_fun <- function(x_1) {
  test_genes_chr <- bplapply(x_1, fun, BPPARAM = mt_param)
  message("bind the ranges")
  x_1 <- bind_ranges(x_1)
  message("make the final tibble")
  test_genes_chr <- test_genes_chr %>%
    bind_rows() %>% 
    bind_cols(as_tibble(x_1)) %>% 
    dplyr::select(name:subjectHits,rnaID)
  test_genes_chr
}
# chr17 ----
mt_param <- MulticoreParam(workers = 14)

test_genes_12_chr17 <- bplapply(test12_chr17, fun, BPPARAM = mt_param)

test12_chr17 <- bind_ranges(test12_chr17)

test_genes_12_chr17 <- test_genes_12_chr17 %>%
  bind_rows() %>% 
  bind_cols(as_tibble(test12_chr17)) %>% 
  dplyr::select(name:subjectHits,rnaID)

chr_21_Y_1_22_18_20_13_14_19_16_15_17 <- bind_rows(test_genes_1,
                                             test_genes_2,
                                             test_genes_3_chr1,
                                             test_genes_4_chr22, 
                                             test_genes_5_chr18,
                                             test_genes_6_chr20,
                                             test_genes_7_chr13,
                                             test_genes_8_chr14,
                                             test_genes_9_chr19,
                                             test_genes_10_chr16,
                                             test_genes_11_chr15,
                                             test_genes_12_chr1)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14_19_16_15_17, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13_14_19_16_15_17.txt")
# chr8 ----
mt_param <- MulticoreParam(workers = 14)

test_genes_13_chr8 <- bplapply(test13_chr8, fun, BPPARAM = mt_param)

test13_chr8 <- bind_ranges(test13_chr8)

test_genes_13_chr8 <- test_genes_13_chr8 %>%
  bind_rows() %>% 
  bind_cols(as_tibble(test13_chr8)) %>% 
  dplyr::select(name:subjectHits,rnaID)

chr_21_Y_1_22_18_20_13_14_19_16_15_17_8 <- bind_rows(test_genes_1,
                                                   test_genes_2,
                                                   test_genes_3_chr1,
                                                   test_genes_4_chr22, 
                                                   test_genes_5_chr18,
                                                   test_genes_6_chr20,
                                                   test_genes_7_chr13,
                                                   test_genes_8_chr14,
                                                   test_genes_9_chr19,
                                                   test_genes_10_chr16,
                                                   test_genes_11_chr15,
                                                   test_genes_12_chr17,
                                                   test_genes_13_chr8)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14_19_16_15_17_8, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13_14_19_16_15_17_8.txt")

# chrX ----
mt_param <- MulticoreParam(workers = 14)

test_genes_14_chrX <- bplapply(test14_chrX, fun, BPPARAM = mt_param)

test14_chrX <- bind_ranges(test14_chrX)

test_genes_14_chrX <- test_genes_14_chrX %>%
  bind_rows() %>% 
  bind_cols(as_tibble(test14_chrX)) %>% 
  dplyr::select(name:subjectHits,rnaID)

chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X <- bind_rows(test_genes_1,
                                                     test_genes_2,
                                                     test_genes_3_chr1,
                                                     test_genes_4_chr22, 
                                                     test_genes_5_chr18,
                                                     test_genes_6_chr20,
                                                     test_genes_7_chr13,
                                                     test_genes_8_chr14,
                                                     test_genes_9_chr19,
                                                     test_genes_10_chr16,
                                                     test_genes_11_chr15,
                                                     test_genes_12_chr17,
                                                     test_genes_13_chr8,
                                                     test_genes_14_chrX)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X.txt")


# chrM ----
test_genes_15_chrM <- annotate_fun(test15_chrM)

chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X <- bind_rows(test_genes_1,
                                                       test_genes_2,
                                                       test_genes_3_chr1,
                                                       test_genes_4_chr22, 
                                                       test_genes_5_chr18,
                                                       test_genes_6_chr20,
                                                       test_genes_7_chr13,
                                                       test_genes_8_chr14,
                                                       test_genes_9_chr19,
                                                       test_genes_10_chr16,
                                                       test_genes_11_chr15,
                                                       test_genes_12_chr17,
                                                       test_genes_13_chr8,
                                                       test_genes_14_chrX)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X.txt")
