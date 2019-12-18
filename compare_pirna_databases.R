## load libraries -----
library(plyranges)
library(tidyverse)
## load pirnadbs -----
piRNA_dbs_files <- list.files("/home/0/piRNA_DBs", full.names = TRUE)
# piRBase
pirbase <- piRNA_dbs_files[2] %>% 
  read_bed() %>% 
  as_tibble() %>% 
  rename(pirbase = "name") %>% 
  select(-score) %>% 
  mutate(pirbase_coor = str_c(.$pirbase, .$seqnames, .$start, .$end, .$strand, sep = "_")) %>% 
  as_granges() %>% 
  keepStandardChromosomes(pruning.mode="coarse")
# piRNADB
pirnadb <- piRNA_dbs_files[5] %>% 
  read_tsv(comment = "#", col_names = c("seqnames", "x2", "x3", "start", "end", "x6", "strand", "x7", "piRNAdb" )) %>% 
  select(-x2, -x3, -x6, -x7) %>% 
  mutate(seqnames = if_else(seqnames == "chrMT", "chrM", as.character(seqnames))) %>% 
  mutate(pirnadb_coor = str_c(.$piRNAdb, .$seqnames, .$start, .$end, .$strand, sep = "_")) %>% 
  as_granges() %>% 
  arrange(start) %>% 
  keepStandardChromosomes(pruning.mode="coarse")
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
  arrange(start) %>% 
  keepStandardChromosomes(pruning.mode="coarse")
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
length(p_DB_U_chr[[24]]) %/% 15

(l_chr <- plyr::round_any(length(p_DB_U_chr[[24]]), 100, ceiling)) 
(l_chr %/% 14)

(chunks_chr <- seq(1,l_chr, by = l_chr %/% 14))

p_DB_U_chr_1 <- p_DB_U_chr[[24]][chunks_chr[1]:(chunks_chr[2]-1)]
p_DB_U_chr_2 <- p_DB_U_chr[[24]][chunks_chr[2]:(chunks_chr[3]-1)]
p_DB_U_chr_3 <- p_DB_U_chr[[24]][chunks_chr[3]:(chunks_chr[4]-1)]
p_DB_U_chr_4 <- p_DB_U_chr[[24]][chunks_chr[4]:(chunks_chr[5]-1)]
p_DB_U_chr_5 <- p_DB_U_chr[[24]][chunks_chr[5]:(chunks_chr[6]-1)]
p_DB_U_chr_6 <- p_DB_U_chr[[24]][chunks_chr[6]:(chunks_chr[7]-1)]
p_DB_U_chr_7 <- p_DB_U_chr[[24]][chunks_chr[7]:(chunks_chr[8]-1)]
p_DB_U_chr_8 <- p_DB_U_chr[[24]][chunks_chr[8]:(chunks_chr[9]-1)]
p_DB_U_chr_9 <- p_DB_U_chr[[24]][chunks_chr[9]:(chunks_chr[10]-1)]
p_DB_U_chr_10 <- p_DB_U_chr[[24]][chunks_chr[10]:(chunks_chr[11]-1)]
p_DB_U_chr_11 <- p_DB_U_chr[[24]][chunks_chr[11]:(chunks_chr[12]-1)]
p_DB_U_chr_12 <- p_DB_U_chr[[24]][chunks_chr[12]:(chunks_chr[13]-1)]
p_DB_U_chr_13 <- p_DB_U_chr[[24]][chunks_chr[13]:(chunks_chr[14]-1)]
p_DB_U_chr_14 <- p_DB_U_chr[[24]][chunks_chr[14]:length(p_DB_U_chr[[24]])]

#p_DB_U_chr_1_13 <- p_DB_U_chr[[1]][300001:325000]
#p_DB_U_chr_1_14 <- p_DB_U_chr[[1]][325001:338735]

test14_chrX <- list(p_DB_U_chr_1, p_DB_U_chr_2, 
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

#
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

test_genes_12_chr17 <- annotate_fun(test12_chr17)

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
                                             test_genes_12_chr17)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14_19_16_15_17, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13_14_19_16_15_17.txt")
# chr8 ----
mt_param <- MulticoreParam(workers = 14)

test_genes_13_chr8 <- annotate_fun(test13_chr8)

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

test_genes_14_chrX <- annotate_fun(test14_chrX)

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

# function to slice the dataset ----
slice_dataframe_fun <- function(n_1,by_n) {
  (l_chr <- plyr::round_any(length(p_DB_U_chr[[n_1]]), by_n, ceiling)) 
  (l_chr %/% 14)
  (chunks_chr <- seq(1,l_chr, by = l_chr %/% 14))
  
  p_DB_U_chr_1 <- p_DB_U_chr[[n_1]][chunks_chr[1]:(chunks_chr[2]-1)]
  p_DB_U_chr_2 <- p_DB_U_chr[[n_1]][chunks_chr[2]:(chunks_chr[3]-1)]
  p_DB_U_chr_3 <- p_DB_U_chr[[n_1]][chunks_chr[3]:(chunks_chr[4]-1)]
  p_DB_U_chr_4 <- p_DB_U_chr[[n_1]][chunks_chr[4]:(chunks_chr[5]-1)]
  p_DB_U_chr_5 <- p_DB_U_chr[[n_1]][chunks_chr[5]:(chunks_chr[6]-1)]
  p_DB_U_chr_6 <- p_DB_U_chr[[n_1]][chunks_chr[6]:(chunks_chr[7]-1)]
  p_DB_U_chr_7 <- p_DB_U_chr[[n_1]][chunks_chr[7]:(chunks_chr[8]-1)]
  p_DB_U_chr_8 <- p_DB_U_chr[[n_1]][chunks_chr[8]:(chunks_chr[9]-1)]
  p_DB_U_chr_9 <- p_DB_U_chr[[n_1]][chunks_chr[9]:(chunks_chr[10]-1)]
  p_DB_U_chr_10 <- p_DB_U_chr[[n_1]][chunks_chr[10]:(chunks_chr[11]-1)]
  p_DB_U_chr_11 <- p_DB_U_chr[[n_1]][chunks_chr[11]:(chunks_chr[12]-1)]
  p_DB_U_chr_12 <- p_DB_U_chr[[n_1]][chunks_chr[12]:(chunks_chr[13]-1)]
  p_DB_U_chr_13 <- p_DB_U_chr[[n_1]][chunks_chr[13]:(chunks_chr[14]-1)]
  p_DB_U_chr_14 <- p_DB_U_chr[[n_1]][chunks_chr[14]:length(p_DB_U_chr[[n_1]])]
  
  test_chr <- list(p_DB_U_chr_1, p_DB_U_chr_2, 
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
  
  test_chr
}


# chrM ----
test15_chrM <- slice_dataframe_fun(23,10)
test_genes_15_chrM <- annotate_fun(test15_chrM)

chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M <- bind_rows(test_genes_1,
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
                                                       test_genes_14_chrX,
                                                       test_genes_15_chrM)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M.txt")

# chr9 ----
test16_chr9 <- slice_dataframe_fun(22,100)
test_genes_16_chr9 <- annotate_fun(test16_chr9)

chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9 <- bind_rows(test_genes_1,
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
                                                       test_genes_14_chrX,
                                                       test_genes_15_chrM,
                                                       test_genes_16_chr9)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9.txt")

# chr11 ----
test17_chr11 <- slice_dataframe_fun(3,100)
identical(length(p_DB_U_chr[[3]]), length(bind_ranges(test17_chr11)))
test_genes_17_chr11 <- annotate_fun(test17_chr11)

chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11 <- bind_rows(test_genes_1,
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
                                                         test_genes_14_chrX,
                                                         test_genes_15_chrM,
                                                         test_genes_16_chr9,
                                                         test_genes_17_chr11)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11.txt")
# chr2 ----
test18_chr2 <- slice_dataframe_fun(12,100)
identical(length(p_DB_U_chr[[12]]), length(bind_ranges(test18_chr2)))
test_genes_18_chr2 <- annotate_fun(test18_chr2)

chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2 <- bind_rows(test_genes_1,
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
                                                              test_genes_14_chrX,
                                                              test_genes_15_chrM,
                                                              test_genes_16_chr9,
                                                              test_genes_17_chr11,
                                                              test_genes_18_chr2)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2.txt")
# chr3 ----
test19_chr3 <- slice_dataframe_fun(16,100)
identical(length(p_DB_U_chr[[16]]), length(bind_ranges(test19_chr3)))
test_genes_19_chr3 <- annotate_fun(test19_chr3)

chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3 <- bind_rows(test_genes_1,
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
                                                                test_genes_14_chrX,
                                                                test_genes_15_chrM,
                                                                test_genes_16_chr9,
                                                                test_genes_17_chr11,
                                                                test_genes_18_chr2,
                                                                test_genes_19_chr3)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3.txt")
# chr7 ----
test20_chr7 <- slice_dataframe_fun(20,100)
identical(length(p_DB_U_chr[[20]]), length(bind_ranges(test20_chr7)))
test_genes_20_chr7 <- annotate_fun(test20_chr7)

chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3_7 <- bind_rows(test_genes_1,
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
                                                                   test_genes_14_chrX,
                                                                   test_genes_15_chrM,
                                                                   test_genes_16_chr9,
                                                                   test_genes_17_chr11,
                                                                   test_genes_18_chr2,
                                                                   test_genes_19_chr3,
                                                                   test_genes_20_chr7)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3_7, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3_7.txt")
# chr4 ----
test21_chr4 <- slice_dataframe_fun(17,100)
identical(length(p_DB_U_chr[[17]]), length(bind_ranges(test21_chr4)))
test_genes_21_chr4 <- annotate_fun(test21_chr4)

chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3_7_4 <- bind_rows(test_genes_1,
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
                                                                    test_genes_14_chrX,
                                                                    test_genes_15_chrM,
                                                                    test_genes_16_chr9,
                                                                    test_genes_17_chr11,
                                                                    test_genes_18_chr2,
                                                                    test_genes_19_chr3,
                                                                    test_genes_20_chr7,
                                                                    test_genes_21_chr4)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3_7_4, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3_7_4.txt")
# chr5 ----
test22_chr5 <- slice_dataframe_fun(18,100)
identical(length(p_DB_U_chr[[18]]), length(bind_ranges(test22_chr5)))
test_genes_22_chr5 <- annotate_fun(test22_chr5)

chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3_7_4_5 <- bind_rows(test_genes_1,
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
                                                                      test_genes_14_chrX,
                                                                      test_genes_15_chrM,
                                                                      test_genes_16_chr9,
                                                                      test_genes_17_chr11,
                                                                      test_genes_18_chr2,
                                                                      test_genes_19_chr3,
                                                                      test_genes_20_chr7,
                                                                      test_genes_21_chr4,
                                                                      test_genes_22_chr5)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3_7_4_5, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3_7_4_5.txt")
# chr6 ----
test23_chr6 <- slice_dataframe_fun(19,100)
identical(length(p_DB_U_chr[[19]]), length(bind_ranges(test23_chr6)))
test_genes_23_chr6 <- annotate_fun(test23_chr6)

chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3_7_4_5_6 <- bind_rows(test_genes_1,
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
                                                                        test_genes_14_chrX,
                                                                        test_genes_15_chrM,
                                                                        test_genes_16_chr9,
                                                                        test_genes_17_chr11,
                                                                        test_genes_18_chr2,
                                                                        test_genes_19_chr3,
                                                                        test_genes_20_chr7,
                                                                        test_genes_21_chr4,
                                                                        test_genes_22_chr5,
                                                                        test_genes_23_chr6)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3_7_4_5_6, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3_7_4_5_6.txt")
# chr12 ----
test24_chr12 <- slice_dataframe_fun(4,100)
identical(length(p_DB_U_chr[[4]]), length(bind_ranges(test24_chr12)))
test_genes_24_chr12 <- annotate_fun(test24_chr12)

chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3_7_4_5_6_12 <- bind_rows(test_genes_1,
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
                                                                          test_genes_14_chrX,
                                                                          test_genes_15_chrM,
                                                                          test_genes_16_chr9,
                                                                          test_genes_17_chr11,
                                                                          test_genes_18_chr2,
                                                                          test_genes_19_chr3,
                                                                          test_genes_20_chr7,
                                                                          test_genes_21_chr4,
                                                                          test_genes_22_chr5,
                                                                          test_genes_23_chr6,
                                                                          test_genes_24_chr12)
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3_7_4_5_6_12, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3_7_4_5_6_12.txt")
# chr10 ----
test25_chr10 <- slice_dataframe_fun(2,100)
identical(length(p_DB_U_chr[[2]]), length(bind_ranges(test25_chr10)))
test_genes_25_chr10 <- annotate_fun(test25_chr10)

chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3_7_4_5_6_12_10 <- bind_rows(
  test_genes_3_chr1,
  test_genes_18_chr2,
  test_genes_19_chr3,
  test_genes_21_chr4,
  test_genes_22_chr5,
  test_genes_23_chr6,
  test_genes_20_chr7,
  test_genes_13_chr8,
  test_genes_16_chr9,
  test_genes_25_chr10,
  test_genes_17_chr11,
  test_genes_24_chr12,
  test_genes_7_chr13,
  test_genes_8_chr14,
  test_genes_11_chr15,
  test_genes_10_chr16,
  test_genes_12_chr17,
  test_genes_5_chr18,
  test_genes_9_chr19,
  test_genes_6_chr20,
  test_genes_1,
  test_genes_4_chr22, 
  test_genes_14_chrX,
  test_genes_2,
  test_genes_15_chrM)
 
p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y_1_22_18_20_13_14_19_16_15_17_8_X_M_9_11_2_3_7_4_5_6_12_10, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chr_all.txt")


map() %>% bind_rows(.id = .x)
# summarize stats for the dbs -----
dbs <- read_tsv("/home/0/IPP/spar_ipp/all_DB_genes_chr_all.txt", col_names = TRUE, 
                col_types = cols(
                  .default = col_character(),
                  start = col_double(),
                  end = col_double(),
                  width = col_double(),
                  #Cluster_Acess = col_logical(),
                  #cl_db = col_logical(),
                  distance = col_double(),
                  insideDistance = col_double(),
                  exonnumber = col_double(),
                  nexons = col_double(),
                  geneL = col_double(),
                  #codingL = col_logical(),
                  Geneid = col_double(),
                  subjectHits = col_double()
                ))


# all entries
dbs %>% nrow()
# all regions 
dbs %>% group_by(rnaID) %>% n_groups()
# length of the regions
dbs$width %>% summary()
# common between the three databases
dbs %>% filter(!is.na(dashr_srna), !is.na(piRNAdb),!is.na(pirbase)) %>% group_by(rnaID) %>% n_groups()
# common between the three databases and 
dbs %>% filter(!is.na(dashr_srna), !is.na(piRNAdb),!is.na(pirbase), !is.na(Cluster_Acess)) %>% group_by(rnaID) %>% n_groups()
# common between the three databases and clusters of piRNADB
dbs %>% filter(!is.na(dashr_srna), !is.na(piRNAdb),!is.na(pirbase), !is.na(cl_db)) %>% group_by(rnaID) %>% n_groups()
# common between the three databases clusterDB and clusters of piRNADB
dbs %>% filter(!is.na(dashr_srna), !is.na(piRNAdb),!is.na(pirbase), !is.na(cl_db), !is.na(Cluster_Acess)) %>% group_by(rnaID) %>% n_groups()
# common between the three databases not in clusters but inside exon
dbs %>% filter(!is.na(dashr_srna), !is.na(piRNAdb),!is.na(pirbase), is.na(Cluster_Acess), is.na(cl_db), subregion == "inside exon") %>% group_by(rnaID) %>% n_groups()
# common between the three databases in clusters and inside exon
dbs %>% filter(!is.na(dashr_srna), !is.na(piRNAdb),!is.na(pirbase), !is.na(Cluster_Acess) | !is.na(cl_db), subregion == "inside exon") %>% group_by(rnaID) %>% n_groups()
# common between the three databases in clusters and inside introns
dbs %>% filter(!is.na(dashr_srna), !is.na(piRNAdb),!is.na(pirbase), !is.na(Cluster_Acess) | !is.na(cl_db), subregion == "inside intron") %>% group_by(rnaID)
# regions found in piRBase
dbs %>% filter(!is.na(pirbase)) %>% group_by(rnaID) %>% n_groups()
# find how many pirnas are inside the regions 
dbs %>% group_by(rnaID) %>% summarise(n = dplyr::n()) %>% 
summarise_at(vars(n) ,list(min = min, Q1=~quantile(., probs = 0.25),
                               median=median, Q3=~quantile(., probs = 0.75),
                               max=max))
summr_rnaID <- dbs %>% group_by(rnaID) %>% 
  summarise(n = dplyr::n()) %>% arrange(desc(n))
dbs %>% filter(rnaID == "RNA3141770") 
dbs %>% filter(rnaID == "RNA3141770") %>% select(dashr_srna) %>% deframe %>% as_factor() %>% levels() %>% length()
dbs %>% filter(rnaID == "RNA3141770") %>% select(piRNAdb) %>% deframe %>% as_factor() %>% levels() %>% length()
dbs %>% filter(rnaID == "RNA3141770") %>% select(pirbase) %>% deframe %>% as_factor() %>% levels() %>% length()
# how many are only in piRBase
dbs %>% filter(!is.na(pirbase),is.na(piRNAdb),is.na(dashr_srna) ) %>% group_by(rnaID) %>% n_groups()
summarized_rnaID <- dbs %>% filter(!is.na(pirbase),is.na(piRNAdb),is.na(dashr_srna) ) %>% 
  group_by(rnaID) %>% summarise(n = dplyr::n()) %>% arrange(desc(n))


# same for long ranges ----
## create a union of all dbs 
dashr_db_red <- dashr_db %>% 
  filter(dashr_type == "piRNA") %>% 
  select(-dashr_type) %>% 
  reduce_ranges_directed()

pirbase_red <- pirbase %>% reduce_ranges_directed()

pirnadb_red <- pirnadb %>% reduce_ranges_directed()
# concat them and reduce
pirna_DB_union <- c(dashr_db_red, pirbase_red, pirnadb_red) %>% 
  reduce_ranges_directed()
## pirbase
pirbase_red %>% 
  as_tibble() %>% 
  group_by(seqnames) %>% 
  summarise_at(vars(width) ,list(min = min, Q1=~quantile(., probs = 0.25),
                                 median=median, Q3=~quantile(., probs = 0.75),
                                 max=max)) %>% 
  arrange(as.character(seqnames)) %>% tail(20)
# evaluating summary statistics we find that most of the regions
# that are around ~34 base pairs for that reason we will remove 
# all regions bigger than 39

pirbase_long <- pirbase_red %>% 
  as_tibble() %>% 
  filter(width >= 40)
#dashr
dashr_db_red %>% 
  as_tibble() %>% 
  group_by(seqnames) %>% 
  summarise_at(vars(width) ,list(min = min, Q1=~quantile(., probs = 0.25),
                                 median=median, Q3=~quantile(., probs = 0.75),
                                 max=max)) %>% 
  arrange(as.character(seqnames)) %>% tail(20)

dashr_db_Long <- dashr_db_red %>% 
  as_tibble() %>% 
  filter(width >= 40)
#pirnadb
pirnadb_red %>% 
  as_tibble() %>% 
  group_by(seqnames) %>% 
  summarise_at(vars(width) ,list(min = min, Q1=~quantile(., probs = 0.25),
                                 median=median, Q3=~quantile(., probs = 0.75),
                                 max=max)) %>% 
  arrange(as.character(seqnames)) %>% tail(20)

pirnadb_long <- pirnadb_red %>% 
  as_tibble() %>% 
  filter(width >= 40)
# concat them and reduce
pirna_DB_union <- c(as_granges(pirbase_long), as_granges(dashr_db_Long), as_granges(pirnadb_red)) %>% 
  reduce_ranges_directed() %>%  
  arrange(start) %>%
  as_tibble() %>%  
  mutate(rnaID = 1:length(.$strand) %>%
  str_c("RNA_cluster", .)) %>% 
  as_granges()

pirna_DB_union %>% 
  as_tibble() %>% 
  group_by(seqnames) %>% 
  summarise_at(vars(width) ,list(min = min, Q1=~quantile(., probs = 0.25),
                                 median=median, Q3=~quantile(., probs = 0.75),
                                 max=max)) %>% 
  arrange(as.character(seqnames)) %>% tail(20)
# make the file for the 3 dbds
#p_DB_U <- pirna_DB_union %>% 
  join_overlap_left_directed(dashr_db %>% 
                               filter(dashr_type == "piRNA") %>% 
                               select(-dashr_type)) %>% 
  join_overlap_left_directed(pirnadb) %>% 
  join_overlap_left_directed(pirbase) %>% 
  join_overlap_left_directed(pirnadb_cl) %>%
  join_overlap_left_directed(pirna_cl_db)
