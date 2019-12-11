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
seq(1,340000,by= 25000)

p_DB_U_chr_1_1 <- p_DB_U_chr[[1]][1:25000]
p_DB_U_chr_1_2 <- p_DB_U_chr[[1]][25001:50000]
p_DB_U_chr_1_3 <- p_DB_U_chr[[1]][50001:75000]
p_DB_U_chr_1_4 <- p_DB_U_chr[[1]][75001:100000]
p_DB_U_chr_1_5 <- p_DB_U_chr[[1]][100001:125000]
p_DB_U_chr_1_6 <- p_DB_U_chr[[1]][125001:150000]
p_DB_U_chr_1_7 <- p_DB_U_chr[[1]][150001:175000]
p_DB_U_chr_1_8 <- p_DB_U_chr[[1]][175001:200000]
p_DB_U_chr_1_9 <- p_DB_U_chr[[1]][200001:225000]
p_DB_U_chr_1_10 <- p_DB_U_chr[[1]][225001:250000]
p_DB_U_chr_1_11 <- p_DB_U_chr[[1]][250001:275000]
p_DB_U_chr_1_12 <- p_DB_U_chr[[1]][275001:300000]
p_DB_U_chr_1_13 <- p_DB_U_chr[[1]][300001:325000]
p_DB_U_chr_1_14 <- p_DB_U_chr[[1]][325001:338735]

test3_chr1 <- list(p_DB_U_chr_1_1, p_DB_U_chr_1_2, 
                   p_DB_U_chr_1_3, p_DB_U_chr_1_4,
                   p_DB_U_chr_1_5, p_DB_U_chr_1_6,
                   p_DB_U_chr_1_7, p_DB_U_chr_1_8,
                   p_DB_U_chr_1_9, p_DB_U_chr_1_10,
                   p_DB_U_chr_1_11, p_DB_U_chr_1_12,
                   p_DB_U_chr_1_13, p_DB_U_chr_1_14)

test1 <- list(p_DB_U_chr_Y_1, p_DB_U_chr_Y_2, 
              p_DB_U_chr_Y_3, p_DB_U_chr_Y_4,
              p_DB_U_chr_Y_5, p_DB_U_chr_Y_6, 
              p_DB_U_chr_Y_7, p_DB_U_chr_Y_8, 
              p_DB_U_chr_Y_9, p_DB_U_chr_Y_10,
              p_DB_U_chr_Y_11, p_DB_U_chr_Y_12)

test2 <- list(p_DB_U_chr_21_1, p_DB_U_chr_21_2, 
              p_DB_U_chr_21_3, p_DB_U_chr_21_4,
              p_DB_U_chr_21_5, p_DB_U_chr_21_6)
test_genes_2 <- bplapply(test2, fun, BPPARAM = mt_param)


test2 <- bind_ranges(test2)
test_genes_2 <- test_genes_2 %>%
  bind_rows() %>% 
  bind_cols(as_tibble(test2)) %>% 
  dplyr::select(name:subjectHits,rnaID)

chr_21_Y <- bind_rows(test_genes_1,test_genes_2)

p_DB_U %>% 
  as_tibble %>% 
  inner_join(chr_21_Y, by = "rnaID") %>% 
  write_tsv("all_DB_genes_chrY_21.txt")

mt_param <- MulticoreParam(workers = 14)
test_genes_3_chr1 <- bplapply(test3_chr1, fun, BPPARAM = mt_param) 
