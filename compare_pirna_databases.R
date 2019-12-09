## load libraries -----
library(plyranges)
library(tidyverse)
## load pirnadb -----
piRNA_dbs_files <- list.files("/home/0/piRNA_DBs", full.names = TRUE)
# piRBase
pirbase <- piRNA_dbs_files[2] %>% read_bed()
# piRNADB
pirnadb <- piRNA_dbs_files[5] %>% 
  read_tsv(comment = "#", col_names = c("seqnames", "X2", "x3", "start", "end", "x6", "strand", "x7", "piRNA" )) %>% 
  select(-X2, -x3, -x6, -x7) %>% 
  mutate(seqnames = if_else(seqnames == "chrMT", "chrM", as.character(seqnames))) %>% 
  as_granges() %>% 
  arrange(start)
# piRNADB cluster
pirnadb_cl <- piRNA_dbs_files[4] %>% 
  read_tsv() %>% select(-Build_Code) %>% mutate(Strand = if_else(Strand == "biDirectional","*", Strand),
                                                Chromosome = str_replace(.$Chromosome,"^", "chr")) %>% 
  as_granges(seqnames = Chromosome,
             start = Start,
             end = End, 
             strand = Strand
             ) 
# dashr DB cluster
dashr_db <- piRNA_dbs_files[1] %>% read_tsv(col_names = F) %>%  as_granges(seqnames = X1,
                                                                           start = X2,
                                                                           end = X3, 
                                                                           strand = X6) 
# cluster db 
pirna_cl_db <- piRNA_dbs_files[3] %>% 
  read_gff() %>% 
  as_tibble %>% 
  mutate(score = 1:length(.$score),
         seqnames = str_replace(.$seqnames,"^", "chr")) %>% 
  unite(type, type:score) %>% 
  select(-group, -phase, -source) %>% 
  as_granges()
# check all databases to find overlaps
pirna_cl_db %>% join_overlap_intersect_directed(pirnadb_cl)
