# Collapse_unannot_peaks.r
## Konstantinos Geles
### Collapsing and creating the matrix for downstream DE  
## load libraries -----
library(plyranges, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(data.table, quietly = TRUE)
## add todate ----
todate <- format(Sys.time(), "%d_%b_%Y")
## create the dir for the analysis -----
#my_exp <- "whatever"
#dir.create(str_glue("./{my_exp}_analysis"))
#dat_path <- str_glue("{my_exp}_analysis")
# import data ------
path <- "/home/0/Project_piRNA/3_IPP_COLO205_SPAR_results"
smallRNA_files <- dir(path, full.names = TRUE,
                      pattern = "peaks_unannot.with_con.+.xls",
                      recursive = TRUE)

main_Grange <- rbindlist(sapply(smallRNA_files,fread,
                       simplify=FALSE,
                       verbose=getOption("datatable.verbose", TRUE)),
                use.names= TRUE,
                idcol="sample_file")  %>%
  rename(peakChr = "#peakChr") %>% 
  as_tibble %>% 
  mutate(sample_file = sample_file %>% 
           str_remove(".+SPAR_results/") %>%
           str_remove(".trimmed_.+") %>%
           str_remove(pattern = "COLO205_IPP_") %>%
           str_remove(pattern = "noAb_")
           ) %>% 
  select(sample_file:peakStrand) %>% 
  select(-peakID) %>% 
  as_granges(seqnames = peakChr,
                         start = peakChrStart,
                         end = peakChrEnd,
                         strand = peakStrand) %>% 
  arrange(start)
##import seq info ----
sInfo <- Seqinfo(genome="hg38")
seqlevels(sInfo) <- seqlevels(main_Grange)
seqinfo(main_Grange)<- sInfo
## reduce the granges of all samples in one global grange -----
common_reduced_ranges <- main_Grange %>% 
  reduce_ranges_directed()
## make a list with samples names ---- 
smpls_lst <- main_Grange %>% 
  as_tibble() %>% 
  distinct(sample_file) %>% 
  arrange(sample_file) %>% 
  deframe() %>% 
  as.list() %>% 
  set_names(.)
## prepare the list with granges samples -----
grSamples_l <- map(smpls_lst, ~ main_Grange  %>% 
      filter(sample_file == .x) %>% 
      select(!!str_c(.x, "_peakExpr") := "peakExpressionValue",everything()))
## join each sample with the reduced Granges of all samples -----
common_reduced_ranges <- common_reduced_ranges %>%
  list(.) %>% 
  c(., grSamples_l) %>% 
  purrr::reduce(join_overlap_left_directed) %>% 
  mutate(unannot_peak = selfmatch(.)) %>%
  dplyr::mutate_at(dplyr::vars(dplyr::ends_with('peakExpr')), 
                   function(.) ifelse(is.na(.), 0, .)) %>% 
  select(-starts_with("sample_file")) %>% 
  group_by(unannot_peak)
# remake the reduced Granges with info about unannot_peaks -----
my_reduced_Granges <- common_reduced_ranges  %>% 
  reduce_ranges_directed
# calculate the median values for all peaks ----
my_medianL <- map(smpls_lst, ~common_reduced_ranges %>% 
      select(!!str_c(.x, "_peakExpr"), unannot_peak) %>% 
      as_tibble() %>% 
      group_by(unannot_peak) %>% 
      summarize(!!str_c(.x, "_peakExpr") := median(.data[[!!str_c(.x, "_peakExpr")]], na.rm = TRUE)))

my_medianL <- my_medianL %>%
  purrr::reduce(inner_join)
## make the final Granges Obect -----
my_unAnnot_GR <- my_reduced_Granges %>% 
  as_tibble() %>% 
  inner_join(my_medianL) %>%
  mutate(unannot_peak = unannot_peak %>% 
           str_c("peak_",.)) %>% 
  as_granges() %>% 
  write_rds("unannot_peak_Ranges.rds")
## create the collapsed matrix of annotated -----
path 
smallRNA_files <- smallRNA_files <- dir(path, full.names = TRUE,
                                        pattern = "smRNA_gene_expression.xls",
                                        recursive = TRUE)
### load the list of files in one table -----
DT <- rbindlist(sapply(smallRNA_files, fread,
                       simplify=FALSE,
                       verbose=getOption("datatable.verbose", TRUE)),
                use.names= TRUE, idcol="sample_file")  %>%
  rename(smallRNA = "#Gene") %>% 
  select(-RPM) %>% 
  as_tibble() %>% 
  mutate(sample_file = sample_file %>% 
           str_remove(".trimmed_.+") %>% 
           str_remove(".+/")
           )
### separate first column to multiple ----
#separate("#Gene",c("chr","start","end","strand","smallRNA","DQ"), sep = ":") %>% 
#as_tibble()
#DT <- DT %>% 
#  str_remove(".+SPAR_results/") %>%
#  str_remove(".trimmed_.+") %>%
#  str_remove(pattern = "COLO205_IPP_") %>%
#  str_remove(pattern = "noAb_")

### make matrix for DE analysis -----
dt <- DT %>% 
  group_by(sample_file) %>% 
  spread(key = "sample_file", value = "ReadCount" ) %>% 
  rename(Peaks_GeneClass = "GeneClass")
# collapse not annot and annot in one matrix ----
my_unAnnot_GR <- my_unAnnot_GR %>% 
  as_tibble() %>% 
  unite(smallRNA,seqnames:strand, sep = "_") %>% 
  rename(Peaks_GeneClass = "unannot_peak")
print("are the name identical?")
identical(names(my_unAnnot_GR),names(dt))
all_peaks_smRNAs <- dt[names(my_unAnnot_GR)] %>% 
  bind_rows(my_unAnnot_GR) %>% 
  write_tsv("all_peaks_smRNAs.txt")
