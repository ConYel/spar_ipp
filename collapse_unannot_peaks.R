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

#my_exp <- "COLO205_IPP"
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
test_GR <- common_reduced_ranges  
##################################
test1 <- common_reduced_ranges %>% 
  join_overlap_left_directed( grSamples_l[[1]]) %>% 
  join_overlap_left_directed( grSamples_l[[2]]) %>% 
  join_overlap_left_directed( grSamples_l[[3]]) %>% 
  mutate(revmap = selfmatch(.)) %>%
  dplyr::mutate_at(dplyr::vars(dplyr::ends_with('peakExpr')), 
                   function(.) ifelse(is.na(.), 0, .)) %>% 
  select(-starts_with("samples_file"))



test2 <- test1 %>% 
  group_by(revmap) %>% 
  reduce_ranges_directed(
    C_Ab_1_peakExpr = median(C_Ab_1_peakExpr, na.rm = TRUE),
    C_Ab_2_peakExpr = median(C_Ab_2_peakExpr, na.rm = TRUE),
    C_Ab_3_peakExpr = median(C_Ab_3_peakExpr, na.rm = TRUE))

#########################################################################
test_GR <- test_GR %>%
  list(.) %>% 
  c(., grSamples_l) %>% 
  purrr::reduce(join_overlap_left_directed) %>% 
  mutate(unannot_peak = selfmatch(.)) %>%
  dplyr::mutate_at(dplyr::vars(dplyr::ends_with('peakExpr')), 
                   function(.) ifelse(is.na(.), 0, .)) %>% 
  select(-starts_with("sample_file")) %>% 
  group_by(unannot_peak)

my_reduced_Granges <- test_GR  %>% 
  reduce_ranges_directed

my_medianL <- map(smpls_lst, ~test_GR %>% 
      select(!!str_c(.x, "_peakExpr"), unannot_peak) %>% 
      as_tibble() %>% 
      group_by(unannot_peak) %>% 
      summarize(!!str_c(.x, "_peakExpr") := median(.data[[!!str_c(.x, "_peakExpr")]], na.rm = TRUE)))

my_medianL <- my_medianL %>%
  purrr::reduce(inner_join)

my_unAnnot_GR <- my_reduced_Granges %>% 
  as_tibble() %>% 
  inner_join(my_medianL) %>% 
  as_granges()
