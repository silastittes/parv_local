suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library(rdmc))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ape))
theme_set(cowplot::theme_cowplot(15))
suppressPackageStartupMessages(library(patchwork))


# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-s", "--sweep_file", type="character", help="File that contains the allele frequencies for the sweep region.")
parser$add_argument("-n", "--neutral_file", help="File that containes the allele frequencies for the neutral loci.")
parser$add_argument("-I", "--pop_ids", help="The indices for populations that experienced the sweep. Hyphen separated. ex. 1-2-5.")
parser$add_argument("-g", "--gmap", help="Genetic map to get recombination rate. Must contained tab delimited and named columns: 'chr', 'cm', 'pos'. Chromosomes will be prefixed with 'chr'.")
#parser$add_argument("-", "--", help="")

args <- parser$parse_args()


FREQ_POPS = c(
    "chrom",
    "start",
    "end",
    "v5--LR--Amatlan_de_Canas",
    "v5--LR--Crucero_Lagunitas",
    "v5--LR--Los_Guajes",
    "v5--LR--random1_Palmar_Chico",
    "v5--LR--San_Lorenzo",
    "v5--Teo--Amatlan_de_Canas",
    "v5--Teo--Crucero_Lagunitas",
    "v5--Teo--El_Rodeo",
    "v5--Teo--Los_Guajes",
    "v5--Teo--random1_Palmar_Chico",
    "v5--Teo--San_Lorenzo"
)


neutral_file <- args$neutral_file
neutral_freqs <- vroom::vroom(file = neutral_file,   
    delim = "\t",
    col_names = FREQ_POPS) %>%
    sample_n(10000)


s_file <- args$sweep_file
str_split(s_file, "pops", simplify = FALSE)

sel_vec <- str_split(args$pop_ids, "-", simplify = TRUE) %>% 
           as_vector() %>% 
           as.numeric()

sweep_freqs <- vroom::vroom(file = s_file,   
    delim = "\t",
    col_names = FREQ_POPS) 

gmap <- args$gmap
gen_map_all_chr <- vroom::vroom(gmap, delim = "\t") %>% 
  drop_na() %>%
  mutate(cm = cm + abs(min(cm))) %>%
  group_by(chr) %>% 
  group_modify(~{
    df1 <- slice(.x, -nrow(.x))
    df2 <- slice(.x, -1)
    to_keep <- df2$cm > df1$cm & df2$pos > df1$pos
    df1 <- df1[to_keep, ]
    df2 <- df2[to_keep, ]
    cm_mb <- tibble(cm_mb = 1e6*(df2$cm - df1$cm)/(df2$pos - df1$pos))
    cm_bp <- tibble(rr = (df2$cm - df1$cm)/(df2$pos - df1$pos)/100)
    bind_cols(df2, cm_mb, cm_bp)
  }) %>% 
  mutate(chr = paste0("chr", chr))

get_rr <- function(genetic_df, sweep_chr, sweep_positions){
  chr_df <- filter(genetic_df, chr == sweep_chr)
  median(approx(x = chr_df$pos, y = chr_df$rr, xout = sweep_positions)$y)
}


pos_vec <- select(sweep_freqs, end) %>% pull(end)

sweep_mat <- sweep_freqs %>% 
    select(-c(chrom, start, end)) %>% 
    t()

neut_mat <- 
    neutral_freqs %>% 
    select(-c(chrom, start, end)) %>% 
    t()


sweep_chr <- sweep_freqs$chrom[1]
rr <- get_rr(gen_map_all_chr, sweep_chr, sweep_freqs$end)

param_list <-
  parameter_barge(
    Ne =  50000,
    rec = rr,
    neutral_freqs = neut_mat,
    selected_freqs = sweep_mat,
    selected_pops = sel_vec,
    positions = pos_vec,
    n_sites = 20,
    sample_sizes = rep(10, nrow(neut_mat)),
    num_bins = 100,
    sels = 10^seq(-4, 0, length.out = 20),
    times = c(1e2, 1e4, 1e6),
    gs = 10^seq(-3, -1, length.out = 3),
    migs = 10^(seq(-4, -2, length.out = 2)),      
    sources = sel_vec,
    locus_name = s_file,
    cholesky = TRUE
  )

mode_wrapper <- function(barge, mode) {
       cle_out <- try(mode_cle(barge, mode))
       if(class(cle_out)[1] == 'try-error'){
           barge$cholesky  <- FALSE
           cle_out <- suppressWarnings(mode_cle(barge, mode))
           barge$cholesky  <- TRUE
       }
    return(cle_out)
}

#fit composite likelihood models
print("neutral")
neut_cle <- mode_wrapper(param_list, mode = "neutral")
print("ind")
ind_cle <- mode_wrapper(param_list, mode = "independent")
print("standing")
sv_cle <- mode_wrapper(param_list, mode = "standing")
print("mig")
mig_cle <- mode_wrapper(param_list, mode = "migration")


all_mods <-
  bind_rows(
    ind_cle,
    mig_cle,
    sv_cle
  ) %>% 
    mutate(sel_pop_ids = paste(FREQ_POPS[sel_vec+3], collapse = "; "),
          neut_cle = unique(neut_cle$cle)) 

out_file <- args$out_file
write_delim(x = all_mods, file = out_file, delim = "\t")

