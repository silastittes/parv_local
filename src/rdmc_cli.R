suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library(rdmc))
suppressPackageStartupMessages(library(tidyverse))


# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-s", "--sweep_file", type="character", help="File that contains the allele frequencies for the sweep region.")
parser$add_argument("-n", "--neutral_file", help="File that contains the allele frequencies for the neutral loci.")
parser$add_argument("-I", "--pop_ids", help="The indices for populations that experienced the sweep. Hyphen separated. ex. 1-2-5.")
parser$add_argument("-g", "--gmap", help="Genetic map to get recombination rate. Must contained tab delimited and named columns: 'chr', 'cm', 'pos'. Chromosomes will be prefixed with 'chr'.")
parser$add_argument("-o", "--out_file", help="Name of the outfile to write the fitted model data frame to.")
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
  ungroup() %>% 
  mutate(chr = paste0("chr", chr))

get_rr <- function(genetic_df, sweep_chr, sweep_positions){
  chr_df <- filter(genetic_df, chr == sweep_chr)
  median(approx(x = chr_df$pos, y = chr_df$rr, xout = sweep_positions)$y)
}

get_cm <- function(genetic_df, sweep_chr, sweep_start, sweep_end){
  chr_df <- filter(genetic_df, chr == sweep_chr)
  cm_start <- approx(x = chr_df$pos, y = chr_df$cm, xout = sweep_start)$y
  cm_end <- approx(x = chr_df$pos, y = chr_df$cm, xout = sweep_end)$y
  cm_end - cm_start
}



MIN_FREQ <- 1/20
DEFAULT_SITES <- 1e4
MAX_SITES <- 1e5
MIN_SITES <- 1e3
#snps per cM to get ~ constant density along different sized sweeps 
SNP_K  <- 250000 


neutral_file <- args$neutral_file
neutral_freqs <- vroom::vroom(file = neutral_file,   
    delim = "\t",
    col_names = FREQ_POPS) %>%
    mutate(varz = apply(select(., -c(chrom, start, end)), 1, max)) %>% 
    filter(varz >= MIN_FREQ) %>%
    sample_n(50000) %>% 
    select(-varz) 


s_file <- args$sweep_file

#get sweep start and end positions from file name
start <- str_split(s_file, "start", simplify = TRUE) %>% 
    `[`(2) %>% 
    str_split("_", simplify = TRUE) %>% 
    `[`(1) %>% 
    as.numeric(c)

end <- str_split(s_file, "end", simplify = TRUE) %>% 
    `[`(2) %>% 
    str_split("_", simplify = TRUE) %>% 
    `[`(1) %>% 
    as.numeric(c)

#get sweep size in cM. get how many sites to randomly sample along sweep
sweep_cM <- get_cm(gen_map_all_chr, "chr1", start, end)
n_snps <- round(SNP_K*sweep_cM)
if(is.na(n_snps)) n_snps <- DEFAULT_SITES

n_sites <- case_when(
    is.na(n_snps) ~ NA_real_,
    n_snps >= MIN_SITES && n_snps <= MAX_SITES ~ n_snps,
    n_snps < MIN_SITES ~ MIN_SITES,
    n_snps > MAX_SITES ~ MAX_SITES,
)



sel_vec <- str_split(args$pop_ids, "-", simplify = TRUE) %>% 
           as_vector() %>% 
           as.numeric()

sweep_file <- vroom::vroom(file = s_file,   
    delim = "\t",
    col_names = FREQ_POPS) 

if(nrow(sweep_file) > 0){

    sweep_freqs <-
    sweep_file %>% 
        mutate(varz = apply(select(., -c(chrom, start, end)), 1, max)) %>% 
        filter(varz >= MIN_FREQ) %>%
        select(-varz) %>% 
        sample_n(min(nrow(.), n_sites)) %>% 
        arrange(start)

    final_snp_count = nrow(sweep_freqs) 

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
        num_bins = 1000,
        sels = 10^seq(-5, -1, length.out = 15),
        times = c(1e2, 1e3, 1e4, 1e5),
        gs = 10^seq(-3, -1, length.out = 3),
        migs = 10^(seq(-3, -1, length.out = 2)),
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

    #combine data, scale cle relative to neutral
    all_mods <-
      bind_rows(
        ind_cle,
        mig_cle,
        sv_cle
      ) %>% 
        mutate(
            sel_pop_ids = paste(FREQ_POPS[sel_vec+3], collapse = "; "),
            neut_cle = unique(neut_cle$cle),
            n_snps = final_snp_count,
            sweepsize_cM = sweep_cM,
            sweep_rr = rr,
            sweep_start_bp = start,
            sweep_end_bp = end,
            sweep_size_bp = end - start 
        )

} else {
    all_mods <- tibble()    
}

#write to file
out_file <- args$out_file
write_delim(all_mods, out_file, delim = "\t")

