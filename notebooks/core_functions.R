#core functions for parv local project. to load/clean data used in multiple analyses 
#assumes the r-environment conda has been loaded before use

library(tidyverse)


pop = c("Los_Guajes",
        "Crucero_Lagunitas",
        "El_Rodeo",
        "Amatlan_de_Canas", 
        "San_Lorenzo",
        "Palmar_Chico",
        "random1_Palmar_Chico",
        "random2_Palmar_Chico",    
        "random")
color = c("mediumaquamarine", 
          "gold",
          "darkorange", 
          "cornflowerblue",
          "purple",
          "violetred",
          "violetred",
          "violetred",
          "grey")

color_df <- 
  tibble(
    pop = paste0(c("LR_", "Teo_"), rep(pop, each = 2)),
    color = rep(color, each = 2)
    )

hexs <- gplots::col2hex(color_df$color)
names(hexs) <- color_df$pop


get_pi <- function(bp = 1000){
    bp_string = str_glue("{bp}BP_theta.txt")
    all_files = list.files(path = "../data/angsd_pi", full.names = TRUE)
    pi_files = all_files[grep(pattern = bp_string, all_files)]

    pi_df <- 
    pi_files[1:2] %>% 
        map_df(~{
           vroom::vroom(.x, delim = "\t",
           skip = 1,
           col_names = c("info", "chr", "WinCenter", "tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")) %>%
           na.omit() %>%
           filter(nSites > 100000 * 0.1) %>% 
            #filter(nSites > 1000 * 0.5) %>% 
           mutate(
                pi = tP/nSites,
                pop = str_remove_all(string = .x, "(../data/angsd_pi/|.100000BP_theta.txt|.1000000BP_theta.txt)"),
                ) 
        }) %>%
        mutate(ssp_pop = str_remove(pop, "v5--")) %>% 
        separate(col = pop, into = c("ref", "subspecies", "pop"), sep = "--") %>% 
        mutate(pop = str_replace(pop, "LR_|Teo_", "")) %>%
        full_join(., mutate(color_df, pop = str_replace(pop, "LR_|Teo_", "")), by = "pop") %>% 
        mutate(pop = str_replace_all(pop, "_", " ")) %>%
        separate(info, sep = "[\\)\\(,]", into = c(letters[1:7], "WinStart", "WinStop", "g2")) %>% 
        dplyr::select(-c(letters[1:7], g2)) %>%
        mutate(start = as.numeric(WinStart), 
           end = as.numeric(WinStop),
           size = end - start)
 

        return(pi_df)
}
