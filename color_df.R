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

col2hex <- function(colors){
  rgb_mat <- col2rgb(colors)
  sapply(1:ncol(rgb_mat), function(x) {
    c_set <- rgb_mat[,x]
    rgb(c_set[1], c_set[2], c_set[3], maxColorValue=255)
  })
}


hexs <- col2hex(color_df$color)
names(hexs) <- color_df$pop

make_legplot <- function(text_size = 2.5){
    color_df %>%
        mutate(pop = str_replace(pop, "LR_|Teo_", "")) %>% 
        filter(!pop %in% c("random1_Palmar_Chico", "random2_Palmar_Chico")) %>%
        mutate(pop = str_replace_all(pop, "_", " ")) %>% 
        mutate(pop = str_replace_all(pop, "Amatlan de Canas", "Amatlán de Cañas")) %>% 
        mutate(pop = str_replace_all(pop, " ", "\n")) %>%
        mutate(pop = str_replace_all(pop, "\nde", " de")) %>%
        mutate(pop = ifelse(pop == "random", "rangewide", pop)) %>% 
        arrange(pop) %>% 
        distinct() %>% 
        mutate(pt = rep(0, n()), idx = 1:n()) %>% 
        ggplot(aes(pop, pt, fill = color)) +
        geom_tile(color = "white", lwd = 1) +
        geom_text(aes(pop, -1.2, label = pop, hjust = "center"), size = text_size, check_overlap = TRUE, inherit.aes = FALSE) +
        scale_fill_identity() +
        theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),legend.position="none",
              panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank()) +
        ylim(-2.3, 1)
}



