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


