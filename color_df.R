pop = c("Los_Guajes",
        "Crucero_Lagunitas",
        "El_Rodeo",
        "Amatlan_de_Canas", 
        "San_Lorenzo",
        "Palmar_Chico")
color = c("mediumaquamarine", 
          "gold",
          "darkorange", 
          "cornflowerblue",
          "purple",
          "violetred")

color_df <- 
  tibble(
    pop = paste0(c("LR_", "Teo_"), rep(pop, each = 2)),
    color = rep(color, each = 2)
    )

hexs <- gplots::col2hex(color_df$color)
names(hexs) <- color_df$pop


