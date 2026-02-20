# Cluster disturbances and characterise disturbance regimes 
# Author: Sofía Miguel Romero
###############################################################################

#### Prepare session #### 
rm(list = ls(all.names = TRUE))

if ("pacman" %in% installed.packages()) {
  pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
}

pacman::p_load(
  conflicted, terra, tidyverse, tidyterra, here, factoextra, MASS, mclust,
  skimr, sf, tmap, maptiles, tmaptools, RColorBrewer, broom, cowplot, dlookr,
  ggbiplot, mapview, ggcorrplot, patchwork, ggspatial, hrbrthemes)

conflicts_prefer(dplyr::filter, dplyr::select, purrr::map)

#### Import functions and variables #### 

list.scales  <- c("_25km", "_50km", "_100km")
list.agents  <- c("fire", "wind_barkbeetle")
n.periods    <- 3   # Options: 1, 2, 3

source(here("lib/00_functions_jae.R"))

#### Loop across grid scales and agents ####

for (grid.scale in list.scales) {

  message("--------------------------------------------------")
  message("Processing scale: ", grid.scale)
  message("--------------------------------------------------")
  
  for (var.agent in list.agents) {
    
    message("Agent: ", var.agent)
    
    source(here("lib/00_variables_jae.R"))
    
    #### 1. Load and prepare regime data ####

    regime.file <- readRDS(here(
      "data",
      "reference",
      "patch_summary",
      paste0("regime_", n.periods, "_period", grid.scale, ".rds")
    ))
    
    regime.data <- regime.file %>%
      pull(period.regime.filter) %>%
      .[[1]] %>%
      map(
        ~ rename(.x, grid_id = grid_id2) %>%
          filter(agent == var.agent) %>%
          select(
            grid_id,
            agent,
            log_size,
            log_vcsize,
            log_sev,
            log_vcsev,
            log_sumfrequency,
            log_frequency,
            log_vcfrequency,
            log_rate,
            log_sumrate,
            log_sdrate,
            log_vcrate,
            log_sumrate,
            log_n_years
          ) %>%
          filter(across(everything(), is.finite))
        # %>% filter(!grid_id %in% grids_filtered)
      )
    
    #### 3. Preliminary analyses and visualisation ####
    glimpse(regime.data)
    distribution.plots <- plot_attribute_distribution(regime.data, write_out = TRUE)
    print(distribution.plots)
    correlations.plot <- plot_correlation(bind_rows(regime.data))
    print(correlations.plot)
    
    
    #### 4. Create principal component analysis #### 
    
    ## Prepare data for PCA
    # Filter NA's
    period.data <- regime.data %>%
      map(~ .x %>%
            filter(if_all(-grid_id, ~ !is.na(.))) %>%
            as_tibble()) %>%
      rename_periods(. , n.periods)
    # Select variables of interest
    period.subset <- period.data %>%
      map(\(period_subset)
          period_subset %>%
            select(all_of(keep.vars)))
    
    ## Fit PCA model
    pca.mod <-  period.subset[[reference.period]] %>%
      prcomp(scale = TRUE)
    
    saveRDS(pca.mod,
            file = here(
              'output',
              'model',
              'cluster_creation',
              paste0('pca_model_', var.agent, '_', grid.scale, ".rds")
            ))
    pca.mod <- readRDS(file = here(
      'output',
      'model',
      'cluster_creation',
      paste0('pca_model_', var.agent, '_', grid.scale, ".rds")
    ))
    
    ## Apply PCA transformation to other periods
    period.pca <- apply_pca_to_periods(n.periods,
                                       pca.mod,
                                       period.subset,
                                       reference.period,
                                       predict.period,
                                       pca_select)
    
    saveRDS(period.pca,
            file = here(
              'output',
              'model',
              'cluster_creation',
              paste0('decades_pca_', var.agent, '_', grid.scale, ".rds")
            ))
    period.pca <- readRDS(file = here(
      'output',
      'model',
      'cluster_creation',
      paste0('decades_pca_', var.agent, '_', grid.scale, ".rds")
    ))
    
    ## Visualize PCA
    plot_pca(pca.mod, write_out = FALSE)

    #### 5. Cluster analysis ####
    
    period.reference <- purrr::pluck(period.pca, paste0('period', reference.period))
    
    set.seed(666)
    mclust <- Mclust(period.reference, #  initialization=list(hcPairs=hcRandomPairs(period.reference,seed=1)),
                     ModelNames = var.modelname, G = var.g)
    summary(mclust)
    saveRDS(mclust,
            here(
              'output',
              'model',
              'cluster_creation',
              paste0('mclust_model_', var.agent, '_', grid.scale, ".rds")
            ))
    mclust <- readRDS(here(
      'output',
      'model',
      'cluster_creation',
      paste0('mclust_model_', var.agent, '_', grid.scale, ".rds")
    ))
    
    # Predict classification for other periods
    mclust.list <- generate_mclusters(mclust,
                                      period.pca,
                                      reference.period,
                                      predict.period,
                                      n.periods)
    mclust.class <- mclust.list$classification
    mclust.prob <- mclust.list$probability
    
    
    #### 6. Characterise resulting clusters #### 
    
    ##  Calculate the bootstrap distribution of mixture parameters
    set.seed(666)
    mclust.bootstrap <- MclustBootstrap(mclust, nboot = 9999, type = "bs")
    
    saveRDS(mclust.bootstrap,
            here(
              'output',
              'model',
              'cluster_creation',
              paste0(
                'MclustBootstrap_9999_',
                var.agent,
                '_',
                grid.scale,
                ".rds"
              )
            ))
    mclust.bootstrap <- readRDS(here(
      'output',
      'model',
      'cluster_creation',
      paste0(
        'MclustBootstrap_9999_',
        var.agent,
        '_',
        grid.scale,
        ".rds"
      )
    ))
    
    mclust.stats <- list(
      mean = summary(mclust.bootstrap, what = "ave", conf.level = 0.95),
      ci = summary(mclust.bootstrap, what = "ci", conf.level = 0.95),
      se = summary(mclust.bootstrap, what = "se", conf.level = 0.95)
    )
    
    print(mclust.stats$mean)
    
    mclust.stats$mean <- t(mclust.stats$mean$mean) %>%
      as_tibble() %>%
      mutate(group = as.factor(row.names(.)))
    
    
    # Plot PCA with clusters 
    
    loadings.pca <- pca.mod$rotation[, 1:2] %>%
      as.data.frame(.)
    loadings.pca$varnames <- rownames(loadings.pca)
    
    component_formulas <- apply(pca.mod$rotation[, 1:3], 2, function(component) {
      terms <- paste(round(component, 3), " * ", rownames(loadings.pca), sep = " ")
      formula <- paste(terms, collapse = " + ")
      return(formula)
    })
    
    component.equations <- paste0(
      "Principal Component ",
      1:length(component_formulas),
      ": ",
      component_formulas,
      collapse = "\n\n"
    )
    component.equations.plot <- ggplot() +
      annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = component.equations,
        size = 5,
        hjust = 0.5,
        vjust = 0.5
      ) +
      theme_void()
    component.equations.plot <- cowplot::ggdraw(component.equations.plot) +
      cowplot::draw_plot(logo.ggdraw,
                         y = 0.03,
                         x = 0.17,
                         scale = 0.7)
    
    plot(component.equations.plot)
    
    cowplot::save_plot(
      here(
        'output',
        'plots',
        'disturbance_analysis',
        paste0("component_formulas_", var.agent, '_', grid.scale, ".svg")
      ),
      component.equations.plot,
      base_width = 10,
      base_height =  5,
      base_asp = 1,
      dpi = 400
    )
    
    pca.classified <- map2(
      pc.combination,
      plot.titles,
      ~ create_ggbiplot(.x, .y, var.agent, var.colors, reference.period, '50km')
    ) %>%
      patchwork::wrap_plots(ncol = 1) +
      patchwork::plot_layout(guides = "collect") &
      theme(legend.position = "boottom", text = element_text(size = 15))
    saveRDS(pca.classified,
            here(
              'output',
              'plots',
              'disturbance_analysis',
              paste0("pca_", var.agent, '_', grid.scale, ".rds")
            ))
    
    pca.classified2 <- cowplot::ggdraw(pca.classified) +
      cowplot::draw_plot(
        logo.ggdraw,
        vjust = -0.13,
        hjust = 0.11,
        scale = 0.7
      )
    
    plot(pca.classified2)
    
    cowplot::save_plot(
      here(
        'output',
        'plots',
        'disturbance_analysis',
        paste0(
          "PCA_space_with_clusters_",
          var.agent,
          '_',
          grid.scale,
          ".svg"
        )
      ),
      pca.classified,
      ncol = 1,
      nrow = 3,
      base_asp = 4,
      dpi = 400
    )
    
    
    #### 7. Create basic figures #### 
    
    if (grid.scale == '_50km') {
      period.data.classified <- process_period_data(
        n.periods,
        period.data,
        period.pca,
        mclust.class,
        mclust.prob,
        decades.list,
        var.agent,
        grid.scale
      ) %>%
        mutate(
          class2 = case_when(
            class == 'B' & agent == 'fire' ~  'Large & frequent',
            class == 'A' & agent == 'fire' ~ 'Moderate',
            class == 'C' & agent == 'fire' ~   'Severe & rare',
            class == 'A' &
              agent == 'wind_barkbeetle' ~    'Mild & Rare',
            class == 'B' &
              agent == 'wind_barkbeetle' ~   'Severe & large',
            class == 'C' & agent == 'wind_barkbeetle' ~  'Moderate',
            class == 'D' &
              agent == 'wind_barkbeetle' ~ 'Frequent & small'
          )
        )
    } else if (grid.scale == '_25km') {
      period.data.classified <- process_period_data(
        n.periods,
        period.data,
        period.pca,
        mclust.class,
        mclust.prob,
        decades.list,
        var.agent,
        grid.scale
      ) %>%
        mutate(
          class2 = case_when(
            class == 'A' & agent == 'fire' ~   'Moderate',
            class == 'B' & agent == 'fire' ~  'Large & frequent' ,
            class == 'C' & agent == 'fire' ~ 'Severe & rare' ,
            
            class == 'A' &
              agent == 'wind_barkbeetle' ~  'Mild & Rare',
            class == 'B' &
              agent == 'wind_barkbeetle' ~ 'Severe & large',
            class == 'C' & agent == 'wind_barkbeetle' ~  'Moderate',
            class == 'D' &
              agent == 'wind_barkbeetle' ~  'Frequent & small'
          )
        )
    } else if (grid.scale == '_100km') {
      period.data.classified <- process_period_data(
        n.periods,
        period.data,
        period.pca,
        mclust.class,
        mclust.prob,
        decades.list,
        var.agent,
        grid.scale
      ) %>%
        mutate(
          class2 = case_when(
            class == 'C' & agent == 'fire' ~   'Moderate',
            class == 'A' & agent == 'fire' ~ 'Large & frequent',
            class == 'B' & agent == 'fire' ~ 'Severe & rare',
            
            class == 'B' &
              agent == 'wind_barkbeetle' ~  'Severe & large',
            class == 'C' &
              agent == 'wind_barkbeetle' ~  'Moderate' ,
            class == 'D' &
              agent == 'wind_barkbeetle' ~ 'Frequent & small',
            class == 'A' &
              agent == 'wind_barkbeetle' ~ 'Mild & Rare'
          )
        )
    }
    
  saveRDS(period.data.classified,
            here(
              'output',
              'model',
              'cluster_creation',
              paste0('clusters_database', var.agent, '_', grid.scale, ".rds")
            ))
  
  period.data.classified <- readRDS(here(
      'output',
      'model',
      'cluster_creation',
      paste0('clusters_database', var.agent, '_', grid.scale, ".rds")
    ))
    
    if (n.periods != 1) {
      pca.plot <- generate_pca_plot_for_periods(regime.data,
                                                period.data.classified,
                                                loadings.pca,
                                                var.agent,
                                                grid.scale)
    }
    
    periods.classified <- classify_periods(
      n.periods,
      europe,
      period.data.classified,
      class.levels,
      var.agent,
      grid.scale,
      decades.list,
      selector.forest,
      visualize.class = 'reclass',
      reclassify = TRUE
    )
    
    plot(periods.classified)
    
    cowplot::save_plot(
      here(
        'output',
        'plots',
        'disturbance_analysis',
        paste0("cluster_maps_", var.agent, '_', grid.scale, ".svg")
      ),
      periods.classified,
      base_width = 10,
      base_height =  5,
      base_asp = 1,
      dpi = 400
    )
    
    europe.classified <- readRDS(here(
      'output',
      'model',
      'cluster_creation',
      paste0('mclust_data_', var.agent, '_', grid.scale, "_perc.rds")
    )) %>% as.tibble()
    
    percentage_reclassed <- europe.classified %>%
      as.tibble() %>%
      group_by(class2, period) %>%
      summarise(
        num_diferentes = sum(class2 != reclass),
        porcentaje = (num_diferentes / n()) * 100
      )
    
    write.csv(percentage_reclassed,
              here(
                'output',
                'plots',
                'disturbance_analysis',
                paste0('grids_reclassified_', var.agent, '_', grid.scale, ".csv")
              ))
    
    europe.classified.p3 <- europe.classified %>%
      as.tibble() %>%
      filter(period == '2011-2023' & reclass != "Residual")  %>%
      left_join(
        period.data.classified %>%
          select(grid_id, period, log_size, log_sev, log_sumfrequency),
        by = c('grid_id', 'period')
      )
    
    period.boxplots <- period_boxplot(
      europe.classified.p3,
      c('Size', 'Severity', 'Frequency'),
      var.colors.boxplot,
      var.agent = var.agent,
      logo.ggdraw
    )
    plot(period.boxplots)
    
    saveRDS(period.boxplots,
            here(
              'output',
              'plots',
              'disturbance_analysis',
              paste0("boxplots_", var.agent, '_', grid.scale, ".rds")
            ))
    
    
    cowplot::save_plot(
      here(
        'output',
        'plots',
        'disturbance_analysis',
        paste0("boxplots_", var.agent, '_', grid.scale, ".svg")
      ),
      period.boxplots,
      base_width = 7,
      base_height =  8,
      base_asp = 1,
      dpi = 600
    )
    
    plot(period.boxplots)
    
    annual.regime <- readRDS(here(
      "data",
      "reference",
      "patch_summary",
      paste0("regime_", n.periods, "_period", grid.scale, ".rds")
    )) %>%
      pull(annual.regime.filter) %>% .[[1]]
    period3.regime <- summarise_period(annual.regime,
                                       vars = c('grid_id2', 'agent'),
                                       decades.list =  decades.list) %>%
      .[[3]] %>%
      ungroup() %>%
      filter(agent == var.agent) %>%
      rename(
        grid_id   = grid_id2,
        Severity  = sev,
        Size      = size,
        Frequency = sumfrequency
      ) %>%
      select(grid_id, Size, Severity, Frequency) %>%
      filter(across(everything(), is.finite))
    
    period3.regime.class <- europe.classified %>%
      filter(period == '2011-2023' & reclass != "Residual")  %>%
      left_join(period3.regime, by = 'grid_id')
    
    period.boxplots.absvalues <-  period_boxplot_absvalues(
      period3.regime.class,
      c('Size', 'Severity', 'Frequency'),
      var.colors.boxplot,
      var.agent = var.agent,
      logo.ggdraw
    )
    
    saveRDS(
      period.boxplots.absvalues,
      here(
        'output',
        'plots',
        'disturbance_analysis',
        paste0("boxplots_absvalues_", var.agent, '_', grid.scale, ".rds")
      )
    )
    
    cowplot::save_plot(
      here(
        'output',
        'plots',
        'disturbance_analysis',
        paste0("boxplots_absvalues_", var.agent, '_', grid.scale, ".svg")
      ),
      period.boxplots.absvalues,
      base_width = 7,
      base_height =  8,
      base_asp = 1,
      dpi = 600
    )
    
    legend_vertical <- ggplot(europe.classified, aes(y = probability, fill = class2)) +
      geom_histogram() +
      scale_fill_manual(values = var.colors.boxplot) +
      theme_minimal()
    
    cowplot::save_plot(
      here(
        'output',
        'plots',
        'disturbance_analysis',
        paste0("legend_", var.agent, '_vertical', ".svg")
      ),
      legend_vertical,
      base_width = 10,
      base_height =  5,
      base_asp = 1,
      dpi = 400
    )
    
    legend_horizontal <- ggplot(europe.classified, aes(y = probability, fill = class2)) +
      geom_histogram() +
      scale_fill_manual(values = var.colors.boxplot) +
      theme_minimal() +
      theme(legend.position = 'bottom')
    
    cowplot::save_plot(
      here(
        'output',
        'plots',
        'disturbance_analysis',
        paste0("legend_", var.agent, '_hor', ".svg")
      ),
      legend_horizontal,
      base_width = 10,
      base_height =  5,
      base_asp = 1,
      dpi = 400
    )
    
    europe.classified <- readRDS(here(
      'output',
      'model',
      'cluster_creation',
      paste0('mclust_data_', var.agent, '_', grid.scale, "_perc.rds")
    ))
    
    period.classified.shp <- map(1:n.periods, function(period_subset) {
      europe.period <- europe.classified %>%
        filter(period %in% paste(decades.list[[period_subset]][1], decades.list[[period_subset]][2], sep = '-') |
                 period == '0') %>%
        filter(!grid_id %in% selector.forest$grid_id &
                 probability >= 0.4)
      return(
        europe %>%
          dplyr::left_join(as.data.frame(vect(
            europe.period
          )), by = 'grid_id') %>%
          filter(!grid_id %in% selector.forest$grid_id) %>%
          mutate(across(
            c('class2', 'reclass'), ~ if_else(is.na(.), as.factor('Residual'), .)
          ))
      )
      
    })
    
    pca.grob <- ggplotGrob(
      pca.classified +
        envalysis::theme_publish() +
        theme(
          legend.position = "none",
          panel.background = element_rect(fill = "#FFFFFFB3", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.grid = element_blank()
        ) +
        ylim(c(-4, 4))
    )
    
    params <- list(
      var = rep('reclass', n.periods),
      var_data = period.classified.shp %>%
        map(
          ~ .x %>%
            st_as_sf() %>%
            group_by(reclass) %>%
            summarise(geometry = st_union(geometry), .groups = "drop")
        ),
      var_period = decades.list[1:n.periods],
      var_legend = rep('none', n.periods)
    )
    params$var_legend[n.periods] <- 'bottom'
    period.classified.plot <- purrr::pmap(params, create_cluster_map_basic)
    
    saveRDS(
      period.classified.plot[[3]],
      here(
        'output',
        'plots',
        'disturbance_analysis',
        paste0(
          "disturbance_map_3rdperiod_",
          var.agent,
          '_',
          grid.scale,
          ".rds"
        )
      )
    )
    
    cowplot::save_plot(
      here(
        'output',
        'plots',
        'disturbance_analysis',
        paste0(
          "map_period3_cluster_pca_",
          var.agent,
          '_',
          grid.scale,
          ".svg"
        )
      ),
      period.classified.plot[[3]],
      base_width = 20,
      base_height =  10,
      base_asp = 1,
      dpi = 400
    )
    
  }
}
