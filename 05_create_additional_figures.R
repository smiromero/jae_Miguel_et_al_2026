###############################################################################
# Additional disturbance regime analysi. Here we calculate cluster transitions,
# trend analysis (Sen slope), statistical tests, and create additional tables 
# and figures. 
# Author: Sofía Miguel Romero
###############################################################################

#### Prepare session  ####
rm(list = ls(all.names = TRUE))
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Load required packages
pacman::p_load( conflicted, terra, tidyverse, tidyterra, here, factoextra,
  MASS, mclust, skimr, sf, tmap, maptiles, tmaptools, RColorBrewer, broom,
  cowplot, dlookr, ggbiplot, mapview, ggcorrplot, patchwork, ggspatial,
  hrbrthemes, envalysis, paletteer, fs, treemapify, gt, gtExtras, patchwork,
  circlize, gt, circlize )

# Resolve conflicts of common functions preferring some from dplyr/purrr
conflicts_prefer(dplyr::filter, dplyr::select, purrr::map)

##### Load functions and variables ###########################################

list.scales <- c("_50km")
list.agents <- c("fire", "wind_barkbeetle")
n.periods <- 3

##### Loop by scale and agent ###########################################

for (grid.scale in list.scales) {
  message("Analizando escala de grid: ", grid.scale)
  
  for (var.agent in list.agents) {
    message(" - Agente: ", var.agent)
    
    source(here("lib/00_functions.R"))
    source(here("lib/00_variables.R"))
    
    #### 1. Load and prepare cluster pattern data ####
    
    pattern <- readRDS(here(
      "output",
      "model",
      "cluster_creation",
      paste0("mclust_data_", var.agent, '_', grid.scale, "_perc.rds")
    ))
    
    # Convert the data from long to wide format to see the class in each period by grid
    pattern.wide <- pattern %>%
      dplyr::select(-probability, -class2) %>%
      mutate(across(c(period, reclass), as.character)) %>%
      pivot_wider(names_from = period, values_from = reclass) %>%
      mutate(change = as.factor(if_else(
        `1985-1997` == `2011-2023`, NA, `2011-2023`
      )), `1985-1997` = `1985-1997`)
    
    #### 2. Classification maps for initial period (P1) and middle period (P2) ####
    
    period1.map <- create_cluster_map_basic(
      var_data   = pattern.wide %>%
        mutate(period1 = as.factor(`1985-1997`)) %>%
        group_by(period1) %>%
        summarise(geometry = st_union(geometry)),
      var_period = c("", ""),
      var_legend = "none",
      var        = "period1"
    )
    
    # Create map of the cluster pattern in the second period (1998-2010)
    period2.map <- create_cluster_map_basic(
      var_data   = pattern.wide %>%
        mutate(period2 = as.factor(`1998-2010`)) %>%
        group_by(period2) %>%
        summarise(geometry = st_union(geometry)),
      var_period = c("", ""),
      var_legend = "none",
      var        = "period2"
    )
    
    # Save the maps of period 1 and period 2 in SVG format
    cowplot::save_plot(
      here(
        "output",
        "plots",
        "disturbance_analysis",
        
        paste0("period1_map_", var.agent, "_", grid.scale, ".svg")
      ),
      period1.map,
      dpi = 600
    )
    cowplot::save_plot(
      here(
        "output",
        "plots",
        "disturbance_analysis",
        
        paste0("period2_map_", var.agent, "_", grid.scale, ".svg")
      ),
      period2.map,
      dpi = 600
    )
    
    #### 3. Prepare cluster transition data between periods ####
    
    pattern.change <- terra::aggregate(terra::vect(pattern.wide %>%
                                                     filter(!is.na(change))),
                                       by = "change") %>%
      rename(grid_id = mean_grid_id) %>%
      dplyr::select(grid_id, change, 'X1985-1997', 'X1998-2010', 'X2011-2023')
    
    pattern.dissolve <- terra::aggregate(
      terra::vect(pattern.wide %>% mutate(`1985-1997` = as.factor(`1985-1997`))),
      by = "1985-1997")
    
    #### 4. Maps of disturbance pattern transitions ####
    
    p1.transitions.map <- cluster_transitions_map(
      pattern = pattern,
      pattern.dissolve   = pattern.dissolve,
      pattern.change     = pattern.change,
      world.boundaries   = world.boundaries,
      europe.boundaries  = europe.boundaries,
      var.colors.boxplot = var.colors.boxplot
    ) + labs(y = 'Cluster changes', fill = "")
    
    cowplot::save_plot(
      here(
        "output",
        "plots",
        "disturbance_analysis",
        paste0(
          "pattern_transitions_map_with_period1_",
          var.agent,
          "_",
          grid.scale,
          ".svg"
        )
      ),
      p1.transitions.map,
      ncol = 1,
      nrow = 1,
      base_asp = 4,
      dpi = 400
    )
    
    
    #### 5. Cluster map for the first period (for reference) ####
    
    pattern.1st.period.map <- create_cluster_map_basic(
      var_data   = pattern.wide %>% mutate(x1st = as.factor(`1985-1997`)),
      var_period = c("", ""),
      var_legend = "none",
      var        = "x1st"
    )
    saveRDS(
      pattern.1st.period.map,
      here(
        "output",
        "plots",
        "disturbance_analysis",
        paste0(
          "pattern_1stperiod_cluster_",
          var.agent,
          "_",
          grid.scale,
          ".rds"
        )
      )
    )
    cowplot::save_plot(
      here(
        "output",
        "plots",
        "disturbance_analysis",
        paste0(
          "pattern_1stperiod_cluster_",
          var.agent,
          "_",
          grid.scale,
          ".svg"
        )
      ),
      pattern.1st.period.map,
      ncol = 1,
      nrow = 1,
      base_asp = 4,
      dpi = 400
    )

    #### 6. Calculation of annual trends (Sen's slope) in disturbance attributes ####
    
    annual.regime <- readRDS(here(
      "data",
      "reference",
      "patch_summary",
      paste0("regime_", n.periods, "_period", grid.scale, ".rds")
    )) %>%
      pull(annual.regime.filter) %>% .[[1]] %>%
      filter(agent == var.agent) %>%
      rename(
        grid_id   = grid_id2,
        Severity  = mean_sev,
        Size      = mean_size,
        Frequency = frequency
      ) %>%
      select(grid_id, year, Size, Severity, Frequency) %>%
      filter(across(everything(), is.finite)) %>%
      mutate(
        period = case_when(
          between(year, 1985, 1997) ~ "1985-1997",
          between(year, 1998, 2010) ~ "1998-2010",
          between(year, 2011, 2023) ~ "2011-2023"
        ) %>% factor(levels = c(
          "1985-1997", "1998-2010", "2011-2023"
        ))
      )
    
    annual.regime.long <- annual.regime %>%
      pivot_longer(
        cols = -c(grid_id, year, period),
        names_to = "variable",
        values_to = "valor"
      ) %>%
      mutate(variable = factor(variable, levels = c("Size", "Severity", "Frequency"))) %>%
      arrange(grid_id, variable, year)
    
    annual.trends <- annual.regime.long %>%
      group_by(grid_id, variable) %>%
      filter(n_distinct(year) >= 10) %>%
      summarise(
        trend_abs = round(trend::sens.slope(valor, conf.level = 0.95)$estimates, 4),
        trend_relative = round((
          trend::sens.slope(valor)$estimate / mean(valor)
        ) * 100, 2),
        p_value = trend::sens.slope(valor, conf.level = 0.95)$p.value,
        .groups = "drop"
      ) %>%
      mutate(trend_abs = if_else(variable == "Severity", trend_abs * 100, trend_abs)) #%>%
    # filter( p_value <= 0.05)
    
    
    annual.trends.params <- list(
      trend_df = list(annual.trends, annual.trends, annual.trends),
      var      = c("Size", "Severity", "Frequency"),
      ylab     = list("Size (ha)", "Severity (0-1)", expression_ha_year)  #ç
    )
    
    annual.trends.maps <- purrr::pmap(annual.trends.params, cluster_trend_map)
    
    trend.combined <- cowplot::plot_grid(
      annual.trends.maps[[1]] + theme(plot.margin = margin(0, 0, 0, 0)),
      annual.trends.maps[[2]] + theme(plot.margin = margin(0, 0, 0, 0)),
      annual.trends.maps[[3]] + theme(plot.margin = margin(0, 0, 0, 0)),
      nrow = 2,
      align = "hv",
      axis = "tblr"
    )
    
    
    trend.combined.transitions <- cowplot::plot_grid(
      p1.transitions.map + theme(plot.margin = margin(0, 0, 0, 0)) ,
      annual.trends.maps[[1]] + theme(plot.margin = margin(0, 0, 0, 0)),
      annual.trends.maps[[2]] + theme(plot.margin = margin(0, 0, 0, 0)),
      annual.trends.maps[[3]] + theme(plot.margin = margin(0, 0, 0, 0)),
      nrow = 2,
      align = "hv",
      axis = "tblr",
      greedy = TRUE
    )
    
    cowplot::save_plot(
      here(
        "output",
        "plots",
        "disturbance_analysis",
        paste0(
          "disturbance_trends_combined_transitions_",
          var.agent,
          "_",
          grid.scale,
          ".svg"
        )
      ),
      trend.combined.transitions,
      base_width = 12,
      base_height = 10,
      dpi = 600
    )
    
    cowplot::save_plot(
      here(
        "output",
        "plots",
        "disturbance_analysis",
        paste0(
          "trend_combined_pvalue",
          var.agent,
          "_",
          grid.scale,
          ".svg"
        )
      ),
      trend.combined,
      base_width = 8,
      base_height = 15,
      dpi = 600
    )
    saveRDS(trend.combined,
            here(
              "output",
              "plots",
              "disturbance_analysis",
              paste0(
                "disturbance_trends_combined_",
                var.agent,
                "_",
                grid.scale,
                ".rds"
              )
            ))
    
    save.names <- c("size", "severity", "frequency")
    purrr::walk2(annual.trends.maps, save.names, function(mapa, name) {
      cowplot::save_plot(
        here(
          "output",
          "plots",
          "disturbance_analysis",
          paste0(
            "annual_trend_",
            name,
            "_",
            var.agent,
            "_",
            grid.scale,
            ".svg"
          )
        ),
        mapa,
        ncol = 1,
        nrow = 1,
        base_asp = 4,
        dpi = 400
      )
    })
    
    #### 8. Histograms of relative trends by biome ####
    
    annual.trends.hist <- annual.trends %>%
      left_join(europe_df, by = "grid_id") %>%
      drop_na(biome3)
    
    annual.trends.biome.hist <- ggplot(annual.trends.hist,
                                       aes(
                                         x = trend_relative,
                                         fill = biome3,
                                         color = biome3
                                       )) +
      geom_histogram(
        aes(y = after_stat(density)),
        position = "identity",
        bins = 30,
        alpha = 0.4,
        size = 0.5
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "dotted",
        size = 1,
        color = "grey30"
      ) +
      facet_wrap(~ variable, scales = "free", ncol = 1) +
      scale_fill_manual(values = biome.colors, name = NULL) +
      scale_color_manual(values = biome.colors, name = NULL) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_x_continuous(expand = c(0, 0)) +
      envalysis::theme_publish(base_size = 12) +
      theme(
        axis.text       = element_text(color = "grey30"),
        axis.title      = element_text(color = "grey30"),
        axis.ticks      = element_line(color = "grey30"),
        axis.line       = element_line(color = "grey30"),
        strip.text      = element_text(size = 12, color = "grey30"),
        legend.text     = element_text(color = "grey30"),
        legend.title    = element_blank(),
        legend.position = "bottom",
        panel.grid.minor    = element_blank(),
        panel.grid.major.y  = element_line(color = "grey90"),
        panel.background    = element_rect(fill = "transparent", color = NA),
        plot.background     = element_rect(fill = "transparent", color = NA),
        legend.background   = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        plot.margin         = margin(0, 10, 0, 0),
        panel.spacing       = unit(0, "pt")
      ) +
      labs(x = "Trend (%)", y = "Density", fill = "")
    
    annual.trends.total.hist <- ggplot(annual.trends.hist, aes(x = trend_relative, 
                                                               fill = 'grey10')) +
      geom_histogram(
        aes(y = after_stat(density)),
        position = "identity",
        bins = 30,
        alpha = 0.4,
        size = 0.5
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "dotted",
        size = 1,
        color = "grey30"
      ) +
      facet_wrap(~ variable, scales = "free_y", ncol = 1) +
      scale_fill_manual(values = biome.colors, name = NULL) +
      scale_color_manual(values = biome.colors, name = NULL) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_x_continuous(expand = c(0, 0), limits = c(-10, 10)) +
      envalysis::theme_publish(base_size = 12) +
      theme(
        axis.text       = element_text(color = "grey30"),
        axis.title      = element_text(color = "grey30"),
        axis.ticks      = element_line(color = "grey30"),
        axis.line       = element_line(color = "grey30"),
        strip.text      = element_text(size = 12, color = "grey30"),
        legend.text     = element_text(color = "grey30"),
        legend.title    = element_blank(),
        legend.position = "bottom",
        panel.grid.minor    = element_blank(),
        panel.grid.major.y  = element_line(color = "grey90"),
        panel.background    = element_rect(fill = "transparent", color = NA),
        plot.background     = element_rect(fill = "transparent", color = NA),
        legend.background   = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        plot.margin         = margin(0, 10, 0, 0),
        panel.spacing       = unit(0, "pt")
      ) +
      labs(x = "Trend (%)", y = "Density", fill = "")
    
    saveRDS(
      annual.trends.biome.hist,
      here(
        "output",
        "plots",
        "disturbance_analysis",
        paste0("biome_trend_hist_", var.agent, "_", grid.scale, ".rds")
      )
    )
    cowplot::save_plot(
      here(
        "output",
        "plots",
        "disturbance_analysis",
        paste0("biome_trend_hist_", var.agent, "_", grid.scale, ".svg")
      ),
      annual.trends.biome.hist,
      base_width = 4,
      base_height = 6,
      dpi = 400
    )
    
    saveRDS(
      annual.trends.total.hist,
      here(
        "output",
        "plots",
        "disturbance_analysis",
        paste0("trend_hist_", var.agent, "_", grid.scale, ".rds")
      )
    )
    cowplot::save_plot(
      here(
        "output",
        "plots",
        "disturbance_analysis",
        paste0("trend_hist_", var.agent, "_", grid.scale, ".svg")
      ),
      annual.trends.total.hist,
      base_width = 4,
      base_height = 6,
      dpi = 400
    )
    
    #### 9. Summary of trends at the European level and by biome ####
    
    # Trend for Europe
    
    annual.trends.europe <- calculate_trends(data = annual.regime.long,
                                             var_agent = var.agent,
                                             nombre_archivo = "europe_trend_") %>%
      mutate(nivel = "Europa",
             periodo = "1985–2023",
             agent = var.agent)
    
    # Period 1–2 Europe
    annual.trends.europe.period1to2 <- calculate_trends(
      data = annual.regime.long,
      periodo = c(1985, 2010),
      var_agent = var.agent,
      nombre_archivo = "europe_trend_1_to_2"
    ) %>%
      mutate(nivel = "Europa",
             periodo = "1985–2010",
             agent = var.agent)
    
    # Period 2–3 Europe
    annual.trends.europe.period2to3 <- calculate_trends(
      data = annual.regime.long,
      periodo = c(1998, 2023),
      var_agent = var.agent,
      nombre_archivo = "europe_trend_2_to_3"
    ) %>%
      mutate(nivel = "Europa",
             periodo = "1998–2023",
             agent = var.agent)
    
    # Trend by biome
    annual.trends.biome <- calculate_trends(
      data = annual.regime.long,
      europe_df = europe_df,
      agrupar_por_bioma = TRUE,
      var_agent = var.agent,
      nombre_archivo = "biome_trend_"
    ) %>%
      mutate(nivel = "Bioma",
             periodo = "1985–2023",
             agent = var.agent)
    
    # Period 1–2 by biome
    annual.trends.biome.period1to2 <- calculate_trends(
      data = annual.regime.long,
      europe_df = europe_df,
      periodo = c(1985, 2010),
      agrupar_por_bioma = TRUE,
      var_agent = var.agent,
      nombre_archivo = "biome_trend_1_to_2"
    ) %>%
      mutate(nivel = "Bioma",
             periodo = "1985–2010",
             agent = var.agent)
    
    # Period 2–3 by biome
    annual.trends.biome.period2to3 <- calculate_trends(
      data = annual.regime.long,
      europe_df = europe_df,
      periodo = c(1998, 2023),
      agrupar_por_bioma = TRUE,
      var_agent = var.agent,
      nombre_archivo = "biome_trend_2_to_3"
    ) %>%
      mutate(nivel = "Bioma",
             periodo = "1998–2023",
             agent = var.agent)
    
    all.trends <- bind_rows(
      annual.trends.europe,
      annual.trends.europe.period1to2,
      annual.trends.europe.period2to3,
      annual.trends.biome,
      annual.trends.biome.period1to2,
      annual.trends.biome.period2to3
    )
    
    readr::write_csv(all.trends,
                     here(
                       "output",
                       "plots",
                       "disturbance_analysis",
                       
                       paste0("all.trends.", var.agent, ".csv")
                     ))
    
    
    
    print(annual.trends.europe)
    print(annual.trends.biome)
    
    summary_biome <- summarise_trends_by_group(
      group_var = "biome3",
      annual_regime_long = annual.regime.long,
      annual_trends_hist = annual.trends.hist,
      group_colors = biome.colors,
      var_agent = var.agent
    )
    summary_reclass <- summarise_trends_by_group(
      group_var = "reclass",
      annual_regime_long = annual.regime.long %>%
        left_join(
          pattern %>% as.data.frame() %>% select(grid_id, period, reclass),
          by = c("grid_id", "period")
        ) %>%
        filter(reclass != "Residual"),
      annual_trends_hist = NULL,
      group_colors = var.colors.boxplot,
      var_agent = var.agent
    )
    
    print(summary_biome)
    print(summary_reclass)
    
    
    #### 12: Summary by period and cluster type ####
    
    annual.regime2 <- readRDS(here(
      "data",
      "reference",
      "patch_summary",
      paste0("regime_", n.periods, '_period', grid.scale, ".rds")
    )) %>%
      pull(annual.regime.filter) %>%
      .[[1]] %>%
      filter(agent == var.agent) %>%
      rename('grid_id' = 'grid_id2') %>%
      filter(across(everything(), is.finite)) %>%
      mutate(
        period = case_when(
          between(year, 1985, 1997) ~ '1985-1997' ,
          between(year, 1998, 2010) ~ '1998-2010' ,
          between(year, 2011, 2023) ~ '2011-2023'
        ) %>%
          factor()
      )
    
    
    period.summary <- summarise_period(annual.regime2,
                                       vars = c('grid_id', 'agent'),
                                       decades.list =  decades.list) %>%
      map2(., 1:n.periods, ~ .x %>%
             ungroup() %>%
             mutate(period = .y)) %>%
      bind_rows() %>%
      mutate(
        period = case_when(
          period == 1 ~ '1985-1997' ,
          period == 2 ~ '1998-2010' ,
          period == 3 ~ '2011-2023'
        ) %>%
          factor()
      ) %>%
      arrange(grid_id, period) %>%
      left_join(as.data.frame(pattern), by = c('grid_id', 'period'))
    
    
    grid_summary_stats_europe <- period.summary %>%
      group_by(period) %>%
      summarise(
        mean_size = mean(size, na.rm = TRUE),
        p95_size = quantile(size, 0.95, na.rm = TRUE),
        mean_sev = mean(sev, na.rm = TRUE),
        mean_frequency = mean(frequency, na.rm = TRUE),
        mean_sumfrequency = mean(sumfrequency, na.rm = TRUE),
        mean_sum_area = mean(sum_area, na.rm = TRUE)
      )
    grid_summary_stats_europe
    
    write.csv(
      grid_summary_stats_europe,
      file = here(
        "output",
        "plots",
        "disturbance_analysis",
        
        paste0("table1_europe_statistics_", var.agent, ".csv")
      ),
      row.names = FALSE
    )
    grid_summary_stats_biome <- period.summary %>%
      left_join(europe_df, by = 'grid_id') %>%
      group_by(period, biome3) %>%
      summarise(
        mean_size = mean(size, na.rm = TRUE),
        p95_size = quantile(size, 0.95, na.rm = TRUE),
        mean_sev = mean(sev, na.rm = TRUE),
        mean_frequency = mean(frequency, na.rm = TRUE),
        mean_sumfrequency = mean(sumfrequency, na.rm = TRUE),
        
        mean_sum_area = mean(sum_area, na.rm = TRUE)
      ) %>%
      arrange(biome3)
    
    grid_summary_stats_biome
    
    write.csv(
      grid_summary_stats_biome,
      file = here(
        "output",
        "plots",
        "disturbance_analysis",
        
        paste0("table1_biome_statistics_", var.agent, ".csv")
      ),
      row.names = FALSE
    )
    
    
    
    grid_summary_stats <- period.summary %>%
      group_by(period, reclass) %>%
      summarise(
        mean_size = mean(size, na.rm = TRUE),
        sd_size = sd(size, na.rm = TRUE),
        p5_size = quantile(size, 0.05, na.rm = TRUE),
        p95_size = quantile(size, 0.95, na.rm = TRUE),
        
        mean_sev = mean(sev, na.rm = TRUE),
        sd_sev = sd(sev, na.rm = TRUE),
        p5_sev = quantile(sev, 0.05, na.rm = TRUE),
        p95_sev = quantile(sev, 0.95, na.rm = TRUE),
        
        mean_sumfrequency = mean(sumfrequency, na.rm = TRUE),
        sd_sumfrequency = sd(sumfrequency, na.rm = TRUE),
        p5_sumfrequency = quantile(sumfrequency, 0.05, na.rm = TRUE),
        p95_sumfrequency = quantile(sumfrequency, 0.95, na.rm = TRUE),
        
        mean_frequency = mean(frequency, na.rm = TRUE),
        sd_frequency = sd(frequency, na.rm = TRUE),
        p5_frequency = quantile(frequency, 0.05, na.rm = TRUE),
        p95_frequency = quantile(frequency, 0.95, na.rm = TRUE),
        
        mean_sum_area = mean(sum_area, na.rm = TRUE),
        percent_grids = (n() * 100) / 2415
      ) %>%
      mutate(
        size_stats = paste0(
          round(mean_size, 2),
          " ± ",
          round(sd_size, 2),
          " (",
          round(p5_size, 2),
          ", ",
          round(p95_size, 2),
          ")"
        ),
        sev_stats = paste0(
          round(mean_sev, 2),
          " ± ",
          round(sd_sev, 2),
          " (",
          round(p5_sev, 2),
          ", ",
          round(p95_sev, 2),
          ")"
        ),
        sumfrequency_stats = paste0(
          round(mean_sumfrequency, 2),
          " ± ",
          round(sd_sumfrequency, 2),
          " (",
          round(p5_sumfrequency, 2),
          ", ",
          round(p95_sumfrequency, 2),
          ")"
        ),
        frequency_stats = paste0(
          round(mean_frequency, 2),
          " ± ",
          round(sd_frequency, 2),
          " (",
          round(p5_frequency, 2),
          ", ",
          round(p95_frequency, 2),
          ")"
        )
      ) %>%
      select(
        period,
        reclass,
        size_stats,
        sev_stats,
        frequency_stats,
        sumfrequency_stats,
        mean_sum_area,
        percent_grids
      )
    
    saveRDS(grid_summary_stats,
            here(
              'output',
              'plots',
              'disturbance_analysis',
              
              paste0("table1_statistics", var.agent, ".rds")
            ))
    
    gt_stats <- grid_summary_stats %>%
      gt() %>%
      cols_label(
        period = "Period",
        reclass = "Cluster",
        size_stats = "Size (ha)",
        sev_stats = "Severity (0-1)",
        frequency_stats = "Frequency (events year-1)",
        sumfrequency_stats = "Sum Frequency (events period-1)",
        mean_sum_area = 'Total area (ha)',
        percent_grids = "Grids (%)"
      ) %>%
      tab_style(
        style = cell_text(
          font = c(google_font(name = "Arial"), default_fonts()),
          size = px(11)
        ),
        locations = cells_body(columns = everything())
      ) %>%
      tab_style(
        style = cell_text(
          font = c(google_font(name = "Arial"), default_fonts()),
          weight = "bold",
          size = px(12)
        ),
        locations = cells_column_labels(columns = everything())
      ) %>%
      tab_style(
        style = cell_text(
          font = c(google_font(name = "Arial"), default_fonts()),
          weight = "bold",
          size = px(12)
        ),
        locations = cells_row_groups(groups = everything())
      ) %>%
      cols_align(align = "center", columns = everything()) %>%
      tab_spanner(
        label = "Disturbance",
        columns = c(
          "size_stats",
          "sev_stats",
          "frequency_stats",
          "sumfrequency_stats",
          "mean_sum_area"
        )
      ) %>%
      tab_spanner(label = "Grids", columns = "percent_grids") %>%
      tab_style(style = cell_text(
        font = c(google_font(name = "Arial"), default_fonts()),
        size = px(12)
      ),
      locations = cells_column_spanners()) %>%
      gt::fmt_number(
        columns = c(
          "size_stats",
          "sev_stats",
          "frequency_stats",
          "sumfrequency_stats",
          "mean_sum_area",
          "percent_grids"
        ),
        decimals = 2
      )
    
    gt_stats %>%
      gtsave(here(
        'output',
        'plots',
        'disturbance_analysis',
        paste0("table1_statistics_", var.agent, ".html")
      ))
    
    write.csv(
      grid_summary_stats,
      file = here(
        "output",
        "plots",
        "disturbance_analysis",
        
        paste0("table1_statistics_", var.agent, ".csv")
      ),
      row.names = FALSE
    )
    
    #### 10.1 Tests for significant differences between biomes (Wilcoxon)) ####
    
    biome.diff <- annual.regime %>%
      left_join(europe_df, by = "grid_id") %>%
      drop_na(biome3) %>%
      group_by(grid_id) %>%
      summarise(
        biome3    = unique(biome3),
        biome7    = unique(biome7),
        Size      = mean(Size, na.rm = TRUE),
        Severity  = mean(Severity, na.rm = TRUE),
        Frequency = mean(Frequency, na.rm = TRUE),
        .groups   = "drop"
      )
    
    wilcoxon.results.biome3 <- wilcoxon_tests(data = biome.diff, group_var = "biome3")
    wilcoxon.results.biome7 <- wilcoxon_tests(data = biome.diff, group_var = "biome7")
    
    wilcoxon_subset <- wilcoxon.results.biome7 %>%
      filter(
        group1 %in% c("broadleafMixed", "mediterranean", "boreal") &
          group2 %in% c("broadleafMixed", "mediterranean", "boreal")
      )
    print(wilcoxon_subset)
    
    write.csv(
      wilcoxon.results.biome3,
      file = here(
        "output",
        "plots",
        "disturbance_analysis",
        
        paste0(var.agent, "_wilcoxon_results_biome3.csv")
      ),
      row.names = FALSE
    )
    write.csv(
      wilcoxon.results.biome7,
      file = here(
        "output",
        "plots",
        "disturbance_analysis",
        
        paste0(var.agent, "_wilcoxon_results_biome7.csv")
      ),
      row.names = FALSE
    )
    
    
    #### 10.2 Tests for significant differences between regimes (Wilcoxon) ####
    period.summary.wilcoxon <- period.summary %>%
      filter(reclass != 'Residual') %>%
      rename(Size = size,
             Severity = sev,
             Frequency = frequency)
    wilcoxon.results.reclass <- wilcoxon_tests(data = period.summary.wilcoxon, group_var = "reclass")
    write.csv(
      wilcoxon.results.reclass,
      file = here(
        "output",
        "plots",
        "disturbance_analysis",
        
        paste0(var.agent, "_wilcoxon_results_reclass.csv")
      ),
      row.names = FALSE
    )
    
    #### 11. Additional exploratory analysis: relationships between variables and 
    cluster/biome ####
    
    annual.reg <- annual.regime %>%
      left_join(europe_df, by = "grid_id") %>%
      drop_na(biome3) %>%
      group_by(grid_id) %>%
      summarise(
        biome3    = unique(biome3),
        Size      = mean(Size, na.rm = TRUE),
        Severity  = mean(Severity, na.rm = TRUE),
        Frequency = mean(Frequency, na.rm = TRUE),
        .groups   = "drop"
      )
    
    annual.regimes.class.wide <- annual.reg %>%
      select(-biome3) %>%
      right_join(pattern.wide %>% as.data.frame() %>% select(-change), by = "grid_id") %>%
      left_join(europe_df, by = "grid_id") %>%
      mutate(
        `1985-1997` = as.factor(`1985-1997`),
        `1998-2010` = as.factor(`1998-2010`),
        `2011-2023` = as.factor(`2011-2023`)
      )
    
    annual.regimes.class <- annual.regimes.class.wide %>%
      pivot_longer(
        cols = c(`1985-1997`, `1998-2010`, `2011-2023`),
        names_to = "period",
        values_to = "cluster"
      ) %>%
      mutate(
        period = factor(period, levels = c(
          "1985-1997", "1998-2010", "2011-2023"
        )),
        period_abbr = recode(
          period,
          "1985-1997" = "P1",
          "1998-2010" = "P2",
          "2011-2023" = "P3"
        )
      )
    
    ggpairs.biome <- GGally::ggpairs(
      data   = annual.reg,
      columns = c("biome3", "Size", "Severity", "Frequency"),
      title   = paste0("Biome summary (", var.agent, ")"),
      mapping = aes(color = biome3)
    ) +
      scale_color_manual(values = biome.colors) +
      scale_fill_manual(values = biome.colors)
    
    cowplot::save_plot(
      here(
        "output",
        "plots",
        "disturbance_analysis",
        
        paste0("ggpairs_biome_", var.agent, ".svg")
      ),
      ggpairs.biome,
      base_width = 12,
      base_height = 10
    )
    
    ggpairs.reclass <- GGally::ggpairs(
      data   = annual.regimes.class,
      columns = c(
        "Size",
        "Severity",
        "Frequency",
        "period_abbr",
        "biome3",
        "cluster"
      ),
      title   = paste0("Cluster summary (", var.agent, ")"),
      mapping = aes(color = cluster)
    )
    
    if (var.agent == 'fire') {
      annual.regimes.class <- annual.regimes.class %>%
        mutate(
          cluster = factor(
            cluster,
            levels = c("Residual", "Moderate", 'Severe & rare', 'Large & frequent')
          ),
          biome3 = factor(biome3, levels = c(
            'Boreal', 'Temperate', 'Mediterranean'
          ))
        ) %>%
        arrange(cluster, biome3, period_abbr)
    } else if (var.agent == 'wind_barkbeetle') {
      annual.regimes.class <- annual.regimes.class %>%
        mutate(
          cluster = factor(
            cluster,
            levels = c(
              "Residual",
              "Moderate",
              "Frequent & small",
              'Mild & Rare',
              'Severe & large'
            )
          ),
          biome3 = factor(biome3, levels = c(
            'Boreal', 'Temperate', 'Mediterranean'
          ))
        ) %>%
        arrange(cluster, biome3, period_abbr)
    }
    
    biome.period.reclass <- ggplot(annual.regimes.class, aes(x = fct_rev(period_abbr), fill = fct_rev(cluster))) +
      geom_bar(position = "fill", color = "grey60") +
      scale_fill_manual(values = var.colors.boxplot) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(x = "Period", y = "Percentage of total", fill = "Cluster") +
      facet_grid(. ~ biome3, scales = "free_x") +
      envalysis::theme_publish(base_size = 18) +
      theme(
        legend.position = "bottom",
        axis.title.y = element_text(size = 12),
        panel.spacing.x = unit(0.5, "lines")
      ) +
      coord_flip()
    
    cowplot::save_plot(
      here(
        "output",
        "plots",
        "disturbance_analysis",
        
        paste0(
          "percentage_grids_per_biome_period1to3_",
          var.agent,
          '_',
          grid.scale,
          ".svg"
        )
      ),
      biome.period.reclass,
      base_width = 12,
      base_height = 6,
      dpi = 400
    )
    
    saveRDS(biome.period.reclass,
            here(
              'output',
              'plots',
              'disturbance_analysis',
              paste0(
                "percentage_grids_per_biome_period1to3_",
                var.agent,
                '_',
                grid.scale,
                ".rds"
              )
            ))
    
    biome.period3.reclass <-  ggplot(
      annual.regimes.class %>%
        filter(period == '2011-2023') ,
      aes(x = fct_rev(biome3), fill = fct_rev(cluster))
    ) +
      geom_bar(
        stat = "count",
        alpha = 0.95,
        color = 'grey60',
        position = "fill"
      ) +
      envalysis::theme_publish() +
      scale_fill_manual(values = c(var.colors.boxplot)) +
      labs(y = "Percentage of grids") +
      coord_flip() +  # ← Gira el gráfico
      envalysis::theme_publish(base_size = 18) +
      scale_y_continuous(
        labels = function(x)
          x * 100,
        breaks = seq(0, 1, by = 0.25)
      ) +
      theme(
        legend.position = 'none',
        axis.title.y = element_blank(),
        axis.text.y = element_text(
          angle = 90,
          hjust = 0.5,
          vjust = 0.5
        )
      )
    saveRDS(
      biome.period3.reclass,
      here(
        'output',
        'plots',
        'disturbance_analysis',
        paste0(
          "percentage_grids_per_biome_period3_",
          var.agent,
          '_',
          grid.scale,
          ".rds"
        )
      )
    )
    

    #### 13. Combined figure: Boxplots by period and PCA of attributes (example for an agent) ####
    
    period.boxplots <- readRDS(here(
      "output",
      "plots",
      "disturbance_analysis",
      
      paste0("boxplots_", var.agent, "_", list.scales[1], ".rds")
    ))
    pca.classified <- readRDS(here(
      "output",
      "plots",
      "disturbance_analysis",
      
      paste0("pca_", var.agent, "_", list.scales[1], ".rds")
    ))
    
    boxplot.pca.plot <- cowplot::plot_grid(
      pca.classified + theme(plot.margin = margin(5, 5, 5, 5)),
      period.boxplots + theme(plot.margin = margin(5, 5, 5, 5)),
      rel_widths = c(1, 1.5),
      nrow = 1,
      align = "hv",
      axis = "tblr"
    )
    
    cowplot::save_plot(
      here(
        "output",
        "plots",
        "disturbance_analysis",
        
        paste0("pca_plus_boxplot_", var.agent, ".svg")
      ),
      boxplot.pca.plot,
      base_width = 10,
      base_height = 5,
      dpi = 600
    )
    
    
    #### 14. Cluster map of the last period with percentage of biomes  ####
    
    period3.classified.map <- readRDS(here(
      "output",
      "plots",
      "disturbance_analysis",
      paste0(
        "disturbance_map_3rdperiod_",
        var.agent,
        "_",
        list.scales[1],
        ".rds"
      )
    ))
    biome.reclass.hist <- readRDS(here(
      "output",
      "plots",
      "disturbance_analysis",
      
      paste0(
        "percentage_grids_per_biome_period3_",
        var.agent,
        "_",
        list.scales[1],
        ".rds"
      )
    ))
    
    biome.reclass.hist <- biome.reclass.hist +
      envalysis::theme_publish(base_size = 18) +
      theme(
        legend.position = "none",
        axis.title.y    = element_blank(),
        axis.text.y     = element_text(
          angle = 45,
          hjust = 1,
          vjust = 0.5
        )
      )
    
    classified.biome.plot <- cowplot::plot_grid(
      period3.classified.map + theme(plot.margin = margin(0, 0, 0, 0)),
      biome.reclass.hist + theme(plot.margin = margin(0, 0, 0, 0)),
      nrow = 2,
      rel_heights = c(2, 0.5)
    )
    
    cowplot::save_plot(
      here(
        "output",
        "plots",
        "disturbance_analysis",
        
        paste0("cluster_map_with_biome_breakdown_", var.agent, ".svg")
      ),
      classified.biome.plot,
      base_width = 10,
      base_height = 12,
      dpi = 600
    )
    
    
    
    #### 15: Percentage change between periods (1985–1997 vs. 2011–2023) ####
    
    changes.by.period <- period.summary %>%
      filter(period %in% c('1985-1997', '2011-2023')) %>%
      arrange(grid_id, period) %>%
      group_by(grid_id) %>%
      mutate(previous_reclass = dplyr::lag(reclass)) %>%
      filter(#!is.na(previous_reclass) &
        previous_reclass != reclass) %>%
      ungroup() %>%
      count(previous_reclass, reclass) %>%
      mutate(percentage = n / sum(n) * 100)
    
    gt(changes.by.period) %>%
      gtsave(here(
        'output',
        'plots',
        'disturbance_analysis',
        paste0("period_1_to3_transitions_", var.agent, ".html")
      ))
    
    chord_data <- changes.by.period[, c("previous_reclass", "reclass", "n")]
    
    svg(
      here(
        'output',
        'plots',
        'disturbance_analysis',
        
        paste0(var.agent, "_chord_diagram_period1_to_3.svg")
      ),
      width = 10,
      height = 10
    )
    
    chordDiagram(
      chord_data,
      grid.col = var.colors2,
      annotationTrack = "grid",
      preAllocateTracks = 1,
      directional = 1,
      direction.type = "arrows",
      link.arr.type = "triangle",
      link.arr.length = 0.1,
      
      link.border = "grey20",
      link.lwd = 0.6,
      link.arr.lwd = 0.2,
      link.arr.col = 'grey20',
    )
    
    circos.trackPlotRegion(
      track.index = 2,
      panel.fun = function(x, y) {
        xlim <- get.cell.meta.data("xlim")
        ylim <- get.cell.meta.data("ylim")
        sector.name <- get.cell.meta.data("sector.index")
        circos.rect(
          xlim[1],
          ylim[1],
          xlim[2],
          ylim[2],
          border = "grey20",
          lwd = 0.8,
          col = NA
        )
        circos.text(
          mean(xlim),
          ylim[1] + 7,
          sector.name,
          facing = "inside",
          niceFacing = TRUE,
          cex = 1.5
        )
        circos.axis(
          h = "top",
          labels.cex = 1.2,
          major.tick = 1e-10,
          minor.tick = 1e-10,
          sector.index = sector.name,
          track.index = 2
        )
      },
      bg.border = NA
    )
    
    dev.off()

  }
  
}


####  Combined graphs and  figures ####


#### 16. Map of the study region with biomes and highlighted areas ####

biome.colors <- c(
  "Boreal Forests/Taiga"                  = "#532E57FF",
  "Mediterranean Forests, Woodlands & Scrub" = "#C1B178FF",
  "Temperate Broadleaf & Mixed Forests"   = "#296656FF",
  "Temperate Conifer Forests"             = "#1B4F72FF",
  "Temperate Grasslands, Savannas & Shrublands" = "#f9f9f9e5",
  "Tundra"                                = "#f9f9f9e5"
)

grids_filtered <- europe %>%
  left_join(readRDS(here(
    "data",
    "reference",
    "patch_summary",
    paste0("selector_forest", grid.scale, ".rds")
  )), by = join_by('grid_id' == 'grid_id2')) %>%
  filter(selector == FALSE)

study.area.plot <- ggplot() +
  
  geom_sf(
    data = biome.boundaries ,
    aes(fill = BIOME_NAME),
    color = NA,
    alpha = 0.8,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = biome.colors) +
  geom_sf(
    data = pattern,
    aes(color = reclass),
    fill = NA,
    size = 0.05,
    alpha = 1,
    show.legend = FALSE
  ) +
  scale_color_manual(values = biome.colors) +
  geom_sf(
    data = grids_filtered,
    alpha = 1,
    fill = "#f9f9f9e5",
    color = "#f9f9f9e5",
    show.legend = FALSE
  ) +
  geom_sf(
    data = world.boundaries %>%
      filter(
        !name %in% c("Russian Federation", "Iceland", "Svalbard and Jan Mayen Islands")
      ),
    color = "grey30",
    fill = NA,
    size = 0.2,
    alpha = 1,
    show.legend = FALSE
  ) +
  coord_sf(crs = st_crs(3035)) +
  annotation_scale(
    location = "bl",
    width_hint = 0.1,
    style = 'ticks',
    text_family = "Arial",
    text_col = "grey50"
  ) +
  annotation_north_arrow(
    location = "br",
    which_north = "true",
    height = unit(1, "cm"),
    width = unit(1, "cm"),
    style = north_arrow_fancy_orienteering(text_family = "Arial", text_col = "grey50")
  ) +
  labs(fill = "") +
  theme_minimal() +
  theme(
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0.5,
      l = 0
    ),
    legend.position = c(0.99, 0.99),
    legend.justification = c("right", "top"),
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    text = element_text(color = "grey50", size = 14),
    axis.title = element_blank(),
    axis.text = element_text(color = "grey50"),
    panel.grid.major = element_line(color = "grey50", size = 0.2),
    panel.grid.minor = element_line(color = "grey50", size = 0.1)
  ) +
  envalysis::theme_publish(base_linewidth = 0.2)

cowplot::save_plot(
  here(
    "output",
    "plots",
    "disturbance_analysis",
    
    paste0("study_region.svg")
  ),
  study.area.plot,
  base_width = 10,
  base_height = 6,
  dpi = 600
)


## 16. Comparison of trends by biome between agents (Fire vs Wind) ##

pinpoint.fire <- readRDS(
  here(
    "output",
    "plots",
    "disturbance_analysis",
    "biome3_trend_pointplot_fire.rds"
  )
)

fire.hist <- readRDS(here(
  "output",
  "plots",
  "disturbance_analysis",
  paste0("trend_hist_fire__50km.rds")
)) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

pinpoint.wind <- readRDS(
  here(
    "output",
    "plots",
    "disturbance_analysis",
    "biome3_trend_pointplot_wind_barkbeetle.rds"
  )
)
wind.hist <- readRDS(here(
  "output",
  "plots",
  "disturbance_analysis",
  paste0("trend_hist_wind_barkbeetle__50km.rds")
)) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )


classified.biome.pinpoint.plot.fire <- cowplot::plot_grid(
  pinpoint.fire + theme(plot.margin = margin(0, 0, 0, 0)),
  fire.hist + theme(plot.margin = margin(0, 0, 0, 0), axis.text.y = element_blank()),
  nrow = 1,
  align = "hv",
  axis = "tblr"
)

cowplot::save_plot(
  here(
    "output",
    "plots",
    "disturbance_analysis",
    "combined_pinpoints_biome_fire.svg"
  ),
  classified.biome.pinpoint.plot.fire,
  base_width = 7,
  base_height = 4,
  dpi = 600
)

classified.biome.pinpoint.plot.wind <- cowplot::plot_grid(
  pinpoint.wind + theme(plot.margin = margin(0, 0, 0, 0)),
  wind.hist + theme(plot.margin = margin(0, 0, 0, 0), axis.text.y = element_blank()),
  nrow = 1,
  align = "hv",
  axis = "tblr"
)

cowplot::save_plot(
  here(
    "output",
    "plots",
    "disturbance_analysis",
    "combined_pinpoints_biome_wind.svg"
  ),
  classified.biome.pinpoint.plot.wind,
  base_width = 7,
  base_height = 4,
  dpi = 600
)


####  17: Combined comparison of agents for the most recent period (2011–2023) ####

grid_summary_windbarkbeetle <- read.csv(
  here(
    'output',
    'plots',
    'disturbance_analysis',
    "table1_statistics_wind_barkbeetle.csv"
  )
) %>%
  as_tibble() %>%
  filter(period == '2011-2023') %>%
  mutate(agent = "Wind & bark beetle")

grid_summary_fire <- read.csv(here(
  'output',
  'plots',
  'disturbance_analysis',
  "table1_statistics_fire.csv"
)) %>%
  as_tibble() %>%
  filter(period == '2011-2023') %>%
  mutate(agent = "Fire")

new_table <- bind_rows(grid_summary_windbarkbeetle, grid_summary_fire) %>%
  mutate(
    agent = factor(agent),
    period = factor(period),
    cluster = factor(reclass)
  ) %>%
  select(-reclass) %>%
  relocate(c(agent, cluster), .after = period)

gt_stats_new <- new_table %>%
  gt() %>%
  cols_label(
    agent = "Agent",
    cluster = "Cluster",
    size_stats = "Size (ha)",
    sev_stats = "Severity (0-1)",
    frequency_stats = "Frequency (events year-1)",
    sumfrequency_stats = "Sum Frequency (events period-1)",
    mean_sum_area = 'Total area (ha)',
    percent_grids = "Grids (%)"
  ) %>%
  cols_hide(columns = period) %>%
  tab_style(style = cell_text(
    font = c(google_font(name = "Arial"), default_fonts()),
    size = px(11)
  ),
  locations = cells_body(columns = everything())) %>%
  tab_style(
    style = cell_text(
      font = c(google_font(name = "Arial"), default_fonts()),
      weight = "bold",
      size = px(12)
    ),
    locations = cells_column_labels(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(
      font = c(google_font(name = "Arial"), default_fonts()),
      weight = "bold",
      size = px(12)
    ),
    locations = cells_row_groups(groups = everything())
  ) %>%
  cols_align(align = "center", columns = everything()) %>%
  tab_spanner(
    label = "Disturbance",
    columns = c(
      "size_stats",
      "sev_stats",
      "frequency_stats",
      "sumfrequency_stats",
      "mean_sum_area"
    )
  ) %>%
  tab_spanner(label = "Grids", columns = "percent_grids") %>%
  tab_style(style = cell_text(
    font = c(google_font(name = "Arial"), default_fonts()),
    size = px(12)
  ), locations = cells_column_spanners()) %>%
  fmt_number(
    columns = c(
      "size_stats",
      "sev_stats",
      "frequency_stats",
      "mean_sum_area",
      "percent_grids"
    ),
    decimals = 2
  )

gt_stats_new %>%
  gtsave(here(
    'output',
    'plots',
    'disturbance_analysis',
    "table1_statistics_period3.png"
  ))

gt_stats_new %>%
  gtsave(here(
    'output',
    'plots',
    'disturbance_analysis',
    "table1_statistics_period3.html"
  ))

write.csv(
  new_table,
  file = here(
    'output',
    'plots',
    'disturbance_analysis',
    "table1_statistics_period3.csv"
  ),
  row.names = FALSE
)
