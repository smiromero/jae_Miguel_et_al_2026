# Functions:


summarise_annual_regime <- function(data, vars, var_forest, var_grid) {
  data %>%
    group_by(across(all_of(vars)), year) %>%
    reframe(
      sum_area = sum(area_ha , na.rm = TRUE),
      q95_size = quantile(area_ha, 0.95, na.rm = TRUE),
      mean_size = mean(area_ha , na.rm = TRUE),
      mean_sev = mean(severity_mean , na.rm = TRUE),
      forest_ha = unique(na.omit(!!sym(var_forest))),
      grid_ha = unique(na.omit(!!sym(var_grid))),
      forest_cover = (forest_ha * 100) / grid_ha,
      n_events = n(),
      n_events = if_else(sum_area == 0, 0, n_events),
      frequency = (n_events / forest_cover) ,
      rate = (sum_area * 100) / forest_ha,

    ) %>%
    arrange(across(all_of(vars)), year) %>%
    ungroup()
  
}


summarise_period <- function(regime.data, vars, decades.list) {
  map(decades.list, \(decade) {
    decade_data <- regime.data %>%
      dplyr::filter(between(year, decade[1], decade[2])) %>%
      period_stats_calculation(., vars = vars)
  })
  
}


period_stats_calculation <- function(data, vars) {
  data %>%
    mutate(across(everything(), ~ if_else(is.nan(.), NA, .))) %>%
    group_by(across(all_of(vars))) %>%
    #filter(if_all(contains('mean'), ~ !is.na(.))) %>%
    reframe(
      size = mean(mean_size , na.rm = TRUE),
      sdsize = sd(mean_size, na.rm = T),
      vcsize = sdsize / size,
      q95_size = mean(q95_size, na.rm = TRUE),
      
      sev = mean(mean_sev, na.rm = TRUE),
      sdsev = sd(mean_sev, na.rm = T),
      vcsev =  sdsev / sev,
      
      ratex = mean(rate, na.rm = TRUE),
      sdrate = sd(rate, na.rm = TRUE),
      vcrate = sdrate / ratex,
      sum_area = sum(sum_area, na.rm = TRUE),
      sumrate = (sum_area * 100) / unique(forest_ha),
      
      frequencyx = mean(frequency, na.rm = TRUE),
      sdfrequency = sd(frequency, na.rm = TRUE),
      vcfrequency = sdfrequency / frequencyx,
      
      
      sum_events = sum(n_events , na.rm = TRUE),
      sumfrequency = sum_events / unique(forest_cover),
      
      
      n_years = sum(if_else(n_events > 0, 1, 0))
    ) %>% # suma de eventos por ha
    mutate(across(-all_of(c(
      vars, 'frequencyx'
    )), ~ round(.x, digits = 3))) %>%
    rename('rate' = 'ratex', 'frequency' = 'frequencyx')  %>%
    group_by(agent) %>%
    mutate(across(
      contains('rate') | contains('frequency'),
      ~ if_else(. == 0, NA, .)
    ))
  #   across(c(size, sdsize, vcsize, sev, sdsev, vcsev, rate, sdrate, vcrate, sumrate,
  #            frequency, sdfrequency, vcfrequency, sumfrequency),
  #          ~log(. ),
  #          .names = "log_{.col}")) %>%
  # filter(across(everything(), is.finite)) %>%
  # group_by(agent) %>%
  # mutate(
  #   across(contains('log') | contains('z'),  ~ scales::rescale(., to = c(0, 1)))) %>%
  # select(grid_id2, agent, contains('log'))
  
}


create_ggbiplot <- function(choices,
                            plot_title,
                            agent,
                            color_palette,
                            reference.period,
                            grid_scale) {
  rownames(pca.mod$rotation) <- ifelse(
    rownames(pca.mod$rotation) == "log_sumfrequency",
    "Frequency",
    ifelse(
      rownames(pca.mod$rotation) == "log_size",
      "Size",
      ifelse(
        rownames(pca.mod$rotation) == "log_sev",
        "Severity",
        rownames(pca.mod$rotation)
      )
    )
  )
  rownames(loadings.pca) <- ifelse(
    rownames(loadings.pca) == "log_sumfrequency",
    "Frequency",
    ifelse(
      rownames(loadings.pca) == "log_size",
      "Size",
      ifelse(
        rownames(loadings.pca) == "log_sev",
        "Severity",
        rownames(loadings.pca)
      )
    )
  )
  
  
  decades.i <- purrr::pluck(period.data, paste0('period', reference.period))
  decades.mclust.class.i <- purrr::pluck(mclust.class, paste0('period', reference.period))
  
  #
  if (agent == 'fire' | agent == 'harvest') {
    color_palette <- var.colors[c(3, 2, 4)]
    labels_custom <-  c("A", "B", "C")
    
    if (grid_scale == '50km') {
      levels_custom <- c(2, 1, 3)
      
    } else if (grid_scale == '25km') {
      levels_custom <- c(3, 1, 2)
    } else if (grid_scale == '100km') {
      levels_custom <- c(1, 3, 2)
    }
    
    decades.mclust.class.i <- factor(decades.mclust.class.i,
                                     levels = levels_custom,
                                     labels = c("A", "B", "C"))
    arrow_size <-  4
    
  } else if (agent == 'wind_barkbeetle') {
    color_palette <- var.colors[c(2, 3, 4, 5)]
    labels_custom <-  c("A", "B", "C", "D")
    
    if (grid_scale == '50km') {
      levels_custom <- c(1, 2, 3, 4)
    } else if (grid_scale == '25km') {
      levels_custom <- c(3, 4, 5, 1)
    } else if (grid_scale == '100km') {
      levels_custom <- c(1, 2, 3, 4)
    }
    
    decades.mclust.class.i <- factor(
      decades.mclust.class.i,
      levels = levels_custom,
      labels = c("A", "B", "C", "D")
    )
    arrow_size <-  4
  }
  
  pca_plot <- ggbiplot(
    pca.mod,
    choices = choices,
    labels = NULL,
    groups = factor(decades.mclust.class.i),
    varname.size = 0,
    ellipse = TRUE,
    ellipse.linewidth = 1,
    var.axes = TRUE,
    var.factor = 3,
    var.scale  = 1,
    arrow.size = arrow_size,
    alpha = 0.3,
    alpha_arrow = 0.7,
    repel = FALSE
  ) +
    scale_color_manual(values = color_palette, guide = "none")  
  pca_plot +
    geom_point(
      data = mclust.stats$mean,
      aes_string(
        x = paste0("PC", choices[1]),
        y = paste0("PC", choices[2]),
        fill = factor(
          mclust.stats$mean$group,
          levels = levels_custom,
          labels = labels_custom
        )
      ),
      size = 4,
      shape = 21,
      stroke = 1,
      color = "grey20"  
    ) +
    ggrepel::geom_text_repel(
      data = loadings.pca,
      aes(
        x = !!sym(paste0("PC", choices[1])) * 5,
        y = !!sym(paste0("PC", choices[2])) * 5,
        label = row.names(loadings.pca)
      ),
      size = 4,
      color = "grey20",
      bg.color = 'grey90',
      max.overlaps = 5
    ) +
    
    scale_fill_manual(values = color_palette, name = "Regime") +
    envalysis::theme_publish() + 
    labs(x = "PCA 1", y = "PCA 2") +
    theme(legend.position = "none",
          plot.title = element_text(size = 14, face = "bold"))   
}

agent_data_decades <- function(regime.data, decades.list, var.agent) {
  map(decades.list, \(decade) {
    decade_data <- regime.data %>%
      dplyr::filter(between(year, decade[1], decade[2]) &
                      agent == var.agent) %>%
      summarise_regime(.,
                       vars = c('grid_id'),
                       vars_landscape = c('grid_id'))
    
    decade_summarised <- calculate_log(decade_data, var.agent) %>%
      dplyr::select(grid_id, contains('log')) %>%
      mutate(period = paste(decade[1], decade[2], sep = '-'))
    
    decade_summarised
  })
}


plot_attribute_distribution <- function(agent.data.decades, write_out = TRUE) {
  purrr::imap(agent.data.decades, \(agent.decade, period) {
    agent.decade.plot <- agent.decade %>%
      select(-grid_id, -agent) %>%
      map2(., names(.), \(var, varname) {
        ggplot() +
          geom_histogram(aes(var)) +
          xlab(varname)
      }) %>%
      cowplot::plot_grid(plotlist = .)
    
    if (write_out) {
      ggsave(
        here::here(
          "output",
          "plots",
          "disturbance_analysis",
          paste0("distribution_disturbance_", period, ".jpg")
        ),
        agent.decade.plot,
        device = "png",
        height = 4,
        dpi = 300
      )
      message("El gráfico ha sido guardado")
    } else {
      message("El gráfico no ha sido guardado")
    }
    agent.decade.plot
  })
}

plot_correlation <- function(data, logo_path = (here(paste0(
  "data/reference/icons/", var.agent, ".png"
)))) {
  corr_matrix <- data %>%
    bind_rows() %>%
    dplyr::select(-c(grid_id, agent)) %>%
    cor(use = "complete.obs")
  
  
  corr_plot <- ggcorrplot::ggcorrplot(
    corr_matrix,
    type = "upper",
    hc.order = TRUE,
    lab = TRUE,
    lab_size = 6,
    colors = c("blue", "white", "red"),
    tl.cex = 50 / .pt,
    tl.col = "black",
    tl.srt = 45
  
  logo <- cowplot::ggdraw() +
    cowplot::draw_image(
      logo_path,
      x = 0.85,
      y = 0.5,
      width = 0.15,
      height = 0.15,
      scale = 0.6
    )
  
  final_plot <- cowplot::ggdraw(corr_plot) +
    cowplot::draw_plot(logo, 0, 0, 1, 1)
  
  
  ggsave(
    here(
      'output',
      'plots',
      'disturbance_analysis',
      paste0("corrplots_", var.agent , "_alldisturbances_1985_2023.png")
    ),
    plot = final_plot,
    device = "tiff",
    dpi = 300
  )
  
  message("El gráfico ha sido guardado")
  print(final_plot)
}

plot_correlation_article <- function(agent_name,
                                     var.agent,
                                     grid.scale,
                                     n.periods) {
  corr_data <- readRDS(here(
    "data",
    "reference",
    "patch_summary",
    paste0("regime_", n.periods, "_period", grid.scale, "_TUM.rds")
  )) %>%
    pull(period.regime.filter) %>%
    .[[1]] %>%
    map( ~ .x %>%
           rename(grid_id = grid_id2) %>%
           filter(across(everything(), is.finite)) %>%
           filter(agent == agent_name)) %>%
    bind_rows() %>%
    dplyr::select(where(is.numeric)) %>%
    dplyr::select(-contains("log")) %>%
    dplyr::select(where( ~ sd(., na.rm = TRUE) != 0)) %>%
    cor(use = "complete.obs")
  
  corr_plot <- ggcorrplot::ggcorrplot(
    corr_data,
    hc.order = TRUE,
    type = "lower",
    lab = TRUE,
    lab_size = 5,
    method = "square",
    colors = c("#6D9EC1", "white", "#E46767"),
    title = paste("Matriz de correlación -", agent_name),
    ggtheme = ggplot2::theme_minimal()
  ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 16,
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 15,
        face = "bold"
      ),
      axis.text.y = element_text(size = 15, face = "bold"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 15)
    )
  cowplot::save_plot(
    here(
      "output",
      "plots",
      "disturbance_analysis",
      "TUM",
      paste0(
        "correlation_plot_",
        var.agent,
        "_",
        grid.scale,
        "_",
        agent_name,
        "_TUM.svg"
      )
    ),
    corr_plot,
    base_width = 12,
    base_height = 12,
    dpi = 600
  )
}


plot_pca <- function(pca, write_out = TRUE) {
  
  eigenvalues <- get_eigenvalue(pca)
  print(eigenvalues)
  
  common_theme <- theme(
    text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18)
  )
  
  eigenvalues.plot <- fviz_eig(pca, main = '') +
    common_theme +
    geom_bar(stat = "identity",
             fill = "grey30",
             color = 'grey30')
  
  pca.var.plot <- fviz_pca_var(pca, col.var = 'contrib', repel = TRUE) +
    scale_color_gradient(low = "#ffba08", high = "#6a040f") +
    common_theme
  
  pca.var.contrib.plot <- fviz_contrib(pca, choice = 'var') +
    geom_bar(stat = "identity",
             fill = "grey30",
             color = 'grey30') +
    common_theme +
    ggtitle("")
  combined_plot <- plot_grid(
    eigenvalues.plot,
    pca.var.contrib.plot,
    labels = c("A", "B"),
    ncol = 1,
    rel_heights = c(1, 1.25),
    label_size = 20
  )
  combined_plot2 <-  plot_grid(pca.var.plot, labels = "C", label_size = 20)
  combined_plot3 <- plot_grid(combined_plot,
                              combined_plot2,
                              ncol = 2,
                              rel_heights = c(2, 1)) 
  plot(combined_plot3)
  
  if (write_out == TRUE) {
    ggsave(
      here(
        'output',
        'plots',
        'disturbance_analysis',
        paste0("pca_", var.agent, ".png")
      ),
      plot = combined_plot3,
      device = "png",
      height = 4,
      dpi = 300
    )
    
    message("El gráfico ha sido guardado")
  } else {
    message("El gráfico no ha sido guardado")
  }
  
    print(combined_plot3)
}



cluster_maps <- function(agent_shp, decades_list, show_legend = FALSE) {
  tm_shape(st_as_sf(agent_shp)) + 
    tm_fill(col = "gray90", alpha = 1) + 
    tm_shape(st_as_sf(world.boundaries)) + 
    tm_fill(col = "gray80", alpha = 1) +  
    tm_shape(st_as_sf(europe.boundaries)) +  
    tm_fill(col = "gray80", alpha = 1) + 
    tm_shape(st_as_sf(agent_shp)) + 
    tm_borders() +
    tm_fill(
      col = "class",
      palette = var.colors[-6],
      alpha = 1,
      legend.show = show_legend
    ) + 
    tm_grid(n.x = 4, n.y = 3, alpha = 0.4) +  
    tm_layout(
      main.title = paste(decades_list[1], decades_list[2], sep = '-'),
      legend.text.size = 4 
    )
}

create_cluster_map <- function(var_data, var_period, var_legend, var = 'reclass') {
  if (var.agent == 'fire') {
    color_scale <- c(
      "Western Mediterranean" = "#CD6463",
      'Large & frequent' = "#CD6463",
      "Central Eastern Mediterranean" = "#F8D67B",
      'Moderate' = "#F8D67B",
      "Temperate Boreal" = "#88B680",
      'Severe & rare' = "#88B680",
      "Residual" = "#f9f9f9e5"
    )
  } else {
    color_scale <- c(
      "Central European" = "#CD6463" ,
      'Frequent & small' = "#CD6463" ,
      
      "Temperate Boreal" = "#F8D67B",
      'Moderate' = "#F8D67B",
      
      "Atlantic" = '#9c89b8',
      
      'Severe & large' = '#9c89b8',
      "Mediterranean" = "#88B680" ,
      'Mild & Rare' = "#88B680" ,
      
      "Residual" = "#f9f9f9e5"
    )
    
  }
  
  ggplot() +
    
    geom_sf(
      data = st_as_sf(world.boundaries %>% filter(
        !name %in% c(
          "Russian Federation",
          'Iceland',
          'Svalbard and Jan Mayen Islands'
        )
      )),
      color = 'grey60',
      fill = '#DDDDDD',
      size = 0.2,
      show.legend = FALSE
    ) +
    
   
    geom_sf(
      data = st_as_sf(var_data),
      aes_string(fill = var),
      color = NA,
     
      size = 0.05,
      alpha = 0.9,
      show.legend = FALSE
    ) +
    geom_sf(
      data = st_as_sf(grids.noforest),
      color = 'grey85',
      fill = 'white',
      size = 0.05,
      alpha = 0.4,
      show.legend = FALSE
    ) +
    geom_sf(
      data = st_as_sf(europe.boundaries),
      color = 'grey60',
      fill = NA,
      size = 0.05,
      alpha = 0.1,
      show.legend = FALSE
    ) +
    labs(fill = "") +
    scale_fill_manual(values = color_scale) +
    
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
      style = north_arrow_fancy_orienteering(text_family = "Arial", text_col = "grey50")  # Color gris para la brújula
    ) +
    
    ggtitle(paste0(var_period[1], '- ', var_period[2])) +
    theme(
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 0.5,
        l = 0
      ),
      legend.position = "none",
      legend.key.size = unit(0.5, "cm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      text = element_text(color = "grey50", size = 14),
      axis.title = element_blank(),
      axis.text = element_text(color = "grey50"),
      panel.grid.major = element_line(color = "grey50", size = 0.2),
      panel.grid.minor = element_line(color = "grey50", size = 0.1)  
    ) +
    envalysis::theme_publish(base_linewidth = 0.2) +
    coord_sf(crs = st_crs(3035))
  
}



create_cluster_map_basic <- function(var_data, var_period, var_legend, var = 'reclass') {
  if (var.agent == 'fire') {
    color_scale <- c(
      "Western Mediterranean" = "#CD6463",
      'Large & frequent' = "#CD6463",
      "Central Eastern Mediterranean" = "#F8D67B",
      'Moderate' = "#F8D67B",
      "Temperate Boreal" = "#88B680",
      'Severe & rare' = "#88B680",
      "Residual" = "#f9f9f9e5"
    )
  } else {
    color_scale <- c(
      "Central European" = "#CD6463" ,
      'Frequent & small' = "#CD6463" ,
      
      "Temperate Boreal" = "#F8D67B",
      'Moderate' = "#F8D67B",
      
      "Atlantic" = '#9c89b8',
      
      'Severe & large' = '#9c89b8',
      "Mediterranean" = "#88B680" ,
      'Mild & Rare' = "#88B680" ,
      
      "Residual" = "#f9f9f9e5"
    )
    
  }

  ggplot() +
    geom_sf(
      data = st_as_sf(europe.boundaries),
      color = NA,
      fill = 'grey30',
      size = 0.05,
      alpha = 0.2,
      show.legend = FALSE
    ) +
    
    geom_sf(
      data = st_as_sf(grids.noforest),
      color = 'grey85',
      fill = NA,
      size = 0.05,
      alpha = 0.4,
      show.legend = FALSE
    ) +
    
    geom_sf(
      data = st_as_sf(var_data),
      aes_string(fill = var),
      color = NA,
      linewidth = 0,
      
      alpha = 0.9,
      show.legend = FALSE
    ) +
    
    geom_sf(
      data = st_as_sf(europe.boundaries),
      color = 'grey20',
      fill = NA,
      size = 0.05,
      alpha = 0.7,
      show.legend = FALSE
    ) +
    
    labs(fill = "") +
    scale_fill_manual(values = color_scale) +
    
    theme(
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      ),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA),
      legend.position = "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      text = element_text(color = "grey50", size = 14)  
      
    ) +
    coord_sf(crs = st_crs(3035))
  
}




agent_icon <- function(agent, icons_path = "data/reference/icons/") {
  logo_agent <- png::readPNG(here::here(paste0(icons_path, agent, ".png")), native = TRUE)
  
  logo_ggdraw <- cowplot::ggdraw() +
    cowplot::draw_image(
      logo_agent,
      x = 0.85,
      y = 0.85,
      width = 0.15,
      height = 0.15,
      scale = 0.6
    )
  
  return(logo_ggdraw)
}



generate_mclusters  <- function(mclust,
                                period.pca,
                                reference.period,
                                predict.period,
                                n.periods) {
  mclust.class <- vector("list", n.periods)
  mclust.prob <- vector("list", n.periods)
  
  mclust.class[[1]] <- mclust$classification
  mclust.prob[[1]] <- as.data.frame(mclust$z) %>%
    mutate(max.prob = pmap_dbl(., max)) %>%
    select(max.prob)
  
  if (n.periods > 1) {
    for (i in 2:n.periods) {
      pred <- predict(mclust, purrr::pluck(period.pca, paste0('period', predict.period[i - 1])))
      mclust.class[[i]] <- pred$classification
      mclust.prob[[i]] <- as.data.frame(pred$z) %>%
        mutate(max.prob = pmap_dbl(., max)) %>%
        select(max.prob)
    }
  }
  
  period.names <- c(paste0('period', reference.period),
                    paste0('period', predict.period[1:(n.periods - 1)]))
  names(mclust.class) <- period.names
  names(mclust.prob) <- period.names
  
  return(list(classification = mclust.class, probability = mclust.prob))
}


remove_outliers <- function(data, quantile_threshold = 0.99) {
  mahalanobis_distance <- mahalanobis(data, colMeans(data), cov(data))
  
  threshold <- quantile(mahalanobis_distance, quantile_threshold)
  
  outliers <- mahalanobis_distance > threshold
  
  filtered_data <- data[!outliers, ]
  
  return(
    list(
      filtered_data = filtered_data,
      outliers = outliers,
      mahalanobis_distance = mahalanobis_distance
    )
  )
}

create_pca_periods_plot <- function(data, x_var, y_var, grid_ids, loadings_pca) {
  color_values <- c("#45648C", "#76A660", "#D95A4E")
  
  ggplot(data %>% filter(grid_id %in% grid_ids),
         aes(
           x = .data[[x_var]],
           y = .data[[y_var]],
           color = factor(period)
         )) +
    geom_point(aes(shape = factor(period)),
               show.legend = FALSE,
               alpha = 0.4) +
    stat_ellipse(
      geom = "polygon",
      level = 0.95,
      alpha = 0.3,
      size = 1,
      aes(fill = factor(period), color = factor(period))
    ) +
    geom_segment(
      data = loadings_pca,
      aes(
        x = 0,
        y = 0,
        xend = .data[[x_var]] * 5,
        yend = .data[[y_var]] * 5
      ),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "black"
    ) +
    ggrepel::geom_text_repel(
      data = loadings_pca,
      aes(
        x = .data[[x_var]] * 5,
        y = .data[[y_var]] * 5,
        label = row.names(loadings_pca)
      ),
      size = 6,
      color = "grey20",
      bg.color = 'grey90',
      max.overlaps = 8
    ) +
    labs(
      x = x_var,
      y = y_var,
      color = "period ",
      fill = "period ",
      size = "Group"
    ) +
    envalysis::theme_publish() +
    scale_color_manual(values = color_values) +
    scale_fill_manual(values = color_values) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 15))
}

process_period_data <- function(n.periods,
                                period.data,
                                period.pca,
                                mclust.class,
                                mclust.prob,
                                decades.list,
                                var.agent,
                                grid.scale) {
  results_list <- list()
  
  classify_data <- function(period_num,
                            period_data,
                            period_pca,
                            mclust_class,
                            mclust_prob,
                            decade_range,
                            var_agent,
                            grid_scale) {
    data.frame(
      period_data,
      period_pca,
      period = paste(decade_range[1], decade_range[2], sep = '-'),
      class = as.factor(mclust_class),
      probability = mclust_prob$max.prob
    ) %>%
      mutate(
        grid_scale = grid_scale,
        var_agent = var.agent,
        class = case_when(
          grid_scale == "_25km" &  var_agent == 'fire' & class == 1 ~ "A",
          grid_scale == "_25km" &
            var_agent == 'fire' & class == 2 ~ "B",
          grid_scale == "_25km" &
            var_agent == 'fire' & class == 3 ~ "C",
          
          grid_scale == "_50km" &
            var_agent == 'fire' & class == 1 ~ "A",
          grid_scale == "_50km" &
            var_agent == 'fire' & class == 2 ~ "B",
          grid_scale == "_50km" &
            var_agent == 'fire' & class == 3 ~ "C",
          
          grid_scale == "_100km" &
            var_agent == 'fire' & class == 1 ~ "A",
          grid_scale == "_100km" &
            var_agent == 'fire' & class == 2 ~ "B",
          grid_scale == "_100km" &
            var_agent == 'fire' & class == 3 ~ "C",
          
          grid_scale == "_25km" &
            var_agent == 'wind_barkbeetle' & class == 1 ~ "A",
          grid_scale == "_25km" &
            var_agent == 'wind_barkbeetle' & class == 2 ~ "B",
          grid_scale == "_25km" &
            var_agent == 'wind_barkbeetle' & class == 3 ~ "C",
          grid_scale == "_25km" &
            var_agent == 'wind_barkbeetle' & class == 4 ~ "D",
          
          grid_scale == "_50km" &
            var_agent == 'wind_barkbeetle' & class == 1 ~ "A",
          grid_scale == "_50km" &
            var_agent == 'wind_barkbeetle' & class == 2 ~ "B",
          grid_scale == "_50km" &
            var_agent == 'wind_barkbeetle' & class == 3 ~ "C",
          grid_scale == "_50km" &
            var_agent == 'wind_barkbeetle' & class == 4 ~ "D",
          
          grid_scale == "_100km" &
            var_agent == 'wind_barkbeetle' & class == 1 ~ "A",
          grid_scale == "_100km" &
            var_agent == 'wind_barkbeetle' & class == 2 ~ "B",
          grid_scale == "_100km" &
            var_agent == 'wind_barkbeetle' & class == 3 ~ "C",
          grid_scale == "_100km" &
            var_agent == 'wind_barkbeetle' & class == 4 ~ "D",
          
          
          
        ),
        class = as.factor(class)
      ) %>%
      select(-var_agent, -grid_scale)
  }
  
  for (i in 1:n.periods) {
    results_list[[i]] <- classify_data(
      i,
      period.data[[paste0('period', i)]],
      period.pca[[paste0('period', i)]],
      mclust.class[[paste0('period', i)]],
      mclust.prob[[paste0('period', i)]],
      decades.list[[i]],
      var.agent,
      grid.scale
    )
  }
  
  period.data.classified <- do.call(rbind, results_list)
  
  return(period.data.classified)
}






generate_pca_plot_for_periods <- function(regime.data,
                                          period.data.classified,
                                          loadings.pca,
                                          var.agent,
                                          grid.scale) {
  get_grid_ids_intersect <- function(regime_data) {
    Reduce(base::intersect,
           list(regime_data[[1]]$grid_id, regime_data[[2]]$grid_id))
  }
  
  grid.id.intersect <- get_grid_ids_intersect(regime.data)
  
  params.pca.period <- list(
    data = list(period.data.classified),
    x_var = list("PC1"),
    y_var = list("PC2"),
    grid_ids = list(grid.id.intersect),
    loadings_pca = list(loadings.pca)
  )
  
  pca.periods.plot <- pmap(params.pca.period, create_pca_periods_plot)
  
  names(pca.periods.plot) <- c("PCA 1&2")
  
  pca.periods.plot.join <- plot_grid(pca.periods.plot$`PCA 1&2`)
  
  pca.periods.plot.join <- cowplot::ggdraw(pca.periods.plot.join) +
    cowplot::draw_plot(logo.ggdraw,
                       y = 0.1,
                       x = 0.1,
                       scale = 0.8)
  
  plot(pca.periods.plot.join)
  
  cowplot::save_plot(
    here(
      'output',
      'plots',
      'disturbance_analysis',
      paste0("PCA_space_by_periods_", var.agent, '_', grid.scale, ".svg")
    ),
    pca.periods.plot.join,
    base_width = 10,
    base_height = 5,
    base_asp = 1,
    dpi = 400
  )
  
  return(pca.periods.plot.join)
}


rename_periods <- function(period.data, n.periods) {
  period.names <- paste0("period", 1:n.periods)
  names(period.data) <- period.names
  return(period.data)
}

apply_pca_to_periods <- function(n_periods,
                                 pca_mod,
                                 period_subset,
                                 reference_period,
                                 predict_period,
                                 pca_select) {
  if (n_periods == 1) {
    period_pca <- list(pca_mod$x) %>%
      map( ~ select(as.data.frame(.x), pca_select))
    
    names(period_pca) <- paste0('period', reference_period)
    
  } else if (n_periods == 2) {
    period_pca <- list(pca_mod$x, predict(pca_mod, newdata = period_subset[[predict_period[1]]])) %>% map( ~ select(as.data.frame(.x), pca_select))
    
    names(period_pca) <- c(paste0('period', reference_period),
                           paste0('period', predict_period[1]))
    
  } else if (n_periods == 3) {
    period_pca <- list(
      pca_mod$x,
      predict(pca_mod, newdata = period_subset[[predict_period[1]]]),
      predict(pca_mod, newdata = period_subset[[predict_period[2]]])
    ) %>% map( ~ select(as.data.frame(.x), pca_select))
    
    names(period_pca) <- c(
      paste0('period', reference_period),
      paste0('period', predict_period[1]),
      paste0('period', predict_period[2])
    )
  }
  
  return(period_pca)
}





reclass_cluster <- function(n_periods,
                            europe_classified,
                            decades_list,
                            selector_forest,
                            europe) {
  period.reclass <- map(1:n_periods, function(period_subset) {
    europe.period <- europe_classified %>%
      filter(period %in% paste(decades_list[[period_subset]][1], decades_list[[period_subset]][2], sep = '-') |
               period == '0') %>%
      filter(!grid_id %in% selector_forest$grid_id)
    
    period.classified.shp <- europe %>%
      dplyr::left_join(europe.period, by = 'grid_id') %>%
      filter(!grid_id %in% selector_forest$grid_id) %>%
      mutate(class2 = if_else(is.na(class2), factor("Residual"), class2),
             reclass = class2) %>%
      st_as_sf()  
    
    return(period.classified.shp)
  })
  
  return(period.reclass)
}

reclass_by_neighbors <- function(n_periods,
                                 period_reclass,
                                 europe_classified) {
  period.reclass <- purrr::map(1:n_periods, function(period) {
    period.reclass.i <- purrr::map2(period_reclass[[period]]$geometry, 1:nrow(period_reclass[[period]]), function(geoms, rows) {
      indices <- st_touches(period_reclass[[period]], geoms, sparse = FALSE)
      
      poligonos_vecinos <- period_reclass[[period]][indices, ]
      
      if (period_reclass[[period]][rows, ]$class2 != "Residual" &
          !period_reclass[[period]][rows, ]$class2 %in% poligonos_vecinos$class2 &
          nrow(poligonos_vecinos) > 0) {
        print(TRUE)
        
        new_class <- poligonos_vecinos %>%
          as_tibble() %>%
          group_by(class2) %>%
          reframe(freq = n(),
                  mean_prob = mean(probability, na.rm =
                                     TRUE)) %>%
          arrange(desc(freq), desc(mean_prob)) %>%
          slice(1) %>%
          pull(class2) %>%
          as.character()
        
        period_reclass[[period]][rows, ]$reclass <- new_class
      } else {
        print(FALSE)
      }
      return(period_reclass[[period]][rows, ])
    }) %>%
      bind_rows()
    
    return(period.reclass.i)
  })
  
  return(period.reclass)
}



classify_periods <- function(n.periods,
                             europe,
                             period.data.classified,
                             class.levels,
                             var.agent,
                             grid.scale,
                             decades.list,
                             selector.forest,
                             visualize.class,
                             reclassify = TRUE) {
  if (reclassify == TRUE) {
    # Definir periodos
    period.names <- list(
      '1' = c('1985-2023'),
      '2' = c('1985-2004', '2005-2023'),
      '3' = c('1985-1997', '1998-2010', '2011-2023')
    )
    
    europe.classified <- europe %>%
      as_tibble() %>%
      select(grid_id) %>%
      mutate(!!!set_names(period.names[[as.character(n.periods)]], paste0("period", seq_len(n.periods)))) %>%
      pivot_longer(cols = contains('period'), names_to = '.x') %>%
      select(-.x) %>%
      rename(period = value) %>%
      left_join(
        period.data.classified %>% select(grid_id, period, class2, probability),
        by = c('grid_id', 'period')
      ) %>%
      mutate_at(vars(-class2), ~ ifelse(is.na(.), 0, .)) %>%
      mutate(
        class2 = factor(if_else(is.na(class2), "Residual", class2), levels = class.levels),
        period = as.factor(period),
        probability = round(probability, 1)
      )
    
    saveRDS(europe.classified,
            here(
              'output',
              'model',
              'cluster_creation',
              paste0(
                'clusters_long_database_',
                var.agent,
                '_',
                grid.scale,
                "_perc_TUM.rds"
              )
            ))
    europe.classified <- readRDS(here(
      'output',
      'model',
      'cluster_creation',
      paste0(
        'clusters_long_database_',
        var.agent,
        '_',
        grid.scale,
        "_perc_TUM.rds"
      )
    ))
    
    # Reclasificación si hay más de un período
    europe.classified <- reclass_cluster(n.periods,
                                         europe.classified,
                                         decades.list,
                                         selector.forest,
                                         europe)  ## eliminar el bind_rows()
    europe.classified <- reclass_by_neighbors(n.periods, europe.classified, europe.classified)  %>%
      map( ~ select(., grid_id, period, class2, reclass, probability)) %>%
      bind_rows()
    
    saveRDS(europe.classified,
            here(
              'output',
              'model',
              'cluster_creation',
              paste0(
                'mclust_data_',
                var.agent,
                '_',
                grid.scale,
                "_perc_TUM.rds"
              )
            ))
    
  }
  europe.classified <- readRDS(here(
    'output',
    'model',
    'cluster_creation',
    paste0('mclust_data_', var.agent, '_', grid.scale, "_perc_TUM.rds")
  ))
  
  period.classified.shp <- map(1:n.periods, function(period_subset) {
    europe.period <- europe.classified %>%
      filter(period %in% paste(decades.list[[period_subset]][1], decades.list[[period_subset]][2], sep = '-') |
               period == '0') %>%
      filter(!grid_id %in% selector.forest$grid_id &
               probability >= 0.4)
    
    return(
      europe %>%
        dplyr::left_join(as.data.frame(vect(europe.period)), by = 'grid_id') %>%
        filter(!grid_id %in% selector.forest$grid_id) %>%
        mutate(across(
          c('class2', 'reclass'), ~ if_else(is.na(.), as.factor('Residual'), .)
        ))
    )
    
  })
  
  params_orig <- list(
    var = rep('class2', n.periods),
    var_data = period.classified.shp,
    var_period = decades.list[1:n.periods],
    var_legend = rep('none', n.periods)
  )
  periods.orig.plot <- purrr::pmap(params_orig, create_cluster_map)
  periods.orig.plot <- plot_grid(
    plotlist = periods.orig.plot,
    ncol = n.periods,
    labels = LETTERS[1:n.periods],
    align = "hv",
    axis = "tblr",
    label_size = 15,
    label_y = 0.8
  )
  cowplot::save_plot(
    here(
      'output',
      'plots',
      'disturbance_analysis',
      paste0("cluster_original_maps_", var.agent, '_', grid.scale, ".svg")
    ),
    periods.orig.plot,
    base_width = 10,
    base_height =  5,
    base_asp = 1,
    dpi = 400
  )
  
  
  params <- list(
    var = rep(visualize.class, n.periods),
    var_data = period.classified.shp,
    var_period = decades.list[1:n.periods],
    var_legend = rep('none', n.periods)
  )
  params$var_legend[n.periods] <- 'bottom'  # La última leyenda se muestra
  period.classified.plot <- purrr::pmap(params, create_cluster_map)
  
  periods.plot <- plot_grid(
    plotlist = period.classified.plot,
    ncol = n.periods,
    labels = LETTERS[1:n.periods],
    align = "hv",
    axis = "tblr",
    label_size = 15,
    label_y = 0.8
  )
  
  return(periods.plot)
}


period_boxplot <- function(period_data_classified,
                           keep_vars,
                           var_colors,
                           logo_ggdraw,
                           var.agent,
                           legend_posotion = 'none') {
  period_class_long <- period_data_classified %>%
    rename('Size' = 'log_size',
           'Severity' = 'log_sev',
           'Frequency' = 'log_sumfrequency') %>%
    pivot_longer(
      cols = all_of(keep_vars),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(
      reclass = factor(reclass, labels = class.levels[class.levels != 'Residual']) %>%
        forcats::fct_rev(.),
      variable = forcats::fct_rev(variable)
    )
  
  if (var.agent == 'wind_barkbeetle') {
    wdth <- 1.2
  } else if (var.agent == 'fire') {
    wdth <- 1.05
  }
  period_classified_boxplot <- ggplot(period_class_long, aes(y = reclass, x = value, fill = reclass)) +
    geom_violin(
      width = wdth,
      trim = FALSE,
      alpha = 0.3,
      show.legend = FALSE
    ) +
    geom_boxplot(
      width = 0.3,
      alpha = 1,
      outlier.shape = NA,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = var_colors) +
    labs(x = 'Cluster', y = NULL) +
    facet_wrap(
      vars(variable),
      ncol = 3,
      scales = 'fixed',
      labeller = as_labeller(function(x)
        rep("", length(x)))  # ← esta línea elimina los títulos
      
    ) +
    envalysis::theme_publish(base_size = 15) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_blank(),
      legend.position = legend_posotion,
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 0.5,
        vjust = 0.5
      ),
      
      panel.grid.major = element_line(color = alpha("grey30", 0.2), size = 0.1),
      # panel.grid.minor = element_line(color = alpha("grey30", 0.1), size = 0.1),
      panel.grid.major.y = element_blank(),
      # Quita las líneas verticales mayores
      panel.spacing = unit(0.1, "lines"),
      # 👈 Reduce separación entre facetas
      title = element_blank()
    )
  
  period_boxplots <- cowplot::ggdraw(period_classified_boxplot)
  
  return(period_boxplots)
}

period_boxplot_absvalues <- function(period_data_classified,
                                     keep_vars,
                                     var_colors,
                                     logo_ggdraw,
                                     var.agent,
                                     legend_posotion = 'none') {
  period_class_long <- period_data_classified %>%
    pivot_longer(
      cols = all_of(keep_vars),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(
      reclass = factor(reclass, labels = class.levels[class.levels != 'Residual']) %>%
        forcats::fct_rev(.),
      variable = forcats::fct_rev(variable)
    )
  
  if (var.agent == 'wind_barkbeetle') {
    wdth <- 0.95
  } else if (var.agent == 'fire') {
    wdth <- 0.8
  }
  period_classified_boxplot <- ggplot(period_class_long, aes(y = reclass, x = value, fill = reclass)) +
    
    geom_violin(
      width = wdth,
      scale = "width",
      trim = TRUE,
      alpha = 0.3,
      show.legend = FALSE
    ) +
    geom_boxplot(
      width = 0.3,
      alpha = 1,
      outlier.shape = NA,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = var_colors) +
    labs(x = 'Cluster', y = NULL) +
    facet_wrap(
      vars(variable),
      ncol = 3,
      scales = 'free',
      labeller = as_labeller(function(x)
        rep("", length(x))) 
      
    ) +
    envalysis::theme_publish(base_size = 15) +
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_blank(),
      legend.position = legend_posotion,
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 0.5,
        vjust = 0.5
      ),
      
      panel.grid.major = element_line(color = alpha("grey30", 0.2), size = 0.1),
      panel.grid.major.y = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      title = element_blank()
    )
  
  period_boxplots <- cowplot::ggdraw(period_classified_boxplot)
  
  return(period_boxplots)
}


extract_legend <- function(period_data_classified,
                           keep_vars,
                           var_colors) {
  period_class_long <- period_data_classified %>%
    rename('Size' = 'log_size',
           'Severity' = 'log_sev',
           'Frequency' = 'log_sumfrequency') %>%
    pivot_longer(
      cols = all_of(keep_vars),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(
      reclass = factor(reclass, labels = class.levels[class.levels != 'Residual']) %>%
        forcats::fct_rev(.),
      variable = forcats::fct_rev(variable)
    )
  
  plot_with_legend <- ggplot(period_class_long, aes(y = reclass, x = value, fill = reclass)) +
    geom_boxplot(outliers = FALSE) +
    scale_fill_manual(values = var_colors) +
    labs(y = 'Cluster', x = NULL) +
    envalysis::theme_publish() +
    theme(legend.position = "bottom")  
  legends <- cowplot::get_legend(plot_with_legend, return_all = TRUE)
  
  legend <- legends[[1]]
  
  return(legend)
}




cluster_transitions_map <- function(pattern,
                                    pattern.dissolve,
                                    pattern.change,
                                    world.boundaries,
                                    europe.boundaries,
                                    var.colors.boxplot) {
  ggplot() +    
    geom_sf(
      data = pattern.dissolve,
      aes(fill = `1985-1997`),
      color = NA,
      # << AQUÍ
      alpha = 0.3,
      show.legend = FALSE
    ) +
    geom_sf(
      data = pattern.change,
      fill = 'grey10',
      color = NA,
      alpha = 0.7,
      show.legend = FALSE
    ) +
    geom_sf(
      data = pattern.change,
      aes(fill = `X2011-2023`),
      color = NA,
      # << AQUÍ
      alpha = 0.7,
      show.legend = FALSE
    ) +
    geom_sf(
      data = grids_filtered_shp,
      color = 'grey60',
      # Colores de las líneas según 'BIOME_NAME'
      fill = 'grey60',
      # Polígonos transparentes
      size = 0.05,
      alpha = 1,
      
      show.legend = FALSE
    ) +
    geom_sf(
      data = st_as_sf(europe.boundaries),
      color = 'grey20',
      fill = NA,
      size = 0.05,
      alpha = 0.7,
      show.legend = FALSE
    ) +
    
    # Escalas
    scale_fill_manual(values = c(var.colors.boxplot, var.colors.boxplot)) +
    
    labs(fill = "") +
    theme(
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      ),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA),
      legend.position = c(0.02, 0.9),
      # Posición arriba a la izquierda
      legend.justification = c(0, 1),
      # Justificación de la esquina
      legend.location = "plot",
      plot.title.position = "plot",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      text = element_text(color = "grey50", size = 14)  # Color gris para el texto general
      
    ) +
    coord_sf(crs = st_crs(3035))
}


cluster_trend_map <- function(trend_df,
                              var,
                              alpha = 0.05,
                              xlab = NULL,
                              ylab = NULL) {
  
  trend_df <- trend_df %>%
    select(grid_id, variable, trend_relative, p_value) %>%
    pivot_wider(
      names_from  = variable,
      values_from = c(trend_relative, p_value),
      names_glue  = "{variable}_{.value}"
    )
  
  annual.trends.shp <- europe %>%
    right_join(trend_df, by = "grid_id")
  
  pvalue.name <- paste0(var, '_p_value' )
  var.name <- paste0(var, '_trend_relative' )
  
  annual.trends.shp <- annual.trends.shp %>%
    dplyr::mutate(significant = !!sym(pvalue.name) <= alpha,
                  significant = if_else(is.na(significant), FALSE, significant ))
  
  trend_sig    <- annual.trends.shp %>% filter(significant)
  trend_nonsig <- annual.trends.shp %>% filter(!significant)
  
  all_values <- annual.trends.shp[[var.name]]
  
  q <- quantile(all_values, probs = c(0.0, 1), na.rm = TRUE)
  
  maxlim <- max(q)
  minlim <- min(q)
  
  color_scale <- c(minlim, maxlim)
  
  limits_val <- c(color_scale[1], 0, color_scale[2])
  
  values <- scales::rescale(c(
    color_scale[1],
    color_scale[1] / 5,
    0,
    color_scale[2] / 5,
    color_scale[2]
  ))
  
  
  # ---- 3. Plot ----
  ggplot() +
    
    # Fondo: no forest
    geom_sf(
      data = st_as_sf(grids.noforest),
      color = "grey80",
      fill = "white",
      size = 0.05,
      alpha = 0.4,
      show.legend = FALSE
    ) +
    
    geom_sf(
      data = st_as_sf(trend_sig),
      aes_string(fill = var.name),
      color = NA,
      linewidth = 0
    ) +
    
    geom_sf(
      data = st_as_sf(trend_nonsig),
      fill = "white",
      color = NA,
      linewidth = 0
    ) +
    
    geom_sf(
      data = suppressWarnings(st_centroid(st_as_sf(trend_nonsig))),
      aes_string(colour = var.name),
      size = 0.4,
      shape = 16,
      #alpha = 0.8,
      show.legend = FALSE
    ) +
    
    geom_sf(
      data = grids_filtered_shp,
      color = "grey60",
      fill = "grey60",
      size = 0.05,
      show.legend = FALSE
    ) +
    geom_sf(
      data = st_as_sf(europe.boundaries),
      color = "grey20",
      fill = NA,
      size = 0.05,
      alpha = 0.3,
      show.legend = FALSE
    ) +
    
    scale_fill_gradientn(
      colours = pals::ocean.balance(5),
      values  = values,
      limits  = color_scale,
      breaks  = limits_val,
      na.value = "white",
      guide = guide_colorbar(barwidth = 0.3, barheight = 4)
    ) +
    scale_colour_gradientn(
      colours = pals::ocean.balance(5),
      values  = values,
      limits  = color_scale,
      na.value = "white",
      guide = "none"
    ) +
    
    theme(
      plot.margin = margin(0, 0, 0, 0),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA),
      legend.position = c(0.02, 0.9),
      legend.justification = c(0, 1),
      legend.text = element_text(size = 10),
      text = element_text(color = "grey50", size = 14)
    ) +
    coord_sf(crs = st_crs(3035)) +
    labs(y = ylab, fill = xlab)
}





summarise_trends_by_group <- function(group_var = "biome",
                                      annual_regime_long,
                                      annual_trends_hist = NULL,
                                      group_colors,
                                      var_agent = "default_agent") {
  group_sym <- rlang::sym(group_var)
  
  trend_summary <- annual_regime_long %>%
    left_join(europe_df, by = "grid_id") %>%
    group_by(!!group_sym, variable, year) %>%
    summarise(valor = mean(valor, na.rm = TRUE), .groups = "drop") %>%
    group_by(!!group_sym, variable) %>%
    filter(n_distinct(year) >= 10 & !is.na(!!group_sym)) %>%
    summarise(
      trend_abs = round(trend::sens.slope(valor, conf.level = 0.95)$estimates, 4),
      trend_relative = round((
        trend::sens.slope(valor)$estimate / mean(valor)
      ) * 100, 2),
      p = trend::sens.slope(valor, conf.level = 0.95)$p.value,
      .groups = "keep"
    ) %>%
    mutate(trend_abs = if_else(variable == "Severity", trend_abs * 100, trend_abs))
  
  trend_summary_wide <- trend_summary %>%
    select(!!group_sym, variable, trend_relative, p) %>%
    pivot_wider(names_from = variable,
                values_from = c(trend_relative, p))
  
  save_plot(
    here(
      "output",
      "plots",
      "disturbance_analysis",
      paste0(group_var, "_trend_", var_agent, ".svg")
    ),
    gridExtra::tableGrob(trend_summary_wide),
    base_width = 30,
    base_height = 10,
    dpi = 400
  )
  
  
  # Crear el gráfico
  point_plot <- ggplot(
    trend_summary %>%
      mutate(variable = forcats::fct_rev(variable)),
    aes(x = trend_relative, y = variable)
  ) +
    geom_point(
      aes(
        fill = ifelse(p < 0.05, as.character(!!group_sym), NA_character_),
        color = !!group_sym
      ),
      shape = 21,
      size = 5,
      stroke = 3,
      show.legend = FALSE
    ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggrepel::geom_text_repel(
      aes(label = round(trend_relative, 2)),
      size = 3,
      color = "black",
      direction = "both",
      segment.color = "grey50",
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.1,
      force = 0.2,
      force_pull = 0.2,
      max.iter = 10000
    ) +
    scale_fill_manual(values = group_colors, na.value = "white") +
    scale_color_manual(values = group_colors) +
    envalysis::theme_publish(base_size = 12) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(angle = 90, hjust = 0.5),
      plot.margin = margin(7, 7, 7, 7)
    ) +
    xlab("Relative trend")
  
  save_plot(
    here(
      "output",
      "plots",
      "disturbance_analysis",
      paste0(group_var, "_trend_pointplot_", var_agent, ".svg")
    ),
    point_plot,
    base_width = 6,
    base_height = 6,
    base_asp = 1,
    dpi = 600
  )
  saveRDS(point_plot,
          here(
            "output",
            "plots",
            "disturbance_analysis",
            paste0(group_var, "_trend_pointplot_", var_agent, ".rds")
          ))
  
  if (!is.null(annual_trends_hist)) {
    boxplot <- ggplot(annual_trends_hist,
                      aes(x = trend_relative, y = .data[[group_var]], fill = .data[[group_var]])) +
      geom_vline(xintercept = 0,
                 linetype = "dashed",
                 size = 1) +
      geom_boxplot(outlier.shape = NA, alpha = 0.8) +
      facet_wrap( ~ variable, scales = "free_x") +
      scale_fill_manual(values = group_colors) +
      envalysis::theme_publish(base_size = 15) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        plot.margin = margin(7, 7, 7, 7)
      ) +
      labs(x = "Trend (%)", fill = group_var)
    
    save_plot(
      here(
        "output",
        "plots",
        "disturbance_analysis",
        paste0(group_var, "_trend_boxplot_", var_agent, ".svg")
      ),
      boxplot,
      base_width = 12,
      base_height = 6,
      base_asp = 1,
      dpi = 600
    )
    
    return(list(
      summary = trend_summary_wide,
      point_plot = point_plot,
      boxplot = boxplot
    ))
  } else {
    return(list(summary = trend_summary_wide, point_plot = point_plot))
  }
}


wilcoxon_tests <- function(data,
                           group_var,
                           variables = c("Size", "Severity", "Frequency")) {
  combinations <- combn(unique(data[[group_var]]), 2, simplify = FALSE)
  
  wilcoxon_test <- function(var) {
    results <- lapply(combinations, function(pair) {
      subset_data <- data[data[[group_var]] %in% pair, ]
      
      test_result <- wilcox.test(formula = as.formula(paste(var, "~", group_var)), data = subset_data)
      
      data.frame(
        group1 = pair[1],
        group2 = pair[2],
        p_value = round(test_result$p.value, 6),
        var = var
      )
    })
    
    do.call(rbind, results)
  }
  
  all_results <- lapply(variables, wilcoxon_test)
  
  return(do.call(rbind, all_results))
}




calculate_trends <- function(data,
                             europe_df = NULL,
                             periodo = NULL,
                             agrupar_por_bioma = FALSE,
                             var_agent,
                             nombre_archivo,
                             guardar_csv = TRUE) {
  if (!is.null(periodo)) {
    data <- data %>% filter(year >= periodo[1], year <= periodo[2])
  }
  
  if (agrupar_por_bioma) {
    data <- data %>%
      left_join(europe_df, by = "grid_id") %>%
      filter(!is.na(biome3))
  }
  
  agrupadores <- c(if (agrupar_por_bioma)
    "biome3", "variable", "year")
  
  resumen <- data %>%
    group_by(across(all_of(agrupadores))) %>%
    summarise(valor = mean(valor, na.rm = TRUE), .groups = "drop") %>%
    group_by(across(all_of(setdiff(
      agrupadores, "year"
    )))) %>%
    summarise(
      trend_abs = round(trend::sens.slope(valor, conf.level = 0.95)$estimates, 4),
      trend_relative = round((
        trend::sens.slope(valor)$estimate / mean(valor)
      ) * 100, 2),
      p_value = trend::sens.slope(valor, conf.level = 0.95)$p.value,
      .groups = "drop"
    ) %>%
    mutate(trend_abs = if_else(variable == "Severity", trend_abs * 100, trend_abs))
  
  dir_out <- here("output", "plots", "disturbance_analysis", "TUM")
  
  cowplot::save_plot(
    file = file.path(dir_out, paste0(nombre_archivo, var_agent, "_TUM.svg")),
    plot = gridExtra::tableGrob(resumen),
    base_width = 10,
    base_height = 5,
    dpi = 400
  )
  
  if (guardar_csv) {
    readr::write_csv(resumen, file.path(dir_out, paste0(
      nombre_archivo, var_agent, "_TUM.csv"
    )))
  }
  
  return(resumen)
}
