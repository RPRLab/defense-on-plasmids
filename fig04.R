library(patchwork)

# READ DATA --------------------------------------------------------------------

merged_metadata_raw <- readr::read_tsv(
  "data/merged_master_table.tsv"
)

to_remove <- c("Hok/Sok", "MazEF")

clean_pair <- function(subtype, type) {
  if (is.na(subtype) || is.na(type)) {
    return(c(Unified_subtype_name = NA_character_, Type = NA_character_))
  }
  
  subtype_vec <- strsplit(subtype, ";", fixed = TRUE)[[1]]
  type_vec    <- strsplit(type, ";", fixed = TRUE)[[1]]
  
  keep <- !(subtype_vec %in% to_remove)
  
  subtype_new <- subtype_vec[keep]
  type_new    <- type_vec[keep]
  
  c(
    Unified_subtype_name =
      if (length(subtype_new) == 0) NA_character_
    else paste(subtype_new, collapse = ";"),
    
    Type =
      if (length(type_new) == 0) NA_character_
    else paste(type_new, collapse = ";")
  )
}

cleaned <- t(mapply(
  clean_pair,
  merged_metadata_raw$Unified_subtype_name,
  merged_metadata_raw$Type
))

merged_metadata <- merged_metadata_raw
merged_metadata$Unified_subtype_name <- cleaned[, "Unified_subtype_name"]
merged_metadata$Type <- cleaned[, "Type"]

plsdb_metadata <- readxl::read_xlsx(
  "data/plsdb_plasmid-host_metadata_with_taxonomy.xlsx"
)

plsdb_metadata_rep <- plsdb_metadata |>
  dplyr::filter(representative == TRUE)

plsdb_filter <- plsdb_metadata |>
  dplyr::filter(representative == TRUE) |>
  dplyr::pull(plasmid_seqid)

imgpr_metadata <- readxl::read_xlsx(
  "data/imgpr_plasmid-host_metadata.xlsx"
) 

imgpr_filter <- imgpr_metadata |>
  dplyr::filter(representative == TRUE) |>
  dplyr::filter(!drop) |>
  dplyr::pull(plasmid_seqid)

merged_metadata |>
  dplyr::filter(stringr::str_detect(id, "IMGPR")) |>
  dplyr::filter(id %in% imgpr_filter) |>
  nrow()

imgpr_metadata |>
  dplyr::filter(representative == TRUE) |>
  dplyr::summarise(n = dplyr::n(), .by = drop) |>
  dplyr::mutate(p = n / sum(n) * 100) |>
  dplyr::summarise(
    total = sum(n),
    filtered = sum(n[!drop]),
    removed_pct = sum(p[drop])
  ) |>
  glue::glue_data(
    "Total (representative) plasmids: {format(total, big.mark = ',')}\n",
    "Final filtered plasmids: {format(filtered, big.mark = ',')}\n",
    "Proportion removed: {format(removed_pct, digits = 3)}%"
  )

plsdb_plasmid_defense_raw <- readxl::read_xlsx(
  "data/plsdb_plasmid_defense.xlsx"
)

plsdb_plasmid_defense <- plsdb_plasmid_defense_raw |>
  # REMOVE HOK/SOK AND MAZEF
  dplyr::filter(!type %in% c("Hok/Sok", "MazEF")) |>
  dplyr::mutate(replicon = "plasmid", feature = "defense") |>
  dplyr::mutate(
    feature_id = glue::glue(
      "{plasmid_seqid}", "_defense_",
      "{stringr::str_pad(dplyr::row_number(), 2, 'left', '0')}"
    ),
    .by = plasmid_seqid
  )

imgpr_plasmid_defense_raw <- readxl::read_xlsx(
  "data/imgpr_plasmid_defense.xlsx"
)

imgpr_plasmid_defense <- imgpr_plasmid_defense_raw |>
  # REMOVE HOK/SOK AND MAZEF
  dplyr::filter(!type %in% c("Hok/Sok", "MazEF")) |>
  dplyr::mutate(replicon = "plasmid", feature = "defense") |>
  dplyr::mutate(
    feature_id = glue::glue(
      "{plasmid_seqid}", "_defense_",
      "{stringr::str_pad(dplyr::row_number(), 2, 'left', '0')}"
    ),
    .by = plasmid_seqid
  )

systems_count_per_source <- imgpr_metadata |>
  dplyr::filter(representative == TRUE) |>
  dplyr::filter(!drop) |>
  dplyr::distinct(plasmid_seqid, source_type) |>
  dplyr::left_join(
    imgpr_plasmid_defense,
    by = dplyr::join_by(plasmid_seqid)
  ) |>
  dplyr::bind_rows(
    plsdb_metadata |>
      dplyr::filter(representative == TRUE) |>
      dplyr::distinct(plasmid_seqid) |>
      dplyr::left_join(
        plsdb_plasmid_defense,
        by = dplyr::join_by(plasmid_seqid)
      )
  ) |>
  dplyr::mutate(
    source_type = dplyr::if_else(is.na(source_type), "Isolate", source_type)
  ) |>
  dplyr::mutate(
    category = dplyr::if_else(
      source_type == "Isolate",
      "Isolate",
      "Metagenomic"
    )
  ) |>
  dplyr::filter(!is.na(type)) |>
  dplyr::summarise(n = dplyr::n(), .by = category)

systems_count_per_source |>
  purrr::pwalk(
    \(category, n) {
      cli::cli_alert_info(
        "{n} plasmid-encoded defense in {category} derived plasmids"
      )
    }
  )

imgpr_plasmid_source <- imgpr_metadata |>
  dplyr::filter(representative == TRUE) |>
  dplyr::filter(!drop) |>
  dplyr::mutate(
    source = dplyr::case_when(
      source_type == "Isolate" ~ "Isolate",
      .default = "Metagenomic"
    )
  ) |>
  dplyr::distinct(plasmid_seqid, drop, source)

imgpr_length <- imgpr_metadata |>
  dplyr::filter(representative == TRUE) |>
  dplyr::filter(!drop) |>
  dplyr::distinct(plasmid_seqid, plasmid_length)

imgpr_has_defense_vs_length <- imgpr_plasmid_source |>
  dplyr::right_join(imgpr_length, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::left_join(
    imgpr_plasmid_defense,
    by = dplyr::join_by(plasmid_seqid)
  ) |>
  dplyr::mutate(
    has_defense = dplyr::if_else(
      any(feature == "defense", na.rm = TRUE),
      "With Defense",
      "No Defense"
    ),
    .by = plasmid_seqid
  ) |>
  dplyr::distinct(plasmid_seqid, plasmid_length, source, has_defense)


plot_S12 <- imgpr_has_defense_vs_length |>
  ggplot2::ggplot(
    ggplot2::aes(x = plasmid_length, fill = has_defense)
  ) +
  ggplot2::geom_density(adjust = 0.5, alpha = 0.5, colour = NA) +
  ggplot2::scale_fill_manual(values = c("#c1c6d4", "#6e788d")) +
  # ggplot2::scale_colour_manual(values = c("#c5cad7", "#6e788d")) +
  ggplot2::scale_x_log10(
    expand = ggplot2::expansion(mult = c(0, 0)),
    labels = scales::label_log(base = 10),
    breaks = c(10, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7),
    limits = c(10^2, 10^7)
  ) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
  ggplot2::labs(x = "Plasmid length (bp)", y = "Density") +
  ggplot2::geom_vline(
    xintercept = c(10000, 50000, 300000),
    colour = "black",
    linewidth = 0.24,
    linetype = "dashed"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    axis.text.x = ggplot2::element_text(size = 7, colour = "black"),
    axis.text.y = ggplot2::element_text(size = 7, colour = "black"),
    axis.ticks.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.ticks.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(linewidth = 0.24),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(
      linewidth = 0.24,
      fill = NA,
      colour = NA
    ),
    strip.text = ggplot2::element_text(
      size = 7,
      colour = "black",
      margin = ggplot2::margin(0.08, 0, 0.08, 0, "cm")
    ),
    legend.key.height = grid::unit(0.25, "cm"),
    legend.key.width = grid::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "top",
    legend.title = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank()
  ) +
  ggplot2::facet_wrap(~source, ncol = 1) +
  ggplot2::coord_cartesian(clip = "off")

plot_S12

ggplot2::ggsave(
  filename = "plots/figS12.pdf",
  width = 60,
  height = 70,
  units = "mm",
  dpi = 300
)


imgpr_filter <- imgpr_metadata |>
  dplyr::filter(representative == TRUE) |>
  dplyr::filter(!drop) |>
  dplyr::select(plasmid_seqid) |>
  dplyr::pull(plasmid_seqid)

all_reps <- merged_metadata |>
  dplyr::filter(representative == "Yes") |>
  dplyr::select(id) |>
  dplyr::rename(plasmid_seqid = id)

plsdb_reps <- all_reps |>
  dplyr::filter(!stringr::str_detect(plasmid_seqid, "^IMGPR_"))

imgpr_reps <- all_reps |>
  dplyr::filter(stringr::str_detect(plasmid_seqid, "^IMGPR_")) |>
  dplyr::filter(plasmid_seqid %in% imgpr_filter)

merged_reps <- dplyr::bind_rows(plsdb_reps, imgpr_reps)

merged_mobility <- merged_metadata |>
  dplyr::mutate(mobility = dplyr::if_else(is.na(predicted_mobility), "Non-mobilizable", stringr::str_to_sentence(predicted_mobility))) |>
  dplyr::select(id, mobility) |>
  dplyr::rename(plasmid_seqid = id)

merged_length <- merged_metadata |>
  dplyr::select(id, length) |>
  dplyr::rename(plasmid_seqid = id, plasmid_length = length)

merged_has_defense <- merged_metadata |>
  dplyr::mutate(has_defense = dplyr::if_else(!is.na(Type), "With Defense", "No Defense")) |>
  dplyr::select(id, has_defense) |>
  dplyr::rename(plasmid_seqid = id)

merged_has_amr <- merged_metadata |>
  dplyr::mutate(has_amr = dplyr::if_else(!is.na(Class), "With AMR", "No AMR")) |>
  dplyr::select(id, has_amr) |>
  dplyr::rename(plasmid_seqid = id)

merged_source <- merged_metadata |>
  dplyr::mutate(
    source = dplyr::case_when(
      source_type == "Isolate" ~ "Isolate",
      is.na(source_type) ~ "Isolate",
      .default = "Metagenomic"
    )
  ) |>
  dplyr::select(id, source) |>
  dplyr::rename(plasmid_seqid = id)

plsdb_is_complete <- plsdb_metadata |>
  dplyr::select(plasmid_seqid, plasmid_complete) |>
  dplyr::rename(is_complete = plasmid_complete)

imgpr_is_complete <- merged_metadata |>
  dplyr::filter(stringr::str_detect(id, "^IMGPR")) |>
  dplyr::select(id, putatively_complete) |>
  dplyr::rename(is_complete = putatively_complete, plasmid_seqid = id) |>
  dplyr::mutate(is_complete = dplyr::if_else(is_complete == "Yes", TRUE, FALSE))

merged_is_complete <- dplyr::bind_rows(plsdb_is_complete, imgpr_is_complete)

all_reps_useful_meta <- all_reps |>
  dplyr::left_join(merged_mobility, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::left_join(merged_length, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::left_join(merged_source, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::left_join(merged_has_defense, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::left_join(merged_has_amr, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::left_join(merged_is_complete, by = dplyr::join_by(plasmid_seqid))

# FIGURE 4A --------------------------------------------------------------------

has_defense_vs_len_mob <- all_reps_useful_meta |>
  dplyr::mutate(
    plasmid_seqid,
    plasmid_length,
    mobility = factor(
      mobility, levels = c("Conjugative", "Mobilizable", "Non-mobilizable")
    ),
    has_defense = has_defense == "With Defense",
    .keep = "none"
  ) |>
  dplyr::mutate(
    log_length = log10(plasmid_length),
    length_kb = plasmid_length / 1000,
    log10_kb = log10(length_kb)
  )

has_defense_vs_len_mob |>
  dplyr::group_by(mobility) |>
  dplyr::summarise(
    n = dplyr::n(),
    defense_positive = sum(has_defense),
    defense_negative = sum(!has_defense),
    defense_prevalence = mean(has_defense),
    .groups = "drop"
  )

support <- has_defense_vs_len_mob |>
  dplyr::summarise(
    n = dplyr::n(),
    lower_01 = stats::quantile(
      plasmid_length, probs = 0.01, na.rm = TRUE
    ),
    upper_99 = stats::quantile(
      plasmid_length, probs = 0.99, na.rm = TRUE
    ),
    .by = mobility
  )

support

common_lower <- max(support$lower_01)
common_upper <- min(support$upper_99)

cli::cli_alert_info(
  "Shared support: {round(common_lower / 1000, 1)} to {round(common_upper / 1000, 1)} kb"
)

has_defense_vs_len_mob_common <- has_defense_vs_len_mob |>
  dplyr::filter(
    plasmid_length >= common_lower,
    plasmid_length <= common_upper
  )

retention_summary <- has_defense_vs_len_mob |>
  dplyr::mutate(
    retained = plasmid_length >= common_lower &
      plasmid_length <= common_upper
  ) |>
  dplyr::summarise(
    original_n = dplyr::n(),
    retained_n = sum(retained),
    retained_percent = mean(retained) * 100,
    .by = mobility
  )

retention_summary

retention_by_defense <- has_defense_vs_len_mob |>
  dplyr::mutate(
    retained = plasmid_length >= common_lower &
      plasmid_length <= common_upper
  ) |>
  dplyr::summarise(
    original_n = dplyr::n(),
    retained_n = sum(retained),
    retained_percent = mean(retained) * 100,
    .by = c(mobility, has_defense)
  )

retention_by_defense

has_defense_vs_len_mob_common |>
  dplyr::summarise(
    n = dplyr::n(),
    min_kb = min(length_kb),
    median_kb = median(length_kb),
    max_kb = max(length_kb),
    defense_prevalence = mean(has_defense),
    .by = mobility
  )

fit_gam <- mgcv::gam(
  has_defense ~ mobility + s(log_length, by = mobility, k = 6),
  data = has_defense_vs_len_mob_common,
  family = stats::binomial(),
  method = "REML"
)

summary(fit_gam)

mgcv::gam.check(fit_gam)

fit_shared_ml <- mgcv::gam(
  has_defense ~ mobility + s(log_length, k = 6),
  data = has_defense_vs_len_mob_common,
  family = stats::binomial(),
  method = "ML"
)

fit_specific_ml <- mgcv::gam(
  has_defense ~ mobility + s(log_length, by = mobility, k = 6),
  data = has_defense_vs_len_mob_common,
  family = stats::binomial(),
  method = "ML"
)

smooth_comparison <- stats::anova(
  fit_shared_ml,
  fit_specific_ml,
  test = "Chisq"
)

smooth_comparison

mobility_levels <- levels(has_defense_vs_len_mob_common$mobility)

length_sequence <- 10^seq(
  log10(common_lower),
  log10(common_upper),
  length.out = 400
)

prediction_grid <- tidyr::expand_grid(
  plasmid_length = length_sequence,
  mobility = mobility_levels
) |>
  dplyr::mutate(
    mobility = factor(mobility, levels = mobility_levels),
    log_length = log10(plasmid_length),
    length_kb = plasmid_length / 1000,
    log10_kb = log10(length_kb)
  )

beta_hat <- stats::coef(fit_gam)

model_vcov <- stats::vcov(fit_gam, unconditional = TRUE)

prediction_matrix <- stats::predict(
  fit_gam, newdata = prediction_grid,  type = "lpmatrix"
)

prediction_eta <- drop(
  prediction_matrix %*% beta_hat
)

prediction_eta_se <- sqrt(
  rowSums((prediction_matrix %*% model_vcov) * prediction_matrix)
)

prediction_data <- prediction_grid |>
  dplyr::mutate(
    eta = prediction_eta,
    eta_se = prediction_eta_se,
    predicted_probability = stats::plogis(eta),
    lower_ci = stats::plogis(eta - 1.96 * eta_se),
    upper_ci = stats::plogis(eta + 1.96 * eta_se)
  )

calculate_gam_contrast <- function(
    comparator,
    fitted_model,
    lengths,
    mobility_levels
) {
  
  newdata_conjugative <- tibble::tibble(
    plasmid_length = lengths,
    mobility = factor(
      "Conjugative",
      levels = mobility_levels
    ),
    log_length = log10(plasmid_length),
    length_kb = plasmid_length / 1000,
    log10_kb = log10(length_kb)
  )
  
  newdata_comparator <- tibble::tibble(
    plasmid_length = lengths,
    mobility = factor(
      comparator,
      levels = mobility_levels
    ),
    log_length = log10(plasmid_length),
    length_kb = plasmid_length / 1000,
    log10_kb = log10(length_kb)
  )
  
  X_conjugative <- stats::predict(
    fitted_model,
    newdata = newdata_conjugative,
    type = "lpmatrix"
  )
  
  X_comparator <- stats::predict(
    fitted_model,
    newdata = newdata_comparator,
    type = "lpmatrix"
  )
  
  beta <- stats::coef(fitted_model)
  
  V_beta <- stats::vcov(
    fitted_model,
    unconditional = TRUE
  )
  
  eta_conjugative <- drop(
    X_conjugative %*% beta
  )
  
  eta_comparator <- drop(
    X_comparator %*% beta
  )
  
  probability_conjugative <- stats::plogis(
    eta_conjugative
  )
  
  probability_comparator <- stats::plogis(
    eta_comparator
  )
  
  gradient_conjugative <-
    X_conjugative *
    as.numeric(
      probability_conjugative *
        (1 - probability_conjugative)
    )
  
  gradient_comparator <-
    X_comparator *
    as.numeric(
      probability_comparator *
        (1 - probability_comparator)
    )
  
  gradient_difference <-
    gradient_conjugative -
    gradient_comparator
  
  contrast_se <- sqrt(
    rowSums(
      (gradient_difference %*% V_beta) *
        gradient_difference
    )
  )
  
  contrast_estimate <-
    probability_conjugative -
    probability_comparator
  
  tibble::tibble(
    plasmid_length = lengths,
    length_kb = lengths / 1000,
    log10_kb = log10(lengths / 1000),
    comparator = comparator,
    conjugative_probability = probability_conjugative,
    comparator_probability = probability_comparator,
    estimate = contrast_estimate,
    standard_error = contrast_se,
    lower_ci = contrast_estimate -
      1.96 * contrast_se,
    upper_ci = contrast_estimate +
      1.96 * contrast_se
  )
}

contrast_data <- dplyr::bind_rows(
  calculate_gam_contrast(
    comparator = "Mobilizable",
    fitted_model = fit_gam,
    lengths = length_sequence,
    mobility_levels = mobility_levels
  ),
  calculate_gam_contrast(
    comparator = "Non-mobilizable",
    fitted_model = fit_gam,
    lengths = length_sequence,
    mobility_levels = mobility_levels
  )
) |>
  dplyr::mutate(
    comparator = factor(
      comparator,
      levels = c(
        "Mobilizable",
        "Non-mobilizable"
      )
    )
  )

contrast_significance <- contrast_data |>
  dplyr::mutate(
    direction = dplyr::case_when(
      lower_ci > 0 ~ "Conjugative higher",
      upper_ci < 0 ~ "Conjugative lower",
      TRUE ~ "No clear difference"
    )
  )

contrast_significance |>
  dplyr::summarise(
    minimum_length_kb = min(length_kb),
    maximum_length_kb = max(length_kb),
    .by = c(comparator, direction)
  )

contrast_significance |>
  dplyr::filter(direction == "Conjugative higher") |>
  dplyr::summarise(
    first_supported_length_kb = min(length_kb),
    .by = comparator
  )

mobility_colours <- c(
  "Conjugative" = "#E69F00",
  "Mobilizable" = "#009E73",
  "Non-mobilizable" = "#0072B2"
)

x_breaks_kb <- c(10, 20, 40, 80, 160)
x_breaks_log10 <- log10(x_breaks_kb)


plot_fig04_A_gam <- ggplot2::ggplot(
  prediction_data,
  ggplot2::aes(
    x = log10_kb,
    y = predicted_probability,
    colour = mobility,
    fill = mobility
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = lower_ci,
      ymax = upper_ci
    ),
    alpha = 0.14,
    colour = NA
  ) +
  ggplot2::geom_line(
    linewidth = 0.48,
  ) +
  ggplot2::scale_x_continuous(
    breaks = x_breaks_log10,
    labels = x_breaks_kb,
    limits = range(
      prediction_data$log10_kb
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0)
    )
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    limits = c(0, NA),
    expand = ggplot2::expansion(
      mult = c(0, 0.04)
    )
  ) +
  ggplot2::scale_colour_manual(
    values = mobility_colours
  ) +
  ggplot2::scale_fill_manual(
    values = mobility_colours
  ) +
  ggplot2::labs(
    x = "Plasmid length (kb; log scale)",
    y = "Predicted defense\nprevalence",
    colour = NULL,
    fill = NULL
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    axis.text.x = ggplot2::element_text(size = 7, colour = "black"),
    axis.text.y = ggplot2::element_text(size = 7, colour = "black"),
    axis.ticks.x = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "black"
    ),
    axis.ticks.y = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "black"
    ),
    axis.line = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    legend.position = "none",
    plot.background = ggplot2::element_blank()
  )

plot_fig04_A_gam


plot_fig04_A_contrast <- ggplot2::ggplot(
  contrast_data,
  ggplot2::aes(
    x = log10_kb,
    y = estimate,
    colour = comparator,
    fill = comparator,
    group = comparator
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.24,
    colour = "black"
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = lower_ci,
      ymax = upper_ci
    ),
    alpha = 0.16,
    colour = NA
  ) +
  ggplot2::geom_line(
    linewidth = 0.48,
  ) +
  ggplot2::scale_colour_manual(
    values = mobility_colours
  ) +
  ggplot2::scale_fill_manual(
    values = mobility_colours
  ) +
  ggplot2::scale_x_continuous(
    breaks = x_breaks_log10,
    labels = x_breaks_kb,
    limits = range(contrast_data$log10_kb),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  ggplot2::labs(
    x = "Plasmid length (kb; log scale)",
    y = "Predicted prevalence difference\n(conjugative - comparator)",
    colour = "Comparator",
    fill = "Comparator"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    axis.text.x = ggplot2::element_text(size = 7, colour = "black"),
    axis.text.y = ggplot2::element_text(size = 7, colour = "black"),
    axis.ticks.x = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "black"
    ),
    axis.ticks.y = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "black"
    ),
    axis.line = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    legend.position = "none",
    plot.background = ggplot2::element_blank()
  )

plot_fig04_A_contrast

density_data <- has_defense_vs_len_mob_common |>
  dplyr::mutate(
    log10_kb = log10(plasmid_length / 1000)
  )

plot_fig04_A_density <- ggplot2::ggplot(
  density_data,
  ggplot2::aes(
    x = log10_kb,
    colour = mobility
  )
) +
  ggplot2::geom_density(
    alpha = 0.22,
    linewidth = 0.48,
    adjust = 1
  ) +
  ggplot2::scale_x_continuous(
    breaks = x_breaks_log10,
    labels = x_breaks_kb,
    limits = range(
      prediction_data$log10_kb
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0)
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0, 0)
    )
  ) +
  ggplot2::scale_colour_manual(
    values = mobility_colours
  ) +
  ggplot2::labs(
    x = "Plasmid length (kb; log scale)",
    y = "Density"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    axis.text.x = ggplot2::element_text(size = 7, colour = "black"),
    axis.text.y = ggplot2::element_text(size = 7, colour = "black"),
    axis.ticks.x = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "black"
    ),
    axis.ticks.y = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "black"
    ),
    axis.line = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    legend.position = "none",
    plot.background = ggplot2::element_blank()
  )

plot_fig04_A_density

plot_fig04_A <-
  plot_fig04_A_density + plot_fig04_A_gam + plot_fig04_A_contrast +
  patchwork::plot_layout(ncol = 1, heights = c(0.5, 1, 1))

plot_fig04_A

plot_fig04_A |>
  ggplot2::ggsave(
  filename = "plots/fig04_A.pdf", width = 62, height = 90, units = "mm",
)



# FIGURE 4B --------------------------------------------------------------------


combined_plasmid_feature_categories <- all_reps_useful_meta |>
  dplyr::mutate(
    has_defense = has_defense == "With Defense",
    has_amr = has_amr == "With AMR"
  ) |>
  dplyr::summarise(
    has_defense = any(has_defense, na.rm = TRUE),
    has_amr = any(has_amr, na.rm = TRUE),
    plasmid_length = dplyr::first(plasmid_length),
    .by = c(plasmid_seqid, source)
  ) |>
  dplyr::mutate(
    category = purrr::map2_chr(
      has_defense,
      has_amr,
      \(defense, amr) {
        sets <- c(
          if (isTRUE(defense)) "DS",
          if (isTRUE(amr)) "AMR"
        )
        
        if (length(sets) == 0) {
          "NA"
        } else {
          paste(sort(sets), collapse = "-")
        }
      }
    ),
    size_bin = dplyr::case_when(
      dplyr::between(plasmid_length, 0, 10000) ~ "0-10",
      dplyr::between(plasmid_length, 10001, 50000) ~ "10-50",
      dplyr::between(plasmid_length, 50001, 300000) ~ "50-300",
      plasmid_length > 300000 ~ ">300",
      .default = NA_character_
    )
  )

venn_counts <- combined_plasmid_feature_categories |>
  dplyr::count(category) |>
  tidyr::complete(
    category = c("DS", "AMR", "AMR-DS", "NA"),
    fill = list(n = 0)
  )

venn_counts |>
  dplyr::mutate(p = n / sum(n)) |>
  purrr::pwalk(
    \(category, n, p) {
      cli::cli_alert_info(
        "{category}: {round(p * 100, 2)}% ({n})"
      )
    }
  )

n_ds <- venn_counts |>
  dplyr::filter(category == "DS") |>
  dplyr::pull(n)

n_amr <- venn_counts |>
  dplyr::filter(category == "AMR") |>
  dplyr::pull(n)

n_overlap <- venn_counts |>
  dplyr::filter(category == "AMR-DS") |>
  dplyr::pull(n)

n_neither <- venn_counts |>
  dplyr::filter(category == "NA") |>
  dplyr::pull(n)

total_population <- sum(venn_counts$n)

fit_04B <- eulerr::euler(
  c(
    DS = n_ds,
    AMR = n_amr,
    "DS&AMR" = n_overlap
  ),
  input = "disjoint",
  shape = "circle",
  complement = n_neither
)

box_mm <- 28.5

circle_diameter_mm <- 2 * box_mm / sqrt(pi)

measured_box_mm <- 26.5

scale_factor <- box_mm / measured_box_mm

euler_box_plot <- plot(
  fit_04B,
  
  fills = list(
    mode = "disjoint",
    fill = c(
      DS       = "#87a9ca",
      AMR      = "#e3c589",
      "DS&AMR" = "#b0a276"
    )
  ),
  
  edges = list(
    col = NA,
    lwd = 0
  ),
  
  complement = list(
    fill = NA,
    col = NA,
    lwd = 0,
    label = ""
  ),
  
  labels = FALSE,
  quantities = FALSE,
  
  margin = grid::unit(0, "mm"),
  padding = grid::unit(0, "mm")
)

grDevices::pdf(
  file = "plots/fig04_B.pdf",
  width = box_mm / 25.4,
  height = box_mm / 25.4,
  useDingbats = FALSE
)

grid::grid.newpage()

grid::grid.circle(
  x = 0.5,
  y = 0.5,
  r = grid::unit(circle_diameter_mm / 2, "mm"),
  gp = grid::gpar(
    fill = "#EFEFEF",
    col = NA
  )
)

grid::pushViewport(
  grid::viewport(
    x = 0.5,
    y = 0.5,
    width = grid::unit(box_mm * scale_factor, "mm"),
    height = grid::unit(box_mm * scale_factor, "mm"),
    clip = "on"
  )
)

print(
  euler_box_plot,
  newpage = FALSE
)

grid::popViewport()
grDevices::dev.off()

# FIGURE 4C --------------------------------------------------------------------

A_x_order_features <- c(
  "AMR",
  "DS",
  "AMR-DS"
)

plot_04_C <- combined_plasmid_feature_categories |>
  dplyr::summarise(n = dplyr::n(), .by = c(source, category)) |>
  dplyr::mutate(p = n / sum(n), .by = c(source)) |>
  dplyr::filter(category != "NA") |>
  ggplot2::ggplot(ggplot2::aes(
    x = factor(category, levels = A_x_order_features),
    y = p,
    fill = source
  )) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::scale_fill_manual(values = c("#e3e5eb", "#a7adb7")) +
  ggupset::axis_combmatrix(sep = "-", clip = "off") +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(),
    labels = scales::percent_format()
  ) +
  ggplot2::labs(x = "", y = "Proportion of plasmids") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    axis.text.x = ggplot2::element_text(size = 7, colour = "black"),
    axis.text.y = ggplot2::element_text(size = 7, colour = "black"),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(
      linewidth = 0.24,
      fill = NA,
      colour = NA
    ),
    strip.text = ggplot2::element_text(
      size = 7,
      colour = "black",
      margin = ggplot2::margin(0.08, 0, 0.08, 0, "cm")
    ),
    legend.key.height = grid::unit(0.25, "cm"),
    legend.key.width = grid::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "top",
    legend.title = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank()
  ) +
  ggupset::theme_combmatrix(
    combmatrix.panel.point.size = 1.5,
    combmatrix.panel.line.size = 1,
    combmatrix.panel.striped_background = FALSE
  )

plot_04_C

plot_04_C |>
  ggplot2::ggsave(
    filename = "plots/fig04_C.pdf",
    width = 28.45,
    height = 57.35,
    units = "mm",
    dpi = 300
  )

defense_vs_amr <- all_reps_useful_meta |>
  dplyr::mutate(
    defense = has_defense == "With Defense",
    amr = has_amr == "With AMR"
  ) |>
  dplyr::group_by(source) |>
  dplyr::group_modify(\(.x, .y) {
    tab <- with(.x, table(defense, amr))
    test <- mcnemar.test(tab)
    
    tibble::tibble(
      n = nrow(.x),
      defense_n = sum(.x$defense),
      defense_percent = 100 * mean(.x$defense),
      amr_n = sum(.x$amr),
      amr_percent = 100 * mean(.x$amr),
      chi_squared = unname(test$statistic),
      df = unname(test$parameter),
      p_value = test$p.value
    )
  })

defense_vs_amr |>
  purrr::pwalk(
    \(source, n, defense_n, defense_percent, amr_n, amr_percent, chi_squared, df, p_value) {
      cli::cli_alert_info(
        "{source}:
      Defense detected in {format(defense_n, big.mark = ',')} plasmids ({round(defense_percent)}%), 
      AMR detected in {format(amr_n, big.mark = ',')} plasmids ({round(amr_percent)}%)
      McNemar's chi-squared = {format(round(chi_squared, 1), big.mark = ',')}
      p = {p_value}\n\n"
      )
    }
  )

# FIGURE 4D --------------------------------------------------------------------

functional_categories_plsdb_raw <- readr::read_tsv(
  "data/plsdb_gene_annotations.tsv"
)

plsdb_plasmid_defense_start <- plsdb_plasmid_defense |>
  dplyr::distinct(plasmid_seqid, start, .keep_all = TRUE) |>
  dplyr::rename(system_end = end)

functional_categories_plsdb_processed <- functional_categories_plsdb_raw |>
  tidyr::separate_longer_delim(merged_operon, ",") |>
  dplyr::select(-type) |>
  dplyr::rename(plasmid_seqid = seqid) |>
  dplyr::left_join(
    plsdb_plasmid_defense_start,
    by = dplyr::join_by(plasmid_seqid, start)
  ) |>
  dplyr::mutate(
    type = if (all(is.na(type))) {type[NA_integer_]} else {max(type, na.rm = TRUE)},
    subtype = if (all(is.na(subtype))) {subtype[NA_integer_]} else {max(subtype, na.rm = TRUE)}, .by = c(plasmid_seqid, merged_operon)
  ) |>
  dplyr::mutate(
    merged_type = type, merged_subtype = subtype
  ) |>
  dplyr::rename(seqid = plasmid_seqid) |>
  dplyr::select(names(functional_categories_plsdb_raw))

functional_categories_plsdb <- functional_categories_plsdb_processed |>
  dplyr::right_join(all_reps_useful_meta, by = dplyr::join_by(seqid == plasmid_seqid)) |>
  dplyr::mutate(
    category = dplyr::case_when(
      !is.na(merged_type) ~ "Defense",
      !is.na(`Element type`) ~ "AMR",
      !is.na(target) ~ "Mobilome",
      COG_category == "-" ~ "Unknown",
      COG_category == "S" ~ "Unknown",
      is.na(COG_category) ~ "Unknown",
      .default = COG_category
    )
  ) |>
  dplyr::select(seqid, prot_id, category)

functional_categories_imgpr_raw <- readr::read_tsv(
  "data/imgpr_gene_annotations.tsv"
)

eggnog <- readr::read_tsv(
  "data/imgpr_eggnog.tsv"
)

functional_categories_imgpr_raw_processed <- functional_categories_imgpr_raw |>
  dplyr::mutate(
    merged_type = dplyr::case_when(
      # REMOVE HOK/SOK AND MAZEF
      merged_type %in% c("Mok_Hok_Sok", "MazEF") ~ NA,
      stringr::str_detect(merged_type, "^Anti_") ~ NA,
      stringr::str_detect(merged_type, "^PDC-") ~ NA,
      .default = merged_type
    )
  ) |>
  dplyr::select(seqid, prot_id, merged_type, `Element type`, target) |>
  dplyr::left_join(eggnog, by = dplyr::join_by(prot_id))

functional_categories_imgpr <- functional_categories_imgpr_raw_processed |>
  dplyr::mutate(seqid = stringr::str_remove_all(seqid, "\\|.*")) |>
  dplyr::right_join(all_reps_useful_meta, by = dplyr::join_by(seqid == plasmid_seqid)) |>
  dplyr::mutate(
    category = dplyr::case_when(
      !is.na(merged_type) ~ "Defense",
      !is.na(`Element type`) ~ "AMR",
      !is.na(target) ~ "Mobilome",
      COG_category == "-" ~ "Unknown",
      COG_category == "S" ~ "Unknown",
      is.na(COG_category) ~ "Unknown",
      .default = COG_category
    )
  ) |>
  dplyr::select(seqid, prot_id, category)

functional_categories <-
  dplyr::bind_rows(
    functional_categories_plsdb, 
    functional_categories_imgpr
  ) |>
  dplyr::filter(
    seqid %in% all_reps$plasmid_seqid
  ) |>
  dplyr::filter(!is.na(prot_id))

special_cats <- c("AMR", "Defense", "Mobilome", "Unknown")

functional_categories_simple <- functional_categories |>
  dplyr::select(seqid, prot_id, category) |>
  dplyr::mutate(category = stringr::str_remove_all(category, ",")) |>
  dplyr::right_join(all_reps_useful_meta, by = dplyr::join_by(seqid == plasmid_seqid)) |>
  dplyr::filter(!is.na(prot_id)) |>
  dplyr::select(seqid, prot_id, category, plasmid_length, source) |>
  dplyr::mutate(
    length_bin = dplyr::case_when(
      plasmid_length <= 1e4  ~ "0-10",
      plasmid_length <= 5e4  ~ "10-50",
      plasmid_length <= 3e5  ~ "50-300",
      plasmid_length > 3e5   ~ ">300"
    ),
    cat_list = dplyr::if_else(
      category %in% special_cats,
      as.list(category),
      stringr::str_split(category, "")
    )
  ) |>
  tidyr::unnest(cat_list) |>
  dplyr::rename(category_expanded = cat_list) |>
  dplyr::filter(category_expanded != "S")

genes_per_plasmid <- functional_categories_simple |>
  dplyr::distinct(seqid, length_bin, source, prot_id) |>
  dplyr::count(seqid, length_bin, source, name = "n_genes")

genes_per_category <- functional_categories_simple |>
  dplyr::distinct(seqid, prot_id, length_bin, source, category_expanded) |>
  dplyr::count(seqid, length_bin, source, category_expanded, name = "n_cat_genes")

resistance_plasmids <- functional_categories_simple |>
  dplyr::filter(category_expanded %in% c("AMR", "Defense")) |>
  dplyr::distinct(seqid)

genes_per_plasmid_res <- genes_per_plasmid |>
  dplyr::semi_join(resistance_plasmids, by = "seqid")

genes_per_category_res <- genes_per_category |>
  dplyr::semi_join(resistance_plasmids, by = "seqid")

props_res <- genes_per_category_res |>
  dplyr::left_join(
    genes_per_plasmid_res,
    by = c("seqid", "length_bin", "source")
  ) |>
  dplyr::mutate(prop = n_cat_genes / n_genes)

result_res <- props_res |>
  dplyr::group_by(length_bin, source, category = category_expanded) |>
  dplyr::summarise(
    mean_prop  = mean(prop),
    sd_prop    = sd(prop),
    n_plasmids = dplyr::n_distinct(seqid),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    se = sd_prop / sqrt(n_plasmids),
    ci_low = pmax(mean_prop - 1.96 * se, 0),
    ci_high = pmin(mean_prop + 1.96 * se, 1)
  )

bin_sizes_res <- genes_per_plasmid_res |>
  dplyr::distinct(seqid, length_bin, source) |>
  dplyr::count(length_bin, source, name = "n_all")

source_map <- c(Isolate = "I", Metagenomic = "M")
pad_width <- 8

facet_labels_res <- bin_sizes_res |>
  dplyr::mutate(src_short = dplyr::recode(source, !!!source_map)) |>
  dplyr::select(length_bin, src_short, n_all) |>
  tidyr::complete(length_bin, src_short, fill = list(n_all = 0L)) |>
  tidyr::pivot_wider(names_from = src_short, values_from = n_all) |>
  dplyr::mutate(
    I_txt = stringr::str_pad(scales::comma(I), width = pad_width, side = "left"),
    M_txt = stringr::str_pad(scales::comma(M), width = pad_width, side = "left"),
    label = paste0(length_bin, " kb\n", I_txt, "\n", M_txt)
  )

facet_labels_res <- setNames(facet_labels_res$label, facet_labels_res$length_bin)

facet_labels_res

facet_order <- c("0-10", "10-50", "50-300", ">300")
cog_order <- c("AMR","Defense","Mobilome","Unknown","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","T","U","V","W","X","Y","Z")

pd <- ggplot2::position_dodge(width = 0.9)

sources_all <- c("Isolate", "Metagenomic")

result_res_plot <- result_res |>
  tidyr::complete(
    length_bin,
    category,
    source = sources_all,
    fill = list(
      mean_prop = 0,
      sd_prop = 0,
      se = 0,
      ci_low = 0,
      ci_high = 0,
      n_plasmids = 0
    )
  )

plot_04_D <- result_res_plot |>
  dplyr::filter(category %in% c("AMR","Defense","Mobilome","Unknown")) |>
  ggplot2::ggplot(ggplot2::aes(
    x = factor(category, levels = cog_order),
    y = mean_prop,
    fill = source,
    group = source
  )) +
  ggplot2::geom_col(linewidth = 0.24, position = pd, width = 0.9) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = ci_low, ymax = ci_high),
    width = 0.25, linewidth = 0.24, position = pd
  ) +
  ggplot2::scale_fill_manual(values = c("#e3e5eb", "#a7adb7")) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(),
    labels = scales::percent_format()
  ) +
  ggplot2::labs(x = "", y = "Proportion of plasmid genes") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    axis.text.x = ggplot2::element_text(
      size = 7,
      colour = "black",
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = ggplot2::element_text(size = 7, colour = "black"),
    axis.ticks.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.ticks.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    legend.key.height = grid::unit(0.25, "cm"),
    legend.key.width = grid::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "top",
    legend.title = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank()
  ) +
  ggplot2::facet_wrap(~factor(length_bin, facet_order), nrow = 1)

plot_04_D

plot_04_D |>
  ggplot2::ggsave(
    filename = "plots/fig04_D_tmp.pdf",
    width = 60.8,
    height = 63.4,
    units = "mm",
    dpi = 300
  )  

plot_04C_sup <- result_res_plot |>
  ggplot2::ggplot(ggplot2::aes(
    x = factor(category, levels = cog_order),
    y = mean_prop,
    fill = source,
    group = source
  )) +
  ggplot2::geom_col(linewidth = 0.24, position = pd, width = 0.9) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = ci_low, ymax = ci_high),
    width = 0.25, linewidth = 0.24, position = pd
  ) +
  ggplot2::scale_fill_manual(values = c("#e3e5eb", "#a7adb7")) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(),
    labels = scales::percent_format()
  ) +
  ggplot2::labs(x = "", y = "Proportion of plasmid genes") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    axis.text.x = ggplot2::element_text(
      size = 7,
      colour = "black",
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = ggplot2::element_text(size = 7, colour = "black"),
    axis.ticks.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.ticks.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    legend.key.height = grid::unit(0.25, "cm"),
    legend.key.width = grid::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "top",
    legend.title = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank()
  ) +
  ggplot2::facet_wrap(~factor(length_bin, facet_order), nrow = 1)

plot_04C_sup

plot_04C_sup |>
  ggplot2::ggsave(
    filename = "plots/figS13.pdf",
    width = 300,
    height = 63.4,
    units = "mm",
    dpi = 300
  )  

within_tests <- props_res |>
  dplyr::filter(category_expanded %in% c("AMR", "Defense")) |>
  dplyr::group_by(source, length_bin) |>
  dplyr::group_modify(\(.x, .y) {
    defense <- .x |>
      dplyr::filter(category_expanded == "Defense") |>
      dplyr::pull(prop)
    amr <- .x |>
      dplyr::filter(category_expanded == "AMR") |>
      dplyr::pull(prop)
    test <- stats::wilcox.test(defense, amr, paired = FALSE, exact = FALSE)
    effect <- effectsize::rank_biserial(defense, amr)
    tibble::tibble(
      n_defense = length(defense),
      n_amr = length(amr),
      mean_defense = mean(defense),
      mean_amr = mean(amr),
      med_defense = stats::median(defense),
      med_amr = stats::median(amr),
      rank_biserial = effect$r_rank_biserial,
      eff_ci_low = effect$CI_low,
      eff_ci_high = effect$CI_high,
      p_value = test$p.value
    )
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    p_adjusted = stats::p.adjust(
      p_value,
      method = "holm"
    )
  )

within_tests |>
  dplyr::select(!p_value)

# FIGURE 4E & 4F ---------------------------------------------------------------

imgpr_ecosystems <- imgpr_metadata |>
  dplyr::select(plasmid_seqid, dplyr::starts_with("ecosystem")) |>
  dplyr::filter(
    plasmid_seqid %in% all_reps_useful_meta$plasmid_seqid
  ) |>
  dplyr::rename(eco_1 = ecosystem,  eco_2 = ecosystem_category)

drop_eco <- combined_plasmid_feature_categories |>
  dplyr::left_join(imgpr_ecosystems, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::filter(!is.na(eco_2)) |>
  dplyr::filter(!stringr::str_detect(eco_1, "Engineered")) |>
  dplyr::mutate(eco_2 = paste0(eco_1, ": ", eco_2)) |>
  dplyr::summarise(n = dplyr::n(), .by = c(category, eco_2)) |>
  dplyr::mutate(p = n / sum(n), .by = category) |>
  dplyr::filter(p < 0.005) |>
  dplyr::distinct(eco_2) |>
  dplyr::pull()

features_by_ecosystem <- combined_plasmid_feature_categories |>
  dplyr::left_join(imgpr_ecosystems, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::filter(!is.na(eco_2)) |>
  dplyr::filter(!stringr::str_detect(eco_1, "Engineered")) |>
  dplyr::mutate(eco_2 = paste0(eco_1, ": ", eco_2)) |>
  dplyr::mutate(eco_2 = stringr::str_replace_all(eco_2, ": .*:", ":")) |>
  dplyr::summarise(n = dplyr::n(), .by = c(category, eco_2)) |>
  dplyr::mutate(p = n / sum(n), .by = category) |>
  dplyr::mutate(eco_2 = dplyr::if_else(p <= 0.01, "Other", eco_2)) |>
  dplyr::mutate(p = sum(p), .by = c(category, eco_2)) |>
  dplyr::distinct(category, eco_2, .keep_all = TRUE)

plot_04_E <- features_by_ecosystem |>
  ggplot2::ggplot(ggplot2::aes(
    x = factor(category, levels = A_x_order_features),
    y = p,
    fill = eco_2
  )) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(
    values = c("#0072B2", "#56B4E9", "#CC79A7", "#009E73", "#F0E442", "#E69F00", "#D55E00", "#C5CAD7")
  ) +
  ggupset::axis_combmatrix(sep = "-", clip = "off") +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(),
    labels = scales::percent_format()
  ) +
  ggplot2::labs(x = "", y = "Proportion of plasmids") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    axis.text.x = ggplot2::element_text(size = 7, colour = "black"),
    axis.text.y = ggplot2::element_text(size = 7, colour = "black"),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(
      linewidth = 0.24,
      fill = NA,
      colour = NA
    ),
    strip.text = ggplot2::element_text(
      size = 7,
      colour = "black",
      margin = ggplot2::margin(0.08, 0, 0.08, 0, "cm")
    ),
    legend.key.height = grid::unit(0.25, "cm"),
    legend.key.width = grid::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "top",
    legend.title = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank()
  ) +
  ggupset::theme_combmatrix(
    combmatrix.panel.point.size = 1.5,
    combmatrix.panel.line.size = 1,
    combmatrix.panel.striped_background = FALSE
  )

plot_04_E

imgpr_category_by_ecosystem <- imgpr_ecosystems |>
  dplyr::left_join(
    combined_plasmid_feature_categories,
    by = dplyr::join_by(plasmid_seqid)
  ) |>
  dplyr::mutate(
    category = dplyr::if_else(category == "NA", NA, category)
  ) |>
  dplyr::filter(
    !is.na(eco_1) & !stringr::str_detect(eco_1, "Mixed|Engineered")
  ) |>
  dplyr::mutate(
    category = dplyr::case_when(
      is.na(category) ~ "Nothing",
      .default = category
    )
  ) |>
  dplyr::mutate(eco = paste0(eco_1, ": ", eco_2)) |>
  dplyr::mutate(eco = stringr::str_replace_all(eco, ": .*:", ":")) |>
  dplyr::summarise(n = dplyr::n(), .by = c(eco, category)) |>
  dplyr::arrange(dplyr::desc(n)) |>
  dplyr::mutate(p = n / sum(n), .by = eco) |>
  dplyr::mutate(eco_n = sum(n), .by = eco) |>
  dplyr::mutate(eco_p = eco_n / sum(n)) |>
  dplyr::mutate(eco = dplyr::if_else(eco_p <= 0.01, "Other", eco)) |>
  dplyr::summarise(n = sum(n), .by = c(eco, category)) |>
  dplyr::mutate(p = n / sum(n), .by = eco) |>
  dplyr::mutate(eco_n = sum(n), .by = eco) |>
  dplyr::mutate(eco_p = eco_n / sum(n)) |>
  dplyr::mutate(eco = stringr::str_remove_all(eco, "^.*: "))

x_order_eco <- c(
  "Aquatic",
  "Terrestrial",
  "Human",
  "Mammals",
  "Insects",
  "Plants",
  "Other"
)

plot_04_F <- imgpr_category_by_ecosystem |>
  dplyr::filter(category != "") |>
  dplyr::filter(category != "Nothing") |>
  ggplot2::ggplot(ggplot2::aes(
    x = factor(eco, x_order_eco),
    y = p,
    fill = category
  )) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(
    values = c("#E69F00", "#009E73", "#0072B2")
  ) +
  ggplot2::scale_y_continuous(
    # limits = c(0, 0.25),
    expand = ggplot2::expansion(),
    labels = scales::percent_format(),
    breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25)
  ) +
  ggplot2::labs(x = "", y = "Proportion") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    axis.text.x = ggplot2::element_text(
      size = 7,
      colour = "black",
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = ggplot2::element_text(size = 7, colour = "black"),
    axis.ticks.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.ticks.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(
      linewidth = 0.24,
      fill = NA,
      colour = NA
    ),
    strip.text = ggplot2::element_text(
      size = 7,
      colour = "black",
      margin = ggplot2::margin(0.08, 0, 0.08, 0, "cm")
    ),
    legend.key.height = grid::unit(0.25, "cm"),
    legend.key.width = grid::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "top",
    legend.title = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank()
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_04_F

layout <- "
AB
"

fig_04_EF <- 
  plot_04_E +
  plot_04_F +
  patchwork::plot_layout(
    design = layout,
    guides = "collect",
    axes = "collect_x", 
    widths = c(1, 1)
  ) &
  ggplot2::theme(legend.position = 'top')

fig_04_EF

fig_04_EF |>
  ggplot2::ggsave(
    filename = "plots/fig04_EF.pdf",
    width = 127,
    height = 66.5,
    units = "mm",
    dpi = 300
  )

# TESTS RELATED TO 4E & 4F

imgpr_category_by_ecosystem |>
  dplyr::summarise(n = sum(n), .by = category) |>
  # dplyr::summarise(n = sum(n))
  purrr::pwalk(
    \(category, n) {
      cli::cli_alert_info(
        "{category}: {n} samples"
      )
    }
  )

amr_by_eco <- imgpr_category_by_ecosystem |>
  dplyr::mutate(
    amr_status = dplyr::if_else(
      category %in% c("AMR", "AMR+DS"),
      "AMR_present",
      "AMR_absent"
    ),
    human_status = dplyr::if_else(
      eco == "Human",
      "Human",
      "Non-human"
    )
  ) |>
  dplyr::group_by(human_status, amr_status) |>
  dplyr::summarise(
    n = sum(n),
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    names_from = amr_status,
    values_from = n,
    values_fill = 0
  )

amr_mat <- amr_by_eco |>
  dplyr::arrange(
    factor(
      human_status,
      levels = c("Human", "Non-human")
    )
  ) |>
  dplyr::select(
    human_status,
    AMR_present,
    AMR_absent
  ) |>
  tibble::column_to_rownames("human_status") |>
  as.matrix()

amr_mat

amr_fisher <- stats::fisher.test(amr_mat)

amr_fisher$p.value
amr_fisher$estimate
amr_fisher$conf.int



category_eco_mat <- features_by_ecosystem |>
  dplyr::filter(
    category %in% c("DS", "AMR", "AMR-DS")
  ) |>
  dplyr::select(category, eco_2, n) |>
  tidyr::pivot_wider(
    names_from = eco_2,
    values_from = n,
    values_fill = 0
  ) |>
  dplyr::arrange(
    factor(
      category,
      levels = c("DS", "AMR", "AMR-DS")
    )
  ) |>
  tibble::column_to_rownames("category") |>
  as.matrix()

pielou_evenness <- function(x) {
  p <- x / sum(x)
  p <- p[p > 0]
  
  -sum(p * log(p)) / log(length(x))
}

observed_evenness <- apply(
  category_eco_mat,
  MARGIN = 1,
  FUN = pielou_evenness
)

observed_evenness

observed_differences <- c(
  "DS_vs_AMR" =
    observed_evenness[["DS"]] -
    observed_evenness[["AMR"]],
  
  "DS_vs_AMR-DS" =
    observed_evenness[["DS"]] -
    observed_evenness[["AMR-DS"]]
)

observed_differences

set.seed(5437)

n_permutations <- 99999

permuted_tables <- stats::r2dtable(
  n = n_permutations,
  r = rowSums(category_eco_mat),
  c = colSums(category_eco_mat)
)

permuted_differences <- vapply(
  permuted_tables,
  FUN = function(x) {
    evenness <- apply(
      x,
      MARGIN = 1,
      FUN = pielou_evenness
    )
    
    c(
      "DS_vs_AMR" =
        evenness[[1]] - evenness[[2]],
      
      "DS_vs_AMR-DS" =
        evenness[[1]] - evenness[[3]]
    )
  },
  FUN.VALUE = c(
    "DS_vs_AMR" = 0,
    "DS_vs_AMR-DS" = 0
  )
)

evenness_p <- c(
  "DS_vs_AMR" = (
    1 + sum(
      permuted_differences["DS_vs_AMR", ] >=
        observed_differences[["DS_vs_AMR"]]
    )
  ) / (
    n_permutations + 1
  ),
  
  "DS_vs_AMR-DS" = (
    1 + sum(
      permuted_differences["DS_vs_AMR-DS", ] >=
        observed_differences[["DS_vs_AMR-DS"]]
    )
  ) / (
    n_permutations + 1
  )
)

evenness_p

evenness_p_adjusted <- stats::p.adjust(
  evenness_p,
  method = "holm"
)

evenness_p_adjusted

evenness_results <- tibble::tibble(
  comparison = names(observed_differences),
  difference = unname(observed_differences),
  p_value = unname(evenness_p),
  p_adjusted = unname(evenness_p_adjusted)
)

evenness_results
