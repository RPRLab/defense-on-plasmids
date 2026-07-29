library(patchwork)

# READ PLSDB DATA --------------------------------------------------------------

defense_unification <- readxl::read_xlsx(
  "data/defense_system_unification.xlsx"
)

plsdb_metadata <- readxl::read_xlsx(
  "data/plsdb_plasmid-host_metadata.xlsx",
  col_types = c(
    "text", "text", "logical", "text", "text", "logical", "numeric", "numeric", 
    "text", "text", "text", "text", "text", "text", "text", "text", "text", 
    "text", "text", "text", "text", "text", "text", "text", "numeric", 
    "numeric", "text", "text", "text", "text", "text", "text", "text", "text", 
    "text", "text", "text", "numeric", "numeric", "text", "text", "numeric", 
    "text", "numeric", "text", "numeric"
  )
)

plsdb_metadata_rep <- plsdb_metadata |>
  dplyr::filter(representative == TRUE)

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

plsdb_host_defense_raw <- readxl::read_xlsx(
  "data/plsdb_host_defense.xlsx"
)

plsdb_host_defense <- plsdb_host_defense_raw |>
# REMOVE HOK/SOK AND MAZEF
dplyr::filter(!type %in% c("Hok/Sok", "MazEF")) |>
  dplyr::mutate(replicon = "host", feature = "defense") |>
  dplyr::mutate(
    feature_id = glue::glue(
      "{host_acc}", "_defense_",
      "{stringr::str_pad(dplyr::row_number(), 2, 'left', '0')}"
    ),
    .by = host_acc
  )

# CUSTOM THEME -----------------------------------------------------------------

theme_custom <- function(grid = c("none", "x", "y")) {
  grid <- match.arg(grid)
  
  line_axis <- ggplot2::element_line(
    colour = "black",
    linewidth = 0.24,
    lineend = "round"
  )
  
  line_grid <- ggplot2::element_line(
    colour = "#EBEBEB",
    linewidth = 0.24,
    lineend = "round"
  )
  
  ggplot2::theme(
    axis.title = ggplot2::element_text(colour = "black", size = 7),
    axis.text  = ggplot2::element_text(colour = "black", size = 7),
    axis.ticks = line_axis,
    axis.line  = line_axis,
    
    legend.position = "none",
    
    panel.background = ggplot2::element_blank(),
    panel.border     = ggplot2::element_blank(),
    plot.background  = ggplot2::element_blank(),
    
    panel.grid.major.x = if (grid == "x") line_grid else ggplot2::element_blank(),
    panel.grid.minor.x = if (grid == "x") line_grid else ggplot2::element_blank(),
    panel.grid.major.y = if (grid == "y") line_grid else ggplot2::element_blank(),
    panel.grid.minor.y = if (grid == "y") line_grid else ggplot2::element_blank()
  )
}

# ANALYSIS ---------------------------------------------------------------------

# NUMBER OF REPLICONS (TOTAL)

stat_n_plasmids_total <- plsdb_metadata |>
  dplyr::distinct(plasmid_seqid) |>
  nrow()

cli::cli_alert_info("Plasmids (total): {scales::comma(stat_n_plasmids_total)}")
  
stat_n_chromosomes_total <- plsdb_metadata |>
  dplyr::distinct(host_acc) |>
  nrow()

cli::cli_alert_info("Chromosomes (total): {scales::comma(stat_n_chromosomes_total)}")

# NUMBER OF REPLICONS (REPRESENTATIVE)

stat_n_plasmids_rep <- plsdb_metadata_rep |>
  dplyr::distinct(plasmid_seqid) |>
  nrow()

cli::cli_alert_info("Plasmids (representative): {scales::comma(stat_n_plasmids_rep)}")

stat_n_chromosomes_rep <- plsdb_metadata_rep |>
  dplyr::distinct(host_acc) |>
  nrow()

cli::cli_alert_info("Chromosomes (representative): {scales::comma(stat_n_chromosomes_rep)}")
  
count_representative_replicons <- tibble::tribble(
 ~"replicon",          ~"n_replicon",
         "P",    stat_n_plasmids_rep,
         "C", stat_n_chromosomes_rep
)

x_order <- c("P", "C")

plot_01B <- count_representative_replicons |>
  ggplot2::ggplot(ggplot2::aes(
    x = factor(replicon, levels = x_order),
    y = n_replicon
  )) +
  ggplot2::geom_col(ggplot2::aes(fill = factor(replicon, levels = x_order))) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.05)),
    labels = scales::label_comma()
  ) +
  ggplot2::scale_fill_manual(values = c("#5496CE", "#c5cad7")) +
  ggplot2::labs(x = "", y = "No. sequences") +
  theme_custom(grid = "y") +
  ggplot2::coord_cartesian(clip = "off")

plot_01B

# NUMBER OF SYSTEMS IN REPLICONS

plsdb_metadata_rep |>
  dplyr::distinct(plasmid_seqid) |>
  dplyr::left_join(plsdb_plasmid_defense_raw, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::filter(!is.na(type)) |>
  dplyr::mutate(category = dplyr::if_else(!type %in% c("MazEF", "Hok/Sok"), "other", type)) |>
  dplyr::summarise(n = dplyr::n(), .by = category) |>
  dplyr::mutate(p = round(n / sum(n) * 100))

defense_per_plasmid <- plsdb_metadata_rep |>
  dplyr::distinct(plasmid_seqid) |>
  dplyr::left_join(plsdb_plasmid_defense, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::mutate(has_defense = dplyr::if_else(is.na(type), FALSE, TRUE)) |>
  dplyr::summarise(n_systems = dplyr::n(), .by = c(plasmid_seqid, has_defense)) |>
  dplyr::mutate(n_systems = dplyr::if_else(has_defense, n_systems, 0)) |>
  dplyr::mutate(replicon = "P") |>
  dplyr::select(-c(plasmid_seqid, has_defense))

defense_per_chromosome <- plsdb_metadata_rep |>
  dplyr::distinct(host_acc) |>
  dplyr::left_join(plsdb_host_defense, by = dplyr::join_by(host_acc)) |>
  dplyr::mutate(has_defense = dplyr::if_else(is.na(type), FALSE, TRUE)) |>
  dplyr::summarise(n_systems = dplyr::n(), .by = c(host_acc, has_defense)) |>
  dplyr::mutate(n_systems = dplyr::if_else(has_defense, n_systems, 0)) |>
  dplyr::mutate(replicon = "C") |>
  dplyr::select(-c(host_acc, has_defense))

defense_per_replicon <-
  dplyr::bind_rows(
    defense_per_plasmid,
    defense_per_chromosome
  )

defense_in_replicons <- defense_per_replicon |>
  dplyr::summarise(n_systems = sum(n_systems), .by = replicon)

stat_n_systems_in_plasmids <- defense_in_replicons |>
  dplyr::filter(replicon == "P") |>
  dplyr::pull(n_systems)
  
cli::cli_alert_info("Defense systems in plasmids: {scales::comma(stat_n_systems_in_plasmids)}")

stat_n_systems_in_chromsomes <- defense_in_replicons |>
  dplyr::filter(replicon == "C") |>
  dplyr::pull(n_systems)

cli::cli_alert_info("Defense systems in chromosomes: {scales::comma(stat_n_systems_in_chromsomes)}")

x_order <- c("P", "C")

plot_01C <- defense_in_replicons |>
  ggplot2::ggplot(ggplot2::aes(
    x = factor(replicon, levels = x_order),
    y = n_systems
  )) +
  ggplot2::geom_col(ggplot2::aes(fill = factor(replicon, levels = x_order))) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.01)),
    labels = scales::label_comma(),
    limits = c(0, 125000)
  ) +
  ggplot2::scale_fill_manual(values = c("#5496CE", "#c5cad7")) +
  ggplot2::labs(x = "", y = "No. defense systems") +
  theme_custom(grid = "y") +
  ggplot2::coord_cartesian(clip = "off")

plot_01C

# NUMBER OF SYSTEM TYPES IN REPLICONS

stat_n_system_types <- defense_unification |>
  dplyr::filter(keep == TRUE) |>
  dplyr::distinct(`unified type`) |>
  nrow()

cli::cli_alert_info("Total detectable system types: {stat_n_system_types}")

stat_n_system_types_in_plasmids <- plsdb_metadata_rep |>
  dplyr::distinct(plasmid_seqid) |>
  dplyr::left_join(plsdb_plasmid_defense, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::filter(!is.na(type)) |>
  dplyr::distinct(type) |>
  nrow()

cli::cli_alert_info("System types in plasmids: {stat_n_system_types_in_plasmids}")

cli::cli_alert_info("Percent of types in plasmids: {round(stat_n_system_types_in_plasmids / stat_n_system_types * 100)}%")

stat_n_system_types_in_chromosomes <- plsdb_metadata_rep |>
  dplyr::distinct(host_acc) |>
  dplyr::left_join(plsdb_host_defense, by = dplyr::join_by(host_acc)) |>
  dplyr::filter(!is.na(type)) |>
  dplyr::distinct(type) |>
  nrow()

cli::cli_alert_info("System types in chromosomes: {stat_n_system_types_in_chromosomes}")

cli::cli_alert_info("Percent of types in chromosomes: {round(stat_n_system_types_in_chromosomes / stat_n_system_types * 100)}%")

count_replicon_defense_types <- tibble::tribble(
  ~"replicon",                         ~"n_types",
          "P",    stat_n_system_types_in_plasmids,
          "C", stat_n_system_types_in_chromosomes
)

x_order <- c("P", "C")

plot_01D <- count_replicon_defense_types |>
  ggplot2::ggplot(ggplot2::aes(x = factor(replicon, levels = x_order), y = n_types)) +
  ggplot2::geom_col(
    ggplot2::aes(fill = factor(replicon, levels = x_order)),
    position = "dodge"
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.05)),
    labels = scales::label_comma()
  ) +
  ggplot2::scale_fill_manual(values = c("#5496CE", "#c5cad7")) +
  ggplot2::labs(x = "", y = "No. system types") +
  theme_custom(grid = "y") +
  ggplot2::coord_cartesian(clip = "off")

plot_01D

# NUMBER OF REPLICONS WITH AT LEAST ONE SYSTEM

stat_n_plasmids_with_defense <- plsdb_metadata_rep |>
  dplyr::distinct(plasmid_seqid) |>
  dplyr::left_join(plsdb_plasmid_defense, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::mutate(has_defense = dplyr::if_else(is.na(type), FALSE, TRUE)) |>
  dplyr::distinct(plasmid_seqid, has_defense) |>
  dplyr::filter(has_defense == TRUE) |>
  nrow()

cli::cli_alert_info(
  "Plasmids with ≥1 system: {scales::comma(stat_n_plasmids_with_defense)}"
)

stat_p_plasmids_with_defense <- stat_n_plasmids_with_defense / stat_n_plasmids_rep

cli::cli_alert_info(
  "Plasmids with ≥1 system (%): {round(stat_p_plasmids_with_defense * 100)}%"
)

stat_n_chromosomes_with_defense <- plsdb_metadata_rep |>
  dplyr::distinct(host_acc) |>
  dplyr::left_join(plsdb_host_defense, by = dplyr::join_by(host_acc)) |>
  dplyr::mutate(has_defense = dplyr::if_else(is.na(type), FALSE, TRUE)) |>
  dplyr::distinct(host_acc, has_defense) |>
  dplyr::filter(has_defense == TRUE) |>
  nrow()

cli::cli_alert_info(
  "Chromosomes with ≥1 system: {scales::comma(stat_n_chromosomes_with_defense)}"
)

stat_p_chromosomes_with_defense <- stat_n_chromosomes_with_defense / stat_n_chromosomes_rep

cli::cli_alert_info(
  "Chromosomes with ≥1 system (%): {scales::comma(stat_p_chromosomes_with_defense *100)}%"
)

proportion_replicons_with_defense <- tibble::tribble(
  ~"replicon",                  ~"p_w_defense",
          "P",    stat_p_plasmids_with_defense,
          "C", stat_p_chromosomes_with_defense
)

x_order <- c("P", "C")

plot_01E <- proportion_replicons_with_defense |>
  ggplot2::ggplot(ggplot2::aes(
    x = factor(replicon, levels = x_order),
    y = p_w_defense
  )) +
  ggplot2::geom_col(ggplot2::aes(fill = factor(replicon, levels = x_order))) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.05)),
    labels = scales::label_percent()
  ) +
  ggplot2::scale_fill_manual(values = c("#5496CE", "#c5cad7")) +
  ggplot2::labs(x = "", y = "Proportion with systems") +
  theme_custom(grid = "y") +
  ggplot2::coord_cartesian(clip = "off")

plot_01E

layout <- "ABCD"

plot_01BCDE <- plot_01B + plot_01C + plot_01D + plot_01E +
  patchwork::plot_layout(design = layout)

plot_01BCDE

plot_01BCDE |>
  ggplot2::ggsave(
    filename = "plots/fig01_BCDE.pdf",
    width = 110,
    height = 35,
    units = "mm",
    dpi = 300
  )

# FIGURE 1D - TOTAL DEFENSE IN REPLICON

defense_per_replicon_filt <- defense_per_replicon |>
  dplyr::filter(dplyr::between(n_systems, 1, 100))

defense_per_replicon_stats <- defense_per_replicon_filt |>
  dplyr::summarise(
    avg = round(mean(n_systems), 2),
    med = median(n_systems),
    q1 = quantile(n_systems, 0.25),
    q3 = quantile(n_systems, 0.75),
    .by = replicon
  )

defense_per_replicon_stats |>
  purrr::pwalk(
    \(replicon, avg, med, q1, q3) {
      label <- dplyr::recode_values(replicon, "P" ~ "Plasmid", "C" ~ "Chromosome")
      cli::cli_alert_info(
        "{label}: mean = {round(avg)}, median = {med}, IQR = {q1}-{q3}"
      )
    }
  )

stat_n_plasmids_w_5_sys <- defense_per_replicon_filt |>
  dplyr::filter(replicon == "P") |>
  dplyr::filter(n_systems >= 5) |>
  nrow()

cli::cli_alert_info("Plasmids with ≥5 systems: {stat_n_plasmids_w_5_sys}")

cli::cli_alert_info("Proportion plasmids with 5+ systems: {round(stat_n_plasmids_w_5_sys / stat_n_plasmids_rep * 100)}%")

x_order <- c("P", "C")

box_width <- 0.75

plot_01F <- defense_per_replicon_filt |>
  ggplot2::ggplot(ggplot2::aes(
    x = factor(replicon, levels = x_order),
    y = n_systems
  )) +
  ggrastr::rasterise(
    dpi = 600,
    gghalves::geom_half_point(
      ggplot2::aes(colour = replicon),
      alpha = 0.5,
      size = 0.2,
      shape = 16,
      transformation = ggplot2::position_jitter(height = 0),
      range_scale = 0.6
    )
  ) +
  gghalves::geom_half_boxplot(
    width = box_width,
    outlier.colour = NA,
    linewidth = 0.24,
    errorbar.length = 0,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = defense_per_replicon_stats,
    ggplot2::aes(y = avg),
    position = ggplot2::position_nudge(x = -box_width / 4),
    colour = "#c5373d",
    size = 1,
    shape = 16
  ) +
  ggplot2::scale_colour_manual(values = c("#c5cad7", "#5496CE")) +
  ggplot2::scale_y_log10(
    expand = ggplot2::expansion(mult = c(0.05, 0)),
    breaks = function(limits) {
      scales::breaks_log(base = 10)(limits) |>
        (\(x) x[x >= 1])()
    },
    labels = scales::label_log(base = 10),
    limits = function(l) c(l[1], 10^ceiling(log10(l[2]))),
    minor_breaks = function(limits) {
      scales::minor_breaks_log()(limits) |>
        (\(x) x[x >= 1])()
    }
  ) +
  ggplot2::labs(x = "", y = "No. systems per sequence (log scale)") +
  theme_custom(grid = "y") +
  ggplot2::theme(
    axis.title.y = ggtext::element_markdown()
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_01F

plasmid_defense_per_length <- plsdb_metadata_rep |>
  dplyr::distinct(plasmid_seqid, plasmid_length) |>
  dplyr::left_join(plsdb_plasmid_defense, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::filter(!is.na(type)) |>
  dplyr::summarise(n_systems = dplyr::n(), .by = c(plasmid_seqid, plasmid_length)) |>
  dplyr::mutate(
    sys_per_mbp = n_systems / (plasmid_length / 1000000),
    sys_per_kbp = n_systems / (plasmid_length / 1000)
  ) |>
  dplyr::select(-c(plasmid_seqid)) |>
  dplyr::rename(length = plasmid_length) |>
  dplyr::mutate(replicon = "P")

chromosome_defense_per_length <- plsdb_metadata_rep |>
  dplyr::distinct(host_acc, host_length) |>
  dplyr::slice_head(n = 1, by = host_acc) |>
  dplyr::left_join(plsdb_host_defense, by = dplyr::join_by(host_acc)) |>
  dplyr::filter(!is.na(type)) |>
  dplyr::summarise(n_systems = dplyr::n(), .by = c(host_acc, host_length)) |>
  dplyr::mutate(
    sys_per_mbp = n_systems / (host_length / 1000000),
    sys_per_kbp = n_systems / (host_length / 1000)
  ) |>
  dplyr::select(-c(host_acc)) |>
  dplyr::rename(length = host_length) |>
  dplyr::mutate(replicon = "C")

replicon_defense_per_length <-
  dplyr::bind_rows(
    plasmid_defense_per_length,
    chromosome_defense_per_length
  ) |>
  dplyr::mutate(
    log10_sys_per_mbp = log10(sys_per_mbp),
    log10_sys_per_kbp = log10(sys_per_kbp),
    log10_n_sys = log10(n_systems)
  )

defense_per_length_stats <- replicon_defense_per_length |>
  dplyr::summarise(
    avg = round(mean(sys_per_mbp), 2),
    med = round(median(sys_per_mbp), 2),
    .by = replicon
  )

defense_per_length_stats |>
  purrr::pwalk(
    \(replicon, avg, med) {
      label <- dplyr::recode_values(replicon, "P" ~ "Plasmid", "C" ~ "Chromosome")
      cli::cli_alert_info(
        "{label}: mean = {avg}, median = {med}"
      )
    }
  )

box_width <- 0.75

plot_01G <- replicon_defense_per_length |>
  ggplot2::ggplot(ggplot2::aes(
    x = factor(replicon, levels = x_order),
    y = sys_per_mbp
  )) +
  ggrastr::rasterise(
    dpi = 600,
    gghalves::geom_half_point(
      ggplot2::aes(colour = replicon),
      alpha = 0.5,
      size = 0.2,
      shape = 16,
      transformation = ggplot2::position_jitter(height = 0),
      range_scale = 0.6
    )
  ) +
  gghalves::geom_half_boxplot(
    width = box_width,
    outlier.colour = NA,
    linewidth = 0.24,
    errorbar.length = 0,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = defense_per_length_stats,
    ggplot2::aes(y = avg),
    position = ggplot2::position_nudge(x = -box_width / 4),
    colour = "#c5373d",
    size = 1,
    shape = 16
  ) +
  ggplot2::scale_colour_manual(values = c("#c5cad7", "#5496CE")) +
  ggplot2::scale_y_log10(
    expand = ggplot2::expansion(mult = c(0.05, 0)),   # lower only
    breaks = scales::log_breaks(base = 10),
    labels = function(x) scales::label_math(10^.x)(log10(x)),
    limits = function(l) c(l[1], 10^ceiling(log10(l[2]))),  # upper -> next major tick
    minor_breaks = function(limits) {
      limits <- as.numeric(limits)
      limits <- limits[is.finite(limits) & limits > 0]
      if (length(limits) < 2) return(NULL)
      lo <- floor(log10(min(limits)))
      hi <- ceiling(log10(max(limits)))
      as.vector(outer(2:9, 10^(lo:hi)))
    }
  ) +
  ggplot2::labs(x = "", y = "No. systems per Mb (log scale)") +
  theme_custom(grid = "y") +
  ggplot2::theme(
    axis.title.y = ggtext::element_markdown()
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_01G

plasmid_defense_encoding_proportion <- plsdb_metadata_rep |>
  dplyr::distinct(plasmid_seqid, plasmid_length) |>
  dplyr::left_join(plsdb_plasmid_defense, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::filter(!is.na(type)) |>
  dplyr::mutate(defense_length = end - start) |>
  dplyr::summarise(
    defense_length = sum(defense_length),
    .by = c(plasmid_seqid, plasmid_length)
  ) |>
  dplyr::mutate(proportion_defense = defense_length / plasmid_length) |>
  dplyr::mutate(replicon = "P")

defense_proportion_stats <- plasmid_defense_encoding_proportion |>
  dplyr::summarise(
    avg = round(mean(proportion_defense) * 100, 2),
    med = round(median(proportion_defense) * 100, 2),
    .by = replicon
  )

cli::cli_alert_info(
  "Plasmid defense-encoding proportion: mean = {defense_proportion_stats$avg}%, median = {defense_proportion_stats$med}%"
)

stat_p_plasmid_gt50 <- plasmid_defense_encoding_proportion |>
  dplyr::filter(proportion_defense >= 0.5) |>
  nrow()

cli::cli_alert_info("Plasmids with ≥50% encoding systems: {stat_p_plasmid_gt50}")

custom_percent <- function(accuracy_below_1 = 0.1, accuracy_above_1 = 1) {
  function(x) {
    accuracy <- ifelse(abs(x) < 0.01, accuracy_below_1, accuracy_above_1)
    scales::percent(x, accuracy = accuracy)
  }
}

box_width <- 0.75

plot_01H <- plasmid_defense_encoding_proportion |>
  ggplot2::ggplot(ggplot2::aes(
    x = factor(replicon, levels = x_order),
    y = proportion_defense
  )) +
  ggrastr::rasterise(
    dpi = 600,
    gghalves::geom_half_point(
      ggplot2::aes(colour = replicon),
      alpha = 0.5,
      size = 0.2,
      shape = 16,
      transformation = ggplot2::position_jitter(height = 0),
      range_scale = 0.6
    )
  ) +
  gghalves::geom_half_boxplot(
    width = box_width,
    outlier.colour = NA,
    linewidth = 0.24,
    errorbar.length = 0,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = defense_proportion_stats,
    ggplot2::aes(y = avg / 100),
    position = ggplot2::position_nudge(x = -box_width / 4),
    colour = "#c5373d",
    size = 1,
    shape = 16
  ) +
  ggplot2::geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    colour = "#c5373d",
    linewidth = 0.24
  ) +
  ggplot2::scale_colour_manual(values = c("#5496CE")) +
  ggplot2::scale_y_log10(
    expand = ggplot2::expansion(mult = c(0.05, 0)),
    breaks = scales::log_breaks(base = 10),
    labels = scales::label_percent(accuracy = 1),
    limits = function(l) c(l[1], 10^ceiling(log10(l[2]))),
    minor_breaks = function(limits) {
      limits <- as.numeric(limits)
      limits <- limits[is.finite(limits) & limits > 0]
      if (length(limits) < 2) return(NULL)
      lo <- floor(log10(min(limits)))
      hi <- ceiling(log10(max(limits)))
      as.vector(outer(2:9, 10^(lo:hi)))
    }
  ) +
  ggplot2::labs(x = "", y = "Defense-encoding proportion (log scale)") +
  theme_custom(grid = "y") +
  ggplot2::theme(
    axis.title.y = ggtext::element_markdown()
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_01H

layout <- "ABC"

plot_01FGH <- plot_01F + plot_01G + plot_01H +
  patchwork::plot_layout(design = layout, widths = c(2, 2, 1.1))

plot_01FGH

plot_01FGH |>
  ggplot2::ggsave(
    filename = "plots/fig01_FGH.pdf",
    width = 80,
    height = 55,
    units = "mm",
    dpi = 300
  )

# FIGURE 1B — HEATMAP OF SYSTEM PREVALENCE BY TAXA -----------------------------

# READ GTDB METADATA

bac120_meta <- readr::read_tsv("data/gtdb/bac120_metadata_r220.tsv")
ar53_meta <- readr::read_tsv("data/gtdb/ar53_metadata_r220.tsv")
gtdb_meta <- dplyr::bind_rows(bac120_meta, ar53_meta)

gtdb_colnames <- purrr::map_chr(
  c("d", "p", "c", "o", "f", "g", "s"),
  ~ paste0("gtdb_", .)
)

gtdb_taxonomy <- gtdb_meta |>
  dplyr::distinct(accession, gtdb_taxonomy) |>
  tidyr::separate_wider_delim(gtdb_taxonomy, ";", names = gtdb_colnames) |>
  dplyr::mutate(id = stringr::str_extract(accession, "[0-9]{9}"))

# RESOLVE TAXONOMY ASSIGNMENT

taxonomy_gtdb <- plsdb_metadata |>
  dplyr::select(host_acc, dplyr::starts_with("taxonomy_")) |>
  dplyr::select(!taxonomy_strain) |>
  dplyr::mutate(id = stringr::str_extract(host_acc, "[0-9]{9}")) |>
  dplyr::left_join(gtdb_taxonomy, by = dplyr::join_by(id)) |>
  dplyr::distinct()

missing_taxa_summary <- taxonomy_gtdb |>
  dplyr::filter(is.na(accession)) |>
  dplyr::summarise(n = dplyr::n(), .by = taxonomy_phylum) |>
  dplyr::arrange(dplyr::desc(n))

missing_taxa_summary

missing_taxa <- missing_taxa_summary |>
  dplyr::pull(taxonomy_phylum)

taxa_majority <- taxonomy_gtdb |>
  dplyr::filter(taxonomy_phylum %in% missing_taxa) |>
  dplyr::summarise(n = dplyr::n(), .by = c(taxonomy_phylum, gtdb_p)) |>
  dplyr::arrange(taxonomy_phylum, dplyr::desc(n)) |>
  dplyr::slice_max(n, by = taxonomy_phylum, with_ties = FALSE)

taxa_majority

taxonomy_missing <- taxonomy_gtdb |>
  dplyr::filter(is.na(accession)) |>
  dplyr::select(!c(dplyr::starts_with("gtdb_"), accession)) |>
  dplyr::left_join(taxa_majority, by = dplyr::join_by(taxonomy_phylum)) |>
  dplyr::select(-n) |>
  dplyr::rename(resolved_p = gtdb_p)

taxonomy_resolved <- taxonomy_gtdb |>
  dplyr::filter(!is.na(accession)) |>
  dplyr::bind_rows(taxonomy_missing) |>
  dplyr::arrange(host_acc) |>
  dplyr::mutate(
    resolved_p = dplyr::if_else(is.na(resolved_p), gtdb_p, resolved_p),
    resolved_d = paste0("d__", taxonomy_superkingdom)
  ) |>
  dplyr::distinct(host_acc, resolved_d, resolved_p)

plsdb_metadata_taxa_out <- plsdb_metadata |>
  dplyr::left_join(taxonomy_resolved, by = dplyr::join_by(host_acc)) |>
  dplyr::rename(gtdb_domain = resolved_d, gtdb_phylum = resolved_p)

plsdb_metadata_taxa_out |> 
  writexl::write_xlsx(
    "data/plsdb_plasmid-host_metadata_with_taxonomy.xlsx"
  )

# SUMMARISE SYSTEM COUNT BY TAXA

system_count_by_taxa <- plsdb_metadata_rep |>
  dplyr::select(plasmid_seqid, host_acc) |>
  dplyr::left_join(taxonomy_resolved, by = dplyr::join_by(host_acc)) |>
  dplyr::mutate(n_plasmid_in_phylum = dplyr::n(), .by = c(resolved_p)) |>
  dplyr::left_join(
    plsdb_plasmid_defense,
    by = dplyr::join_by(plasmid_seqid)
  ) |>
  dplyr::filter(type != "NA") |>
  dplyr::distinct(
    plasmid_seqid,
    resolved_d,
    resolved_p,
    n_plasmid_in_phylum,
    type
  ) |>
  dplyr::summarise(
    n_system_in_phylum = dplyr::n(),
    .by = c(resolved_d, resolved_p, n_plasmid_in_phylum, type)
  ) |>
  dplyr::mutate(
    prevalence = n_system_in_phylum / n_plasmid_in_phylum
  ) |>
  dplyr::mutate(n_system_in_total = sum(n_system_in_phylum), .by = type) |>
  tidyr::complete(resolved_p, type, fill = list(prevalence = NA)) |>
  dplyr::mutate(
    resolved_d = max(resolved_d, na.rm = TRUE),
    n_plasmid_in_phylum = max(n_plasmid_in_phylum, na.rm = TRUE),
    n_system_in_phylum = max(n_system_in_phylum, na.rm = TRUE),
    .by = resolved_p
  ) |>
  dplyr::mutate(
    n_system_in_total = max(n_system_in_total, na.rm = TRUE),
    .by = type
  )

# FILTER

system_prevalence_by_taxa_filt <- system_count_by_taxa |>
  dplyr::filter(n_system_in_total >= 100) |>
  dplyr::filter(n_plasmid_in_phylum >= 20) |>
  dplyr::mutate(prevalence = dplyr::if_else(prevalence > 0.5, 0.5, prevalence))

# PREPARE GTDB TREES

bac120_p_tree <- treeio::read.tree("data/gtdb/bac120_r220_phylum.tree")

bac120_p_tip_labels <- bac120_p_tree[["tip.label"]]

tip_labels_keep_bac120 <- system_prevalence_by_taxa_filt |>
  dplyr::filter(resolved_d == "d__Bacteria") |>
  dplyr::distinct(resolved_p) |>
  dplyr::pull()

tip_labels_drop_bac120 <- bac120_p_tip_labels |>
  setdiff(tip_labels_keep_bac120)

bac120_p_tree_reduced <- bac120_p_tree |>
  tidytree::drop.tip(tip_labels_drop_bac120)

bac120_p_tree_reduced_midpoint <- bac120_p_tree_reduced |>
  phytools::midpoint_root()

bac120_p_tree_clean <- bac120_p_tree_reduced_midpoint |>
  dplyr::as_tibble() |>
  dplyr::mutate(label = stringr::str_remove(label, "p__")) |>
  tidytree::as.phylo()

plot_bac_tree <- bac120_p_tree_clean |>
  ggtree::ggtree(
    linewidth = 0.24,
    color = "black"
  ) +
  ggtree::geom_treescale(
    fontsize = 7,
    linesize = 0.24,
    color = "black",
    width = 0.2,
    y = 2
  ) +
  ggplot2::theme(
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")
  )

plot_bac_tree

ar53_p_tree <- treeio::read.tree("data/gtdb/ar53_r220_phylum.tree")

ar53_p_tip_labels <- ar53_p_tree[["tip.label"]]

tip_labels_keep_ar53 <- system_prevalence_by_taxa_filt |>
  dplyr::filter(resolved_d == "d__Archaea") |>
  dplyr::distinct(resolved_p) |>
  dplyr::pull()

tip_labels_drop_ar53 <- ar53_p_tip_labels |>
  setdiff(tip_labels_keep_ar53)

ar53_p_tree_reduced <- ar53_p_tree |>
  tidytree::drop.tip(tip_labels_drop_ar53)

if (length(tip_labels_keep_ar53) == 1) {
  ar53_p_tree_reduced <- ar53_p_tree_reduced |>
    phytools::bind.tip(
      tip.label = "_dummy",
      where = ape::Ntip(ar53_p_tree_reduced) + 1,
      edge.length = 1
    )
}

ar53_p_tree_reduced_midpoint <- ar53_p_tree_reduced |>
  phytools::midpoint_root()

ar53_p_tree_clean <- ar53_p_tree_reduced_midpoint |>
  dplyr::as_tibble() |>
  dplyr::mutate(label = stringr::str_remove(label, "p__")) |>
  tidytree::as.phylo()

plot_ar_tree <- ar53_p_tree_clean |>
  ggtree::ggtree(
    linewidth = 0.24,
    color = "black"
  ) +
  ggtree::geom_treescale(
    fontsize = 7,
    linesize = 0.24,
    color = "black",
    width = 0.2,
    y = 2
  ) +
  ggplot2::theme(
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")
  )

plot_ar_tree

y_order <- dplyr::bind_rows(
  plot_bac_tree$data |> dplyr::mutate(order = 2),
  plot_ar_tree$data |> dplyr::mutate(order = 1)
) |>
  dplyr::as_tibble() |>
  dplyr::filter(isTip) |>
  dplyr::arrange(order, y) |>
  dplyr::pull(label)

x_order <- system_prevalence_by_taxa_filt |>
  dplyr::mutate(resolved_p = stringr::str_remove(resolved_p, "p__")) |>
  dplyr::arrange(dplyr::desc(n_system_in_total)) |>
  dplyr::distinct(type) |>
  dplyr::pull()

plot_heatmap <- system_prevalence_by_taxa_filt |>
  dplyr::mutate(resolved_p = stringr::str_remove(resolved_p, "p__")) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(type, x_order),
      y = factor(resolved_p, y_order),
      fill = prevalence
    )
  ) +
  ggplot2::geom_tile(colour = "#DEDEDE", linewidth = 0.24) +
  ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::scale_fill_gradient2(
    low = "#c5e5fb",
    mid = "#026eae",
    high = "black",
    na.value = "white",
    limits = c(0, 0.3),
    midpoint = 0.15,
    breaks = c(0, 0.1, 0.2, 0.3),
    labels = scales::label_percent()
  ) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(ticks.colour = c("black", NA))
  ) +
  ggplot2::labs(fill = "") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      size = 5,
      colour = "black",
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = ggplot2::element_text(size = 5, colour = "black"),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(linewidth = 0.24),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    legend.key.height = ggplot2::unit(0.1, "cm"),
    legend.key.width = ggplot2::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "bottom",
    legend.title.position = "left",
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(2, 2, 0, 0)
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_heatmap

count_sys <- system_prevalence_by_taxa_filt |>
  dplyr::filter(!is.na(n_system_in_total)) |>
  dplyr::distinct(type, n_system_in_total)

plot_count_sys <- count_sys |>
  ggplot2::ggplot(ggplot2::aes(
    x = factor(type, x_order),
    y = n_system_in_total
  )) +
  ggplot2::geom_col(fill = "#c5cad7", linewidth = 0.24, width = 0.7) +
  ggplot2::scale_y_log10(
    name = "n",
    expand = ggplot2::expansion(mult = c(0, 0.1)),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    rect = ggplot2::element_rect(linewidth = 0.24),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 5, colour = "black"),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.border = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(
      color = "black",
      linewidth = 0.24,
      lineend = "round"
    ),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "#ebebeb"
    ),
    panel.grid.minor.y = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "#ebebeb"
    ),
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_count_sys

count_tax <- system_prevalence_by_taxa_filt |>
  dplyr::distinct(resolved_p, n_plasmid_in_phylum)

plot_count_tax <- count_tax |>
  dplyr::mutate(resolved_p = stringr::str_remove(resolved_p, "p__")) |>
  dplyr::filter(resolved_p %in% y_order) |>
  ggplot2::ggplot(
    ggplot2::aes(y = factor(resolved_p, y_order), x = n_plasmid_in_phylum)
  ) +
  ggplot2::geom_col(fill = "#c5cad7", linewidth = 0.24, width = 0.7) +
  ggplot2::scale_x_log10(
    name = "n",
    expand = ggplot2::expansion(mult = c(0, 0.1)),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    rect = ggplot2::element_rect(linewidth = 0.24),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      size = 5,
      colour = "black",
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    axis.ticks.y = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.border = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(
      color = "black",
      linewidth = 0.24,
      lineend = "round"
    ),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "#ebebeb"
    ),
    panel.grid.minor.x = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "#ebebeb"
    ),
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_count_tax

layout <- "
#DDDDDDDDDDDDDDDDDDDDDD##
#DDDDDDDDDDDDDDDDDDDDDD##
ACCCCCCCCCCCCCCCCCCCCCCEE
ACCCCCCCCCCCCCCCCCCCCCCEE
ACCCCCCCCCCCCCCCCCCCCCCEE
ACCCCCCCCCCCCCCCCCCCCCCEE
BCCCCCCCCCCCCCCCCCCCCCCEE
BCCCCCCCCCCCCCCCCCCCCCCEE
BCCCCCCCCCCCCCCCCCCCCCCEE
BCCCCCCCCCCCCCCCCCCCCCCEE
"

plot_01I <-
  plot_bac_tree +
  plot_ar_tree +
  plot_heatmap +
  plot_count_sys +
  plot_count_tax +
  patchwork::plot_layout(design = layout)

plot_01I

plot_01I |>
  ggplot2::ggsave(
    filename = "plots/fig01_I.pdf",
    width = 100,
    height = 67,
    units = "mm",
    dpi = 300
  )

# FIGURE S02 — UNFILTERED HEATMAP ----------------------------------------------

# FILTER

system_prevalence_by_taxa_filt <- system_count_by_taxa |>
  dplyr::mutate(prevalence = dplyr::if_else(prevalence > 0.5, 0.5, prevalence))

# PREPARE GTDB TREES

bac120_p_tree <- treeio::read.tree("data/gtdb/bac120_r220_phylum.tree")

bac120_p_tip_labels <- bac120_p_tree[["tip.label"]]

tip_labels_keep_bac120 <- system_prevalence_by_taxa_filt |>
  dplyr::filter(resolved_d == "d__Bacteria") |>
  dplyr::distinct(resolved_p) |>
  dplyr::pull()

tip_labels_drop_bac120 <- bac120_p_tip_labels |>
  setdiff(tip_labels_keep_bac120)

bac120_p_tree_reduced <- bac120_p_tree |>
  tidytree::drop.tip(tip_labels_drop_bac120)

bac120_p_tree_reduced_midpoint <- bac120_p_tree_reduced |>
  phytools::midpoint_root()

bac120_p_tree_clean <- bac120_p_tree_reduced_midpoint |>
  dplyr::as_tibble() |>
  dplyr::mutate(label = stringr::str_remove(label, "p__")) |>
  tidytree::as.phylo()

plot_bac_tree <- bac120_p_tree_clean |>
  ggtree::ggtree(
    linewidth = 0.24,
    color = "black"
  ) +
  ggtree::geom_treescale(
    fontsize = 7,
    linesize = 0.24,
    color = "black",
    width = 0.2,
    y = 2
  ) +
  ggplot2::theme(
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")
  )

plot_bac_tree

ar53_p_tree <- treeio::read.tree("data/gtdb/ar53_r220_phylum.tree")

ar53_p_tip_labels <- ar53_p_tree[["tip.label"]]

tip_labels_keep_ar53 <- system_prevalence_by_taxa_filt |>
  dplyr::filter(resolved_d == "d__Archaea") |>
  dplyr::distinct(resolved_p) |>
  dplyr::pull()

tip_labels_drop_ar53 <- ar53_p_tip_labels |>
  setdiff(tip_labels_keep_ar53)

ar53_p_tree_reduced <- ar53_p_tree |>
  tidytree::drop.tip(tip_labels_drop_ar53)

ar53_p_tree_reduced_midpoint <- ar53_p_tree_reduced |>
  phytools::midpoint_root()

ar53_p_tree_clean <- ar53_p_tree_reduced_midpoint |>
  dplyr::as_tibble() |>
  dplyr::mutate(label = stringr::str_remove(label, "p__")) |>
  tidytree::as.phylo()

plot_ar_tree <- ar53_p_tree_clean |>
  ggtree::ggtree(
    linewidth = 0.24,
    color = "black"
  ) +
  ggtree::geom_treescale(
    fontsize = 7,
    linesize = 0.24,
    color = "black",
    width = 0.2,
    y = 2
  ) +
  ggplot2::theme(
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")
  )

plot_ar_tree

y_order <- dplyr::bind_rows(
  plot_bac_tree$data |> dplyr::mutate(order = 2),
  plot_ar_tree$data |> dplyr::mutate(order = 1)
) |>
  dplyr::as_tibble() |>
  dplyr::filter(isTip) |>
  dplyr::arrange(order, y) |>
  dplyr::pull(label)

x_order <- system_prevalence_by_taxa_filt |>
  dplyr::mutate(resolved_p = stringr::str_remove(resolved_p, "p__")) |>
  dplyr::arrange(dplyr::desc(n_system_in_total)) |>
  dplyr::distinct(type) |>
  dplyr::pull()

plot_heatmap <- system_prevalence_by_taxa_filt |>
  dplyr::mutate(resolved_p = stringr::str_remove(resolved_p, "p__")) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(type, x_order),
      y = factor(resolved_p, y_order),
      fill = prevalence
    )
  ) +
  ggplot2::geom_tile(colour = "#DEDEDE", linewidth = 0.24) +
  ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::scale_fill_gradient2(
    low = "#c5e5fb",
    mid = "#026eae",
    high = "black",
    na.value = "white",
    limits = c(0, 0.5),
    midpoint = 0.25,
    breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
    labels = scales::label_percent()
  ) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(ticks.colour = c("black", NA))
  ) +
  ggplot2::labs(fill = "") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      size = 5,
      colour = "black",
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = ggplot2::element_text(size = 5, colour = "black"),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(linewidth = 0.24),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    legend.key.height = ggplot2::unit(0.1, "cm"),
    legend.key.width = ggplot2::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "bottom",
    legend.title.position = "left",
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(2, 2, 0, 0)
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_heatmap

count_sys <- system_prevalence_by_taxa_filt |>
  dplyr::filter(!is.na(n_system_in_total)) |>
  dplyr::distinct(type, n_system_in_total)

plot_count_sys <- count_sys |>
  ggplot2::ggplot(ggplot2::aes(
    x = factor(type, x_order),
    y = n_system_in_total
  )) +
  ggplot2::geom_col(fill = "#c5cad7", linewidth = 0.24, width = 0.7) +
  ggplot2::scale_y_log10(
    name = "n",
    expand = ggplot2::expansion(mult = c(0, 0.1)),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    rect = ggplot2::element_rect(linewidth = 0.24),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 5, colour = "black"),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.border = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(
      color = "black",
      linewidth = 0.24,
      lineend = "round"
    ),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "#ebebeb"
    ),
    panel.grid.minor.y = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "#ebebeb"
    ),
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_count_sys

count_tax <- system_prevalence_by_taxa_filt |>
  dplyr::distinct(resolved_p, n_plasmid_in_phylum)

plot_count_tax <- count_tax |>
  dplyr::mutate(resolved_p = stringr::str_remove(resolved_p, "p__")) |>
  dplyr::filter(resolved_p %in% y_order) |>
  ggplot2::ggplot(
    ggplot2::aes(y = factor(resolved_p, y_order), x = n_plasmid_in_phylum)
  ) +
  ggplot2::geom_col(fill = "#c5cad7", linewidth = 0.24, width = 0.7) +
  ggplot2::scale_x_log10(
    name = "n",
    expand = ggplot2::expansion(mult = c(0, 0.1)),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    rect = ggplot2::element_rect(linewidth = 0.24),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      size = 5,
      colour = "black",
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    axis.ticks.y = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.border = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(
      color = "black",
      linewidth = 0.24,
      lineend = "round"
    ),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "#ebebeb"
    ),
    panel.grid.minor.x = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "#ebebeb"
    ),
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_count_tax

layout <- "
#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD##
#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD##
ACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCEE
ACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCEE
ACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCEE
ACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCEE
BCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCEE
BCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCEE
BCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCEE
BCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCEE
"

plot_S02 <-
  plot_bac_tree +
  plot_ar_tree +
  plot_heatmap +
  plot_count_sys +
  plot_count_tax +
  patchwork::plot_layout(design = layout)

plot_S02

plot_S02 |>
  ggplot2::ggsave(
    filename = "plots/figS02.pdf",
    width = 364.8,
    height = 120,
    units = "mm",
    dpi = 300
  )

# FIG S04 — FOLLOW-UP EXPLORATION OF FIG 01 ------------------------------------

gtdb_pd <- readr::read_tsv("data/gtdb/gtdb_220_pd.tsv", show_col_types = FALSE) |>
  dplyr::mutate(pd = pd / 100)

system_types_by_taxa <- system_count_by_taxa |>
  dplyr::filter(!is.na(prevalence)) |>
  dplyr::distinct(resolved_p, n_plasmid_in_phylum, type) |>
  dplyr::summarise(
    system_types = dplyr::n(),
    .by = c(resolved_p, n_plasmid_in_phylum)
  ) |>
  dplyr::arrange(dplyr::desc(system_types)) |>
  dplyr::mutate(phylum = stringr::str_remove(resolved_p, "p__")) |>
  dplyr::select(phylum, n_plasmid_in_phylum, system_types)

log_model <- lm(
  system_types ~ log(n_plasmid_in_phylum),
  data = system_types_by_taxa
)

summary(log_model)

system_types_by_taxa_vs_model <- system_types_by_taxa |>
  dplyr::mutate(
    predicted_value = predict(log_model),
    residual = residuals(log_model),
    abs_residual = abs(residual)
  ) |>
  dplyr::left_join(gtdb_pd, by = dplyr::join_by(phylum))

residual_sd <- sd(system_types_by_taxa_vs_model$residual)
threshold <- 2 * residual_sd

system_types_by_taxa_vs_model <- system_types_by_taxa_vs_model |>
  dplyr::mutate(
    outlier_label = dplyr::if_else(abs_residual > threshold, phylum, "")
  )

system_types_by_taxa_vs_model |>
  dplyr::filter(abs_residual > threshold & residual < 0) |>
  dplyr::select(phylum, abs_residual) |>
  purrr::pwalk(
    \(phylum, abs_residual) {
      cli::cli_alert_info(
        "Plasmids of {phylum} encode {round(abs_residual)} less system types than expected"
      )
    }
  )

plot_S04 <- system_types_by_taxa_vs_model |>
  ggplot2::ggplot(ggplot2::aes(x = predicted_value, y = residual)) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = outlier_label),
    size = 5 / 2.845,
    color = "black",
    segment.size = 0.24,
    segment.colour = "black",
    min.segment.length = 0,
    nudge_y = -10,
    direction = "x",
    max.overlaps = Inf
  ) +
  ggplot2::geom_point(ggplot2::aes(colour = pd), size = 0.5) +
  ggplot2::scale_colour_viridis_c(
    name = "Contribution to phylogenetic diversity\nin the GTDB reference tree",
    option = "rocket",
    direction = -1,
    limits = c(0, 0.2),
    labels = scales::percent_format()
  ) +
  ggplot2::scale_x_continuous(
    name = "No. predicted system types",
    limits = c(-15, 80),
    expand = ggplot2::expansion()
  ) +
  ggplot2::scale_y_continuous(
    name = "Residual (observed - predicted)",
    limits = c(-50, 50),
    expand = ggplot2::expansion()
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.24
  ) +
  ggplot2::geom_hline(
    yintercept = threshold,
    linetype = "dashed",
    color = "#026eae",
    linewidth = 0.24
  ) +
  ggplot2::geom_hline(
    yintercept = -threshold,
    linetype = "dashed",
    color = "#c5373d",
    linewidth = 0.24
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 7, colour = "black"),
    axis.text.y = ggplot2::element_text(size = 7, colour = "black"),
    axis.ticks.x = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      "black"
    ),
    axis.ticks.y = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      "black"
    ),
    axis.title.x = ggplot2::element_text(size = 7, colour = "black"),
    axis.title.y = ggplot2::element_text(size = 7, colour = "black"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round"
    ),
    panel.grid.minor.y = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round"
    ),
    panel.border = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(
      color = "black",
      linewidth = 0.24,
      lineend = "round"
    ),
    legend.key.width = ggplot2::unit(0.1, "cm"),
    legend.key.height = ggplot2::unit(0.4, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.title = ggplot2::element_text(
      size = 7,
      colour = "black",
      angle = 90,
      hjust = 0.5
    ),
    legend.title.position = "left",
    legend.background = ggplot2::element_blank(),
    legend.margin = ggplot2::margin(),
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(2, 0, 0, 0)
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_S04

plot_S04 |>
  ggplot2::ggsave(
    filename = "plots/figS04.pdf",
    width = 88.2,
    height = 40,
    units = "mm"
  )

# FIGURE S5 — DEFENSE SYSTEM CO-OCCURRENCE -------------------------------------

plsdb_metadata |>
  dplyr::select(plasmid_seqid) |>
  dplyr::left_join(plsdb_plasmid_defense, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::filter(!is.na(type)) |>
  dplyr::summarise(n = dplyr::n(), .by = plasmid_seqid) |>
  plyr::summarise(
    pct_gt_1 = round(mean(n > 1) * 100)
  )

plasmid_defense_types <- plsdb_metadata_rep |>
  dplyr::select(plasmid_seqid) |>
  dplyr::left_join(
    plsdb_plasmid_defense,
    by = dplyr::join_by(plasmid_seqid)
  ) |>
  dplyr::filter(!is.na(feature)) |>
  dplyr::distinct(plasmid_seqid, type) |>
  dplyr::mutate(present = 1L) |>
  dplyr::rename(entity = type)

def_type_matrix <- plasmid_defense_types |>
  tidyr::pivot_wider(
    id_cols = plasmid_seqid,
    names_from = entity,
    values_from = present,
    values_fill = 0L
  ) |>
  tibble::column_to_rownames("plasmid_seqid") |>
  as.matrix()

if (!file.exists("data/plsdb_defense_type_affinity.xlsx")) {
  def_type_affinity <- def_type_matrix |>
    CooccurrenceAffinity::affinity(row.or.col = "col", squarematrix = c("all"))
  
  writexl::write_xlsx(
    def_type_affinity$all,
    "data/plsdb_defense_type_affinity.xlsx"
  )
}

def_type_affinity_all <- readxl::read_xlsx(
  "data/plsdb_defense_type_affinity.xlsx"
)

stat_def_type_affinity <- def_type_affinity_all |>
  tibble::as_tibble() |>
  dplyr::mutate(
    category = dplyr::case_when(
      p_value <= 0.05 & alpha_mle >= 0 ~ "Statistically significant +ve affinity",
      p_value <= 0.05 & alpha_mle < 0 ~ "Statistically significant -ve affinity",
      .default = "Other"
    )
  ) |>
  dplyr::summarise(n = dplyr::n(), .by = category) |>
  dplyr::mutate(p = round(n / sum(n) * 100)) |>
  dplyr::filter(!category == "Other")

purrr::pwalk(
  stat_def_type_affinity,
  \(category, n, p) {
    cli::cli_alert_info(
      "{category}: {n} ({p}%)"
    )
  }
)

plasmid_defense_subtypes <- plsdb_metadata_rep |>
  dplyr::select(plasmid_seqid) |>
  dplyr::left_join(
    plsdb_plasmid_defense,
    by = dplyr::join_by(plasmid_seqid)
  ) |>
  dplyr::filter(!is.na(feature)) |>
  dplyr::distinct(plasmid_seqid, subtype) |>
  dplyr::mutate(present = 1L) |>
  dplyr::rename(entity = subtype)

def_subtype_matrix <- plasmid_defense_subtypes |>
  tidyr::pivot_wider(
    id_cols = plasmid_seqid,
    names_from = entity,
    values_from = present,
    values_fill = 0L
  ) |>
  tibble::column_to_rownames("plasmid_seqid") |>
  as.matrix()

if (!file.exists("data/plsdb_defense_subtype_affinity.xlsx")) {
  def_subtype_affinity <- def_subtype_matrix |>
    CooccurrenceAffinity::affinity(row.or.col = "col", squarematrix = c("all"))
  
  writexl::write_xlsx(
    def_subtype_affinity$all,
    "data/plsdb_defense_subtype_affinity.xlsx"
  )
}

def_subtype_affinity_all <- readxl::read_xlsx(
  "data/plsdb_defense_subtype_affinity.xlsx"
)

stat_def_subtype_affinity <- def_subtype_affinity_all |>
  tibble::as_tibble() |>
  dplyr::mutate(
    category = dplyr::case_when(
      p_value <= 0.05 & alpha_mle >= 0 ~ "Statistically significant +ve affinity",
      p_value <= 0.05 & alpha_mle < 0 ~ "Statistically significant -ve affinity",
      .default = "Other"
    )
  ) |>
  dplyr::summarise(n = dplyr::n(), .by = category) |>
  dplyr::mutate(p = round(n / sum(n) * 100)) |>
  dplyr::filter(!category == "Other")

purrr::pwalk(
  stat_def_subtype_affinity,
  \(category, n, p) {
    cli::cli_alert_info(
      "{category}: {n} ({p}%)"
    )
  }
)

def_type_affinity_data <- def_type_affinity$all |>
  dplyr::mutate(
    p_value = as.double(p_value)
  ) |>
  dplyr::filter(entity_1_count_mA >= 50 & entity_2_count_mB >= 50) |>
  dplyr::mutate(
    alpha_mle = dplyr::case_when(
      obs_cooccur_X == 0 ~ NA,
      .default = alpha_mle
    )
  ) |>
  dplyr::mutate(
    label = dplyr::case_when(
      dplyr::between(p_value, 0.05, 1) ~ "",
      dplyr::between(p_value, 0.01, 0.05) ~ "*",
      dplyr::between(p_value, 0.001, 0.01) ~ "**",
      dplyr::between(p_value, 0, 0.001) ~ "***"
    )
  ) |>
  dplyr::mutate(
    label = dplyr::case_when(
      dplyr::between(p_value, 0, 0.001) ~ "*",
      .default = ""
    )
  ) |>
  dplyr::select(
    entity_1,
    entity_2,
    entity_1_count_mA,
    entity_2_count_mB,
    p_value,
    alpha_mle,
    label
  )


fill_limit <- max(
  abs(max(def_type_affinity_data$alpha_mle, na.rm = TRUE)),
  abs(min(def_type_affinity_data$alpha_mle, na.rm = TRUE))
)

def_type_affinity_y_axis <- rev(colnames(def_type_affinity$occur_mat[-1])) |>
  intersect(def_type_affinity_data$entity_2)

def_type_affinity_x_axis <- colnames(
  def_type_affinity$occur_mat
)[-length(colnames(def_type_affinity$occur_mat))] |>
  intersect(def_type_affinity_data$entity_1)

plot_S05 <- def_type_affinity_data |>
  ggplot2::ggplot(ggplot2::aes(x = entity_1, y = entity_2, fill = alpha_mle)) +
  ggplot2::geom_tile(colour = "#575653") +
  ggplot2::geom_point(
    ggplot2::aes(
      size = dplyr::if_else(dplyr::between(p_value, 0, 0.05), "dot", "nodot")
    ),
    colour = "#575653"
  ) +
  ggplot2::scale_size_manual(values = c(dot = 0.2, nodot = -1), guide = "none") +
  ggplot2::coord_fixed(clip = "off") +
  ggplot2::ylim(def_type_affinity_y_axis) +
  ggplot2::xlim(def_type_affinity_x_axis) +
  ggplot2::scale_fill_gradientn(
    name = "Alpha ",
    colours = c(
      "#026eae",
      "#026eae",
      "#5496ce",
      "#9bcae9",
      "#c5e5fb",
      "#ffffff",
      "#ffffff",
      "#f6ceca",
      "#e9a0a5",
      "#dc6465",
      "#c5373d",
      "#c5373d"
    ),
    na.value = "#E9E0D4",
    limits = c(-fill_limit, fill_limit),
    breaks = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      size = 5,
      colour = "black",
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = ggplot2::element_text(size = 5, colour = "black"),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    legend.key.height = ggplot2::unit(0.1, "cm"),
    legend.key.width = ggplot2::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "bottom",
    legend.title.position = "right",
    legend.ticks.length = ggplot2::unit(-0.2, 'cm'),
    legend.ticks = ggplot2::element_line(colour = "black", linewidth = 0.24),
    legend.frame = ggplot2::element_rect(colour = "black", linewidth = 0.24),
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(2, 2, 0, 0)
  )

plot_S05

plot_S05 |>
  ggplot2::ggsave(
    filename = "plots/figS05.pdf",
    width = 182.4,
    height = 182.4,
    units = "mm"
  )

top_10_def_subtype_affinity <- def_subtype_affinity$all |>
  tibble::as_tibble() |>
  dplyr::filter(entity_1_count_mA >= 50 & entity_2_count_mB >= 50) |>
  dplyr::arrange(dplyr::desc(alpha_mle)) |>
  dplyr::slice_head(n = 15) |>
  dplyr::select(
    entity_1,
    entity_2,
    entity_1_count_mA,
    entity_2_count_mB,
    obs_cooccur_X,
    alpha_mle,
    alpha_medianInt,
    p_value
  ) |>
  dplyr::mutate(
    alpha_medianInt = stringr::str_remove_all(alpha_medianInt, "\\[|\\]")
  ) |>
  tidyr::separate_wider_delim(
    alpha_medianInt,
    delim = ", ",
    names = c("lci", "uci")
  ) |>
  dplyr::mutate(dplyr::across(c(lci, uci), ~ as.double(.x))) |>
  dplyr::mutate(
    p1 = obs_cooccur_X / entity_1_count_mA * 100,
    p2 = obs_cooccur_X / entity_2_count_mB * 100
  )

plot_cooccurring_text <- top_10_def_subtype_affinity |>
  dplyr::mutate(
    y = factor(
      seq_len(nrow(top_10_def_subtype_affinity)),
      levels = rev(seq_len(nrow(top_10_def_subtype_affinity)))
    )
  ) |>
  ggplot2::ggplot() +
  ggtext::geom_richtext(
    ggplot2::aes(
      x = 1,
      y = y,
      label = paste0(
        stringr::str_remove(entity_1, "-Cas") |>
          stringr::str_replace("_", "-"),
        "<span style = 'color:#b3b3b3;'> (",
        round(p1),
        "%)</span>"
      )
    ),
    hjust = 0,
    size = 7 / 2.845,
    fill = NA,
    label.color = NA
  ) +
  ggtext::geom_richtext(
    ggplot2::aes(
      x = 2,
      y = y,
      label = paste0(
        stringr::str_remove(entity_2, "-Cas") |>
          stringr::str_replace("_", "-"),
        "<span style = 'color:#b3b3b3;'> (",
        round(p2),
        "%)</span>"
      )
    ),
    hjust = 0,
    size = 7 / 2.845,
    fill = NA,
    label.color = NA
  ) +
  ggplot2::scale_x_continuous(limits = c(1, 2.8), breaks = NULL) +
  ggplot2::scale_y_discrete() +
  ggplot2::theme_void() +
  ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank()
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_cooccurring_text

y_order <- top_10_def_subtype_affinity |>
  dplyr::mutate(y = paste0(entity_1, " :: ", entity_2)) |>
  dplyr::pull(y) |>
  rev()

plot_cooccurring_alpha <- top_10_def_subtype_affinity |>
  ggplot2::ggplot(ggplot2::aes(
    x = alpha_mle,
    y = factor(paste0(entity_1, " :: ", entity_2), y_order)
  )) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 2,
      xend = alpha_mle,
      yend = factor(paste0(entity_1, " :: ", entity_2), y_order)
    ),
    linewidth = 0.24,
    lineend = "round"
  ) +
  ggplot2::geom_point(size = 0.24) +
  ggplot2::scale_x_continuous(
    name = "Alpha MLE",
    limits = c(2, 5),
    expand = ggplot2::expansion(mult = c(0, 0)),
    breaks = c(2, 3, 4, 5)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 7, colour = "black"),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.ticks.y = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.line.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round"
    ),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round"
    ),
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
  ggplot2::guides(x = ggplot2::guide_axis(minor.ticks = TRUE)) +
  ggplot2::coord_cartesian(clip = "off")

plot_cooccurring_alpha

plot_01J <-
  plot_cooccurring_text + plot_cooccurring_alpha + 
  patchwork::plot_layout(widths = c(4, 1))

plot_01J

plot_01J |>
  ggplot2::ggsave(
    filename = "plots/fig01_J.pdf",
    width = 80,
    height = 65,
    units = "mm",
    dpi = 300
  )

# ENRICHMENT BY TAXA -----------------------------------------------------------

plsdb_metadata <- readxl::read_xlsx(
  "data/plsdb_plasmid-host_metadata_with_taxonomy.xlsx"
)

plsdb_metadata_rep <- plsdb_metadata |>
  dplyr::filter(representative == TRUE)

phylum_stats <- plsdb_metadata_rep |>
  dplyr::distinct(gtdb_phylum, host_acc, plasmid_seqid) |>
  dplyr::arrange(gtdb_phylum, host_acc, plasmid_seqid) |>
  dplyr::mutate(total_hosts = length(unique(host_acc)), .by = gtdb_phylum) |>
  dplyr::mutate(total_plasmids = dplyr::n(), .by = gtdb_phylum) |>
  dplyr::distinct(gtdb_phylum, total_hosts, total_plasmids)

plasmid_system_count_by_taxa <- plsdb_metadata_rep |>
  dplyr::distinct(gtdb_phylum, plasmid_seqid) |>
  dplyr::left_join(
    plsdb_plasmid_defense,
    by = dplyr::join_by(plasmid_seqid)
  ) |>
  dplyr::left_join(phylum_stats, by = dplyr::join_by(gtdb_phylum)) |>
  dplyr::filter(!is.na(type)) |>
  dplyr::distinct(plasmid_seqid, gtdb_phylum, total_hosts, total_plasmids, subtype) |>
  dplyr::summarise(plasmids_w_sys = dplyr::n(), .by = c(gtdb_phylum, total_hosts, total_plasmids, subtype)) |>
  dplyr::arrange(gtdb_phylum, subtype)

chromosome_system_count_by_taxa <- plsdb_metadata |>
  dplyr::filter(representative == TRUE) |>
  dplyr::distinct(gtdb_phylum, host_acc) |>
  dplyr::left_join(
    plsdb_host_defense,
    by = dplyr::join_by(host_acc)
  ) |>
  dplyr::left_join(phylum_stats, by = dplyr::join_by(gtdb_phylum)) |>
  dplyr::filter(!is.na(type)) |>
  dplyr::distinct(host_acc, gtdb_phylum, total_hosts, total_plasmids, subtype) |>
  dplyr::summarise(hosts_w_sys = dplyr::n(), .by = c(gtdb_phylum, total_hosts, total_plasmids, subtype)) |>
  dplyr::arrange(gtdb_phylum, subtype)

system_count_by_taxa_comparison <- plasmid_system_count_by_taxa |>
  dplyr::full_join(
    chromosome_system_count_by_taxa, 
    by = dplyr::join_by(gtdb_phylum, total_hosts, total_plasmids, subtype)
  ) |>
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ tidyr::replace_na(., 0))) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    mat = list(matrix(c(
      plasmids_w_sys, (total_plasmids - plasmids_w_sys),
      hosts_w_sys,    (total_hosts - hosts_w_sys)
    ), nrow = 2)),
    ft = list(fisher.test(mat)),
    p_value = ft$p.value,
    odds_ratio = unname(ft$estimate),
    log_odds_cc = log(((plasmids_w_sys + 0.5) * (total_hosts - hosts_w_sys + 0.5)) /
                        ((total_plasmids - plasmids_w_sys + 0.5) * (hosts_w_sys + 0.5))),
    prevalence_ratio = (plasmids_w_sys / total_plasmids) / (hosts_w_sys / total_hosts),
    total_with_system = plasmids_w_sys + hosts_w_sys
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) |>
  dplyr::mutate(
    significance = dplyr::case_when(
      p_adj > 0.05 ~ "Not Significant",
      p_adj < 0.05 ~ "Significant"
    ),
    category = dplyr::case_when(
      prevalence_ratio > 1 & p_adj < 0.05 ~ "Enriched on Plasmids",
      prevalence_ratio < 1 & p_adj < 0.05 ~ "Enriched on Chromosomes",
      prevalence_ratio > 1 & p_adj > 0.05 ~ "Higher proportion on Plasmids",
      prevalence_ratio < 1 & p_adj > 0.05 ~ "Higher proportion on Chromosomes",
      prevalence_ratio == 1 ~ "No observed difference"
    )
  ) |>
  dplyr::mutate(
    category = dplyr::case_when(
      hosts_w_sys == 0 ~ "Only detected on Plasmids",
      plasmids_w_sys == 0 ~ "Only detected on Chromosomes",
      .default = category
    )
  ) |>
  dplyr::select(gtdb_phylum, subtype, total_hosts, hosts_w_sys, total_plasmids, plasmids_w_sys, total_with_system, prevalence_ratio, odds_ratio, log_odds_cc, p_value, p_adj, significance, category) |>
  dplyr::mutate(
    prevalence_ratio = dplyr::if_else(prevalence_ratio == Inf, NA, prevalence_ratio)
  )

system_count_by_taxa_comparison |>
  dplyr::mutate(phylum = stringr::str_remove(gtdb_phylum, "p__")) |>
  dplyr::select(phylum, dplyr::everything()) |>
  dplyr::select(!gtdb_phylum) |>
  dplyr::rename(system = subtype, hosts_with_system = hosts_w_sys, plasmids_with_system = plasmids_w_sys, p_adjusted = p_adj) |>
  writexl::write_xlsx(
    "data/plsdb_defense_enrichment_by_phylum.xlsx"
  )

system_count_by_taxa_comparison_filt <- system_count_by_taxa_comparison |>
  dplyr::filter(total_hosts >= 20) |>
  dplyr::filter(total_plasmids >= 20) |>
  dplyr::mutate(
    dot = dplyr::case_when(
      stringr::str_detect(category, "Only") ~ "dot",
      significance == "Significant" ~ "star",
      .default = "nodot"
    )
  ) |>
  tidyr::complete(gtdb_phylum, subtype, fill = list(prevalence = NA))

# PREPARE GTDB TREES

bac120_p_tree <- treeio::read.tree("data/gtdb/bac120_r220_phylum.tree")

bac120_p_tip_labels <- bac120_p_tree[["tip.label"]]

tip_labels_keep_bac120 <- system_count_by_taxa_comparison_filt |>
  dplyr::distinct(gtdb_phylum) |>
  dplyr::pull()

tip_labels_drop_bac120 <- bac120_p_tip_labels |>
  setdiff(tip_labels_keep_bac120)

bac120_p_tree_reduced <- bac120_p_tree |>
  tidytree::drop.tip(tip_labels_drop_bac120)

bac120_p_tree_reduced_midpoint <- bac120_p_tree_reduced |>
  phytools::midpoint_root()

bac120_p_tree_clean <- bac120_p_tree_reduced_midpoint |>
  dplyr::as_tibble() |>
  dplyr::mutate(label = stringr::str_remove(label, "p__")) |>
  tidytree::as.phylo()

plot_bac_tree <- bac120_p_tree_clean |>
  ggtree::ggtree(
    linewidth = 0.24,
    color = "black"
  ) +
  ggtree::geom_treescale(
    fontsize = 7,
    linesize = 0.24,
    color = "black",
    width = 0.2,
    y = 2
  ) +
  ggplot2::theme(
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")
  )

plot_bac_tree

ar53_p_tree <- treeio::read.tree("data/gtdb/ar53_r220_phylum.tree")

ar53_p_tip_labels <- ar53_p_tree[["tip.label"]]

gtdb_archaea <- gtdb_taxonomy |> 
  dplyr::filter(gtdb_d == "d__Archaea") |>
  dplyr::distinct(gtdb_p) |> 
  dplyr::pull()

tip_labels_keep_ar53 <- system_count_by_taxa_comparison_filt |>
  dplyr::distinct(gtdb_phylum) |>
  dplyr::filter(gtdb_phylum %in% gtdb_archaea) |>
  dplyr::pull()

tip_labels_drop_ar53 <- ar53_p_tip_labels |>
  setdiff(tip_labels_keep_ar53)

ar53_p_tree_reduced <- ar53_p_tree |>
  tidytree::drop.tip(tip_labels_drop_ar53)

if (length(tip_labels_keep_ar53) == 1) {
  ar53_p_tree_reduced <- ar53_p_tree_reduced |>
    phytools::bind.tip(
      tip.label = "_dummy",
      where = ape::Ntip(ar53_p_tree_reduced) + 1,
      edge.length = 1
    )
}

ar53_p_tree_reduced_midpoint <- ar53_p_tree_reduced |>
  phytools::midpoint_root()

ar53_p_tree_clean <- ar53_p_tree_reduced_midpoint |>
  dplyr::as_tibble() |>
  dplyr::mutate(label = stringr::str_remove(label, "p__")) |>
  tidytree::as.phylo()

plot_ar_tree <- ar53_p_tree_clean |>
  ggtree::ggtree(
    linewidth = 0.24,
    color = "black"
  ) +
  ggtree::geom_treescale(
    fontsize = 7,
    linesize = 0.24,
    color = "black",
    width = 0.2,
    y = 2
  ) +
  ggplot2::theme(
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")
  )

plot_ar_tree

y_order <- dplyr::bind_rows(
  plot_bac_tree$data |> dplyr::mutate(order = 2),
  plot_ar_tree$data |> dplyr::mutate(order = 1)
) |>
  dplyr::as_tibble() |>
  dplyr::filter(isTip) |>
  dplyr::arrange(order, y) |>
  dplyr::pull(label)

x_order <- system_count_by_taxa_comparison_filt |>
  dplyr::filter(total_plasmids >= 50) |>
  dplyr::filter(total_hosts >= 50) |>
  dplyr::summarise(
    score_1 = sum(log_odds_cc > 0, na.rm = TRUE),
    score_2 = sum(log_odds_cc[log_odds_cc > 0], na.rm = TRUE),
    .by = subtype
  ) |>
  dplyr::arrange(desc(score_1), desc(score_2)) |>
  dplyr::slice_head(n = 40) |>
  dplyr::pull(subtype)

clip <- 5

plot_heatmap <- system_count_by_taxa_comparison_filt |>
  dplyr::mutate(
    dot = dplyr::if_else(is.na(dot), "nodot", dot)
  ) |>
  dplyr::filter(subtype %in% x_order) |>
  dplyr::mutate(
    gtdb_phylum = stringr::str_remove(gtdb_phylum, "p__")
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(subtype, x_order),
      y = factor(gtdb_phylum, y_order),
      fill = log_odds_cc
    )
  ) +
  ggplot2::geom_tile(colour = "#DEDEDE", linewidth = 0.24) +
  ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::scale_fill_gradient2(
    low = "#6E788D",
    mid = "white",
    high = "#5496ce",
    na.value = "white",
    name = "Odds Ratio",
    limits = c(-log(clip), log(clip)),
    oob = scales::squish,
    midpoint = 0,
    breaks = log(c(0.2, 0.5, 1, 2, 5)),
    labels = function(z) format(exp(z), digits = 2)
  ) +
  ggnewscale::new_scale(new_aes = "fill") +
  ggplot2::geom_point(
    ggplot2::aes(shape = dot, size = dot, colour = dot, fill = dot),
    stroke = 0.2
  ) +
  ggplot2::scale_shape_manual(
    values = c(dot = 43, star = 21, nodot = 1),
    guide = "none"
  ) +
  ggplot2::scale_fill_manual(
    values = c(dot = "white", star = "white", nodot = "white"),
    guide = "none"
  ) +
  ggplot2::scale_colour_manual(
    values = c(dot = "white", star = "white", nodot = "white"),
    guide = "none"
  ) +
  ggplot2::scale_size_manual(
    values = c(dot = 1, star = 0.3, nodot = -1),
    guide = "none"
  ) +
  ggplot2::scale_linewidth_manual(
    values = c(dot = 2, star = 0.2, nodot = -1),
    guide = "none"
  ) +
  ggplot2::labs(fill = "") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      size = 5,
      colour = "black",
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = ggplot2::element_text(size = 5, colour = "black"),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(linewidth = 0.24),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    legend.key.height = ggplot2::unit(0.1, "cm"),
    legend.key.width = ggplot2::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "bottom",
    legend.title.position = "left",
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(2, 2, 0, 0)
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_heatmap

plot_odds_sum_sys <- system_count_by_taxa_comparison_filt |>
  dplyr::filter(subtype %in% x_order) |>
  dplyr::summarise(
    score_2 = sum(log_odds_cc[log_odds_cc > 0], na.rm = TRUE),
    .by = subtype
  ) |>
  ggplot2::ggplot(ggplot2::aes(
    x = factor(subtype, x_order),
    y = score_2
  )) +
  ggplot2::geom_col(fill = "#c5cad7", linewidth = 0.24, width = 0.7) +
  ggplot2::scale_y_continuous(
    name = "n",
    expand = ggplot2::expansion(mult = c(0, 0.1)),
    breaks = c(0, 4, 8)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    rect = ggplot2::element_rect(linewidth = 0.24),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 5, colour = "black"),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.border = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(
      color = "black",
      linewidth = 0.24,
      lineend = "round"
    ),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "#ebebeb"
    ),
    panel.grid.minor.y = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "#ebebeb"
    ),
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_odds_sum_sys

plot_01K <- plot_odds_sum_sys + plot_heatmap +
  patchwork::plot_layout(heights = c(1, 4))

plot_01K

plot_01K |>
  ggplot2::ggsave(
    filename = "plots/fig01_K.pdf",
    width = 90.2,
    height = 74.2,
    units = "mm",
    dpi = 300
  )

# SUPPLEMENTARY

system_count_by_taxa_comparison <- system_count_by_taxa_comparison |>
  dplyr::mutate(
    dot = dplyr::case_when(
      stringr::str_detect(category, "Only") ~ "dot",
      significance == "Significant" ~ "star",
      .default = "nodot"
    )
  ) |>
  tidyr::complete(gtdb_phylum, subtype, fill = list(prevalence = NA))

# PREPARE GTDB TREES

bac120_p_tree <- treeio::read.tree("data/gtdb/bac120_r220_phylum.tree")

bac120_p_tip_labels <- bac120_p_tree[["tip.label"]]

tip_labels_keep_bac120 <- system_count_by_taxa_comparison |>
  dplyr::distinct(gtdb_phylum) |>
  dplyr::pull()

tip_labels_drop_bac120 <- bac120_p_tip_labels |>
  setdiff(tip_labels_keep_bac120)

bac120_p_tree_reduced <- bac120_p_tree |>
  tidytree::drop.tip(tip_labels_drop_bac120)

bac120_p_tree_reduced_midpoint <- bac120_p_tree_reduced |>
  phytools::midpoint_root()

bac120_p_tree_clean <- bac120_p_tree_reduced_midpoint |>
  dplyr::as_tibble() |>
  dplyr::mutate(label = stringr::str_remove(label, "p__")) |>
  tidytree::as.phylo()

plot_bac_tree <- bac120_p_tree_clean |>
  ggtree::ggtree(
    linewidth = 0.24,
    color = "black"
  ) +
  ggtree::geom_treescale(
    fontsize = 7,
    linesize = 0.24,
    color = "black",
    width = 0.2,
    y = 2
  ) +
  ggplot2::theme(
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")
  )

plot_bac_tree

ar53_p_tree <- treeio::read.tree("data/gtdb/ar53_r220_phylum.tree")

ar53_p_tip_labels <- ar53_p_tree[["tip.label"]]

tip_labels_keep_ar53 <- system_count_by_taxa_comparison |>
  dplyr::distinct(gtdb_phylum) |>
  dplyr::pull()

tip_labels_drop_ar53 <- ar53_p_tip_labels |>
  setdiff(tip_labels_keep_ar53)

ar53_p_tree_reduced <- ar53_p_tree |>
  tidytree::drop.tip(tip_labels_drop_ar53)

ar53_p_tree_reduced_midpoint <- ar53_p_tree_reduced |>
  phytools::midpoint_root()

ar53_p_tree_clean <- ar53_p_tree_reduced_midpoint |>
  dplyr::as_tibble() |>
  dplyr::mutate(label = stringr::str_remove(label, "p__")) |>
  tidytree::as.phylo()

plot_ar_tree <- ar53_p_tree_clean |>
  ggtree::ggtree(
    linewidth = 0.24,
    color = "black"
  ) +
  ggtree::geom_treescale(
    fontsize = 7,
    linesize = 0.24,
    color = "black",
    width = 0.2,
    y = 2
  ) +
  ggplot2::theme(
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")
  )

plot_ar_tree

y_order <- dplyr::bind_rows(
  plot_bac_tree$data |> dplyr::mutate(order = 2),
  plot_ar_tree$data |> dplyr::mutate(order = 1)
) |>
  dplyr::as_tibble() |>
  dplyr::filter(isTip) |>
  dplyr::arrange(order, y) |>
  dplyr::pull(label)

x_order_pos <- system_count_by_taxa_comparison |>
  dplyr::summarise(
    score_1 = sum(log_odds_cc > 0, na.rm = TRUE),
    score_2 = sum(log_odds_cc[log_odds_cc > 0], na.rm = TRUE),
    .by = subtype
  ) |>
  dplyr::arrange(desc(score_1), desc(score_2)) |>
  dplyr::filter(score_1 != 0) |>
  dplyr::pull(subtype)

x_order_neg <- system_count_by_taxa_comparison |>
  dplyr::summarise(
    score_1 = sum(log_odds_cc < 0, na.rm = TRUE),
    score_2 = sum(log_odds_cc[log_odds_cc < 0], na.rm = TRUE),
    .by = subtype
  ) |>
  dplyr::arrange(score_1, score_2) |>
  dplyr::filter(!subtype %in% x_order_pos) |>
  dplyr::pull(subtype)

x_order <- c(x_order_pos, x_order_neg)

clip <- 5

plot_heatmap <- system_count_by_taxa_comparison |>
  dplyr::mutate(
    dot = dplyr::if_else(is.na(dot), "nodot", dot)
  ) |>
  dplyr::filter(subtype %in% x_order) |>
  dplyr::mutate(
    gtdb_phylum = stringr::str_remove(gtdb_phylum, "p__")
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(subtype, x_order),
      y = factor(gtdb_phylum, y_order),
      fill = log_odds_cc
    )
  ) +
  ggplot2::geom_tile(colour = "#DEDEDE", linewidth = 0.24) +
  ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::scale_fill_gradient2(
    low = "#6E788D",
    mid = "white",
    high = "#5496ce",
    na.value = "white",
    name = "Odds Ratio",
    limits = c(-log(clip), log(clip)),
    oob = scales::squish,
    midpoint = 0,
    breaks = log(c(0.2, 0.5, 1, 2, 5)),
    labels = function(z) format(exp(z), digits = 2)
  ) +
  ggnewscale::new_scale(new_aes = "fill") +
  ggplot2::geom_point(
    ggplot2::aes(shape = dot, size = dot, colour = dot, fill = dot),
    stroke = 0.2
  ) +
  ggplot2::scale_shape_manual(
    values = c(dot = 43, star = 21, nodot = 1),
    guide = "none"
  ) +
  ggplot2::scale_fill_manual(
    values = c(dot = "white", star = "white", nodot = "white"),
    guide = "none"
  ) +
  ggplot2::scale_colour_manual(
    values = c(dot = "white", star = "white", nodot = "white"),
    guide = "none"
  ) +
  ggplot2::scale_size_manual(
    values = c(dot = 1, star = 0.3, nodot = -1),
    guide = "none"
  ) +
  ggplot2::scale_linewidth_manual(
    values = c(dot = 2, star = 0.2, nodot = -1),
    guide = "none"
  ) +
  ggplot2::labs(fill = "") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      size = 5,
      colour = "black",
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = ggplot2::element_text(size = 5, colour = "black"),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(linewidth = 0.24),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    legend.key.height = ggplot2::unit(0.1, "cm"),
    legend.key.width = ggplot2::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "bottom",
    legend.title.position = "left",
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(2, 2, 0, 0)
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_heatmap

plot_odds_sum_sys <- system_count_by_taxa_comparison |>
  dplyr::filter(subtype %in% x_order) |>
  dplyr::summarise(
    score_2 = sum(log_odds_cc[log_odds_cc > 0], na.rm = TRUE),
    .by = subtype
  ) |>
  ggplot2::ggplot(ggplot2::aes(
    x = factor(subtype, x_order),
    y = score_2
  )) +
  ggplot2::geom_col(fill = "#c5cad7", linewidth = 0.24, width = 0.7) +
  ggplot2::scale_y_continuous(
    name = "n",
    expand = ggplot2::expansion(mult = c(0, 0.1)),
    breaks = c(0, 2, 4, 6, 8, 10)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    rect = ggplot2::element_rect(linewidth = 0.24),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 5, colour = "black"),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.border = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(
      color = "black",
      linewidth = 0.24,
      lineend = "round"
    ),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "#ebebeb"
    ),
    panel.grid.minor.y = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "#ebebeb"
    ),
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_odds_sum_sys

count_tax <- system_count_by_taxa |>
  dplyr::distinct(resolved_p, n_plasmid_in_phylum)

plot_count_tax <- count_tax |>
  dplyr::mutate(resolved_p = stringr::str_remove(resolved_p, "p__")) |>
  dplyr::filter(resolved_p %in% y_order) |>
  ggplot2::ggplot(
    ggplot2::aes(y = factor(resolved_p, y_order), x = n_plasmid_in_phylum)
  ) +
  ggplot2::geom_col(fill = "#c5cad7", linewidth = 0.24, width = 0.7) +
  ggplot2::scale_x_log10(
    name = "n",
    expand = ggplot2::expansion(mult = c(0, 0.1)),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    rect = ggplot2::element_rect(linewidth = 0.24),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      size = 5,
      colour = "black",
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    axis.ticks.y = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.border = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(
      color = "black",
      linewidth = 0.24,
      lineend = "round"
    ),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "#ebebeb"
    ),
    panel.grid.minor.x = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round",
      colour = "#ebebeb"
    ),
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_count_tax

plot_S06 <- plot_odds_sum_sys + plot_heatmap +
  patchwork::plot_layout(heights = c(1, 4))

plot_S06

plot_S06 |>
  ggplot2::ggsave(
    filename = "plots/figS06.pdf",
    width = 450,
    height = 140,
    units = "mm",
    dpi = 300
  )
