library(patchwork)

# READ PLSDB DATA --------------------------------------------------------------

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


# FIGURE 2A — HAS DEFENSE BY LENGTH --------------------------------------------

has_defense_vs_length <- plsdb_metadata_rep |>
  dplyr::left_join(
    plsdb_plasmid_defense,
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
  dplyr::distinct(plasmid_seqid, plasmid_length, mobility, has_defense) |>
  dplyr::mutate(
    mobility = stringr::str_to_sentence(mobility)
  )

has_defense_vs_length_combined <- has_defense_vs_length |>
  dplyr::mutate(mobility = "All") |>
  dplyr::bind_rows(has_defense_vs_length)

plot_02A <- has_defense_vs_length_combined |>
  ggplot2::ggplot(
    ggplot2::aes(x = plasmid_length, fill = has_defense)
  ) +
  ggplot2::geom_density(adjust = 0.5, alpha = 0.5, colour = NA) +
  ggplot2::scale_fill_manual(values = c("#c1c6d4", "#6e788d")) +
  ggplot2::scale_x_log10(
    expand = ggplot2::expansion(mult = c(0, 0)),
    labels = scales::label_log(base = 10),
    breaks = c(10, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7),
    limits = c(10^2, 10^7)
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.05))
  ) +
  ggplot2::labs(x = "Plasmid length (bp)", y = "Density") +
  ggplot2::geom_vline(
    xintercept = c(10000, 50000, 300000),
    colour = "black",
    linewidth = 0.24,
    linetype = "dashed"
  ) +
  ggplot2::facet_wrap(
    ~mobility,
    ncol = 1,
    axes = "all_x",
    axis.labels = "margins"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    
    axis.text.x = ggplot2::element_text(size = 7, colour = "black"),
    axis.text.y = ggplot2::element_text(size = 7, colour = "black"),
    
    axis.ticks.x = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round"
    ),
    axis.ticks.y = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round"
    ),
    panel.border = ggplot2::element_blank(),
    axis.line.x.bottom = ggplot2::element_line(linewidth = 0.24),
    axis.line.y.left = ggplot2::element_line(linewidth = 0.24),
    axis.line.x.top = ggplot2::element_blank(),
    axis.line.y.right = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      linewidth = 0.24,
      colour = "#E6E6E6"
    ),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    
    strip.background = ggplot2::element_rect(
      linewidth = 0.24,
      fill = NA,
      colour = NA
    ),
    strip.text = ggplot2::element_text(
      size = 5,
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

plot_02A

plot_02A |>
  ggplot2::ggsave(
    filename = "plots/fig02_A.pdf",
    width = 60,
    height = 90,
    units = "mm",
    dpi = 300
  )

stat_p_gt50kb_w_defense <- plsdb_metadata_rep |>
  dplyr::left_join(
    plsdb_plasmid_defense,
    by = dplyr::join_by(plasmid_seqid)
  ) |>
  dplyr::mutate(
    has_defense = dplyr::case_when(!is.na(type) ~ TRUE, .default = FALSE)
  ) |>
  dplyr::mutate(
    size_bin = dplyr::case_when(
      dplyr::between(plasmid_length, 0, 49999) ~ "<50",
      plasmid_length > 49999 ~ "≥50"
    )
  ) |>
  dplyr::summarise(n = dplyr::n(), .by = c(size_bin, has_defense)) |>
  dplyr::arrange(size_bin, has_defense) |>
  dplyr::mutate(p = n / sum(n) * 100, .by = size_bin) |>
  dplyr::filter(has_defense == TRUE)

stat_p_gt50kb_w_defense

functional_categories_plsdb_raw <- readr::read_tsv(
  "data/plsdb_gene_annotations.tsv"
)

plsdb_plasmid_defense_start <- plsdb_plasmid_defense |>
  dplyr::distinct(plasmid_seqid, start, .keep_all = TRUE) |>
  dplyr::rename(system_end = end)

functional_categories_plsdb <- functional_categories_plsdb_raw |>
  tidyr::separate_longer_delim(merged_operon, ",") |>
  dplyr::select(-type) |>
  dplyr::rename(plasmid_seqid = seqid) |>
  dplyr::filter(!is.na(merged_type)) |>
  dplyr::left_join(
    plsdb_plasmid_defense_start,
    by = dplyr::join_by(plasmid_seqid, start)
  ) |>
  dplyr::mutate(
    type = if (all(is.na(type))) {type[NA_integer_]} else {max(type, na.rm = TRUE)},
    subtype = if (all(is.na(subtype))) {subtype[NA_integer_]} else {max(subtype, na.rm = TRUE)}, .by = c(plasmid_seqid, merged_operon)
  ) |>
  dplyr::select(
    plasmid_seqid, start, end, merged_type, merged_subtype, merged_operon, type, subtype
  )

keep <- functional_categories_plsdb_raw |>
  dplyr::distinct(seqid) |>
  dplyr::pull()

all_gene_counts <- functional_categories_plsdb_raw |>
  dplyr::rename(plasmid_seqid = seqid) |>
  dplyr::summarise(total_genes = dplyr::n(), .by = plasmid_seqid)

defense_gene_counts <- functional_categories_plsdb |>
  dplyr::filter(!is.na(type)) |>
  dplyr::summarise(defense_genes = dplyr::n(), .by = plasmid_seqid)

gene_counts <- all_gene_counts |>
  dplyr::left_join(defense_gene_counts, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::mutate(defense_genes = dplyr::if_else(is.na(defense_genes), 0, defense_genes))

analysis_df <- has_defense_vs_length |>
  dplyr::filter(plasmid_seqid %in% keep) |>
  dplyr::left_join(gene_counts, by = "plasmid_seqid") |>
  dplyr::mutate(
    defense_present = defense_genes > 0,
    other_genes = total_genes - defense_genes,
    above_50kb = factor(
      plasmid_length >= 50000,
      levels = c(FALSE, TRUE),
      labels = c("<50 kb", "≥50 kb")
    )
  )

model_50kb_adjusted <- glm(
  defense_present ~ above_50kb + log1p(other_genes) + mobility,
  family = binomial(),
  data = analysis_df
)

broom::tidy(
  model_50kb_adjusted,
  exponentiate = TRUE,
  conf.int = TRUE
) |>
  dplyr::filter(term == "above_50kb≥50 kb")

# FIGURE 2B — HAS DEFENSE BY LENGTH BINNED -------------------------------------

x_order <- c("0-10", "10-50", "50-300", ">300")

has_defense_vs_length_combined_binned <- has_defense_vs_length_combined |>
  dplyr::mutate(n_mob = dplyr::n(), .by = mobility) |>
  dplyr::mutate(
    size_bin = dplyr::case_when(
      dplyr::between(plasmid_length, 0, 10000) ~ "0-10",
      dplyr::between(plasmid_length, 10001, 50000) ~ "10-50",
      dplyr::between(plasmid_length, 50001, 300000) ~ "50-300",
      plasmid_length > 300001 ~ ">300"
    )
  ) |>
  dplyr::mutate(n_size_bin = dplyr::n(), .by = size_bin) |>
  dplyr::mutate(n = dplyr::n(), .by = c(size_bin, mobility, has_defense)) |>
  dplyr::mutate(p = n / n_size_bin) |>
  dplyr::distinct(size_bin, has_defense, mobility, n, p) |>
  dplyr::mutate(category_x = paste0(size_bin, " - ", has_defense)) |>
  dplyr::mutate(mobility = stringr::str_to_sentence(mobility))

has_defense_vs_length_combined_binned |>
  dplyr::filter(size_bin %in% c("50-300", ">300")) |>
  dplyr::filter(mobility != "All") |>
  dplyr::summarise(n = sum(n), .by = c(mobility, has_defense)) |>
  dplyr::arrange(dplyr::desc(n))

has_defense_vs_length |>
  dplyr::filter(has_defense == "With Defense") |>
  dplyr::summarise(n = dplyr::n(), .by = mobility) |>
  dplyr::mutate(p = round(n / sum(n) * 100))

has_defense_vs_length |>
  dplyr::filter(mobility == "Conjugative") |>
  dplyr::summarise(n = dplyr::n(), .by = has_defense) |>
  dplyr::mutate(p = round(n / sum(n) * 100))
  
plot_02B <- has_defense_vs_length_combined_binned |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(size_bin, x_order),
      y = n,
      fill = has_defense
    ),
    colour = NA
  ) +
  ggplot2::scale_fill_manual(values = c("#e3e5eb", "#a7adb7")) +
  ggplot2::geom_col(
    linewidth = 0.24,
    position = ggplot2::position_dodge(width = 0.75),
    width = 0.7
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0)),
    breaks = c(0, 1000, 2000, 3000, 4000),
    limits = c(0, 4000),
    labels = scales::label_comma()
  ) +
  ggplot2::labs(x = "Plasmid length (kb)", y = "No. of plasmids") +
  ggplot2::geom_vline(
    xintercept = c(10000, 50000, 300000),
    colour = "black",
    linewidth = 0.24,
    linetype = "dashed"
  ) +
  ggplot2::facet_wrap(
    ~mobility,
    ncol = 1,
    axes = "all_x",
    axis.labels = "margins",
    scales = "free_y"
  ) +
  ggh4x::scale_y_facet(
    mobility == "All",
    expand = ggplot2::expansion(mult = c(0, 0)),
    breaks = c(0, 2000, 4000, 6000, 8000),
    limits = c(0, 8000),
    labels = scales::label_comma()
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    
    axis.text.x = ggplot2::element_text(size = 7, colour = "black"),
    axis.text.y = ggplot2::element_text(size = 7, colour = "black"),
    
    axis.ticks.x = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round"
    ),
    axis.ticks.y = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round"
    ),
    panel.border = ggplot2::element_blank(),
    axis.line.x.bottom = ggplot2::element_line(linewidth = 0.24),
    axis.line.y.left = ggplot2::element_line(linewidth = 0.24),
    axis.line.x.top = ggplot2::element_blank(),
    axis.line.y.right = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      linewidth = 0.24,
      colour = "#E6E6E6"
    ),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    
    strip.background = ggplot2::element_rect(
      linewidth = 0.24,
      fill = NA,
      colour = NA
    ),
    strip.text = ggplot2::element_text(
      size = 5,
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

plot_02B

plot_02B |>
  ggplot2::ggsave(
    filename = "plots/fig02_B.pdf",
    width = 60,
    height = 86.5,
    units = "mm",
    dpi = 300
  )

plsdb_metadata_rep |>
  dplyr::distinct(plasmid_seqid, mobility) |>
  writexl::write_xlsx("data/plsdb_plasmid_mob.xlsx")

# FIGURE 2C — NUMBER OF DEFENSE BY LENGTH --------------------------------------

nsys_per_length <- plsdb_metadata_rep |>
  dplyr::distinct(plasmid_seqid, plasmid_length) |>
  dplyr::left_join(
    plsdb_plasmid_defense,
    by = dplyr::join_by(plasmid_seqid)
  ) |>
  dplyr::mutate(n_sys = dplyr::n(), .by = plasmid_seqid) |>
  dplyr::mutate(n_sys = dplyr::case_when(is.na(type) ~ 0, .default = n_sys)) |>
  dplyr::mutate(system_length = end - start) |>
  dplyr::mutate(
    total_system_length = sum(system_length),
    .by = plasmid_seqid
  ) |>
  dplyr::distinct(plasmid_seqid, plasmid_length, n_sys, total_system_length) |>
  dplyr::mutate(
    plasmid_length_minus_system = plasmid_length - total_system_length
  ) |>
  dplyr::mutate(
    plasmid_length_minus_system = dplyr::case_when(
      is.na(plasmid_length_minus_system) ~ plasmid_length,
      .default = plasmid_length_minus_system
    )
  )

nsys_per_length_filt <- nsys_per_length |>
  dplyr::mutate(n = dplyr::n(), .by = n_sys) |>
  dplyr::filter(n > 10)

nsys_per_length_count <- nsys_per_length |>
  dplyr::summarise(n = dplyr::n(), .by = n_sys) |>
  dplyr::arrange(dplyr::desc(n_sys))

nsys_per_length_count_filt <- nsys_per_length_count |>
  dplyr::filter(n > 10)

fig_02C <- nsys_per_length_filt |>
  ggplot2::ggplot(ggplot2::aes(
    x = plasmid_length_minus_system,
    y = n_sys,
    group = n_sys
  )) +
  ggplot2::geom_violin(
    width = 0.8,
    adjust = 0.8,
    scale = "width",
    fill = "#e3e5eb",
    colour = NA,
    linewidth = 0.24
  ) +
  ggplot2::geom_boxplot(
    width = 0.2,
    outliers = FALSE,
    linewidth = 0.24,
    fill = "white",
    colour = "#6e788d"
  ) +
  ggtext::geom_richtext(
    data = nsys_per_length_count_filt,
    ggplot2::aes(
      y = n_sys,
      label = paste0("<i>n</i> = ", prettyNum(n, ","))
    ),
    x = Inf,
    hjust = 0,
    size = 7 / 2.845,
    fill = NA,
    label.color = NA
  ) +
  ggplot2::scale_x_log10(
    expand = ggplot2::expansion(),
    labels = scales::label_log(base = 10),
    breaks = c(1, 10, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7)
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.01, 0.01)),
    breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
  ) +
  ggplot2::geom_vline(
    xintercept = c(10000, 50000, 300000),
    colour = "black",
    linewidth = 0.24,
    linetype = "dashed"
  ) +
  ggplot2::labs(x = "Plasmid length (bp)", y = "No. of defense systems") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    axis.text.x = ggplot2::element_text(size = 7, colour = "black"),
    axis.text.y = ggplot2::element_text(size = 7, colour = "black"),
    axis.ticks.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    axis.ticks.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.border = ggplot2::element_rect(linewidth = 0.24),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    legend.key.height = grid::unit(0.25, "cm"),
    legend.key.width = grid::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "top",
    legend.title = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(5.5, 44, 5.5, 5.5)
  ) +
  ggplot2::coord_cartesian(clip = "off", xlim = c(15^2, 12.5^6)) +
  ggplot2::guides(colour = "none")

fig_02C

fig_02C |>
  ggplot2::ggsave(
    filename = "plots/fig02_C.pdf",
    width = 60,
    height = 76,
    units = "mm",
    dpi = 300
  )

# FIGURE 2D - PCN ENRICHMENT ---------------------------------------------------

stat_has_pcn_p <- plsdb_metadata |>
  dplyr::distinct(plasmid_seqid, pcn) |>
  dplyr::mutate(has_pcn = dplyr::if_else(!is.na(pcn), TRUE, FALSE)) |>
  dplyr::summarise(n = dplyr::n(), .by = has_pcn) |>
  dplyr::mutate(p = round(n / sum(n) * 100, 1)) |>
  dplyr::filter(has_pcn == TRUE) |>
  dplyr::pull(p)

cli::cli_alert_info("Plasmids with PCN data: {scales::comma(stat_has_pcn_p)}%")

pcn_stats <- plsdb_metadata |>
  dplyr::distinct(plasmid_seqid, pcn) |>
  dplyr::filter(!is.na(pcn)) |>
  dplyr::left_join(
    plsdb_plasmid_defense,
    by = dplyr::join_by(plasmid_seqid)
  ) |>
  dplyr::mutate(has_defense = dplyr::if_else(!is.na(type), TRUE, FALSE)) |>
  dplyr::distinct(plasmid_seqid, pcn, has_defense) |>
  dplyr::summarise(
    avg = round(mean(pcn), 2),
    med = round(median(pcn), 2),
    q1 = round(quantile(pcn, 0.25), 2),
    q3 = round(quantile(pcn, 0.75), 2),
    .by = has_defense
  ) |>
  dplyr::arrange(dplyr::desc(has_defense))

pcn_stats |>
  purrr::pwalk(
    \(has_defense, avg, med, q1, q3) {
      label <- dplyr::recode_values(has_defense, TRUE ~ "With Defense PCN", FALSE ~ "Without Defense PCN")
      cli::cli_alert_info(
        "{label}: mean = {avg}, median = {med}, IQR = {q1}-{q3}"
      )
    }
  )

pcn_split <- pcn_stats |>
  dplyr::filter(has_defense == TRUE) |>
  dplyr::pull(med) |>
  round(2)

pcn_plasmids <- plsdb_metadata |>
  dplyr::distinct(plasmid_seqid, pcn) |>
  dplyr::filter(!is.na(pcn)) |>
  dplyr::mutate(
    pcn_bin = dplyr::if_else(
      pcn >= pcn_split,
      "high",
      "low"
    )
  )

subtype_presence <- plsdb_plasmid_defense |>
  dplyr::filter(!is.na(subtype)) |>
  dplyr::distinct(plasmid_seqid, subtype)

subtype_presence_absence <- tidyr::crossing(
  plasmid_seqid = pcn_plasmids$plasmid_seqid,
  subtype = unique(subtype_presence$subtype)
) |>
  dplyr::left_join(pcn_plasmids, by = "plasmid_seqid") |>
  dplyr::left_join(
    subtype_presence |>
      dplyr::mutate(present = TRUE),
    by = c("plasmid_seqid", "subtype")
  ) |>
  dplyr::mutate(present = tidyr::replace_na(present, FALSE))

pcn_enrichment <- subtype_presence_absence |>
  dplyr::summarise(
    test = list({
      tab <- table(
        factor(pcn_bin, levels = c("high", "low")),
        factor(present, levels = c(TRUE, FALSE))
      )
      fisher.test(tab)
    }),
    .by = subtype
  ) |>
  dplyr::mutate(
    tidied = purrr::map(test, broom::tidy)
  ) |>
  tidyr::unnest(tidied) |>
  dplyr::select(
    subtype,
    odds_ratio = estimate,
    conf_low = conf.low,
    conf_high = conf.high,
    p_value = p.value
  ) |>
  dplyr::mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    direction = dplyr::case_when(
      odds_ratio > 1 ~ paste0("Enriched at PCN > ", round(pcn_split, 2)),
      odds_ratio < 1 ~ paste0("Enriched at PCN < ", round(pcn_split, 2)),
      TRUE ~ "No direction"
    ),
    significant = p_adjusted < 0.05
  )

pcn_counts <- subtype_presence_absence |>
  dplyr::summarise(
    n_present = sum(present),
    n_total = dplyr::n(),
    prevalence = mean(present),
    .by = c(subtype, pcn_bin)
  ) |>
  tidyr::pivot_wider(
    names_from = pcn_bin,
    values_from = c(n_present, n_total, prevalence),
    names_glue = "{pcn_bin}_{.value}"
  )

pcn_enrichment <- pcn_enrichment |>
  dplyr::left_join(pcn_counts, by = "subtype") |>
  dplyr::mutate(
    prevalence_ratio = high_prevalence / low_prevalence,
    log2_prevalence_ratio = log2(prevalence_ratio)
  ) |>
  dplyr::arrange(p_adjusted)

pcn_enrichment

n <- 31
n_each <- floor(n / 2)
min_system_count <- 10
fdr_cutoff <- 0.05

pcn_enrichment_data <- pcn_enrichment |>
  dplyr::mutate(
    system_count = high_n_present + low_n_present,
    high_n_absent = high_n_total - high_n_present,
    low_n_absent = low_n_total - low_n_present,
    odds_ratio_corrected =
      ((high_n_present + 0.5) * (low_n_absent + 0.5)) /
      ((high_n_absent + 0.5) * (low_n_present + 0.5)),
    log2_or = log2(odds_ratio),
    log2_or_plot = dplyr::if_else(
      is.finite(log2_or),
      log2_or,
      log2(odds_ratio_corrected)
    ),
    log2_conf_low = log2(conf_low),
    log2_conf_high = log2(conf_high),
    significant = p_adjusted < fdr_cutoff,
    pcn_group = dplyr::case_when(
      log2_or_plot > 0 ~ paste0("PCN > ", round(pcn_split, 2)),
      log2_or_plot < 0 ~ paste0("PCN < ", round(pcn_split, 2)),
      TRUE ~ "No difference"
    )
  )

eligible <- pcn_enrichment_data |>
  dplyr::filter(
    system_count >= 10,
    log2_or_plot != 0
  )

selected <- dplyr::bind_rows(
  eligible |>
    dplyr::filter(log2_or_plot > 0) |>
    dplyr::slice_max(log2_or_plot, n = n_each, with_ties = FALSE),
  eligible |>
    dplyr::filter(log2_or_plot < 0) |>
    dplyr::slice_min(log2_or_plot, n = n_each, with_ties = FALSE)
)

pcn_enrichment_sample <- selected |>
  dplyr::bind_rows(
    eligible |>
      dplyr::anti_join(selected, by = "subtype") |>
      dplyr::slice_max(
        abs(log2_or_plot),
        n = n - nrow(selected),
        with_ties = FALSE
      )
  ) |>
  dplyr::arrange(dplyr::desc(log2_or_plot), subtype)

x_order <- pcn_enrichment_sample |>
  dplyr::pull(subtype)

group_levels <- c(
  paste0("PCN > ", round(pcn_split, 2)),
  paste0("PCN < ", round(pcn_split, 2))
)

pcn_enrichment_sample <- pcn_enrichment_sample |>
  dplyr::mutate(
    pcn_group = factor(pcn_group, levels = group_levels)
  )

plot_segment <- pcn_enrichment_sample |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(subtype, levels = x_order),
      y = log2_or
    )
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    colour = "black",
    linewidth = 0.24
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      xend = factor(subtype, levels = x_order),
      y = 0,
      yend = log2_or,
      colour = pcn_group
    ),
    linewidth = 0.4
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      colour = pcn_group,
      shape = significant
    ),
    fill = "white",
    size = 1,
    stroke = 0.4
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      `TRUE` = 16,
      `FALSE` = 21
    ),
    labels = c(
      `TRUE` = "FDR < 0.05",
      `FALSE` = "FDR > 0.05"
    ),
    name = NULL
  ) +
  ggplot2::scale_colour_manual(
    values = c("#5496ce", "#c5373d")
  ) +
  ggplot2::labs(
    y = expression(log[2] * "(odds ratio)"),
    colour = NULL,
    shape = NULL
  ) +
  ggplot2::scale_y_continuous(
    breaks = c(-4, -2, 0, 2, 4),
    minor_breaks = c(-5, -3, -1, 1, 3, 5),
    guide = ggplot2::guide_axis(minor.ticks = TRUE)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 7, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = ggplot2::element_text(size = 7, colour = "black"),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_text(size = 7, colour = "black"),
    axis.ticks = ggplot2::element_line(linewidth = 0.24, lineend = "round", colour = "black"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(colour = "#E6E6E6", linewidth = 0.24),
    panel.grid.minor.y = ggplot2::element_line(colour = "#E6E6E6", linewidth = 0.24),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black", linewidth = 0.24, lineend = "round"),
    legend.key.height = grid::unit(0.3, "cm"),
    legend.key.width = grid::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "bottom",
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(2, 0, 0, 0)
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_segment

plot_count <- pcn_enrichment_sample |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(subtype, levels = x_order),
      y = system_count,
      fill = pcn_group
    )
  ) +
  ggplot2::geom_col(
    linewidth = 0.24,
    width = 0.7
  ) +
  ggplot2::scale_fill_manual(
    values = c("#5496ce", "#c5373d")
  ) +
  ggplot2::scale_y_log10(
    name = "No. systems",
    breaks = scales::breaks_log(n = 4),
    labels = scales::label_log(base = 10),
    expand = ggplot2::expansion(mult = c(0, 0.05))
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 7, colour = "black"),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_text(size = 7, colour = "black"),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.border = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black", linewidth = 0.24, lineend = "round"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(colour = "#E6E6E6", linewidth = 0.24),
    panel.grid.minor = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 5, 0),
    legend.position = "none"
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_count

layout <- "
A
B
"
plot_02D <- plot_count + plot_segment +
  patchwork::plot_layout(design = layout, heights = c(1, 2))

plot_02D

plot_02D |>
  ggplot2::ggsave(
    filename = "plots/fig02_D.pdf",
    width = 100,
    height = 75,
    units = "mm",
    dpi = 300
  )

# FIGURE S08 - FULL PCN DATA ---------------------------------------------------

pcn_summary <- plsdb_metadata |>
  dplyr::distinct(plasmid_seqid, pcn) |>
  dplyr::filter(!is.na(pcn)) |>
  dplyr::left_join(
    plsdb_plasmid_defense,
    by = dplyr::join_by(plasmid_seqid)
  ) |>
  dplyr::filter(!is.na(type)) |>
  dplyr::distinct(plasmid_seqid, subtype, pcn) |>
  dplyr::mutate(
    pcn_average = mean(pcn),
    pcn_median = median(pcn),
    pcn_breadth = max(pcn) - min(pcn),
    n = dplyr::n(),
    .by = subtype
  )

pcn_summary |>
  dplyr::arrange(dplyr::desc(pcn_breadth)) |>
  dplyr::filter(pcn_breadth == max(pcn_breadth)) |>
  dplyr::summarise(
    min_pcn = min(pcn), max_pcn = max(pcn), .by = c(subtype, pcn_breadth)
  )

y_order <- pcn_summary |>
  dplyr::distinct(subtype, pcn_breadth) |>
  dplyr::arrange(pcn_breadth) |>
  dplyr::pull(subtype)

plot_pcn <- pcn_summary |>
  dplyr::filter(n >= 1) |>
  ggplot2::ggplot(
    ggplot2::aes(x = pcn, y = factor(subtype, levels = y_order))
  ) +
  ggplot2::geom_boxplot(
    width = 0.5,
    linewidth = 0.24,
    colour = "black",
    outlier.size = 0.24
  ) +
  ggplot2::scale_x_log10(
    name = "PCN",
    expand = ggplot2::expansion(mult = c(0, 0)),
    breaks = c(2, 8, 32, 128, 512),
    minor_breaks = c(4, 16, 64, 256, 1024),
    limits = c(1, 1024)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "NA"
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
    panel.border = ggplot2::element_rect(linewidth = 0.24),
    panel.grid.major.x = ggplot2::element_line(linewidth = 0.24),
    panel.grid.minor.x = ggplot2::element_line(linewidth = 0.24),
    panel.grid.major.y = ggplot2::element_line(linewidth = 0.24),
    panel.grid.minor.y = ggplot2::element_blank(),
    legend.key.height = grid::unit(0.25, "cm"),
    legend.key.width = grid::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "top",
    legend.title = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank()
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::guides(colour = "none")

plot_pcn

plot_breadth <- pcn_summary |>
  dplyr::filter(n >= 1) |>
  dplyr::distinct(subtype, pcn_breadth) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = pcn_breadth + 1,
      y = factor(subtype, levels = y_order),
      fill = category
    )
  ) +
  ggplot2::geom_col(fill = "#c5cad7") +
  ggplot2::scale_x_log10(
    name = "PCN range",
    expand = ggplot2::expansion(mult = c(0, 0)),
    breaks = c(2, 8, 32, 128, 512),
    minor_breaks = c(4, 16, 64, 256, 1024),
    limits = c(1, 1024)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    axis.text.x = ggplot2::element_text(size = 7, colour = "black"),
    axis.ticks.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.border = ggplot2::element_rect(linewidth = 0.24),
    panel.grid.major.x = ggplot2::element_line(linewidth = 0.24),
    panel.grid.minor.x = ggplot2::element_line(linewidth = 0.24),
    panel.grid.major.y = ggplot2::element_line(linewidth = 0.24),
    panel.grid.minor.y = ggplot2::element_blank(),
    legend.key.height = grid::unit(0.25, "cm"),
    legend.key.width = grid::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "top",
    legend.title = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank()
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_breadth

plot_count <- pcn_summary |>
  dplyr::filter(n >= 1) |>
  dplyr::distinct(subtype, n) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = n,
      y = factor(subtype, levels = y_order),
      fill = category
    )
  ) +
  ggplot2::geom_col(fill = "#c5cad7") +
  ggplot2::scale_x_log10(
    name = "Count",
    expand = ggplot2::expansion(mult = c(0, 0)),
    breaks = c(2, 8, 32, 128, 512),
    minor_breaks = c(4, 16, 64, 256, 1024),
    limits = c(1, 1024)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    line = ggplot2::element_line(linewidth = 0.24),
    axis.text.x = ggplot2::element_text(size = 7, colour = "black"),
    axis.ticks.x = ggplot2::element_line(linewidth = 0.24, lineend = "round"),
    panel.border = ggplot2::element_rect(linewidth = 0.24),
    panel.grid.major.x = ggplot2::element_line(linewidth = 0.24),
    panel.grid.minor.x = ggplot2::element_line(linewidth = 0.24),
    panel.grid.major.y = ggplot2::element_line(linewidth = 0.24),
    panel.grid.minor.y = ggplot2::element_blank(),
    legend.key.height = grid::unit(0.25, "cm"),
    legend.key.width = grid::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "top",
    legend.title = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank()
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_count

plot_S08 <- plot_pcn + plot_breadth + plot_count +
  patchwork::plot_layout(nrow = 1, widths = c(3, 1, 1))

plot_S08

plot_S08 |>
  ggplot2::ggsave(
    filename = "plots/figS08.pdf",
    width = 182.4,
    height = 300,
    units = "mm"
  )

# LENGTH / MOBILITY / PCN CORRELATION ------------------------------------------

correlation_data <- plsdb_metadata |>
  dplyr::left_join(
    plsdb_plasmid_defense,
    by = dplyr::join_by(plasmid_seqid)
  ) |>
  dplyr::mutate(has_defense = dplyr::if_else(is.na(type), 0, 1)) |>
  dplyr::distinct(plasmid_seqid, plasmid_length, has_defense, pcn) |>
  dplyr::filter(!is.na(pcn))

cor_len_def <- cor.test(
  correlation_data$plasmid_length,
  correlation_data$has_defense,
  method = "spearman",
  exact = FALSE
)

cor_len_def

cor_len_pcn <- cor.test(
  correlation_data$plasmid_length,
  correlation_data$pcn,
  method = "spearman",
  exact = FALSE
)

cor_len_pcn

cor_pcn_def <- cor.test(
  correlation_data$pcn,
  correlation_data$has_defense,
  method = "spearman",
  exact = FALSE
)

cor_pcn_def

logistic_model <- glm(
  has_defense ~ pcn + plasmid_length + pcn:plasmid_length,
  data = correlation_data,
  family = binomial()
)

summary(logistic_model)

correlation_data_no_pcn <- plsdb_metadata_rep |>
  dplyr::left_join(
    plsdb_plasmid_defense,
    by = dplyr::join_by(plasmid_seqid)
  ) |>
  dplyr::mutate(has_defense = dplyr::if_else(is.na(type), 0, 1)) |>
  dplyr::distinct(plasmid_seqid, plasmid_length, has_defense)

cor_len_def_no_pcn <- cor.test(
  correlation_data_no_pcn$plasmid_length,
  correlation_data_no_pcn$has_defense,
  method = "spearman",
  exact = FALSE
)

cor_len_def_no_pcn

logistic_model <- glm(
  has_defense ~ plasmid_length,
  data = correlation_data_no_pcn,
  family = binomial()
)

summary(logistic_model)

# FIGURE 2E — INC GROUPS -------------------------------------------------------

inc_data <- plsdb_metadata_rep |>
  dplyr::left_join(
    plsdb_plasmid_defense,
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
  dplyr::distinct(plasmid_seqid, type, rep) |>
  tidyr::separate_longer_delim(rep, ",") |>
  dplyr::filter(!is.na(rep)) |>
  dplyr::filter(!stringr::str_detect(rep, "rep_cluster_")) |>
  dplyr::mutate(n_rep = dplyr::n_distinct(plasmid_seqid), .by = rep) |>
  dplyr::filter(n_rep >= 100) |>
  dplyr::summarise(n = dplyr::n(), .by = c(rep, n_rep, type)) |>
  dplyr::mutate(p = n / sum(n), .by = rep) |>
  dplyr::filter(!is.na(type)) |>
  dplyr::mutate(type_p_sum = sum(p), .by = type) |>
  dplyr::mutate(rep_p_sum = sum(p), .by = rep) |>
  dplyr::arrange(dplyr::desc(type_p_sum)) |>
  dplyr::mutate(top_n = dplyr::cur_group_id(), .by = type) |>
  dplyr::mutate(category = dplyr::if_else(top_n > 7, "Other", type)) |>
  dplyr::mutate(system_diversity = dplyr::n_distinct(type), .by = rep) |>
  dplyr::mutate(diversity_norm = system_diversity / n_rep)

fill_order <- inc_data |>
  dplyr::arrange(dplyr::desc(category)) |>
  dplyr::distinct(category) |>
  dplyr::filter(category != "Other") |>
  dplyr::pull(category) |>
  (\(x) c("Other", x))()

y_order <- inc_data |>
  dplyr::arrange(rep_p_sum) |>
  dplyr::distinct(rep) |>
  dplyr::pull(rep)

y_order

plot_proportion <- inc_data |>
  dplyr::summarise(p = sum(p), .by = c(rep, category)) |>
  ggplot2::ggplot(ggplot2::aes(
    x = p,
    y = factor(rep, y_order),
    fill = factor(category, fill_order)
  )) +
  ggplot2::geom_col(
    width = 0.8
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "#C5CAD7",
      "#E69F00",
      "#56B4E9",
      "#009E73",
      "#F0E442",
      "#0072B2",
      "#D55E00",
      "#CC79A7"
    )
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1),
    expand = ggplot2::expansion(),
    labels = scales::label_percent()
  ) +
  ggplot2::labs(x = "Proportion with defense", y = "Inc group") +
  theme_custom(grid = "x") +
  ggplot2::theme(
    legend.position = "bottom"
  )

plot_proportion

plot_plasmid_count <- inc_data |>
  dplyr::distinct(rep, n_rep) |>
  ggplot2::ggplot(ggplot2::aes(
    x = n_rep,
    y = factor(rep, y_order)
  )) +
  ggplot2::geom_col(
    fill = "#C5CAD7",
    width = 0.8
  ) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(),
    breaks = c(0, 3000),
    minor_breaks = c(1000, 2000)
  ) +
  ggplot2::labs(x = "No.\nplasmids", y = "Inc group") +
  theme_custom(grid = "x") +
  ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )

plot_plasmid_count

plot_system_count <- inc_data |>
  dplyr::distinct(rep, system_diversity) |>
  ggplot2::ggplot(ggplot2::aes(
    x = system_diversity,
    y = factor(rep, y_order)
  )) +
  ggplot2::geom_col(
    fill = "#C5CAD7",
    width = 0.8
  ) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(),
    breaks = c(0, 60),
    minor_breaks = c(20, 40)
  ) +
  ggplot2::labs(x = "No. system\ntypes", y = "Inc group") +
  theme_custom(grid = "x") +
  ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )

plot_system_count

plot_02E <- plot_proportion + plot_plasmid_count + plot_system_count + 
  patchwork::plot_layout(widths = c(5, 1.5, 1.5))

plot_02E

plot_02E |>
  ggplot2::ggsave(
    filename = "plots/fig02_E.pdf",
    width = 90,
    height = 83,
    units = "mm",
    dpi = 300
  )

keep <- inc_data |> 
  dplyr::filter(n_rep >= 100) |>
  dplyr::distinct(rep) |>
  dplyr::pull()

plsdb_metadata_rep |>
  dplyr::distinct(plasmid_seqid, rep) |>
  tidyr::separate_longer_delim(rep, ",") |>
  dplyr::filter(rep %in% keep) |>
  dplyr::filter(!is.na(rep)) |>
  dplyr::filter(!stringr::str_detect(rep, "rep_cluster_")) |>
  writexl::write_xlsx("data/plsdb_plasmid_inc.xlsx")


inc_data_model <- inc_data |>
  dplyr::distinct(rep, n_rep, system_diversity)

inc_log_model <- lm(system_diversity ~ log(n_rep), data = inc_data_model)

summary(inc_log_model)

inc_data_res <- inc_data_model |>
  dplyr::mutate(
    predicted_value = predict(inc_log_model),
    residual = residuals(inc_log_model),
    abs_residual = abs(residual)
  )

residual_sd <- sd(inc_data_res$residual)
threshold <- 1 * residual_sd

inc_data_res_w_out <- inc_data_res |>
  dplyr::mutate(
    outlier_label = dplyr::if_else(abs_residual > threshold, rep, "")
  )

inc_data_res_w_out

inc_data_res_w_out |>
  dplyr::filter(outlier_label != "") |>
  dplyr::arrange(residual) |>
  dplyr::select(rep, residual) |>
  purrr::pwalk(
    \(rep, residual) {
      cli::cli_alert_info(
        "{rep} residual: {round(residual)}"
      )
    }
  )

inc_data_res_w_out |>
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
  ggplot2::geom_point(size = 0.5) +
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
    plot.margin = ggplot2::margin(2, 0, 0, 0),
  ) +
  ggplot2::coord_cartesian(clip = "off")

inc_data |>
  dplyr::slice_max(order_by = p, by = rep) |>
  dplyr::slice_max(order_by = p, by = type) |>
  dplyr::select(rep, type, p) |>
  dplyr::arrange(dplyr::desc(p)) |>
  purrr::pwalk(
    \(rep, type, p) {
      cli::cli_alert_info(
        "{rep} main type: {type} ({round(p * 100)}%)"
      )
    }
  )

