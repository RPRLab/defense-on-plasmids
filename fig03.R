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

plsdb_plasmid_antidefense <- readxl::read_xlsx(
  "data/plsdb_plasmid_antidefense.xlsx"
) |>
  dplyr::mutate(
    replicon = "plasmid",
    feature = "antidefense"
  ) |>
  dplyr::mutate(
    feature_id = glue::glue(
      "{plasmid_seqid}",
      "_antidefense_",
      "{stringr::str_pad(dplyr::row_number(), 2, 'left', '0')}"
    ),
    .by = plasmid_seqid
  )

plsdb_plasmid_amr <- readxl::read_xlsx(
  "data/plsdb_plasmid_amr.xlsx"
) |>
  dplyr::mutate(
    replicon = "plasmid",
    feature = "amr"
  ) |>
  dplyr::mutate(
    feature_id = glue::glue(
      "{plasmid_seqid}",
      "_amr_",
      "{stringr::str_pad(dplyr::row_number(), 2, 'left', '0')}"
    ),
    .by = plasmid_seqid
  )

# READ GTDB METADATA -----------------------------------------------------------

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

# ANALYSIS ---------------------------------------------------------------------

plasmid_features_vs_mobility <- plsdb_plasmid_defense |>
  dplyr::bind_rows(
    plsdb_plasmid_antidefense,
    plsdb_plasmid_amr
  ) |>
  dplyr::left_join(plsdb_metadata_rep, by = dplyr::join_by(plasmid_seqid)) |>
  # Remove non-representatives
  dplyr::filter(!is.na(mobility)) |>
  dplyr::mutate(
    has_defense = dplyr::if_else(
      any(feature == "defense", na.rm = TRUE),
      TRUE,
      FALSE
    ),
    has_antidefense = dplyr::if_else(
      any(feature == "antidefense", na.rm = TRUE),
      TRUE,
      FALSE
    ),
    has_amr = dplyr::if_else(
      any(feature == "amr", na.rm = TRUE),
      TRUE,
      FALSE
    ),
    .by = plasmid_seqid
  ) |>
  dplyr::distinct(
    plasmid_seqid,
    has_defense,
    has_antidefense,
    has_amr,
    mobility
  )

plasmid_feature_categories <- plasmid_features_vs_mobility |>
  dplyr::mutate(
    category = purrr::pmap(
      list(has_defense, has_antidefense, has_amr),
      ~ c(
        if (isTRUE(..1)) "DS" else NULL,
        if (isTRUE(..2)) "ADS" else NULL,
        if (isTRUE(..3)) "AMR" else NULL
      )
    )
  ) |>
  dplyr::mutate(
    category = paste0(sort(unlist(category)), collapse = "-"),
    .by = plasmid_seqid
  )

total_population <- plsdb_metadata_rep |> nrow()
observed_population <- nrow(plasmid_feature_categories)

plasmid_feature_categories |>
  dplyr::summarise(n = dplyr::n(), .by = category) |>
  dplyr::bind_rows(
    list(category = "NA", n = total_population - observed_population)
  ) |>
  dplyr::mutate(p = round(n / total_population * 100)) |>
  dplyr::arrange(category) |>
  dplyr::select(!n) |>
  purrr::pwalk(
    \(category, p) {
      cli::cli_alert_info(
        "{category}: {p}%"
      )
    }
  )

for_venn <- plasmid_features_vs_mobility |>
  dplyr::mutate(
    category = dplyr::case_when(
      has_defense & !has_antidefense & !has_amr ~ "DS",
      !has_defense & has_antidefense & !has_amr ~ "ADS",
      !has_defense & !has_antidefense & has_amr ~ "AMR",
      has_defense & has_antidefense & !has_amr ~ "DS&ADS",
      has_defense & !has_antidefense & has_amr ~ "DS&AMR",
      !has_defense & has_antidefense & has_amr ~ "ADS&AMR",
      has_defense & has_antidefense & has_amr ~ "DS&ADS&AMR",
      TRUE ~ "NA"
    )
  ) |>
  dplyr::distinct(plasmid_seqid, category)  

region_names <- c(
  "DS",
  "ADS",
  "AMR",
  "DS&ADS",
  "DS&AMR",
  "ADS&AMR",
  "DS&ADS&AMR"
)

region_counts <- for_venn |>
  dplyr::filter(category != "NA") |>
  dplyr::count(category) |>
  tidyr::complete(
    category = region_names,
    fill = list(n = 0)
  ) |>
  dplyr::arrange(match(category, region_names))

total_population <- nrow(plsdb_metadata_rep)

observed_population <- sum(region_counts$n)

n_neither <- total_population - observed_population

disjoint_counts <- stats::setNames(
  region_counts$n,
  region_counts$category
)

fit_03A <- eulerr::euler(
  disjoint_counts,
  input = "disjoint",
  shape = "circle",
  complement = n_neither
)

region_colours <- c(
  DS           = "#87a9ca",
  ADS          = "#92c0ac",
  AMR          = "#e3c589",
  "DS&ADS"     = "#5e9d9a",
  "DS&AMR"     = "#b0a276",
  "ADS&AMR"    = "#8cac79",
  "DS&ADS&AMR" = "#729a70"
)

fit_region_order <- names(fit_03A$original.values)

ordered_region_colours <- unname(
  region_colours[fit_region_order]
)

euler_box_plot_03A <- plot(
  fit_03A,
  
  fills = list(
    mode = "disjoint",
    fill = ordered_region_colours,
    alpha = 1
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

box_mm <- 28.5
circle_diameter_mm <- 2 * box_mm / sqrt(pi)

measured_box_mm <- 26.5
scale_factor <- box_mm / measured_box_mm

grDevices::pdf(
  file = "plots/fig03_A.pdf",
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
    col = NA,
    lwd = 0
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
  euler_box_plot_03A,
  newpage = FALSE
)

grid::popViewport()
grDevices::dev.off()



x_order_features <- c(
  "ADS",
  "AMR",
  "DS",
  "ADS-AMR",
  "AMR-DS",
  "ADS-DS",
  "ADS-AMR-DS"
)

plot_03B <- plasmid_feature_categories |>
  ggplot2::ggplot(ggplot2::aes(
    x = factor(category, levels = x_order_features),
    fill = mobility
  )) +
  ggplot2::geom_bar(linewidth = 0.24) +
  ggplot2::scale_fill_manual(values = c("#E69F00", "#009E73", "#0072B2")) +
  ggupset::axis_combmatrix(sep = "-") +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.05)),
    labels = scales::label_comma()
  ) +
  ggplot2::labs(x = "", y = "No. plasmids") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_text(size = 7, colour = "black"),
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
      lineend = "round"
    ),
    panel.grid.minor.y = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round"
    ),
    legend.key.height = grid::unit(0.25, "cm"),
    legend.key.width = grid::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "none",
    legend.title = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank()
  ) +
  ggupset::theme_combmatrix(
    combmatrix.panel.point.size = 1.5,
    combmatrix.panel.line.size = 1,
    combmatrix.panel.striped_background = FALSE
  )

plot_03B

plasmid_feature_cat_stats <- plasmid_feature_categories |>
  dplyr::summarise(n = dplyr::n(), .by = category) |>
  dplyr::mutate(p = n / sum(n)) |>
  dplyr::mutate(p_total = n / total_population)

stat_main_cat <- plasmid_feature_cat_stats |>
  dplyr::slice_max(order_by = p)

stat_main_cat |>
  dplyr::select(!n) |>
  purrr::pwalk(
    \(category, p, p_total) {
      cli::cli_alert_info(
        "Most common category: {category} at {round(p*100)}% ({round(p_total*100)}% total)"
      )
    }
  )

stat_all_feat <- plasmid_feature_cat_stats |>
  dplyr::filter(category == "ADS-AMR-DS")
  
stat_all_feat |>
  dplyr::select(!n) |>
  purrr::pwalk(
    \(category, p, p_total) {
      cli::cli_alert_info(
        "{category} at {round(p*100)}% ({round(p_total*100)}% total)"
      )
    }
  )

stat_all_feat_mob <- plasmid_feature_categories |>
  dplyr::summarise(n = dplyr::n(), .by = c(category, mobility)) |>
  dplyr::filter(category == "ADS-AMR-DS") |>
  dplyr::mutate(p = n / sum(n)) |>
  dplyr::slice_max(order_by = p) |>
  dplyr::select(mobility, p)

stat_all_feat_mob |>
  purrr::pwalk(
    \(mobility, p) {
      cli::cli_alert_info(
        "Main ADS-AMR-DS mobility: {mobility} at {round(p*100)}%"
      )
    }
  )

human_samples <- plsdb_metadata_rep |>
  dplyr::filter(biosample_host_processed == "Homo sapiens")

human_samples_status <- human_samples |>
  dplyr::mutate(
    clinical_status = dplyr::case_when(
      is.na(biosample_host_disease) ~ "non-clinical",
      stringr::str_detect(
        stringr::str_to_lower(stringr::str_squish(biosample_host_disease)),
        paste0(
          "^(",
          paste(
            c(
              "na",
              "null",
              "missing",
              "no data",
              "not applicable",
              "not available",
              "not collected",
              "not defined",
              "not reported",
              "unknown",
              "unkonwn",
              "unkown",
              "uncertain",
              "healthy",
              "no disease",
              "asymptomatic",
              "asymptomatic carriage",
              "carriage",
              "carrier",
              "human carrier",
              "disease carrier",
              "colonisation",
              "colonization",
              "colonized",
              "gut colonization",
              "intestinal colonization",
              "rectal carriage",
              "cre colonization",
              "screen",
              "screening",
              "screening sample",
              "none/screening",
              "verification for resistant bacteria",
              "contaminant"
            ),
            collapse = "|"
          ),
          ")$"
        )
      ) ~ "non-clinical",
      TRUE ~ "clinical"
    )
  ) |>
  dplyr::select(
    plasmid_seqid,
    biosample_host_processed,
    biosample_host_disease,
    clinical_status
  )

clinical_samples <- plasmid_feature_categories |>
  dplyr::left_join(human_samples_status, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::filter(clinical_status == "clinical")

clinical_categories <- clinical_samples |>
  dplyr::mutate(
    disease_processed = stringr::str_to_sentence(biosample_host_disease)
  ) |>
  dplyr::summarise(n = dplyr::n(), .by = disease_processed) |>
  dplyr::mutate(
    disease_category = dplyr::case_when(
      stringr::str_detect(
        disease_processed,
        stringr::regex(
          "diarrhea|gastroenteritis|enteritis|Salmonella|Shigellosis|Enterococcus faecium|Diverticular disease of intestine|Gastric ulcers|Salmonellosis",
          ignore_case = TRUE
        )
      ) ~
        "Gastrointestinal",
      stringr::str_detect(
        disease_processed,
        stringr::regex(
          "UTI|urinary tract infection|urosepsis|cystitis|renal|kidney|Uropathogenic escherichia coli|Urine with mucoviscosity phenotype",
          ignore_case = TRUE
        )
      ) ~
        "Urinary tract",
      stringr::str_detect(
        disease_processed,
        stringr::regex(
          "bacteremia|sepsis|septicemia|bloodstream infection|Bsi|Bacterimia|Septicaemia|cardiovascular|heart|endocarditis|aortic",
          ignore_case = TRUE
        )
      ) ~
        "Bloodstream",
      stringr::str_detect(
        disease_processed,
        stringr::regex(
          "pneumonia|respiratory|lung|bronchial|pulmonary|Typhoid|Cough|Bronchiectasis|Trachaeal secretion",
          ignore_case = TRUE
        )
      ) ~
        "Respiratory",
      stringr::str_detect(
        disease_processed,
        stringr::regex(
          "meningitis|cerebral|intracranial|brain|Hemiplegia",
          ignore_case = TRUE
        )
      ) ~
        "Other condition",
      stringr::str_detect(
        disease_processed,
        stringr::regex(
          "peritonitis|intra-abdominal|biliary|cholecystitis|pancreatitis|appendicitis|intestinal|Bile duct obstruction|hepatitis|liver|hepatocellular|Ascites|Acute pyelonephritis|Hepatolithiasis",
          ignore_case = TRUE
        )
      ) ~
        "Intra-abdominal",
      stringr::str_detect(
        disease_processed,
        stringr::regex(
          "cancer|leukemia|carcinoma|myelodysplastic|Aml|Meningiomas",
          ignore_case = TRUE
        )
      ) ~
        "Other condition",
      stringr::str_detect(
        disease_processed,
        stringr::regex(
          "Fasciitis, necrotizing|infection|bacterial infection|carrier|colonization|screen|surveillance|Infenction|Vre|Leg ulcer|Escherichia coli \\[e. coli\\] as the cause of diseases classified elsewhere|Empyema thoracis|Colonized|Carriage|Infectious disease|Abscess|Colonisation",
          ignore_case = TRUE
        )
      ) ~
        "Other infection",
      .default = "Other condition"
    )
  ) |>
  dplyr::summarise(n = sum(n), .by = disease_category) |>
  dplyr::arrange(dplyr::desc(n)) |>
  dplyr::mutate(p = n / sum(n))

clinical_categories |>
  dplyr::select(!n) |>
  purrr::pwalk(
    \(disease_category, p) {
      cli::cli_alert_info(
        "{disease_category}: {round(p*100)}%"
      )
    }
  )

plot_03C <- clinical_samples |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(category, levels = x_order_features),
      fill = mobility
    )
  ) +
  ggplot2::geom_bar(linewidth = 0.24) +
  ggplot2::scale_fill_manual(values = c("#E69F00", "#009E73", "#0072B2")) +
  ggupset::axis_combmatrix(sep = "-") +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.05)),
    labels = scales::label_comma()
  ) +
  ggplot2::labs(x = "", y = "No. plasmids") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_text(size = 7, colour = "black"),
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
      lineend = "round"
    ),
    panel.grid.minor.y = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round"
    ),
    legend.key.height = grid::unit(0.25, "cm"),
    legend.key.width = grid::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "none",
    legend.title = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank()
  ) +
  ggupset::theme_combmatrix(
    combmatrix.panel.point.size = 1.5,
    combmatrix.panel.line.size = 1,
    combmatrix.panel.striped_background = FALSE
  )

plot_03C

stat_main_group <- clinical_samples |>
  dplyr::summarise(n = dplyr::n(), .by = c(category, mobility)) |>
  dplyr::mutate(p = n / sum(n)) |>
  dplyr::arrange(dplyr::desc(p)) |>
  dplyr::select(!n) |>
  dplyr::slice_head()

stat_main_group |>
  purrr::pwalk(
    \(category, mobility, p) {
      cli::cli_alert_info(
        "Most common category: {category} on {mobility} plasmids at {round(p*100)}%"
      )
    }
  )

clinical_conj_hatrick_species <-  plasmid_feature_categories |>
    dplyr::left_join(human_samples_status, by = dplyr::join_by(plasmid_seqid)) |>
    dplyr::filter(clinical_status == "clinical" & mobility == "conjugative" & category == "ADS-AMR-DS") |>
    dplyr::left_join(plsdb_metadata[c("plasmid_seqid", "taxonomy_species")], by = dplyr::join_by(plasmid_seqid)) |>
    dplyr::mutate(
      species = stringr::str_replace_all(taxonomy_species, "_", " ")
    ) |>
    dplyr::summarise(n = dplyr::n(), .by = species) |>
    dplyr::arrange(dplyr::desc(n), species) |>
    dplyr::mutate(
      species = dplyr::if_else(dplyr::row_number() >= 7, "Other species", species)
    ) |>
    dplyr::summarise(n = sum(n), .by = species) |>
    dplyr::mutate(p = n / sum(n))
    
clinical_conj_hatrick_species

clinical_conj_hatrick_species |>
  dplyr::select(!n) |>
  purrr::pwalk(
    \(species, p) {
      cli::cli_alert_info(
        "{species}: {round(p*100)}%"
      )
    }
  )

all_plasmid_species <- plsdb_metadata_rep |>
  dplyr::mutate(
    species = stringr::str_replace_all(taxonomy_species, "_", " ")
  ) |>
  dplyr::summarise(n = dplyr::n(), .by = species) |>
  dplyr::arrange(dplyr::desc(n), species) |>
  dplyr::summarise(n = sum(n), .by = species) |>
  dplyr::mutate(p = n / sum(n)) |>
  dplyr::filter(species %in% clinical_conj_hatrick_species$species) |>
  dplyr::bind_rows(c(species = "Other species")) |>
  dplyr::mutate(p = dplyr::if_else(is.na(p), 1 - sum(p, na.rm = TRUE), p))

all_plasmid_species

species_comparison <- clinical_conj_hatrick_species |>
  dplyr::mutate(category = "Conjugative + ADS-AMR-DS + clinical") |>
  dplyr::bind_rows(all_plasmid_species) |>
  dplyr::mutate(
    category = dplyr::if_else(is.na(category), "All plasmids", category)
  ) |>
  dplyr::mutate(
    species = species |>
      stringr::str_replace(
        pattern = "^([A-Za-z])[a-zA-Z]*\\s+",
        replacement = "\\1. "
      ) |>
      stringr::str_replace("O. species", "Other species")
  )

x_order_species <- species_comparison |>
  dplyr::filter(category == "Conjugative + ADS-AMR-DS + clinical") |>
  dplyr::pull(species)

plot_03D <- species_comparison |>
  dplyr::filter(species != "Other species") |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(species, levels = x_order_species),
      y = p,
      fill = category
    )
  ) +
  ggplot2::geom_col(linewidth = 0.24, position = "dodge") +
  ggplot2::scale_fill_manual(values = c("#c5cad7", "#E69F00")) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.05)),
    labels = scales::percent_format()
  ) +
  ggplot2::labs(x = "", y = "Proportion of plasmids") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_text(size = 7, colour = "black"),
    axis.text.x = ggplot2::element_text(
      size = 7,
      colour = "black",
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
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
      lineend = "round"
    ),
    panel.grid.minor.y = ggplot2::element_line(
      linewidth = 0.24,
      lineend = "round"
    ),
    legend.key.height = grid::unit(0.25, "cm"),
    legend.key.width = grid::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "none",
    legend.title = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank()
  )

plot_03D

plot_fig03BCD <- plot_03B +
  plot_03B +
  plot_03C +
  plot_03D +
  patchwork::plot_layout(nrow = 1, widths = c(1, 1, 1, 1))

plot_fig03BCD

plot_fig03BCD |>
  ggplot2::ggsave(
    filename = "plots/fig03_BCD.pdf",
    width = 182.4,
    height = 50,
    units = "mm",
    dpi = 300
  )

# PLSDB DEFENSE/AMR CO-OCCURRENCE ----------------------------------------------

drug_class_mechanisms <- tibble::tibble(
  class = c(
    "AMINOGLYCOSIDE",
    "AVILAMYCIN",
    "FUSIDIC ACID",
    "LINCOSAMIDE",
    "MACROLIDE",
    "MUPIROCIN",
    "OXAZOLIDINONE",
    "PHENICOL",
    "PLEUROMUTILIN",
    "STREPTOGRAMIN",
    "TETRACYCLINE",
    "THIOSTREPTON",
    "BETA-LACTAM",
    "GLYCOPEPTIDE",
    "FOSFOMYCIN",
    "TUBERACTINOMYCIN",
    "COLISTIN",
    "STREPTOTHRICIN",
    "QUINOLONE",
    "NITROIMIDAZOLE",
    "RIFAMYCIN",
    "SULFONAMIDE",
    "TRIMETHOPRIM",
    "BLEOMYCIN"
  ),
  mechanism = c(
    "Inhibitor of Protein Synthesis",
    "Inhibitor of Protein Synthesis",
    "Inhibitor of Protein Synthesis",
    "Inhibitor of Protein Synthesis",
    "Inhibitor of Protein Synthesis",
    "Inhibitor of Protein Synthesis",
    "Inhibitor of Protein Synthesis",
    "Inhibitor of Protein Synthesis",
    "Inhibitor of Protein Synthesis",
    "Inhibitor of Protein Synthesis",
    "Inhibitor of Protein Synthesis",
    "Inhibitor of Protein Synthesis",
    "Inhibitor of Cell Wall Synthesis",
    "Inhibitor of Cell Wall Synthesis",
    "Inhibitor of Cell Wall Synthesis",
    "Inhibitor of Cell Wall Synthesis",
    "Disruptor of Cell Membrane Integrity",
    "Disruptor of Cell Membrane Integrity",
    "Inhibitor of DNA Synthesis or Function",
    "Inhibitor of DNA Synthesis or Function",
    "Inhibitor of RNA Synthesis",
    "Inhibitor of Folate Synthesis",
    "Inhibitor of Folate Synthesis",
    "Inducer of DNA Strand Breaks"
  )
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

plasmid_amr_types <- plsdb_metadata_rep |>
  dplyr::select(plasmid_seqid) |>
  dplyr::left_join(plsdb_plasmid_amr, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::filter(!is.na(feature)) |>
  dplyr::distinct(plasmid_seqid, type) |>
  dplyr::mutate(present = 1L) |>
  dplyr::rename(entity = type)

plasmid_amr_subtypes <- plsdb_metadata_rep |>
  dplyr::select(plasmid_seqid) |>
  dplyr::left_join(plsdb_plasmid_amr, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::filter(!is.na(feature)) |>
  dplyr::distinct(plasmid_seqid, subtype) |>
  dplyr::mutate(present = 1L) |>
  dplyr::rename(entity = subtype)

plasmid_amr_classes <- plasmid_amr_types |>
  dplyr::mutate(entity = stringr::str_remove(entity, ":.*")) |>
  dplyr::distinct(plasmid_seqid, entity, present)

amr_present <- plsdb_metadata_rep |>
  dplyr::left_join(plsdb_plasmid_amr, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::mutate(amr_present = dplyr::if_else(!is.na(feature), 1L, 0L)) |>
  dplyr::distinct(plasmid_seqid, amr_present)

defense_present <- plsdb_metadata_rep |>
  dplyr::left_join(plsdb_plasmid_defense, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::mutate(defense_present = dplyr::if_else(!is.na(feature), 1L, 0L)) |>
  dplyr::distinct(plasmid_seqid, defense_present)

amr_vs_defense <- amr_present |>
  dplyr::left_join(defense_present, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::count(amr_present, defense_present) |>
  tidyr::pivot_wider(
    names_from = defense_present,
    values_from = n,
    values_fill = 0
  ) |>
  tibble::column_to_rownames("amr_present") |>
  as.matrix()

amr_vs_defense_chisq <- chisq.test(amr_vs_defense)

amr_vs_defense_phi <- sqrt(unname(amr_vs_defense_chisq$statistic) / sum(amr_vs_defense))

cli::cli_alert_info(
  "χ²: {round(amr_vs_defense_chisq$statistic, 1)}; p: {amr_vs_defense_chisq$p.value}; φ: {round(amr_vs_defense_phi, 3)}"
)

# AFFINITY: DEFENSE TYPE / AMR CLASS -------------------------------------------

def_type_amr_class_matrix <- plasmid_defense_types |>
  dplyr::bind_rows(plasmid_amr_classes) |>
  tidyr::pivot_wider(
    id_cols = plasmid_seqid,
    names_from = entity,
    values_from = present,
    values_fill = 0L
  ) |>
  tibble::column_to_rownames("plasmid_seqid") |>
  as.matrix()

if (!file.exists("data/plsdb_defense_type_amr_class_affinity.xlsx")) {
  def_type_amr_class_affinity <- def_type_amr_class_matrix |>
    CooccurrenceAffinity::affinity(row.or.col = "col", squarematrix = c("all"))
  
  def_type_amr_class_affinity_filt <- def_type_amr_class_affinity$all |>
    tibble::as_tibble() |>
    dplyr::filter(
      entity_1 %in%
        unique(plasmid_defense_types$entity) &
        entity_2 %in%
        unique(
          stringr::str_remove(plasmid_amr_classes$entity, "^.*:")
        ) |
        entity_2 %in%
        unique(plasmid_defense_types$entity) &
        entity_1 %in%
        unique(
          stringr::str_remove(plasmid_amr_classes$entity, "^.*:")
        )
    )
  
  writexl::write_xlsx(
    def_type_amr_class_affinity_filt,
    "data/plsdb_defense_type_amr_class_affinity.xlsx"
  )
}

def_type_amr_class_affinity_filt <- readxl::read_xlsx(
  "data/plsdb_defense_type_amr_class_affinity.xlsx"
)

stat_def_type_amr_class_affinity_filt <- def_type_amr_class_affinity_filt |>
  dplyr::mutate(
    is_correlated = dplyr::if_else(alpha_mle > 0, TRUE, FALSE),
    is_significant = dplyr::if_else(p_value <= 0.05, TRUE, FALSE)
  ) |>
  dplyr::summarise(n = dplyr::n(), .by = c(is_correlated, is_significant)) |>
  dplyr::mutate(p = n / sum(n))

def_type_amr_class_affinity_for_plot <- def_type_amr_class_affinity_filt |>
  dplyr::mutate(p_value = as.double(p_value)) |>
  dplyr::filter(entity_1_count_mA >= 50 & entity_2_count_mB >= 50) |>
  dplyr::mutate(
    alpha_mle = dplyr::case_when(
      obs_cooccur_X == 0 ~ NA,
      .default = alpha_mle
    )
  )

fill_limit <- max(
  abs(max(def_type_amr_class_affinity_for_plot$alpha_mle, na.rm = TRUE)),
  abs(min(def_type_amr_class_affinity_for_plot$alpha_mle, na.rm = TRUE))
)

y_order <- rev(colnames(def_type_amr_class_affinity$occur_mat[-1])) |>
  intersect(def_type_amr_class_affinity_for_plot$entity_2)

x_order <- colnames(
  def_type_amr_class_affinity$occur_mat
)[
  -length(colnames(def_type_amr_class_affinity$occur_mat))
] |>
  intersect(def_type_amr_class_affinity_for_plot$entity_1)

plot_def_type_amr_class_affinity <- def_type_amr_class_affinity_for_plot |>
  ggplot2::ggplot(ggplot2::aes(x = entity_1, y = entity_2, fill = alpha_mle)) +
  ggplot2::geom_tile(colour = "#575653") +
  ggplot2::geom_point(
    ggplot2::aes(
      size = dplyr::if_else(dplyr::between(p_value, 0, 0.05), "dot", "nodot")
    ),
    colour = "#575653"
  ) +
  ggplot2::scale_size_manual(
    values = c(dot = 0.24, nodot = -1),
    guide = "none"
  ) +
  ggplot2::coord_fixed() +
  ggplot2::scale_fill_gradientn(
    colours = c(
      "#026eae",
      "#026eae",
      "#5496ce",
      "#9bcae9",
      "#c5e5fb",
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

plot_def_type_amr_class_affinity

plot_def_type_amr_class_affinity |>
  ggplot2::ggsave(
    filename = "plots/figS10_A.pdf",
    width = 182.4,
    height = 182.4,
    units = "mm"
  )

# AFFINITY: DEFENSE SUBTYPE / AMR CLASS ----------------------------------------

def_subtype_amr_class_matrix <- plasmid_defense_subtypes |>
  dplyr::bind_rows(plasmid_amr_classes) |>
  tidyr::pivot_wider(
    id_cols = plasmid_seqid,
    names_from = entity,
    values_from = present,
    values_fill = 0L
  )|>
  tibble::column_to_rownames("plasmid_seqid") |>
  as.matrix()

if (!file.exists("data/plsdb_defense_subtype_amr_class_affinity.xlsx")) {
  def_subtype_amr_class_affinity <- def_subtype_amr_class_matrix |>
    CooccurrenceAffinity::affinity(row.or.col = "col", squarematrix = c("all"))
  
  def_subtype_amr_class_affinity_filt <- def_subtype_amr_class_affinity$all |>
    tibble::as_tibble() |>
    dplyr::filter(
      entity_1 %in%
        unique(plasmid_defense_subtypes$entity) &
        entity_2 %in% unique(plasmid_amr_classes$entity) |
        entity_2 %in%
        unique(plasmid_defense_subtypes$entity) &
        entity_1 %in% unique(plasmid_amr_classes$entity)
    )
  
  writexl::write_xlsx(
    def_subtype_amr_class_affinity_filt,
    "data/plsdb_defense_subtype_amr_class_affinity.xlsx"
  )
}

def_subtype_amr_class_affinity_filt <- readxl::read_xlsx(
  "data/plsdb_defense_subtype_amr_class_affinity.xlsx"
)

stat_def_subtype_amr_class_affinity_filt <- def_subtype_amr_class_affinity_filt |>
  dplyr::mutate(
    is_correlated = dplyr::if_else(alpha_mle > 0, TRUE, FALSE),
    is_significant = dplyr::if_else(p_value <= 0.05, TRUE, FALSE)
  ) |>
  dplyr::summarise(n = dplyr::n(), .by = c(is_correlated, is_significant)) |>
  dplyr::mutate(p = n / sum(n)) |>
  dplyr::arrange(dplyr::desc(is_significant), dplyr::desc(is_correlated))

def_subtype_amr_class_affinity_for_plot <- def_subtype_amr_class_affinity_filt |>
  dplyr::mutate(p_value = as.double(p_value)) |>
  dplyr::filter(entity_1_count_mA >= 50 & entity_2_count_mB >= 50) |>
  dplyr::mutate(
    alpha_mle = dplyr::case_when(
      obs_cooccur_X == 0 ~ NA,
      .default = alpha_mle
    )
  )

fill_limit <- max(
  abs(max(def_subtype_amr_class_affinity_for_plot$alpha_mle, na.rm = TRUE)),
  abs(min(def_subtype_amr_class_affinity_for_plot$alpha_mle, na.rm = TRUE))
)

y_order <- rev(colnames(def_subtype_amr_class_affinity$occur_mat[-1])) |>
  intersect(def_subtype_amr_class_affinity_for_plot$entity_2)

x_order <- colnames(
  def_subtype_amr_class_affinity$occur_mat
)[
  -length(colnames(def_subtype_amr_class_affinity$occur_mat))
] |>
  intersect(def_subtype_amr_class_affinity_for_plot$entity_1)

plot_def_subtype_amr_class_affinity <- def_subtype_amr_class_affinity_for_plot |>
  ggplot2::ggplot(ggplot2::aes(x = entity_1, y = entity_2, fill = alpha_mle)) +
  ggplot2::geom_tile(colour = "#575653") +
  ggplot2::geom_point(
    ggplot2::aes(
      size = dplyr::if_else(dplyr::between(p_value, 0, 0.05), "dot", "nodot")
    ),
    colour = "#575653"
  ) +
  ggplot2::scale_size_manual(
    values = c(dot = 0.24, nodot = -1),
    guide = "none"
  ) +
  ggplot2::coord_fixed() +
  ggplot2::scale_fill_gradientn(
    colours = c(
      "#026eae",
      "#026eae",
      "#5496ce",
      "#9bcae9",
      "#c5e5fb",
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

plot_def_subtype_amr_class_affinity

plot_def_subtype_amr_class_affinity |>
  ggplot2::ggsave(
    filename = "plots/figS10_B.pdf",
    width = 182.4,
    height = 182.4,
    units = "mm"
  )

# AFFINITY: DEFENSE TYPE / AMR TYPE --------------------------------------------

def_type_amr_type_matrix <- plasmid_defense_types |>
  dplyr::bind_rows(plasmid_amr_types) |>
  dplyr::mutate(entity = stringr::str_remove(entity, "^.*:")) |>
  tidyr::pivot_wider(
    id_cols = plasmid_seqid,
    names_from = entity,
    values_from = present,
    values_fill = 0L
  ) |>
  tibble::column_to_rownames("plasmid_seqid") |>
  as.matrix()

if (!file.exists("data/plsdb_defense_type_amr_type_affinity.xlsx")) {
  def_type_amr_type_affinity <- def_type_amr_type_matrix |>
    CooccurrenceAffinity::affinity(row.or.col = "col", squarematrix = c("all"))
  
  def_type_amr_type_affinity_filt <- def_type_amr_type_affinity$all |>
    tibble::as_tibble() |>
    dplyr::filter(
      entity_1 %in%
        unique(plasmid_defense_types$entity) &
        entity_2 %in%
        unique(
          stringr::str_remove(plasmid_amr_types$entity, "^.*:")
        ) |
        entity_2 %in%
        unique(plasmid_defense_types$entity) &
        entity_1 %in%
        unique(
          stringr::str_remove(plasmid_amr_types$entity, "^.*:")
        )
    )
  
  writexl::write_xlsx(
    def_type_amr_type_affinity_filt,
    "data/plsdb_defense_type_amr_type_affinity.xlsx"
  )
}

def_type_amr_type_affinity_filt <- readxl::read_xlsx(
  "data/plsdb_defense_type_amr_type_affinity.xlsx"
)

stat_def_type_amr_type_affinity_filt <- def_type_amr_type_affinity_filt |>
  dplyr::mutate(
    is_correlated = dplyr::if_else(alpha_mle > 0, TRUE, FALSE),
    is_significant = dplyr::if_else(p_value <= 0.05, TRUE, FALSE)
  ) |>
  dplyr::summarise(n = dplyr::n(), .by = c(is_correlated, is_significant)) |>
  dplyr::mutate(p = n / sum(n)) |>
  dplyr::arrange(dplyr::desc(is_significant), dplyr::desc(is_correlated))

def_type_amr_type_affinity_for_plot <- def_type_amr_type_affinity_filt |>
  dplyr::mutate(p_value = as.double(p_value)) |>
  dplyr::filter(entity_1_count_mA >= 50 & entity_2_count_mB >= 50) |>
  dplyr::mutate(
    alpha_mle = dplyr::case_when(
      obs_cooccur_X == 0 ~ NA,
      .default = alpha_mle
    )
  )

fill_limit <- max(
  abs(max(def_type_amr_type_affinity_for_plot$alpha_mle, na.rm = TRUE)),
  abs(min(def_type_amr_type_affinity_for_plot$alpha_mle, na.rm = TRUE))
)

y_order <- rev(colnames(def_type_amr_type_affinity$occur_mat[-1])) |>
  intersect(def_type_amr_type_affinity_for_plot$entity_2)

x_order <- colnames(
  def_type_amr_type_affinity$occur_mat
)[
  -length(colnames(def_type_amr_type_affinity$occur_mat))
] |>
  intersect(def_type_amr_type_affinity_for_plot$entity_1)

plot_def_type_amr_type_affinity <- def_type_amr_type_affinity_for_plot |>
  ggplot2::ggplot(ggplot2::aes(x = entity_1, y = entity_2, fill = alpha_mle)) +
  ggplot2::geom_tile(colour = "#575653") +
  ggplot2::geom_point(
    ggplot2::aes(
      size = dplyr::if_else(dplyr::between(p_value, 0, 0.05), "dot", "nodot")
    ),
    colour = "#575653"
  ) +
  ggplot2::scale_size_manual(
    values = c(dot = 0.24, nodot = -1),
    guide = "none"
  ) +
  ggplot2::coord_fixed() +
  ggplot2::scale_fill_gradientn(
    colours = c(
      "#026eae",
      "#026eae",
      "#5496ce",
      "#9bcae9",
      "#c5e5fb",
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

plot_def_type_amr_type_affinity

plot_def_type_amr_type_affinity |>
  ggplot2::ggsave(
    filename = "plots/figS10_C.pdf",
    width = 182.4,
    height = 182.4,
    units = "mm"
  )

# AFFINITY: DEFENSE SUBTYPE / AMR TYPE -----------------------------------------

def_subtype_amr_type_matrix <- plasmid_defense_subtypes |>
  dplyr::bind_rows(plasmid_amr_types) |>
  dplyr::mutate(entity = stringr::str_remove(entity, "^.*:")) |>
  tidyr::pivot_wider(
    id_cols = plasmid_seqid,
    names_from = entity,
    values_from = present,
    values_fill = 0L
  ) |>
  tibble::column_to_rownames("plasmid_seqid") |>
  as.matrix()

if (!file.exists("data/plsdb_defense_subtype_amr_type_affinity.xlsx")) {
  def_subtype_amr_type_affinity <- def_subtype_amr_type_matrix |>
    CooccurrenceAffinity::affinity(row.or.col = "col", squarematrix = c("all"))
  
  def_subtype_amr_type_affinity_filt <- def_subtype_amr_type_affinity$all |>
    tibble::as_tibble() |>
    dplyr::filter(
      entity_1 %in%
        unique(plasmid_defense_subtypes$entity) &
        entity_2 %in%
        unique(
          stringr::str_remove(plasmid_amr_types$entity, "^.*:")
        ) |
        entity_2 %in%
        unique(plasmid_defense_subtypes$entity) &
        entity_1 %in%
        unique(
          stringr::str_remove(plasmid_amr_types$entity, "^.*:")
        )
    )
  
  writexl::write_xlsx(
    def_subtype_amr_type_affinity_filt,
    "data/plsdb_defense_subtype_amr_type_affinity.xlsx"
  )
}

def_subtype_amr_type_affinity_filt <- readxl::read_xlsx(
  "data/plsdb_defense_subtype_amr_type_affinity.xlsx"
)

stat_def_subtype_amr_type_affinity_filt <- def_subtype_amr_type_affinity_filt |>
  dplyr::mutate(
    is_correlated = dplyr::if_else(alpha_mle > 0, TRUE, FALSE),
    is_significant = dplyr::if_else(p_value <= 0.05, TRUE, FALSE)
  ) |>
  dplyr::summarise(n = dplyr::n(), .by = c(is_correlated, is_significant)) |>
  dplyr::mutate(p = n / sum(n)) |>
  dplyr::arrange(dplyr::desc(is_significant), dplyr::desc(is_correlated))

def_subtype_amr_type_affinity_for_plot <- def_subtype_amr_type_affinity_filt |>
  dplyr::mutate(p_value = as.double(p_value)) |>
  dplyr::filter(entity_1_count_mA >= 50 & entity_2_count_mB >= 50) |>
  dplyr::mutate(
    alpha_mle = dplyr::case_when(
      obs_cooccur_X == 0 ~ NA,
      .default = alpha_mle
    )
  )

fill_limit <- max(
  abs(max(def_subtype_amr_type_affinity_for_plot$alpha_mle, na.rm = TRUE)),
  abs(min(def_subtype_amr_type_affinity_for_plot$alpha_mle, na.rm = TRUE))
)

y_order <- rev(colnames(def_subtype_amr_type_affinity$occur_mat[-1])) |>
  intersect(def_subtype_amr_type_affinity_for_plot$entity_2)

x_order <- colnames(
  def_subtype_amr_type_affinity$occur_mat
)[
  -length(colnames(def_subtype_amr_type_affinity$occur_mat))
] |>
  intersect(def_subtype_amr_type_affinity_for_plot$entity_1)

plot_def_subtype_amr_type_affinity <- def_subtype_amr_type_affinity_for_plot |>
  ggplot2::ggplot(ggplot2::aes(x = entity_1, y = entity_2, fill = alpha_mle)) +
  ggplot2::geom_tile(colour = "#575653") +
  ggplot2::geom_point(
    ggplot2::aes(
      size = dplyr::if_else(dplyr::between(p_value, 0, 0.05), "dot", "nodot")
    ),
    colour = "#575653"
  ) +
  ggplot2::scale_size_manual(
    values = c(dot = 0.24, nodot = -1),
    guide = "none"
  ) +
  ggplot2::coord_fixed() +
  ggplot2::scale_fill_gradientn(
    colours = c(
      "#026eae",
      "#026eae",
      "#5496ce",
      "#9bcae9",
      "#c5e5fb",
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

plot_def_subtype_amr_type_affinity

plot_def_subtype_amr_type_affinity |>
  ggplot2::ggsave(
    filename = "plots/figS10_D.pdf",
    width = 182.4,
    height = 182.4,
    units = "mm"
  )

# ASSOCIATION SUMMARY ----------------------------------------------------------

stat_affinity <- dplyr::bind_rows(
  stat_def_type_amr_class_affinity_filt |> dplyr::mutate(category = "Defense type vs AMR class"),
  stat_def_subtype_amr_class_affinity_filt |> dplyr::mutate(category = "Defense subtype vs AMR class"),
  stat_def_type_amr_type_affinity_filt |> dplyr::mutate(category = "Defense type vs AMR type"),
  stat_def_subtype_amr_type_affinity_filt |> dplyr::mutate(category = "Defense subtype vs AMR type")
)

stat_combinations <- stat_affinity |>
  dplyr::summarise(n = sum(n), .by = category)

stat_affinity_significant <- stat_affinity |>
  dplyr::filter(is_significant == TRUE) |>
  dplyr::arrange(is_correlated, dplyr::desc(n)) |>
  dplyr::select(!is_significant)

stat_most_positive_cat <- stat_affinity_significant |>
  dplyr::filter(is_correlated == TRUE) |>
  dplyr::slice_head()

stat_most_positive_comb <- stat_combinations |>
  dplyr::filter(category == stat_most_positive_cat[["category"]]) |>
  dplyr::pull(n)

stat_most_positive_cat |>
  purrr::pwalk(
    \(is_correlated, n, p, category) {
      cli::cli_alert_info(
        paste("Up to {n} significant positive associations",
        "({round(n / stat_most_positive_comb*100)}% of all {stat_most_positive_cat[['category']]})"
        )
      )
    }
  )

stat_most_negative_cat <- stat_affinity_significant |>
  dplyr::filter(is_correlated == FALSE) |>
  dplyr::slice_head()

stat_most_negative_comb <- stat_combinations |>
  dplyr::filter(category == stat_most_negative_cat[["category"]]) |>
  dplyr::pull(n)

stat_most_negative_cat |>
  purrr::pwalk(
    \(is_correlated, n, p, category) {
      cli::cli_alert_info(
        paste("Up to {n} significant negative associations",
              "({round(n / stat_most_negative_comb*100)}% of all {stat_most_negative_cat[['category']]})"
        )
      )
    }
  )

class_all_neg <- def_subtype_amr_class_affinity_filt |>
  dplyr::filter(entity_1_count_mA >= 50) |>
  dplyr::mutate(p_value = as.double(p_value)) |> 
  dplyr::filter(p_value < 0.05) |>
  dplyr::mutate(
    all_neg = dplyr::case_when(
      all(alpha_mle < 0) ~ TRUE, .default = FALSE
    ), .by = entity_1
  ) |>
  dplyr::filter(all_neg == TRUE) |>
  dplyr::distinct(entity_1) |>
  dplyr::pull()

type_all_neg <- def_subtype_amr_type_affinity_filt |>
  dplyr::filter(entity_1_count_mA >= 50) |>
  dplyr::mutate(p_value = as.double(p_value)) |> 
  dplyr::filter(p_value < 0.05) |>
  dplyr::mutate(
    all_neg = dplyr::case_when(
      all(alpha_mle < 0) ~ TRUE, .default = FALSE
    ), .by = entity_1
  ) |>
  dplyr::filter(all_neg == TRUE) |>
  dplyr::distinct(entity_1) |>
  dplyr::pull()

intersect_all_neg <- intersect(class_all_neg, type_all_neg)

intersect_all_neg |> length()

def_subtype_amr_class_affinity_filt |>
  dplyr::filter(entity_1 %in% intersect_all_neg) |>
  dplyr::summarise(sum_alpha_mle = sum(alpha_mle), .by = entity_1) |>
  dplyr::arrange(sum_alpha_mle)

# pNDM-MAR CASE STUDY ----------------------------------------------------------

pNDM_MAR <- plsdb_metadata_rep |>
  dplyr::filter(stringr::str_detect(plasmidfinder, "pNDM-MAR"))

pNDM_MAR_accessions <- pNDM_MAR |>
  dplyr::distinct(plasmid_seqid)

vector <- pNDM_MAR_accessions |> dplyr::pull()
chunk_size <- 20
splits <- seq(1, length(vector), by = chunk_size)
chunks <- purrr::map(
  .x = splits,
  .f = ~ vector[.x:min(.x + chunk_size - 1, length(vector))]
)

seqs_list <- purrr::map(
  .x = chunks,
  .f = function(.x) {
    post <- rentrez::entrez_post(db = "nuccore", id = .x, rettype = "fasta")
    fetch <- rentrez::entrez_fetch(
      db = "nuccore",
      rettype = "fasta",
      web_history = post
    )
    fetch
  },
  .progress = TRUE
)

seqs <- unlist(seqs_list) |>
  stringr::str_split(">") |>
  unlist() |>
  as.list() |>
  purrr::discard(~ . == "") |>
  purrr::map(~ paste0(">", .))

names <- stringr::str_extract(seqs, "(?<=>)[^\\s]+")

names(seqs) <- names

outdir <- "data/pNDM-Mar/"

fs::dir_create(outdir)

purrr::iwalk(seqs, ~ write(.x, file = paste0(outdir, .y, ".fna")))



read_skani_matrix <- function(path) {
  readr::read_delim(
    path,
    delim = "\t",
    skip = 1,
    col_names = FALSE,
    show_col_types = FALSE
  ) |>
    tibble::column_to_rownames("X1") |>
    as.matrix()
}

pNDM_MAR_skani_matrix <- read_skani_matrix(
  "data/pNDM-Mar_skani_matrix.txt"
)

pNDM_MAR_skani_af_matrix <- read_skani_matrix(
  "data/pNDM-Mar_skani_matrix.txt.af"
)

colnames(pNDM_MAR_skani_matrix) <- rownames(pNDM_MAR_skani_matrix)
colnames(pNDM_MAR_skani_af_matrix) <- rownames(pNDM_MAR_skani_af_matrix)

pNDM_MAR_dist <- stats::as.dist(1 - pNDM_MAR_skani_matrix / 100)
pNDM_MAR_hc <- hclust(pNDM_MAR_dist)
pNDM_MAR_tree <- tidytree::as.phylo(pNDM_MAR_hc)
pNDM_MAR_tree_rooted <- phytools::midpoint_root(pNDM_MAR_tree)

plot_pNDM_MAR_tree <- ggtree::ggtree(pNDM_MAR_tree_rooted, branch.length = "none")

pNDM_MAR_y_order <- plot_pNDM_MAR_tree$data |>
  dplyr::arrange(y) |>
  dplyr::filter(!is.na(label)) |>
  dplyr::pull(label) |>
  stringr::str_remove_all(".fna")

pNDM_MAR_ds <- pNDM_MAR |>
  dplyr::left_join(
    plsdb_plasmid_defense,
    by = dplyr::join_by(plasmid_seqid)
  ) |>
  dplyr::select(plasmid_seqid, subtype) |>
  dplyr::mutate(present = TRUE) |>
  tidyr::complete(plasmid_seqid, subtype, fill = list(present = FALSE)) |>
  dplyr::filter(!is.na(subtype))

pNDM_MAR_ads <- pNDM_MAR |>
  dplyr::left_join(
    plsdb_plasmid_antidefense,
    by = dplyr::join_by(plasmid_seqid)
  ) |>
  dplyr::select(plasmid_seqid, subtype) |>
  dplyr::mutate(present = TRUE) |>
  tidyr::complete(plasmid_seqid, subtype, fill = list(present = FALSE)) |>
  dplyr::filter(!is.na(subtype)) |>
  dplyr::mutate(
    subtype = dplyr::case_when(
      subtype == "ardb_klca_acric11" ~ "ArdB/KlcA/AcrIC11",
      subtype == "psiab" ~ "PsiB",
      subtype == "arda_ardu" ~ "ArdA/ArdU",
      subtype == "acrie9" ~ "AcrIE9",
      subtype == "ardc" ~ "ArdC",
      subtype == "ardk" ~ "ArdK",
      subtype == "ardr" ~ "ArdR",
      subtype == "vcrx089_090" ~ "VCRx 089/090",
      subtype == "vcrx091_093" ~ "VCRx 091/093",
      subtype == "ardc" ~ "ArdC"
    )
  )

pNDM_MAR_amr <- pNDM_MAR |>
  dplyr::left_join(plsdb_plasmid_amr, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::mutate(class = stringr::str_remove(type, ":.*")) |>
  dplyr::mutate(class = stringr::str_to_title(class)) |>
  dplyr::select(plasmid_seqid, class) |>
  dplyr::mutate(present = TRUE) |>
  tidyr::complete(plasmid_seqid, class, fill = list(present = FALSE)) |>
  dplyr::filter(!is.na(class)) |>
  dplyr::mutate(
    class = dplyr::case_when(
      class == "Aminoglycoside" ~ "AG",
      class == "Beta-Lactam" ~ "BL",
      class == "Sulfonamide" ~ "SUL",
      class == "Tetracycline" ~ "TET",
      class == "Phenicol" ~ "CHL",
      class == "Macrolide" ~ "MAC",
      class == "Quinolone" ~ "QNL",
      class == "Aminoglycoside/Quinolone" ~ "AG/QNL",
      class == "Rifamycin" ~ "RIF",
      class == "Trimethoprim" ~ "TMP",
      class == "Bleomycin" ~ "BLE",
      class == "Colistin" ~ "COL",
      class == "Macrolide/Streptogramin" ~ "MAC/SGM",
      class == "Phenicol/Quinolone" ~ "CHL/QNL",
      class == "Fosfomycin" ~ "FOS",
      class == "Lincosamide" ~ "LIN",
      class == "Streptothricin" ~ "STN"
    )
  )

pNDM_MAR_amr_subtype <- pNDM_MAR |>
  dplyr::left_join(plsdb_plasmid_amr, by = dplyr::join_by(plasmid_seqid)) |>
  dplyr::select(plasmid_seqid, subtype) |>
  dplyr::mutate(present = TRUE) |>
  tidyr::complete(plasmid_seqid, subtype, fill = list(present = FALSE))|>
  dplyr::filter(!is.na(subtype))

pNDM_MAR_position <- pNDM_MAR_y_order |>
  tibble::as_tibble() |>
  dplyr::mutate(position = dplyr::row_number())

pNDM_MAR_host <- pNDM_MAR |>
  dplyr::select(plasmid_seqid, host_acc, dplyr::starts_with("taxonomy_")) |>
  dplyr::mutate(id = stringr::str_extract(host_acc, "[0-9]{9}")) |>
  dplyr::left_join(gtdb_taxonomy, by = dplyr::join_by(id))

pNDM_MAR_location <- pNDM_MAR |>
  dplyr::select(plasmid_seqid, loc_lat, loc_lng) |>
  dplyr::left_join(
    pNDM_MAR_position,
    by = dplyr::join_by(plasmid_seqid == value)
  ) |>
  dplyr::mutate(dplyr::across(c(loc_lat, loc_lng), ~ as.double(.))) |>
  bdc::bdc_country_from_coordinates(lat = "loc_lat", lon = "loc_lng") |>
  dplyr::mutate(country = dplyr::if_else(is.na(country), "—", country)) |>
  dplyr::left_join(
    dplyr::bind_rows(
      dplyr::mutate(pNDM_MAR_ds, feature = "defense") |>
        dplyr::filter(present == TRUE),
      dplyr::mutate(pNDM_MAR_ads, feature = "antidefense") |>
        dplyr::filter(present == TRUE),
      dplyr::mutate(pNDM_MAR_amr, feature = "amr") |>
        dplyr::filter(present == TRUE)
    ),
    by = dplyr::join_by(plasmid_seqid)
  ) |>
  dplyr::summarise(
    n = dplyr::n(),
    .by = c(
      plasmid_seqid,
      loc_lat,
      loc_lng,
      position,
      country,
      decimalLatitude,
      decimalLongitude,
      feature
    )
  ) |>
  tidyr::pivot_wider(names_from = feature, values_from = n) |>
  dplyr::select(-`NA`) |>
  dplyr::mutate(
    sum = sum(defense, antidefense, amr, na.rm = TRUE),
    .by = plasmid_seqid
  )

plot_pNDM_MAR_name <- pNDM_MAR_location |>
  ggplot2::ggplot() +
  ggplot2::geom_text(
    ggplot2::aes(
      x = 0,
      y = factor(plasmid_seqid, pNDM_MAR_y_order),
      label = plasmid_seqid
    ),
    hjust = 0,
    size = 5 / 2.845,
    colour = "black"
  ) +
  ggplot2::scale_x_continuous(limits = c(0, 1)) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    legend.position = "NA",
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )

plot_pNDM_MAR_country <- pNDM_MAR_location |>
  ggplot2::ggplot() +
  ggplot2::geom_text(
    ggplot2::aes(
      x = 0,
      y = factor(plasmid_seqid, pNDM_MAR_y_order),
      label = country
    ),
    hjust = 0,
    size = 5 / 2.845,
    colour = "black"
  ) +
  ggplot2::scale_x_continuous(limits = c(0, 1)) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    legend.position = "NA",
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )


plot_pNDM_MAR_ds_x_order <- pNDM_MAR_ds |>
  dplyr::filter(present == TRUE) |>
  dplyr::summarise(n = dplyr::n(), .by = c(subtype)) |>
  dplyr::arrange(dplyr::desc(n), subtype) |>
  dplyr::pull(subtype)

plot_pNDM_MAR_ds <- pNDM_MAR_ds |>
  ggplot2::ggplot(ggplot2::aes(
    x = factor(subtype, plot_pNDM_MAR_ds_x_order),
    y = factor(plasmid_seqid, pNDM_MAR_y_order),
    fill = present
  )) +
  ggplot2::geom_tile(colour = "#e5e5e9") +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::scale_fill_manual(values = c("white", "#595854")) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      size = 5,
      colour = "black",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_line(colour = "black", linewidth = 0.24),
    axis.ticks.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(linewidth = 0.24),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    legend.key.height = ggplot2::unit(0.1, "cm"),
    legend.key.width = ggplot2::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "NA",
    legend.title.position = "right",
    legend.ticks.length = ggplot2::unit(-0.2, 'cm'),
    legend.ticks = ggplot2::element_line(colour = "black", linewidth = 0.24),
    legend.frame = ggplot2::element_rect(colour = "black", linewidth = 0.24),
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )

plot_pNDM_MAR_amr_x_order <- pNDM_MAR_amr |>
  dplyr::filter(present == TRUE) |>
  dplyr::summarise(n = dplyr::n(), .by = c(class)) |>
  dplyr::arrange(dplyr::desc(n), class) |>
  dplyr::pull(class)

plot_pNDM_MAR_amr <- pNDM_MAR_amr |>
  ggplot2::ggplot(ggplot2::aes(
    x = factor(class, plot_pNDM_MAR_amr_x_order),
    y = factor(plasmid_seqid, pNDM_MAR_y_order),
    fill = present
  )) +
  ggplot2::geom_tile(colour = "#e5e5e9") +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::scale_fill_manual(values = c("white", "#595854")) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      size = 5,
      colour = "black",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_line(colour = "black", linewidth = 0.24),
    axis.ticks.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(linewidth = 0.24),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    legend.key.height = ggplot2::unit(0.1, "cm"),
    legend.key.width = ggplot2::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "NA",
    legend.title.position = "right",
    legend.ticks.length = ggplot2::unit(-0.2, 'cm'),
    legend.ticks = ggplot2::element_line(colour = "black", linewidth = 0.24),
    legend.frame = ggplot2::element_rect(colour = "black", linewidth = 0.24),
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )

plot_pNDM_MAR_ads_x_order <- pNDM_MAR_ads |>
  dplyr::filter(present == TRUE) |>
  dplyr::summarise(n = dplyr::n(), .by = c(subtype)) |>
  dplyr::arrange(dplyr::desc(n), subtype) |>
  dplyr::pull(subtype)

plot_pNDM_MAR_ads <- pNDM_MAR_ads |>
  ggplot2::ggplot(ggplot2::aes(
    x = factor(subtype, plot_pNDM_MAR_ads_x_order),
    y = factor(plasmid_seqid, pNDM_MAR_y_order),
    fill = present
  )) +
  ggplot2::geom_tile(colour = "#e5e5e9") +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::scale_fill_manual(values = c("white", "#595854")) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      size = 5,
      colour = "black",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_line(colour = "black", linewidth = 0.24),
    axis.ticks.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(linewidth = 0.24),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    legend.key.height = ggplot2::unit(0.1, "cm"),
    legend.key.width = ggplot2::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "NA",
    legend.title.position = "right",
    legend.ticks.length = ggplot2::unit(-0.2, 'cm'),
    legend.ticks = ggplot2::element_line(colour = "black", linewidth = 0.24),
    legend.frame = ggplot2::element_rect(colour = "black", linewidth = 0.24),
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )

plot_pNDM_MAR_tree +
  plot_pNDM_MAR_name +
  plot_pNDM_MAR_amr +
  plot_pNDM_MAR_ds +
  plot_pNDM_MAR_ads +
  plot_pNDM_MAR_country +
  patchwork::plot_layout(nrow = 1, widths = c(0.02, 0.05, NA, NA, NA, 0.2))

ggplot2::ggsave(
  "plots/figS11.pdf",
  width = 20,
  height = 30,
  limitsize = FALSE
)

# QUANTITATIVE ANALYSIS --------------------------------------------------------

pNDM_MAR_ds |> 
  dplyr::filter(present == TRUE) |> 
  nrow() |>
  (\(n) cli::cli_alert_info("{n} defense systems in pNDM-MAR-like plasmids"))()

pNDM_MAR_ds |>
  dplyr::filter(present == TRUE) |>
  dplyr::filter(!is.na(subtype)) |>
  dplyr::distinct(subtype) |>
  nrow() |>
  (\(n) cli::cli_alert_info("{n} defense subtypes in pNDM-MAR-like plasmids"))()

pNDM_MAR_amr |> 
  dplyr::filter(present == TRUE) |> 
  nrow() |>
  (\(n) cli::cli_alert_info("{n} AMR systems in pNDM-MAR-like plasmids"))()

pNDM_MAR_amr_subtype |>
  dplyr::filter(present == TRUE) |>
  dplyr::filter(!is.na(subtype)) |>
  dplyr::distinct(subtype) |>
  nrow() |>
  (\(n) cli::cli_alert_info("{n} AMR subtypes in pNDM-MAR-like plasmids"))()

pNDM_MAR_ads |> 
  dplyr::filter(present == TRUE) |> 
  nrow() |>
  (\(n) cli::cli_alert_info("{n} Anti-defense systems in pNDM-MAR-like plasmids"))()

pNDM_MAR_ads |>
  dplyr::filter(present == TRUE) |>
  dplyr::filter(!is.na(subtype)) |>
  dplyr::distinct(subtype) |>
  nrow() |>
  (\(n) cli::cli_alert_info("{n} Anti-defense subtypes in pNDM-MAR-like plasmids"))()



pNDM_MAR_ani_cutoff <- 95
pNDM_MAR_reciprocal_af_cutoff <- 70

pair_indices <- which(upper.tri(pNDM_MAR_skani_matrix), arr.ind = TRUE)

pair_similarity_df <- tibble::tibble(
  plasmid1 = rownames(pNDM_MAR_skani_matrix)[pair_indices[, "row"]],
  plasmid2 = colnames(pNDM_MAR_skani_matrix)[pair_indices[, "col"]],
  ani = pNDM_MAR_skani_matrix[pair_indices],
  af1 = pNDM_MAR_skani_af_matrix[pair_indices],
  af2 = pNDM_MAR_skani_af_matrix[cbind(pair_indices[, "col"], pair_indices[, "row"])]
) |>
  dplyr::mutate(
    reciprocal_af = pmin(af1, af2),
    ani_dist = 1 - ani / 100
  )

binary_jaccard_distance <- function(plasmid1, plasmid2, repertoire_matrix) {
  repertoire1 <- as.logical(repertoire_matrix[plasmid1, , drop = TRUE])
  repertoire2 <- as.logical(repertoire_matrix[plasmid2, , drop = TRUE])
  union_size <- sum(repertoire1 | repertoire2)
  if (union_size == 0) {
    return(NA_real_)
  }
  1 - sum(repertoire1 & repertoire2) / union_size
}

ds_mat <- pNDM_MAR_ds |>
  dplyr::mutate(plasmid_seqid = paste0(plasmid_seqid, ".fna")) |>
  dplyr::distinct() |>
  tidyr::pivot_wider(
    id_cols = plasmid_seqid,
    names_from = subtype,
    values_from = present,
    values_fill = FALSE
  ) |>
  tibble::column_to_rownames("plasmid_seqid") |>
  as.matrix()

ds_mat <- ds_mat[rownames(pNDM_MAR_skani_matrix), , drop = FALSE]

pair_df_ds <- pair_similarity_df |>
  dplyr::mutate(
    defense_n1 = rowSums(ds_mat[plasmid1, , drop = FALSE]),
    defense_n2 = rowSums(ds_mat[plasmid2, , drop = FALSE]),
    jaccard = purrr::map2_dbl(
      plasmid1,
      plasmid2,
      \(p1, p2) binary_jaccard_distance(p1, p2, ds_mat)
    )
  ) |>
  dplyr::filter(defense_n1 + defense_n2 > 0)

pair_df_ds_gt_99 <- pair_df_ds |>
  dplyr::filter(
    ani >= pNDM_MAR_ani_cutoff,
    reciprocal_af >= pNDM_MAR_reciprocal_af_cutoff
  )

stat_pairs_gt_99_w_ds <- nrow(pair_df_ds_gt_99)

cli::cli_alert_info(
  paste0(
    "{stat_pairs_gt_99_w_ds} pairs of pNDM-Mar-like plasmids with ANI ≥ ",
    "{pNDM_MAR_ani_cutoff}% and reciprocal AF ≥ ",
    "{pNDM_MAR_reciprocal_af_cutoff}%, where at least one has defense"
  )
)

jacc_ds_summary <- summary(pair_df_ds_gt_99$jaccard)

cli::cli_alert_info(
  paste0(
    "μ Jaccard = {round(jacc_ds_summary[['Mean']], 2)}; ",
    "μ1/2 Jaccard = {round(jacc_ds_summary[['Median']], 2)}; ",
    "IQR = {round(jacc_ds_summary[['1st Qu.']], 2)}–",
    "{round(jacc_ds_summary[['3rd Qu.']], 2)}"
  )
)

stat_diff_ds <- sum(pair_df_ds_gt_99$jaccard == 1)
stat_diff_ds_p <- mean(pair_df_ds_gt_99$jaccard == 1)

cli::cli_alert_info(
  paste0(
    "{stat_diff_ds} pairs had completely non-overlapping defense repertoires ",
    "({round(stat_diff_ds_p * 100)}%)"
  )
)


amr_mat <- pNDM_MAR_amr_subtype |>
  dplyr::mutate(plasmid_seqid = paste0(plasmid_seqid, ".fna")) |>
  dplyr::distinct() |>
  tidyr::pivot_wider(
    id_cols = plasmid_seqid,
    names_from = subtype,
    values_from = present,
    values_fill = FALSE
  ) |>
  tibble::column_to_rownames("plasmid_seqid") |>
  as.matrix()

amr_mat <- amr_mat[rownames(pNDM_MAR_skani_matrix), , drop = FALSE]

pair_df_amr <- pair_similarity_df |>
  dplyr::mutate(
    amr_n1 = rowSums(amr_mat[plasmid1, , drop = FALSE]),
    amr_n2 = rowSums(amr_mat[plasmid2, , drop = FALSE]),
    jaccard = purrr::map2_dbl(
      plasmid1,
      plasmid2,
      \(p1, p2) binary_jaccard_distance(p1, p2, amr_mat)
    )
  ) |>
  dplyr::filter(amr_n1 + amr_n2 > 0)

pair_df_amr_gt_99 <- pair_df_amr |>
  dplyr::filter(
    ani >= pNDM_MAR_ani_cutoff,
    reciprocal_af >= pNDM_MAR_reciprocal_af_cutoff
  )

stat_pairs_gt_99_w_amr <- nrow(pair_df_amr_gt_99)

cli::cli_alert_info(
  paste0(
    "{stat_pairs_gt_99_w_amr} pairs of pNDM-Mar-like plasmids with ANI ≥ ",
    "{pNDM_MAR_ani_cutoff}% and reciprocal AF ≥ ",
    "{pNDM_MAR_reciprocal_af_cutoff}%, where at least one has AMR"
  )
)

jacc_amr_summary <- summary(pair_df_amr_gt_99$jaccard)

cli::cli_alert_info(
  paste0(
    "μ Jaccard = {round(jacc_amr_summary[['Mean']], 2)}; ",
    "μ1/2 Jaccard = {round(jacc_amr_summary[['Median']], 2)}; ",
    "IQR = {round(jacc_amr_summary[['1st Qu.']], 2)}–",
    "{round(jacc_amr_summary[['3rd Qu.']], 2)}"
  )
)

stat_diff_amr <- sum(pair_df_amr_gt_99$jaccard == 1)
stat_diff_amr_p <- mean(pair_df_amr_gt_99$jaccard == 1)

cli::cli_alert_info(
  paste0(
    "{stat_diff_amr} pairs had completely non-overlapping AMR repertoires ",
    "({round(stat_diff_amr_p * 100)}%)"
  )
)

# WORLD MAP ---

pNDM_MAR_sample <- c(
  "CP110748.1",
  "NZ_CP055163.1",
  "NZ_CP091327.1",
  "NZ_CP020848.1",
  "NZ_CP050068.1",
  "NZ_CP072906.1",
  "NZ_CP084985.1",
  "NZ_CP132632.1",
  "CP132283.1",
  "NZ_CP113193.1",
  "NZ_CP055011.1",
  "CP113221.1",
  "LR890179.1",
  "NZ_CP103666.1",
  "NZ_CP083933.1",
  "NZ_CP115151.1",
  "NZ_CP097229.1",
  "NZ_CP089098.1",
  "NZ_CP067932.1",
  "NZ_CP100078.1"
)

data <- pNDM_MAR_location |>
  dplyr::select(plasmid_seqid, loc_lat, loc_lng)

missing <- data |>
  dplyr::filter(is.na(loc_lat)) |>
  dplyr::pull(plasmid_seqid)

tree <- pNDM_MAR_tree_rooted

tree$tip.label <- stringr::str_remove(tree$tip.label, "\\.fna$")

tree_pruned <- tree |>
  tidytree::drop.tip(missing)

cladogram_layout <- ggtree::fortify(
  tree_pruned,
  branch.length = "none"
)

node_x <- setNames(
  cladogram_layout$x,
  cladogram_layout$node
)

tree_pruned$edge.length <- unname(
  node_x[as.character(tree_pruned$edge[, 2])] -
    node_x[as.character(tree_pruned$edge[, 1])]
)

data_pruned <- data |>
  dplyr::filter(!is.na(loc_lat)) |>
  tibble::column_to_rownames(var = "plasmid_seqid") |>
  dplyr::select(loc_lat, loc_lng) |>
  as.matrix()

cols <- setNames(
  ifelse(
    tree_pruned$tip.label %in% pNDM_MAR_sample,
    "#c5373d",
    grDevices::adjustcolor("#8791a0", alpha.f = 0.18)
  ),
  tree_pruned$tip.label
)

phy_map <- phytools::phylo.to.map(
  tree_pruned,
  data_pruned,
  plot = FALSE,
  direction = "downwards"
)

pdf(
  "plots/fig03_E_map.pdf",
  width = 3.06,
  height = 3.06
)

par(mar = c(0.2, 0.2, 0.2, 0.2))

plot(
  phy_map,
  direction = "downwards",
  ftype = "off",
  split = c(0.20, 0.80),
  colors = cols,
  lty = "solid",
  lwd = 0.55,
  from.tip = TRUE,
  psize = 0.7,
  map.bg = "#BEC3D1",
  xlim = c(-200, 200),
  ylim = c(-100, 100)
)

dev.off()

ape::plot.phylo(
  phy_map$tree,
  direction = "downwards",
  show.tip.label = FALSE,
  plot = FALSE
)

lp <- get(
  "last_plot.phylo",
  envir = ape::.PlotPhyloEnv
)

n_tip <- ape::Ntip(phy_map$tree)

tip_order <- phy_map$tree$tip.label[
  order(lp$xx[seq_len(n_tip)])
]

tip_order


pNDM_MAR_ds_filt <- pNDM_MAR_ds |>
  dplyr::filter(plasmid_seqid %in% pNDM_MAR_sample)

pNDM_MAR_amr_filt <- pNDM_MAR_amr |>
  dplyr::filter(plasmid_seqid %in% pNDM_MAR_sample)

pNDM_MAR_ads_filt <- pNDM_MAR_ads |>
  dplyr::filter(plasmid_seqid %in% pNDM_MAR_sample)

pNDM_MAR_tree_filt_y_order <- tip_order |>
  purrr::keep(~ .x %in% pNDM_MAR_sample) |>
  rev()

plot_pNDM_MAR_ds_filt_x_order <- pNDM_MAR_ds_filt |>
  dplyr::filter(present == TRUE) |>
  dplyr::summarise(n = dplyr::n(), .by = c(subtype)) |>
  dplyr::arrange(dplyr::desc(n), subtype) |>
  dplyr::pull(subtype)

plot_pNDM_MAR_ds_filt <- pNDM_MAR_ds_filt |>
  dplyr::filter(subtype %in% plot_pNDM_MAR_ds_filt_x_order) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(subtype, plot_pNDM_MAR_ds_filt_x_order),
      y = factor(plasmid_seqid, pNDM_MAR_tree_filt_y_order),
      fill = present
    )
  ) +
  ggplot2::geom_tile(colour = "#e5e5e9") +
  ggplot2::coord_fixed(ratio = 1, expand = FALSE, clip = "off") +
  ggplot2::scale_fill_manual(values = c("white", "#595854")) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      size = 5,
      colour = "black",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_line(colour = "black", linewidth = 0.24),
    axis.ticks.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(linewidth = 0.24),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    legend.key.height = ggplot2::unit(0.1, "cm"),
    legend.key.width = ggplot2::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "NA",
    legend.title.position = "right",
    legend.ticks.length = ggplot2::unit(-0.2, 'cm'),
    legend.ticks = ggplot2::element_line(colour = "black", linewidth = 0.24),
    legend.frame = ggplot2::element_rect(colour = "black", linewidth = 0.24),
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )

plot_pNDM_MAR_ds_filt

plot_pNDM_MAR_amr_filt_x_order <- pNDM_MAR_amr_filt |>
  dplyr::filter(present == TRUE) |>
  dplyr::summarise(n = dplyr::n(), .by = c(class)) |>
  dplyr::arrange(dplyr::desc(n), class) |>
  dplyr::pull(class)

plot_pNDM_MAR_amr_filt <- pNDM_MAR_amr_filt |>
  dplyr::filter(class %in% plot_pNDM_MAR_amr_filt_x_order) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(class, plot_pNDM_MAR_amr_filt_x_order),
      y = factor(plasmid_seqid, pNDM_MAR_tree_filt_y_order),
      fill = present
    )
  ) +
  ggplot2::geom_tile(colour = "#e5e5e9") +
  ggplot2::coord_fixed(ratio = 1, expand = FALSE, clip = "off") +
  ggplot2::scale_fill_manual(values = c("white", "#595854")) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      size = 5,
      colour = "black",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_line(colour = "black", linewidth = 0.24),
    axis.ticks.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(linewidth = 0.24),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    legend.key.height = ggplot2::unit(0.1, "cm"),
    legend.key.width = ggplot2::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "NA",
    legend.title.position = "right",
    legend.ticks.length = ggplot2::unit(-0.2, 'cm'),
    legend.ticks = ggplot2::element_line(colour = "black", linewidth = 0.24),
    legend.frame = ggplot2::element_rect(colour = "black", linewidth = 0.24),
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )

plot_pNDM_MAR_amr_filt

plot_pNDM_MAR_ads_filt_x_order <- pNDM_MAR_ads_filt |>
  dplyr::filter(present == TRUE) |>
  dplyr::summarise(n = dplyr::n(), .by = c(subtype)) |>
  dplyr::arrange(dplyr::desc(n), subtype) |>
  dplyr::pull(subtype)

plot_pNDM_MAR_ads_filt <- pNDM_MAR_ads_filt |>
  dplyr::filter(subtype %in% plot_pNDM_MAR_ads_filt_x_order) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(subtype, plot_pNDM_MAR_ads_filt_x_order),
      y = factor(plasmid_seqid, pNDM_MAR_tree_filt_y_order),
      fill = present
    )
  ) +
  ggplot2::geom_tile(colour = "#e5e5e9") +
  ggplot2::coord_fixed(ratio = 1, expand = FALSE, clip = "off") +
  ggplot2::scale_fill_manual(values = c("white", "#595854")) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      size = 5,
      colour = "black",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_line(colour = "black", linewidth = 0.24),
    axis.ticks.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(linewidth = 0.24),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    legend.key.height = ggplot2::unit(0.1, "cm"),
    legend.key.width = ggplot2::unit(0.5, "cm"),
    legend.text = ggplot2::element_text(size = 7, colour = "black"),
    legend.position = "NA",
    legend.title.position = "right",
    legend.ticks.length = ggplot2::unit(-0.2, 'cm'),
    legend.ticks = ggplot2::element_line(colour = "black", linewidth = 0.24),
    legend.frame = ggplot2::element_rect(colour = "black", linewidth = 0.24),
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )

plot_pNDM_MAR_ads_filt

pNDM_MAR_location_filt <- pNDM_MAR |>
  dplyr::filter(plasmid_seqid %in% pNDM_MAR_sample) |>
  dplyr::select(plasmid_seqid, loc_lat, loc_lng) |>
  dplyr::mutate(dplyr::across(c(loc_lat, loc_lng), ~ as.double(.))) |>
  bdc::bdc_country_from_coordinates(lat = "loc_lat", lon = "loc_lng") |>
  dplyr::mutate(country = dplyr::if_else(plasmid_seqid == "NZ_CP084985.1", "Italy", country))

plot_pNDM_MAR_country_filt <- pNDM_MAR_location_filt |>
  ggplot2::ggplot() +
  ggplot2::geom_text(
    ggplot2::aes(
      x = 0,
      y = factor(plasmid_seqid, pNDM_MAR_tree_filt_y_order),
      label = country
    ),
    hjust = 0,
    size = 5 / 2.845,
    colour = "black"
  ) +
  ggplot2::scale_x_continuous(limits = c(0, 1)) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    legend.position = "NA",
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )

plot_pNDM_MAR_country_filt

plot_pNDM_MAR_name_filt <- pNDM_MAR_location_filt |>
  ggplot2::ggplot() +
  ggplot2::geom_text(
    ggplot2::aes(
      x = 0,
      y = factor(plasmid_seqid, pNDM_MAR_tree_filt_y_order),
      label = plasmid_seqid
    ),
    hjust = 0,
    size = 5 / 2.845,
    colour = "black"
  ) +
  ggplot2::scale_x_continuous(limits = c(0, 1)) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    legend.position = "NA",
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )

pNDM_MAR_host_filt <- pNDM_MAR_host |>
  dplyr::filter(plasmid_seqid %in% pNDM_MAR_sample)

plot_pNDM_MAR_host_filt <- pNDM_MAR_host_filt |>
  dplyr::mutate(
    label = stringr::str_replace_all(taxonomy_species, "_", " ")
  ) |>
  dplyr::mutate(
    label = label |>
      stringr::str_replace(
        pattern = "^([A-Za-z]{1})[a-zA-Z]*\\s+",
        replacement = "\\1. "
      )
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = "",
      y = factor(plasmid_seqid, pNDM_MAR_tree_filt_y_order),
      fill = label
    )
  ) +
  ggplot2::geom_tile(colour = "#e5e5e9") +
  ggplot2::scale_fill_manual(
    values = c(
      "#e69f00",
      "#57b4e9",
      "#009e73",
      "#e5e5e9",
      "#f0e443",
      "#d55e00"
    )
  ) +
  ggplot2::coord_fixed(ratio = 1, expand = FALSE, clip = "off") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      size = 5,
      colour = "black",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = ggplot2::element_blank(),
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
    legend.title.position = "right",
    legend.ticks.length = ggplot2::unit(-0.2, 'cm'),
    legend.ticks = ggplot2::element_line(colour = "black", linewidth = 0.24),
    legend.frame = ggplot2::element_rect(colour = "black", linewidth = 0.24),
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )

plot_pNDM_MAR_host_filt <- plot_pNDM_MAR_host_filt +
  ggplot2::theme(
    legend.position = "none"
  )

plot_pNDM_MAR_source_filt <- pNDM_MAR |>
  dplyr::filter(plasmid_seqid %in% pNDM_MAR_sample) |>
  dplyr::mutate(
    isolation_source = stringr::str_to_sentence(biosample_isolation_source)
  ) |>
  dplyr::mutate(
    isolation_source = dplyr::case_when(
      stringr::str_detect(isolation_source, "effluent") ~ "Effluent",
      stringr::str_detect(isolation_source, "liver abscess") ~ "Liver abcess",
      stringr::str_detect(isolation_source, "Intraabdominal abscess") ~
        "Intra-abdominal",
      stringr::str_detect(isolation_source, "Chicken intestinal contents") ~
        "Chicken",
      stringr::str_detect(isolation_source, "Conjugation") ~ "—",
      stringr::str_detect(isolation_source, "Surveillance") ~ "—",
      stringr::str_detect(isolation_source, "Freshwater") ~ "Freshwater",
      stringr::str_detect(isolation_source, "Human fecal") ~ "Feces",
      stringr::str_detect(isolation_source, "Hospital") ~ "Hospital",
      stringr::str_detect(isolation_source, "Manis") ~ "Pangolin",
      stringr::str_detect(isolation_source, "Sewage sludge") ~ "Sewage",
      stringr::str_detect(isolation_source, "Throat swab") ~ "Throat",
      stringr::str_detect(isolation_source, "Unknown") ~ "—",
      .default = isolation_source
    )
  ) |>
  dplyr::mutate(
    isolation_source = dplyr::if_else(
      is.na(isolation_source),
      "—",
      isolation_source
    )
  ) |>
  ggplot2::ggplot() +
  ggplot2::geom_text(
    ggplot2::aes(
      x = 0,
      y = factor(plasmid_seqid, pNDM_MAR_tree_filt_y_order),
      label = isolation_source
    ),
    hjust = 0,
    size = 5 / 2.845,
    colour = "black"
  ) +
  ggplot2::scale_x_continuous(limits = c(0, 1)) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    legend.position = "NA",
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )

plot_fig03_E_heatmap <-
  plot_pNDM_MAR_name_filt +
  plot_pNDM_MAR_host_filt +
  plot_pNDM_MAR_amr_filt +
  plot_pNDM_MAR_ds_filt +
  plot_pNDM_MAR_ads_filt +
  plot_pNDM_MAR_country_filt +
  plot_pNDM_MAR_source_filt +
  patchwork::plot_layout(nrow = 1, widths = c(4.5, NA, NA, NA, NA, 4, 3))

plot_fig03_E_heatmap

plot_fig03_E_heatmap |>
  ggplot2::ggsave(
    filename = "plots/fig03_E_heatmap.pdf",
    width = 140,
    height = 130,
    units = "mm",
    dpi = 300
  )
