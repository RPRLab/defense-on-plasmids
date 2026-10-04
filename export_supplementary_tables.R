dir.create("supplementary_tables", showWarnings = FALSE)

# Supplementary Table 1
plsdb_columns <- names(readxl::read_xlsx("data/plsdb_plasmid-host_metadata_with_taxonomy.xlsx", n_max = 0))
readxl::read_xlsx(
  "data/plsdb_plasmid-host_metadata_with_taxonomy.xlsx",
  col_types = ifelse(plsdb_columns == "pcn", "numeric", "guess")
) |>
  dplyr::select(
    plasmid_seqid, host_acc, representative, infomap_representative,
    plasmid_length, host_length,
    gtdb_domain, gtdb_phylum, taxonomy_genus, taxonomy_species,
    mobility, plasmidfinder, rep, pcn,
    loc_lat, loc_lng, biosample_isolation_source,
    biosample_host_processed, biosample_host_disease
  ) |>
  writexl::write_xlsx("supplementary_tables/table_S01_plsdb_plasmid-host_metadata_with_taxonomy.xlsx")

# Supplementary Table 2
list(Sheet7 = readxl::read_xlsx("data/defense_system_unification.xlsx", sheet = "Sheet7") |>
  dplyr::distinct()) |>
  writexl::write_xlsx("supplementary_tables/table_S02_defense_system_unification.xlsx")

# Supplementary Tables 3, 4, 9, 10, 15 and 16
c(
  table_S03_plsdb_plasmid_defense.xlsx = "plsdb_plasmid_defense.xlsx",
  table_S04_plsdb_host_defense.xlsx = "plsdb_host_defense.xlsx",
  table_S09_plsdb_plasmid_amr.xlsx = "plsdb_plasmid_amr.xlsx",
  table_S10_plsdb_plasmid_antidefense.xlsx = "plsdb_plasmid_antidefense.xlsx",
  table_S15_imgpr_plasmid_defense.xlsx = "imgpr_plasmid_defense.xlsx",
  table_S16_imgpr_plasmid_amr.xlsx = "imgpr_plasmid_amr.xlsx"
) |>
  purrr::iwalk(\(file, table) file.copy(
    file.path("data", file), file.path("supplementary_tables", table), overwrite = TRUE
  ))

# Supplementary Tables 5, 6 and 11–14
c(
  table_S05_plsdb_defense_type_affinity.xlsx = "plsdb_defense_type_affinity.xlsx",
  table_S06_plsdb_defense_subtype_affinity.xlsx = "plsdb_defense_subtype_affinity.xlsx",
  table_S11_plsdb_defense_type_amr_class_affinity.xlsx = "plsdb_defense_type_amr_class_affinity.xlsx",
  table_S12_plsdb_defense_type_amr_type_affinity.xlsx = "plsdb_defense_type_amr_type_affinity.xlsx",
  table_S13_plsdb_defense_subtype_amr_class_affinity.xlsx = "plsdb_defense_subtype_amr_class_affinity.xlsx",
  table_S14_plsdb_defense_subtype_amr_type_affinity.xlsx = "plsdb_defense_subtype_amr_type_affinity.xlsx"
) |>
  purrr::iwalk(\(file, table) readxl::read_xlsx(file.path("data", file)) |>
    dplyr::select(
      entity_1, entity_2, entity_1_count_mA, entity_2_count_mB,
      obs_cooccur_X, total_N, alpha_mle, p_value
    ) |>
    writexl::write_xlsx(file.path("supplementary_tables", table))
  )

# Supplementary Table 7
readxl::read_xlsx("data/plsdb_defense_enrichment_by_phylum.xlsx") |>
  dplyr::select(
    phylum, system, total_hosts, hosts_with_system,
    total_plasmids, plasmids_with_system,
    odds_ratio, log_odds_cc, p_value, p_adjusted
  ) |>
  dplyr::mutate(odds_ratio = as.numeric(odds_ratio)) |>
  writexl::write_xlsx("supplementary_tables/table_S07_plsdb_defense_enrichment_by_phylum.xlsx")

# Supplementary Table 8 is written by fig02.R.

# Supplementary Table 17
readxl::read_xlsx("data/imgpr_plasmid-host_metadata.xlsx") |>
  dplyr::select(
    plasmid_seqid, host_acc, representative, plasmid_complete,
    plasmid_length, mobility, taxonomy_superkingdom,
    taxonomy_phylum, taxonomy_class, taxonomy_order, taxonomy_family,
    taxonomy_genus, taxonomy_species, source_type, ecosystem, ecosystem_category,
    copla_ptu, drop
  ) |>
  writexl::write_xlsx("supplementary_tables/table_S17_imgpr_plasmid-host_metadata.xlsx")

# Supplementary Table 18
readxl::read_xlsx("data/sorensen_dice.xlsx", guess_max = 10000) |>
  dplyr::select(
    plasmid_id, plasmid_complete, plasmid_length, ptu,
    nn_dice_all, nn_dice_complete
  ) |>
  writexl::write_xlsx("supplementary_tables/table_S18_sorensen_dice.xlsx")

# Combined supplementary tables
table_index <- tibble::tribble(
  ~Name, ~Description,
  "PLSDB + RefSeq Archaea plasmid and host metadata",
  "Plasmid and host accessions, lengths, taxonomy, mobility, replicon types, copy numbers and sample metadata. Representative plasmids are marked TRUE in the representative column.",
  "Defense system naming unification",
  "Mapping of DefenseFinder, PADLOC and CRISPRCasTyper names to unified defense system types and subtypes, with inclusion flags.",
  "PLSDB + RefSeq Archaea plasmid-encoded defenses",
  "Coordinates, types and subtypes of defense systems identified on plasmids.",
  "PLSDB + RefSeq Archaea chromosome-encoded defenses",
  "Coordinates, types and subtypes of defense systems identified on host chromosomes.",
  "PLSDB + RefSeq Archaea plasmid defense co-occurrences (types)",
  "Pairwise co-occurrence of defense system types among representative plasmids carrying defenses. Includes plasmid counts, observed co-occurrences, affinity estimates (α) and unadjusted P values.",
  "PLSDB + RefSeq Archaea plasmid defense co-occurrences (subtypes)",
  "Pairwise co-occurrence of defense system subtypes among representative plasmids carrying defenses. Includes plasmid counts, observed co-occurrences, affinity estimates (α) and unadjusted P values.",
  "PLSDB + RefSeq Archaea defense enrichment on plasmids vs chromosomes",
  "Defense enrichment by host phylum, with counts of plasmids and hosts carrying each subtype, odds ratios, continuity-corrected log odds ratios, and unadjusted and adjusted P values.",
  "Analysis of relationship between plasmid features and encoded defenses",
  "Spearman correlations and logistic regressions examining plasmid length, copy number and defense presence. Includes sample sizes, effect estimates, odds ratios and P values.",
  "PLSDB + RefSeq Archaea plasmid AMR",
  "Classifications of AMR genes identified on plasmids.",
  "PLSDB + RefSeq Archaea plasmid anti-defenses",
  "Coordinates and classifications of anti-defense genes identified on plasmids.",
  "PLSDB + RefSeq Archaea plasmid defense type and AMR class co-occurrences",
  "Pairwise co-occurrence of defense system types and AMR classes. Includes feature counts, observed co-occurrences, total plasmids, affinity estimates (α) and unadjusted P values.",
  "PLSDB + RefSeq Archaea plasmid defense type and AMR type co-occurrences",
  "Pairwise co-occurrence of defense system types and AMR types. Includes feature counts, observed co-occurrences, total plasmids, affinity estimates (α) and unadjusted P values.",
  "PLSDB + RefSeq Archaea plasmid defense subtype and AMR class co-occurrences",
  "Pairwise co-occurrence of defense system subtypes and AMR classes. Includes feature counts, observed co-occurrences, total plasmids, affinity estimates (α) and unadjusted P values.",
  "PLSDB + RefSeq Archaea plasmid defense subtype and AMR type co-occurrences",
  "Pairwise co-occurrence of defense system subtypes and AMR types. Includes feature counts, observed co-occurrences, total plasmids, affinity estimates (α) and unadjusted P values.",
  "IMG/PR plasmid-encoded defenses",
  "Coordinates, types and subtypes of defense systems identified on IMG/PR plasmids.",
  "IMG/PR plasmid-encoded AMR",
  "Classifications of AMR genes identified on IMG/PR plasmids.",
  "IMG/PR plasmid and host metadata",
  "Plasmid and host identifiers, length, completeness, mobility, taxonomy, ecological context and PTU assignments. Filtered representative plasmids are marked representative = TRUE and drop = FALSE.",
  "Sørensen–Dice dissimilarity of defense repertoires",
  "Defense-repertoire dissimilarity to neighbouring plasmids, calculated for the full network (nn_dice_all) and the complete-plasmid analysis (nn_dice_complete). Plasmid length, completeness and PTU assignments are also provided."
) |>
  dplyr::mutate(Table = paste0("Table S", dplyr::row_number()), .before = 1)

combined <- openxlsx2::wb_workbook()$add_worksheet("Index", grid_lines = TRUE)
combined$add_data("Index", table_index)
combined$add_font("Index", dims = "A1:C19", name = "Calibri", size = 11)
combined$add_font("Index", dims = "A1:C1", name = "Calibri", bold = TRUE, update = TRUE)
combined$add_cell_style("Index", dims = "A1:C19", vertical = "center", wrap_text = TRUE)
combined$set_col_widths("Index", cols = 1:3, widths = c(14, 54, 100))
combined$set_row_heights("Index", rows = 2:19, heights = 32)
combined$freeze_pane("Index", first_row = TRUE)
purrr::iwalk(table_index$Table, \(table, i) combined$add_hyperlink(
  "Index", dims = paste0("A", i + 1), target = paste0("'", table, "'!A1"), is_external = FALSE
))
combined$add_font("Index", dims = "A2:A19", name = "Calibri", color = openxlsx2::wb_color("0563C1"), underline = "single", update = TRUE)
list.files("supplementary_tables", pattern = "^table_S[0-9]+_.*\\.xlsx$", full.names = TRUE) |>
  purrr::walk(\(file) combined$clone_worksheet(
    old = 1,
    new = paste0("Table S", readr::parse_number(basename(file))),
    from = openxlsx2::wb_load(file)
  ))
combined$set_active_sheet("Index")
combined$set_selected("Index")
combined$save("supplementary_tables.xlsx", overwrite = TRUE)
