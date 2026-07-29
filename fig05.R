path <- "data/merged_master_table.tsv"

merged_metadata_raw <- readr::read_tsv(path, progress = FALSE)

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

# SMALL/LOW-DEFENSE OUTLIERS ---------------------------------------------------

ptu_summary <- merged_metadata |>
  dplyr::filter(!is.na(new_PTU)) |>
  dplyr::summarise(
    n = dplyr::n(),
    defense_prev = mean(!is.na(num_def_sys) & num_def_sys > 0),
    mean_length = mean(length),
    mean_PCN = mean(PCN, na.rm = TRUE),
    non_mobilizable =
      mean(predicted_mobility == "non-mobilizable", na.rm = TRUE),
    mobilizable =
      mean(predicted_mobility == "mobilizable", na.rm = TRUE),
    conjugative =
      mean(predicted_mobility == "conjugative", na.rm = TRUE),
    .by = new_PTU
  )

ptu_summary |>
  dplyr::arrange(mean_length) |>
  dplyr::filter(n >= 100) |>
  dplyr::filter(defense_prev >= 0.25) |>
  dplyr::slice_head() |>
  purrr::pwalk(
    \(new_PTU, n, defense_prev, mean_length, mean_PCN, non_mobilizable, mobilizable, conjugative) {
      cli::cli_alert_info(
        paste0(
          "{new_PTU} (n = {n}) has:\n",
          "- Average length {round(mean_length / 1000, 1)} kb (including putatively incomplete)\n",
          "- Average PCN {round(mean_PCN, 2)}\n",
          "- Defense prevalence of {round(defense_prev * 100)}%\n",
          "- Non mobilizable = {round(non_mobilizable * 100)}%\n",
          "- Mobilizable = {round(mobilizable * 100)}%\n",
          "- Non mobilizable = {round(non_mobilizable * 100)}%\n"
        )
      )
    }
  )

merged_metadata |>
  dplyr::filter(new_PTU == "PTU-E58") |>
  dplyr::filter(putatively_complete == "Yes") |>
  dplyr::summarise(
    avg = round(mean(length), 2),
    med = median(length),
    q1 = quantile(length, 0.25),
    q3 = quantile(length, 0.75)
  ) |>
  purrr::pwalk(
    \(avg, med, q1, q3) {
      cli::cli_alert_info(
        "PTU-E58 length (complete): mean = {round(avg / 1000, 1)} kb, median = {med} bp, IQR = {round(q1 / 1000, 1)}-{round(q3 / 1000, 1)} kb"
      )
    }
  )

merged_metadata |>
  dplyr::filter(new_PTU == "PTU-E58") |>
  tidyr::separate_longer_delim(Unified_subtype_name, ";") |>
  dplyr::distinct(id, Unified_subtype_name) |>
  dplyr::summarise(n = dplyr::n(), .by = Unified_subtype_name) |>
  dplyr::mutate(p = round(n / sum(n) * 100)) |>
  dplyr::filter(!is.na(Unified_subtype_name)) |>
  dplyr::arrange(dplyr::desc(n)) |>
  purrr::pwalk(
    \(Unified_subtype_name, p, n) {
      cli::cli_alert_info(
        "{Unified_subtype_name}: {p}%"
      )
    }
  )

# PUTATIVE PHAGE PLASMIDS ------------------------------------------------------

merged_metadata |>
  dplyr::filter(is_representative == "yes") |>
  dplyr::mutate(putative_phage_plasmid = dplyr::if_else(putative_phage_plasmid == "Yes", TRUE, FALSE)) |>
  dplyr::summarise(n = dplyr::n(), .by = putative_phage_plasmid) |>
  dplyr::mutate(p = n / sum(n) * 100)

merged_metadata |>
  dplyr::filter(is_representative == "yes") |>
  dplyr::mutate(has_defense = dplyr::if_else(!is.na(Type), TRUE, FALSE)) |>
  dplyr::filter(putative_phage_plasmid == "Yes") |>
  dplyr::summarise(n = dplyr::n(), .by = has_defense) |>
  dplyr::mutate(p = n / sum(n) * 100) |>
  dplyr::arrange(dplyr::desc(n))
