
gtdb_r226_taxonomy <- readr::read_tsv(
  "data/gtdb/bac120_taxonomy_r226.tsv", 
  col_names = c("accession", "taxonomy"), col_types = c("c", "c")
)

gtdb_r226_assemblies <- gtdb_r226_taxonomy |>
  dplyr::mutate(accession = stringr::str_remove(accession, "^RS_|^GB_")) |>
  dplyr::pull(accession)

# CHECK WHAT NEEDS TAXONOMY UPDATED

plasmid_host_raw <- readr::read_tsv("data/transfers_plasmid_host.tsv") 
plasmid_plasmid_raw <- readr::read_tsv("data/transfers_plasmid_plasmid.tsv")
whole_plasmid_chromosome_raw <- readr::read_tsv("data/transfers_whole_plasmid_chromosome.tsv")
whole_plasmid_plasmid_raw <- readr::read_tsv("data/transfers_whole_plasmid_plasmid.tsv")
whole_plasmid_raw <- readr::read_tsv("data/transfers_whole_plasmid.tsv")

plasmid_host_assemblies <- dplyr::bind_rows(
  plasmid_host_raw |> dplyr::distinct(assembly_query) |> dplyr::rename(assembly = assembly_query),
  plasmid_host_raw |> dplyr::distinct(assembly_target) |> dplyr::rename(assembly = assembly_target)
) |> dplyr::distinct()

plasmid_plasmid_assemblies <- dplyr::bind_rows(
  plasmid_plasmid_raw |> dplyr::distinct(assembly_query) |> dplyr::rename(assembly = assembly_query),
  plasmid_plasmid_raw |> dplyr::distinct(assembly_target) |> dplyr::rename(assembly = assembly_target)
) |> dplyr::distinct()

whole_plasmid_chromosome_assemblies <- dplyr::bind_rows(
  whole_plasmid_chromosome_raw |> dplyr::distinct(assembly_host) |> dplyr::rename(assembly = assembly_host),
  whole_plasmid_chromosome_raw |> dplyr::distinct(assembly_plasmid) |> dplyr::rename(assembly = assembly_plasmid)
) |> dplyr::distinct()

whole_plasmid_plasmid_assemblies <- dplyr::bind_rows(
  whole_plasmid_plasmid_raw |> dplyr::distinct(assembly_query) |> dplyr::rename(assembly = assembly_query),
  whole_plasmid_plasmid_raw |> dplyr::distinct(assembly_target) |> dplyr::rename(assembly = assembly_target)
) |> dplyr::distinct()

whole_plasmid_assemblies <- dplyr::bind_rows(
  whole_plasmid_raw |> dplyr::distinct(assembly),
  whole_plasmid_raw |> dplyr::distinct(assembly_rep) |> dplyr::rename(assembly = assembly_rep)
) |> dplyr::distinct()

assemblies <- dplyr::bind_rows(
  plasmid_host_assemblies,
  plasmid_plasmid_assemblies,
  whole_plasmid_chromosome_assemblies,
  whole_plasmid_plasmid_assemblies,
  whole_plasmid_assemblies
) |>
  dplyr::distinct() |>
  dplyr::arrange(assembly) |>
  dplyr::filter(!is.na(assembly))

assemblies_todo <- assemblies |>
  dplyr::filter(! assembly %in% gtdb_r226_assemblies)

assemblies_todo |>
  readr::write_tsv("data/assemblies_not_in_r226.txt", col_names = FALSE)


# RUN GTDB-TK ON ASSEMBLIES TO GENERATE SUMMARIES

gtdbtk_bac_summary <- readr::read_tsv(
  "data/gtdbtk.bac120.summary.tsv", 
  na = "N/A", 
  col_types = list("c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c")
)

gtdbtk_ar_summary <- readr::read_tsv(
  "data/gtdbtk.ar53.summary.tsv", 
  na = "N/A", 
  col_types = list("c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c")
)

gtdbtk_taxonomy <- dplyr::bind_rows(
  gtdbtk_bac_summary,
  gtdbtk_ar_summary
) |>
  dplyr::filter(is.na(warnings)) |>
  dplyr::mutate(accession = stringr::str_extract(user_genome, "GC[AF]_[0-9]*\\.[0-9]*")) |>
  dplyr::select(accession, classification) |>
  dplyr::rename(taxonomy = classification)

taxa_levels <- c("domain", "phylum", "class", "order", "family", "genus", "species")

all_taxonomy <- dplyr::bind_rows(
  gtdb_r226_taxonomy,
  gtdbtk_taxonomy
) |>
  dplyr::mutate(accession = stringr::str_remove(accession, "^RS_|^GB_")) |>
  tidyr::separate_wider_delim(taxonomy, ";", names = taxa_levels) |>
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ stringr::str_remove_all(., "^[a-zA-Z]__")))

classify_transfer <- function(df) {
  out <- df |> dplyr::mutate(
    transfer_type = dplyr::case_when(
      query_assembly == target_assembly ~ "Within assembly",
      query_species  == target_species  ~ "Between\nassemblies\n(same species)",
      query_genus    == target_genus    ~ "Between\nspecies\n(same genus)",
      query_family   == target_family   ~ "Between\ngenera\n(same family)",
      query_order    == target_order    ~ "Between\nfamilies\n(same order)",
      query_class    == target_class    ~ "Between\norders\n(same class)",
      query_phylum   == target_phylum   ~ "Between\nclasses\n(same phylum)",
      query_phylum   != target_phylum   ~ "Between\nphyla\n(same domain)"
    ),
    transfer_level = dplyr::case_when(
      query_assembly == target_assembly ~ 1,
      query_species  == target_species  ~ 2,
      query_genus    == target_genus    ~ 3,
      query_family   == target_family   ~ 4,
      query_order    == target_order    ~ 5,
      query_class    == target_class    ~ 6,
      query_phylum   == target_phylum   ~ 7,
      query_phylum   != target_phylum   ~ 8
    )
  )
  out
}

attach_taxonomy <- function(df) {
  taxa_levels <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  out <- df |>
    dplyr::left_join(all_taxonomy, by = dplyr::join_by(query_assembly == accession)) |>
    dplyr::rename_with(~ paste0("query_", .), dplyr::all_of(taxa_levels)) |>
    dplyr::left_join(all_taxonomy, by = dplyr::join_by(target_assembly == accession)) |>
    dplyr::rename_with(~ paste0("target_", .), dplyr::all_of(taxa_levels))
  out
}



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

# PLASMID-CHROMOSOME -----------------------------------------------------------

plasmid_host <- plasmid_host_raw |> 
  dplyr::filter(!type %in% c("Hok/Sok", "MazEF"))

plasmid_host_taxa <- plasmid_host |>
  dplyr::select(assembly_query, assembly_target, query_plasmid, target) |>
  dplyr::rename(
    query_assembly = assembly_query, target_assembly = assembly_target,
    query_id = query_plasmid, target_id = target
  ) |>
  attach_taxonomy() |>
  classify_transfer() |>
  dplyr::mutate(transfer_category = "plasmid-chromosome")

plasmid_host_clean <- plasmid_host_taxa |>
  dplyr::filter(!is.na(query_domain) & !is.na(target_domain))

plasmid_host_missing <- plasmid_host_taxa |>
  dplyr::filter(is.na(query_domain) | is.na(target_domain))

# PLASMID-PLASMID --------------------------------------------------------------

plasmid_plasmid <- plasmid_plasmid_raw |> 
  dplyr::filter(!type %in% c("Hok/Sok", "MazEF"))

plasmid_plasmid_taxa <- plasmid_plasmid |>
  dplyr::select(assembly_query, assembly_target, query_plasmid, target) |>
  dplyr::rename(
    query_assembly = assembly_query, target_assembly = assembly_target,
    query_id = query_plasmid, target_id = target
  ) |>
  attach_taxonomy() |>
  classify_transfer() |>
  dplyr::mutate(transfer_category = "plasmid-plasmid")

plasmid_plasmid_clean <- plasmid_plasmid_taxa |>
  dplyr::filter(!is.na(query_domain) & !is.na(target_domain))

plasmid_plasmid_missing <- plasmid_plasmid_taxa |>
  dplyr::filter(is.na(query_domain) | is.na(target_domain))

# WHOLE-PLASMID-CHROMOSOME -----------------------------------------------------

whole_plasmid_chromosome <- whole_plasmid_chromosome_raw |>
  dplyr::left_join(
    plsdb_plasmid_defense,
    by = dplyr::join_by(Query_name == plasmid_seqid), relationship = "many-to-many"
  ) |>
  dplyr::filter(!is.na(type)) |>
  dplyr::distinct(dplyr::across(dplyr::all_of(names(whole_plasmid_chromosome_raw))))

whole_plasmid_chromosome_taxa <- whole_plasmid_chromosome |>
  dplyr::select(assembly_plasmid, assembly_host, Query_name, Ref_name) |>
  dplyr::rename(
    query_assembly = assembly_plasmid, target_assembly = assembly_host, 
    query_id = Query_name, target_id = Ref_name
  ) |>
  attach_taxonomy() |>
  classify_transfer() |>
  dplyr::mutate(transfer_category = "whole-plasmid-chromosome")

whole_plasmid_chromosome_clean <- whole_plasmid_chromosome_taxa |>
  dplyr::filter(!is.na(query_domain) & !is.na(target_domain))

whole_plasmid_chromosome_missing <- whole_plasmid_chromosome_taxa |>
  dplyr::filter(is.na(query_domain) | is.na(target_domain))

# WHOLE-PLASMID-PLASMID --------------------------------------------------------

whole_plasmid_plasmid <- whole_plasmid_plasmid_raw |>
  dplyr::left_join(
    plsdb_plasmid_defense |> dplyr::filter(!subtype %in% c("Hok/Sok", "MazEF")),
    by = dplyr::join_by(Query_name == plasmid_seqid), relationship = "many-to-many"
  ) |>
  dplyr::filter(!is.na(type)) |>
  dplyr::distinct(dplyr::across(dplyr::all_of(names(whole_plasmid_plasmid_raw)))) |>
  dplyr::left_join(
    plsdb_plasmid_defense |> dplyr::filter(!subtype %in% c("Hok/Sok", "MazEF")),
    by = dplyr::join_by(Ref_name == plasmid_seqid), relationship = "many-to-many"
  ) |>
  dplyr::filter(!is.na(type)) |>
  dplyr::distinct(dplyr::across(dplyr::all_of(names(whole_plasmid_plasmid_raw))))

whole_plasmid_plasmid_taxa <- whole_plasmid_plasmid |>
  dplyr::select(assembly_query, assembly_target, Query_name, Ref_name) |>
  dplyr::rename(
    query_assembly = assembly_query, target_assembly = assembly_target,
    query_id = Query_name, target_id = Ref_name
  ) |>
  attach_taxonomy() |>
  classify_transfer() |>
  dplyr::mutate(transfer_category = "whole-plasmid-plasmid")

whole_plasmid_plasmid_clean <- whole_plasmid_plasmid_taxa |>
  dplyr::filter(!is.na(query_domain) & !is.na(target_domain))

whole_plasmid_plasmid_missing <- whole_plasmid_plasmid_taxa |>
  dplyr::filter(is.na(query_domain) | is.na(target_domain))

# WHOLE-PLASMID ----------------------------------------------------------------

whole_plasmid <- whole_plasmid_raw |>
  dplyr::mutate(
    Type = purrr::map_chr(stringr::str_split(Type, ";"), \(x) {
      x <- x[!x %in% c("Hok/Sok", "MazEF")]
      if (length(x) == 0) NA_character_ else stringr::str_c(x, collapse = ";")
    })
  ) |>
  dplyr::filter(!is.na(Type))

whole_plasmid_taxa <- whole_plasmid |>
  dplyr::filter(!is.na(assembly_rep)) |>
  dplyr::select(assembly, assembly_rep, Member, Representative) |>
  dplyr::rename(
    query_assembly = assembly, target_assembly = assembly_rep, 
    query_id = Member, target_id = Representative
    ) |>
  attach_taxonomy() |>
  classify_transfer() |>
  dplyr::mutate(transfer_category = "whole-plasmid") |>
  dplyr::filter(!(query_assembly == target_assembly & query_id == target_id))

whole_plasmid_clean <- whole_plasmid_taxa |>
  dplyr::filter(!is.na(query_domain) & !is.na(target_domain)) |>
  dplyr::filter(query_assembly != target_assembly)

whole_plasmid_missing <- whole_plasmid_taxa |>
  dplyr::filter(is.na(query_domain) | is.na(target_domain))

# ------------------------------------------------------------------------------

all_missing <- dplyr::bind_rows(
  plasmid_host_missing,
  plasmid_plasmid_missing,
  whole_plasmid_chromosome_missing,
  whole_plasmid_plasmid_missing,
  whole_plasmid_missing
) |>
  dplyr::distinct()

nrow(all_missing)

all_clean <- dplyr::bind_rows(
  plasmid_host_clean,
  plasmid_plasmid_clean,
  whole_plasmid_chromosome_clean,
  whole_plasmid_plasmid_clean,
  whole_plasmid_clean
) |>
  dplyr::distinct()


plasmid_host_coord <- plasmid_host |>
  dplyr::mutate(
    query_coord = stringr::str_remove_all(query, type),
    query_coord = stringr::str_remove_all(query_coord, "\\.0"),
    query_coord = stringr::str_replace_all(query_coord, "_-", ":"),
    target_coord = paste0(target, ":", tstart, "-", tend)
  ) |>
  dplyr::select(query_plasmid, target, query_coord, target_coord) |>
  dplyr::mutate(transfer_category = "plasmid-chromosome")

plasmid_plasmid_coord <- plasmid_plasmid |>
  dplyr::mutate(
    query_coord = stringr::str_remove_all(query, type),
    query_coord = stringr::str_remove_all(query_coord, "\\.0"),
    query_coord = stringr::str_replace_all(query_coord, "_-", ":"),
    target_coord = paste0(target, ":", tstart, "-", tend)
  ) |>
  dplyr::select(query_plasmid, target, query_coord, target_coord) |>
  dplyr::mutate(transfer_category = "plasmid-plasmid")

all_coord <- dplyr::bind_rows(
  plasmid_host_coord, plasmid_plasmid_coord
) |>
  dplyr::distinct() |>
  dplyr::rename(query_id = query_plasmid, target_id = target)

all_clean_w_coord <- all_clean |>
  dplyr::left_join(all_coord, by = dplyr::join_by(query_id, target_id, transfer_category)) |>
  dplyr::arrange(dplyr::desc(transfer_type), transfer_category) |>
  dplyr::mutate(
    query_coord = dplyr::if_else(is.na(query_coord), query_id, query_coord),
    target_coord = dplyr::if_else(is.na(target_coord), target_id, target_coord)
  ) |>
  dplyr::arrange(dplyr::desc(transfer_level), query_coord, target_coord) |>
  dplyr::mutate(transfer_id = paste0(stringr::str_pad(dplyr::row_number(), width = 5, side = "left", pad = "0"), "-", transfer_level)) |>
  dplyr::select(transfer_id, query_coord, target_coord, dplyr::everything())

nrow(all_missing) / (nrow(all_clean) + nrow(all_missing)) * 100

category <- c(
  `plasmid-chromosome`        = "#fecc67",
  `whole-plasmid-chromosome`  = "#83ba59",
  `plasmid-plasmid`           = "#ffa35f",
  `whole-plasmid-plasmid`     = "#6ec2c6",
  `whole-plasmid`             = "#cc95c1"
)

type <- c(
  "Within assembly",
  "Between\nassemblies\n(same species)",
  "Between\nspecies\n(same genus)",
  "Between\ngenera\n(same family)",
  "Between\nfamilies\n(same order)",
  "Between\norders\n(same class)",
  "Between\nclasses\n(same phylum)",
  "Between\nphyla\n(same domain)"
)

manual_filter <- readxl::read_xlsx("data/transfers_manual_filter.xlsx")

drop_ids <- manual_filter |>
  dplyr::pull(transfer_id)

pc <- 0.1 # pseudocount when n == 1

plot <- all_clean_w_coord |>
  dplyr::filter(!transfer_id %in% drop_ids) |>
  dplyr::summarise(n = dplyr::n(), .by = c(transfer_type, transfer_category)) |>
  dplyr::mutate(
    transfer_type     = factor(transfer_type, levels = type),
    transfer_category = factor(transfer_category, levels = names(category))
  ) |>
  tidyr::complete(transfer_type, transfer_category, fill = list(n = NA_integer_)) |>
  dplyr::mutate(n_plot = dplyr::if_else(n == 1L, n + pc, as.numeric(n))) |>
  ggplot2::ggplot(ggplot2::aes(x = transfer_type, y = n_plot, fill = transfer_category)) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.7),
    width = 0.7,
    na.rm = TRUE
  ) +
  ggplot2::scale_x_discrete(labels = type) +
  ggplot2::scale_y_log10(
    breaks = 10^(0:6),
    labels = function(x) paste0("10^", as.integer(log10(x))),
    expand = ggplot2::expansion(mult = c(0, 0.05))
  ) +
  ggplot2::scale_fill_manual(values = category, name = NULL, drop = FALSE) +
  ggplot2::labs(x = NULL, y = "No. transfers (log scale)") +
  ggplot2::theme_bw(base_size = 6) +
  ggplot2::theme(
    legend.position = "top",
    legend.justification = "center",
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 6, lineheight = 0.95, hjust = 0.5)
  )

plot

plot |> 
  ggplot2::ggsave(
    filename = "plots/fig06_F.pdf", 
    height = 40, width = 160, 
    units = "mm"
  )

all_clean_w_coord |> nrow()

all_clean_w_coord |>
  dplyr::filter(transfer_type == "Within assembly") |>
  nrow()

all_clean_w_coord |>
  dplyr::summarise(n = dplyr::n(), .by = transfer_category) |>
  dplyr::arrange(dplyr::desc(n))
