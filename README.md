# Defense on Plasmids

This repository contains code and data related to the manuscript: _**Plasmids
are major vectors of defense systems**_ by _**Payne, Mestre, Zhang, et al.**_

## Set up

This project uses a mixture of R and python for analyses.

To set up the same version of R (`4.4.3`),
[install manually](https://cran.r-project.org/), or use
[rig](https://github.com/r-lib/rig):

```shell
rig add 4.4.3
rig switch 4.4.3
```

To install the required R packages, use `renv`:

```shell
Rscript -e "renv::restore()"
```

To set up the same version of python (`3.14`), and install the required python
packages, use [uv](https://docs.astral.sh/uv/):

```shell
uv sync --locked
```

## Required data

Required data are distributed on Zenodo (https://doi.org/10.5281/zenodo.21671301).

To restore the data, download `data.tar.gz` and `data.tar.gz.sha256` from the
Zenodo record into the repository root, verify the checksum, and extract the
archive:

```shell
shasum -a 256 -c data.tar.gz.sha256
tar -xzf data.tar.gz
```

The archive includes a copy of [`DATA_RIGHTS.md`](DATA_RIGHTS.md), which
documents licensing and attribution for redistributed source material.

## Code

> Code is roughly divided by the main figures and sections it relates to. Some
> generated figures received minor aesthetic edits for their final published
> versions.

Run all code from the repository root. Run the R scripts in order:

```shell
mkdir -p plots
Rscript fig01.R
Rscript fig02.R
Rscript fig03.R
Rscript fig04.R
Rscript fig05.R
Rscript fig06.R
```

`fig01.R` prepares the taxonomy-enriched metadata used by later analyses.
`fig02.R` prepares inputs used by `figED04.ipynb` and `figED05.ipynb`, and
writes Supplementary Table 8. Then run all cells in `fig05.ipynb`,
`fig06.ipynb`, `figED01.ipynb`, `figED04.ipynb`, `figED05.ipynb` and
`figS09.ipynb`, selecting `.venv/bin/python` as the notebook kernel and keeping
the repository root as the working directory. Finally, export the
supplementary tables:

```shell
Rscript export_supplementary_tables.R
```

### `fig01.R`

Produces plots and statistics related to **Figure 1**, **Extended Data Figures
2** and **3**, **Supplementary Figures 2** and **3**, and the results
sections: _**Plasmids are reservoirs of defense systems**_, _**Plasmid-encoded
defenses are diverse and unevenly distributed**_, _**Many pairs of
plasmid-encoded defenses co-occur more often than expected**_, and _**The
distribution of plasmid-encoded defenses is distinct from chromosomes**_.

### `fig02.R`

Produces plots and statistics related to **Figure 2**, **Supplementary Figure
4**, **Supplementary Table 8**, and the results section: _**Defense systems are most
often found on large conjugative plasmids**_.

### `fig03.R`

Produces plots and statistics related to **Figure 3**, **Supplementary Figures
5** and **6**, **Supplementary Tables 11**, **12**, **13**, and **14**, and
the results section: _**Plasmids arm their hosts against phages and
antibiotics**_.

### `fig04.R`

Produces plots and statistics related to **Figure 4**, **Extended Data Figure
6**, **Supplementary Figure 7**, and the results section: _**Defense is more
common than AMR in plasmids across ecosystems**_.

### `fig05.R`

Produces statistics related to the results section: _**An interactive network
for exploring the global plasmidome**_.

### `fig05.ipynb`

Produces plots and statistics related to **Figure 5** and the results section:
_**An interactive network for exploring the global plasmidome**_.

### `fig06.R`

Produces plots and statistics related to **Figure 6** and the results section:
_**Plasmids drive defense system mobilization across broad taxonomic scales**_.

### `fig06.ipynb`

Produces plots and statistics related to **Figure 6**, **Extended Data Figure
7**, and the results section: _**Plasmids drive defense system mobilization
across broad taxonomic scales**_.

### `figED01.ipynb`

Produces plots and statistics related to **Extended Data Figure 1** and the
results section: _**Plasmid-encoded defenses are diverse and unevenly
distributed**_.

### `figED04.ipynb`

Produces plots and statistics related to **Extended Data Figure 4** and the
results section: _**Defense systems are most often found on large conjugative
plasmids**_.

### `figED05.ipynb`

Produces plots and statistics related to **Extended Data Figure 5** and the
results section: _**Defense systems are most often found on large conjugative
plasmids**_.

### `figS09.ipynb`

Produces plots and statistics related to **Supplementary Figure 9** and the
results section: _**An interactive network for exploring the global
plasmidome**_.

### `export_supplementary_tables.R`

Writes the publication versions of **Supplementary Tables 1–7** and **9–18** to
`supplementary_tables/`, using the full input and analysis files in `data/`.
**Supplementary Table 8** is written by `fig02.R` to
`supplementary_tables/table_S08_correlation_between_pcn_length_and_defense.xlsx`.
The export script also combines all 18 tables into `supplementary_tables.xlsx`,
with an index of table names, descriptions and links.

## Bundled data

### Figure 1, Extended Data Figures 1–3, and Supplementary Figures 2–3

- `data/defense_system_unification.xlsx` contains the cross-tool defense-system
  name mapping used for **Supplementary Table 2**.
- `data/plsdb_plasmid-host_metadata.xlsx` contains PLSDB plasmid and host
  metadata before GTDB taxonomy is added. It is the starting version of
  **Supplementary Table 1**.
- `data/plsdb_plasmid_defense.xlsx` and `data/plsdb_host_defense.xlsx` contain
  the defense calls on plasmids and host chromosomes, respectively. They are
  equivalent to **Supplementary Tables 3** and **4**.
- `data/gtdb/bac120_metadata_r220.tsv` and `data/gtdb/ar53_metadata_r220.tsv`
  provide GTDB release 220 taxonomy. `data/gtdb/bac120_r220_phylum.tree` and
  `data/gtdb/ar53_r220_phylum.tree` provide the corresponding phylum-level
  trees. `data/gtdb/gtdb_220_pd.tsv` contains phylum phylogenetic diversity,
  manually recorded from the
  [GTDB release 220 statistics](https://gtdb.ecogenomic.org/stats/r220).
- `data/gtdb/bac120_metadata_r207.tsv` and `data/gtdb/ar53_metadata_r207.tsv`
  are archived release 207 metadata. They are bundled for reference but are not
  read by the current figure scripts.
- `data/plsdb_plasmid-host_metadata_with_taxonomy.xlsx` contains the full data
  for **Supplementary Table 1** after GTDB taxonomy is added by `fig01.R`.
  It is reused by `fig01.R`, `fig04.R`, `fig06.ipynb`, and `figED01.ipynb`.
- `data/plsdb_defense_type_affinity.xlsx` and
  `data/plsdb_defense_subtype_affinity.xlsx` contain defense co-occurrence
  affinity results generated by `fig01.R` and used for
  **Supplementary Tables 5** and **6**.
- `data/plsdb_defense_enrichment_by_phylum.xlsx` contains the plasmid-versus-
  chromosome defense enrichment tests generated by `fig01.R` and used for
  **Supplementary Table 7**.

The GTDB trees and metadata were retrieved from the
[GTDB data repository](https://data.gtdb.aau.ecogenomic.org/).

### Figure 2, Extended Data Figures 4–5, and Supplementary Figure 4

- `data/plsdb_gene_annotations.tsv` contains the functional annotations for
  PLSDB plasmid genes described in the Methods. It is used by `fig02.R`,
  `fig04.R`, and `fig06.ipynb`.
- `data/plsdb_plasmid_mob.xlsx` contains the PLSDB plasmid mobility
  classifications generated by `fig02.R` and then used by `figED04.ipynb`.
- `data/plsdb_plasmid_inc.xlsx` contains the Inc classifications generated by
  `fig02.R` for the plasmids in **Figure 2e** and is then used by
  `figED05.ipynb`.

### Figure 3 and Supplementary Figures 5–6

- `data/plsdb_plasmid_amr.xlsx` and `data/plsdb_plasmid_antidefense.xlsx`
  contain PLSDB AMR and anti-defense calls and are equivalent to **Supplementary
  Tables 9** and **10**.
- `data/plsdb_defense_type_amr_class_affinity.xlsx`,
  `data/plsdb_defense_type_amr_type_affinity.xlsx`,
  `data/plsdb_defense_subtype_amr_class_affinity.xlsx`, and
  `data/plsdb_defense_subtype_amr_type_affinity.xlsx` contain the defense–AMR
  co-occurrence affinity results generated by `fig03.R` and used for
  **Supplementary Tables 11**, **12**, **13**, and **14**.
- `data/pNDM-Mar_skani_matrix.txt` and `data/pNDM-Mar_skani_matrix.txt.af`
  contain the ANI and alignment-fraction matrices used for **Figure 3e**.

### Figures 4–5, Extended Data Figure 6, and Supplementary Figures 7–9

- `data/merged_master_table.tsv` combines PLSDB and IMG/PR plasmid metadata and
  annotations. It incorporates information represented in **Supplementary Tables
  1**, **3**, **4**, **9**, **10**, **15**, **16**, and **17**, and is used by
  `fig04.R` and `fig05.R`.
- `data/imgpr_plasmid_defense.xlsx`, `data/imgpr_plasmid_amr.xlsx`, and
  `data/imgpr_plasmid-host_metadata.xlsx` contain the IMG/PR defense calls, AMR
  calls, and full plasmid/host metadata used for **Supplementary
  Tables 15**, **16**, and **17**.
- `data/sorensen_dice.xlsx` contains the full nearest-neighbour similarity
  results used for **Supplementary Table 18**.
- `data/imgpr_gene_annotations.tsv` contains the IMG/PR plasmid-gene annotations
  used by `fig04.R` and `fig06.ipynb`. `data/imgpr_eggnog.tsv` supplies the
  EggNOG COG categories joined to those annotations by `fig04.R`.
- `data/plsdb_imgpr_plasmid_ptu.xlsx`, `data/plsdb_imgpr_plasmid_has_amr.xlsx`,
  and `data/plsdb_imgpr_plasmid_has_antidef.xlsx` are combined PLSDB/IMG/PR PTU
  assignments and per-plasmid feature indicators used by `figS09.ipynb`.

### Figure 6 and Extended Data Figure 7

- `data/transfers_plasmid_host.tsv`, `data/transfers_plasmid_plasmid.tsv`,
  `data/transfers_whole_plasmid_chromosome.tsv`,
  `data/transfers_whole_plasmid_plasmid.tsv`, and
  `data/transfers_whole_plasmid.tsv` contain the five classes of precomputed
  candidate transfer calls analysed by `fig06.R`, as described in the Methods.
- `data/gtdb/bac120_taxonomy_r226.tsv` supplies GTDB release 226 taxonomy for
  assemblies in the transfer tables.
- `data/assemblies_not_in_r226.txt` is generated by `fig06.R` and lists
  transfer-associated assemblies needing de novo GTDB-Tk classification.
- `data/gtdbtk.bac120.summary.tsv` and `data/gtdbtk.ar53.summary.tsv` are the
  GTDB-Tk bacterial and archaeal summaries for those assemblies. `fig06.R`
  combines them with the release 226 taxonomy before assigning transfer
  distances.
- `data/transfers_manual_filter.xlsx` records manually excluded transfer
  candidates by `transfer_id`.
