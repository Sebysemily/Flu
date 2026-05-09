# NP+MP+NS RTT sensitivity analysis

## Pipeline changes

- Added an NP+MP+NS mini-concatenation in `rules/01_ml_trees.smk`.
  - Input segments: `NP`, `MP`, `NS`.
  - Partitioning: codon-aware `NP` (`NP_cp12`, `NP_cp3`) plus simple `MP` and `NS` partitions.
  - Alignment output: `data/phylogeny/aligned/H5N1_NP_MP_NS.mafft`.
  - Partition output: `data/phylogeny/H5N1_NP_MP_NS.partitions`.

- Added an ML tree for the NP+MP+NS mini-concat in `rules/01_ml_trees.smk`.
  - RAxML-NG output prefix: `results/phylogeny/raxml/np_mp_ns/H5N1_NP_MP_NS`.
  - Main tree output: `results/phylogeny/raxml/np_mp_ns/H5N1_NP_MP_NS.raxml.supportTBE`.

- Added a group-filtered RTT sensitivity panel in `rules/02_pre_beast.smk`.
  - No MRCA/clade-based pruning is used.
  - Excludes `__usa_distal`.
  - Excludes `__usa_neighbor`.
  - Excludes non-American anchors, such as `__eurasian_anchor_*`.
  - Keeps only 6 `__american_anchor` tips: forced `OQ968009` plus 5 additional anchors.
  - Keeps all other Ecuador and regional-context samples.

- Added `code/02_Beast/build_np_mp_ns_sensitivity_panel.py`.
  - Builds `data/pre_beast/sensitivity_np_mp_ns/panel_taxa.tsv`.
  - Writes `data/pre_beast/sensitivity_np_mp_ns/panel_audit.tsv`.
  - Uses group labels in taxon names rather than hardcoded clade membership.

- Added the NP+MP+NS sensitivity outputs to `rule all` in `snakefile`.
  - Running `snakemake` now also builds this sensitivity analysis.
  - The standalone target is `np_mp_ns_rtt_sensitivity`.

## RTT input panel

The NP+MP+NS mini-concat contains 254 complete taxa.

Initial group filtering produced a 190-taxon RTT panel:

| Metric | Count |
|---|---:|
| Complete NP+MP+NS taxa | 254 |
| Included in RTT panel | 190 |
| Excluded `__usa_distal` | 30 |
| Excluded `__usa_neighbor` | 20 |
| Excluded non-American anchors | 12 |
| Excluded extra American anchors | 2 |
| Kept American anchors | 6 |

Panel composition before RTT:

| Role | Count |
|---|---:|
| Ecuador | 68 |
| Regional context | 116 |
| American anchor | 6 |

Sanity checks:

- `__usa_distal`: 0 in panel.
- `__usa_neighbor`: 0 in panel.
- Non-American anchors: 0 in panel.
- American anchors: 6 in panel, including `22-011334-002_OQ968009__american_anchor/Usa/2022-04-13`.

## TreeTime RTT result

TreeTime root-to-tip regression:

| Metric | Value |
|---|---:|
| Rate | `6.837e-04` |
| R-squared | `0.03` |
| Estimated root date | `1831.57` |
| RTT outliers | 12 |
| Panel kept after RTT outlier removal | 178 |

Interpretation:

- The temporal signal is weak in this NP+MP+NS group-filtered panel.
- `r^2 = 0.03` means sampling dates explain very little of the root-to-tip distance variation.
- The inferred root date, `1831.57`, is biologically implausible for this H5N1 context and is another sign that the current panel/rooting is not clock-like.
- The retained American anchors do not behave well in this RTT: all 6 were flagged as TreeTime outliers.

## TreeTime outliers

TreeTime flagged these 12 tips:

| Taxon | Given date | Apparent date | Residual |
|---|---:|---:|---:|
| `Pantanal1403_PQ135843.1__regional_context/Brazil/2023-10-01` | 2023.7493 | 1985.3226 | -15.1979 |
| `1355-N-2025_PV796407.1__regional_context/Brazil/2025-05-28` | 2025.4041 | 2049.5505 | 9.5500 |
| `1245-N3-2025_PV660446.1__regional_context/Brazil/2025-05-12` | 2025.3603 | 2049.0533 | 9.3707 |
| `1246-N7-2025_PV660438.1__regional_context/Brazil/2025-05-12` | 2025.3603 | 2049.6238 | 9.5963 |
| `1246-N5-2025_PV659830.1__regional_context/Brazil/2025-05-12` | 2025.3603 | 2049.6238 | 9.5963 |
| `Choco3501_OQ683502.1__regional_context/Colombia/2022-10-09` | 2022.7712 | 2041.5719 | 7.4357 |
| `22-010159-001_OQ968065__american_anchor/Usa/2022-04-04` | 2022.2562 | 2035.9812 | 5.4283 |
| `22-011334-002_OQ968009__american_anchor/Usa/2022-04-13` | 2022.2808 | 2038.8254 | 6.5435 |
| `22-010936-001_OQ968073__american_anchor/Usa/2022-04-07` | 2022.2644 | 2038.2876 | 6.3372 |
| `22-012373-002_OQ968049__american_anchor/Usa/2022-04-21` | 2022.3027 | 2039.0455 | 6.6218 |
| `22-011138-001_OQ968033__american_anchor/Usa/2022-04-11` | 2022.2753 | 2039.0455 | 6.6327 |
| `22-011138-002_OQ968041__american_anchor/Usa/2022-04-11` | 2022.2753 | 2039.0455 | 6.6327 |

After applying this one-round RTT outlier filter:

- `data/pre_beast/sensitivity_np_mp_ns/panel_taxa.rtt_filtered.tsv` has 178 taxa.
- All 6 American anchors are removed by RTT.
- No `__usa_distal`, `__usa_neighbor`, or non-American anchors are present before or after RTT filtering.

## Recommendation

The current NP+MP+NS sensitivity panel does not show good temporal signal in a single RTT pass. If this analysis is used as an exploratory sensitivity check, the first samples to remove are exactly the 12 TreeTime outliers listed above. After removing them, a second diagnostic RTT would be useful to determine whether the remaining 178-tip panel improves enough to support downstream time-scaled inference.
