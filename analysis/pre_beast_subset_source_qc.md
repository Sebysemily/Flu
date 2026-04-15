# Source-first QC for the BEAST subset

## Inputs

- Final FASTA: `data/final/H5N1_final.fasta`
- Panel taxa: `data/beast_pre/panels/panel_main_taxa.tsv`
- Subset alignment: `data/beast_pre/panels/panel_main_concat.subset.fasta`
- Post-MAFFT alignment: `data/beast_pre/panels/panel_main.subset.mafft`
- Ecuador summary: `data/input/H5N1_EC_summary.csv`
- Context summary: `data/input/H5N1_context_summary.csv`
- Ecuador audit: `data/assembled/ecuador_intermediate_audit.csv`

## Overview

- Taxa in BEAST panel: 62
- Candidate outliers with >=1 flag: 11
- Ecuador taxa in panel: 28
- Context taxa in panel: 34

## Flag Counts

- subset_missing_segments_flag: 0
- high_raw_n_fraction: 2
- high_raw_ambiguous_fraction: 0
- multiple_short_segments: 0
- low_min_segment_length_ratio: 1
- context_segment_count_mismatch: 0
- ecuador_mira_incomplete_history: 0
- high_subset_consensus_distance: 4
- high_post_consensus_distance: 4
- post_consensus_distance_increase: 0
- high_post_singleton_burden: 4
- atypical_post_gap_shift: 0
- post_ungapped_sequence_changed: 0
- any_outlier_flag: 11
- n_taxa_total: 62

## Top Candidate Outliers

| taxon | origin | role | score | flags | raw_N | subset_dist | post_dist |
| --- | --- | --- | ---: | --- | ---: | ---: | ---: |
| 24-032160-004_PQ797990__usa_distal/Usa/2024-10-24 | context | usa_distal | 2 | high_subset_consensus_distance,high_post_consensus_distance | 0.0000 | 0.0747 | 0.0747 |
| 24-032809-001_PQ798038__usa_distal/Usa/2024-11-05 | context | usa_distal | 2 | high_subset_consensus_distance,high_post_consensus_distance | 0.0000 | 0.0745 | 0.0745 |
| 22-011334-002_OQ968009__american_anchor/Usa/2022-04-13 | context | american_anchor | 2 | high_subset_consensus_distance,high_post_consensus_distance | 0.0000 | 0.0416 | 0.0416 |
| 22-011138-001_OQ968033__american_anchor/Usa/2022-04-11 | context | american_anchor | 2 | high_subset_consensus_distance,high_post_consensus_distance | 0.0000 | 0.0409 | 0.0409 |
| PQ318384.1__regional_context/Bolivia/2023-06-06 | context | regional_context | 1 | high_raw_n_fraction | 0.1584 | 0.0048 | 0.0048 |
| PQ318387.1__regional_context/Bolivia/2023-03-03 | context | regional_context | 1 | high_raw_n_fraction | 0.0274 | 0.0042 | 0.0042 |
| 1455-N_PP386944.1__regional_context/Brazil/2023-05-27 | context | regional_context | 1 | low_min_segment_length_ratio | 0.0014 | 0.0063 | 0.0065 |
| 2259-N_PP768254.1__regional_context/Brazil/2023-11-09 | context | regional_context | 1 | high_post_singleton_burden | 0.0000 | 0.0060 | 0.0070 |
| 1941-N_OR860365.1__regional_context/Brazil/2023-09-25 | context | regional_context | 1 | high_post_singleton_burden | 0.0000 | 0.0076 | 0.0067 |
| 2165-SO_OR852419.1__regional_context/Brazil/2023-10-11 | context | regional_context | 1 | high_post_singleton_burden | 0.0000 | 0.0061 | 0.0061 |
| Magdalena3503_OQ683494.1__regional_context/Colombia/2022-11-18 | context | regional_context | 1 | high_post_singleton_burden | 0.0000 | 0.0061 | 0.0055 |
