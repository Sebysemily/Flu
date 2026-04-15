# Observational QC: subset before vs after MAFFT

## Inputs

- Before: `data/beast_pre/panels/panel_main_concat.subset.fasta`
- After: `data/beast_pre/panels/panel_main.subset.mafft`
- Taxa annotations: `data/beast_pre/panels/panel_main_taxa.tsv`

## Overview

- Taxa compared: 62
- Alignment length before: 14933
- Alignment length after: 14223
- Taxa with ungapped sequence changed: 0
- Candidate outliers with >=1 flag: 10

## Flag Counts

- high_n_fraction: 2
- high_ambiguous_fraction: 0
- short_ungapped_length: 0
- high_consensus_distance_after: 4
- consensus_distance_increase: 0
- high_private_insertion_burden_after: 0
- private_insertion_burden_increase: 0
- high_singleton_substitution_burden_after: 4
- atypical_internal_gap_shift: 0
- ungapped_sequence_changed: 0
- any_outlier_flag: 10
- n_taxa_total: 62
- before_alignment_length: 14933
- after_alignment_length: 14223

## Top Candidate Outliers

| taxon | role | score | flags | after_consensus_distance | after_private_insertion | delta_consensus |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| 24-032160-004_PQ797990__usa_distal/Usa/2024-10-24 | usa_distal | 1 | high_consensus_distance_after | 0.0747 | 0.0000 | 0.0000 |
| 24-032809-001_PQ798038__usa_distal/Usa/2024-11-05 | usa_distal | 1 | high_consensus_distance_after | 0.0745 | 0.0000 | 0.0000 |
| 22-011334-002_OQ968009__american_anchor/Usa/2022-04-13 | american_anchor | 1 | high_consensus_distance_after | 0.0416 | 0.0000 | 0.0000 |
| 22-011138-001_OQ968033__american_anchor/Usa/2022-04-11 | american_anchor | 1 | high_consensus_distance_after | 0.0409 | 0.0000 | 0.0000 |
| 2259-N_PP768254.1__regional_context/Brazil/2023-11-09 | regional_context | 1 | high_singleton_substitution_burden_after | 0.0070 | 0.0061 | 0.0010 |
| 1941-N_OR860365.1__regional_context/Brazil/2023-09-25 | regional_context | 1 | high_singleton_substitution_burden_after | 0.0067 | 0.0299 | -0.0009 |
| 2165-SO_OR852419.1__regional_context/Brazil/2023-10-11 | regional_context | 1 | high_singleton_substitution_burden_after | 0.0061 | 0.0076 | 0.0000 |
| Magdalena3503_OQ683494.1__regional_context/Colombia/2022-11-18 | regional_context | 1 | high_singleton_substitution_burden_after | 0.0055 | 0.0015 | -0.0007 |
| PQ318384.1__regional_context/Bolivia/2023-06-06 | regional_context | 1 | high_n_fraction | 0.0048 | 0.0000 | 0.0000 |
| PQ318387.1__regional_context/Bolivia/2023-03-03 | regional_context | 1 | high_n_fraction | 0.0042 | 0.0000 | 0.0000 |
