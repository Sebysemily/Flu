<!--toc:start-->
- [Requirements
----------](#requirements)
- [Useful Commands
---------------](#useful-commands)
- [Current Flow
------------](#current-flow)
- [Main Outputs
-------------------](#main-outputs)
- [Global Parameters
-----------------](#global-parameters)
  - [Note about ignored segments in ML zoom](#note-about-ignored-segments-in-ml-zoom)
- [Tools and Core Parameters
-------------------------](#tools-and-core-parameters)
  - [Panel Subsampling (Augur Filter Logic)](#panel-subsampling-augur-filter-logic)
  - [Bayesian Phylodynamics (BEAST 1.10.4)](#bayesian-phylodynamics-beast-1104)
<!--toc:end-->

H5N1 Ecuador Pipeline

Workflow for the phylogenetic analysis of H5N1 in Ecuador. Reproducibility materials—including the custom Python wrapper scripts and exact parameters used for the augur filter subsampling, complete Snakemake pipeline commands, conda environment files, BEAST XML configurations, inferred tree files, and visualization R scripts—are openly available in this repository. The complete repository is permanently archived and publicly available via Zenodo: **[DOI: TO BE ADDED]**. Locally generated consensus sequences and raw sequencing reads from the Ecuadorian samples have been deposited in **[PLACEHOLDER FOR GENBANK/SRA ACCESSIONS]**. All contextual genomes retrieved from GISAID remain subject to their respective access and acknowledgment terms (Shu and McCauley, 2017).

Requirements
----------

- `snakemake` with `conda` or `mamba`.
- Internet connection to install the environments.

Useful Commands
---------------

```bash

# Complete main flow (skips the GSS).
snakemake --cores all --use-conda

# To re-run model selection (GSS), use the GSS target...
snakemake --cores all --use-conda gss_model_selection

# To build the local Ecuador FASTA inputs from MIRA (instead of epi_set fasta)
snakemake --cores all --use-conda build_inputs
```

Current Flow
------------

1. **Inputs:** An `epi_set` file in `data/context` (which must have at least the columns `Virus name`, `Isolate ID`, `Collection date`, and `Segment`). For the Ecuador input: FASTA in `data/input/` (GISAID) or constructed from MIRA (`mira_base_dir`) → `data/input/H5N1_EC_gisaid_from_mira.fasta`. Extra metadata of Ecuadorian samples in `metadata/flu_filtrado.csv`.
2. **Unified Ingestion:** `process_raw_to_segments` (rule in `rules/build_inputs.smk`) writes `data/phylogeny/by_segment/H5N1_{segment}.fasta` (headers `EPI_ISL`) and `metadata/H5N1_context.csv` (Ecuador + regional context, deduplicated by taxon).
3. **ML Phylogeny:** Nextclade by segment and QC across the 8 segments; robust Maximum Likelihood trees with 1000 ultrafast bootstraps (e.g., `iqtree_fast_codon_segment` in `rules/01_ml_trees.smk`) on all samples for exploratory and structural analyses.
4. **Panel:** Subsets sequences based on geographic/temporal metadata roles (Ecuador core, American anchors, regional context) using `build_panel_context_taxa.py` (rule `augur_filter_context` in `rules/01_ml_trees.smk`) → `data/phylogeny/main_panel/`.
5. **Subset + RTT:** IQ-TREE on the final panel; TreeTime with dates from `H5N1_context.csv` (rules in `rules/02_pre_beast.smk`).
6. **Figures:** `plot_8_segments.R` + `tree_aesthetics.R` produce `figures/main_panel_8_segments_collapsed.png` via the rule `plot_8_segments_composite` (in `rules/01_ml_trees.smk`).
7. **BEAST:** XML templates in `template_beast/`, evaluated and run in `results/beast/` (rules in `rules/03_beast.smk`).

Main Outputs
-------------------

- Processed segments: `data/phylogeny/by_segment/H5N1_{segment}.fasta`
- Final metadata: `metadata/H5N1_context.csv`, Ecuador aditional: `metadata/flu_filtrado.csv`
- ML Trees (IQ-TREE): `results/phylogeny/iq-tree/`
- Model Selection (GSS): `results/beast/GSS/`
- BEAST MCC tree and final results: `results/beast/final_run/`
- Main figures folder : `figures/`
- Main BEAST figure: `figures/main_panel_HA_beast_mcc.png`
- QC (Nextclade HA Reports): `results/qc_metrics/nextclade/HA/nextclade_report.csv`
- QC (Nextclade Discarded Sequences): `results/qc_metrics/nextclade/HA_eliminated.csv`
- QC (MAFFT Trim/Gap Discarded): `results/qc_metrics/trim_gap_n/{segment}_eliminated.csv`
- RTT Panel: `data/phylogeny/rtt/panel_main_concat.filtered.nwk` (tip list = pruned panel)
- Ecological Groups and Reassortments: Classification files in `results/clasifications_ecological_groups.csv` and outlier details in `results/possible_reassortants_ignored.csv`.

Global Parameters
-----------------

Files: `config/config.yml`, `config/paths.yml`.

To avoid cluttering the configuration files, not all parameters used across the pipeline are exposed as global variables (many are set directly in their respective rules or scripts). The parameters in the config file are:

- **Reproducibility & Resources:** `random_seed: 39809473`, `max_threads: 16`.
- **Panel Selection:** Maximum thresholds for context taxa (e.g., `PANEL_REGIONAL_CONTEXT_MAX = 250`) are defined natively in the `snakefile`.
- **TreeTime:** `treetime_parameters.clock_filter` (set to 4.0).
- **Outgroup:** `outgroup_root_sample` (set to EPI_ISL_18133029).

### Note about ignored segments in ML zoom

During the process, samples acting as recombinants or endemic strains are identified and excluded from the main H5N1 visualization due to massive divergence in certain segments (e.g., NS or MP). The metadata extracted from GenoFLU for these key samples (condors, boobies, otters, etc.), where the discordance of their genotypes against typical HPAI lineages is observed, is consolidated in the `results/possible_reassortants_ignored.csv` file.

Tools and Core Parameters
-------------------------

The pipeline relies on several standard bioinformatics tools for evolutionary analysis.  Some parameters can be tweaked in `config/config.yml` others in their rule call or rule file:

- **Alignment (MAFFT):** Uses MAFFT's `--auto` flag to align raw segments efficiently (rule `mafft_align_segment` in `rules/01_ml_trees.smk`).

  ```bash
  mafft --auto --thread 16 <input_segment.fasta> > <aligned_segment.mafft>
  ```

- **Quality Control & Clade Assignment (Nextclade):** Used for robust QC and determining specific viral clades, filtering out bad alignments or overly divergent segments (rule `run_nextclade_alignment_all` in `rules/01_ml_trees.smk`).

  ```bash
  nextclade run -D <dataset_dir> -j 16 --output-tsv <output_report.tsv> <input.fasta>
  ```

- **Mutational Risk Profiling (FluMut):** Screens complete and partial genomes for curated adaptive markers of public health concern (mammalian adaptation, altered virulence, antiviral resistance). We explicitly queried a predefined panel of 10 critical molecular signatures (e.g., PB2 E627K, PB2 D701N, NA H275Y, etc.) to evaluate the pandemic potential of the Ecuadorian outbreak. The full list of evaluated markers is exported to `results/phylogeny/flu_mut/critical_markers.txt`.

  ```bash
  flumut -m <output_markers.tsv> -M <output_mutations.tsv> --relaxed <input.fasta>
  ```

### Panel Subsampling (Augur Filter Logic)

Selects context sequences to build the final phylogenetic panel while keeping the dataset computationally manageable for BEAST. In plain English, the selection logic works by:

  1. *Mandatory Inclusion:* Forcing the inclusion of all local Ecuador sequences and a mandatory outgroup to root the tree.
  2. *Regional Context (South America):* Randomly subsampling other South American sequences up to a strict limit of 250 (`regional_context_max`), distributing them evenly across time (by year and month) so no single heavily-sampled outbreak dominates the tree.
  3. *American Anchors (North America):* Subsampling older North American precursor sequences up to a limit of 100 (`american_anchor_max`) to properly anchor the root of the lineage.
  (rule `augur_filter_context` in `rules/01_ml_trees.smk`).

  ```bash
  python code/01_ml_trees/build_panel_context_taxa.py --metadata <metadata.csv> --ha-fasta <HA.fasta> --out-context-taxa <output_taxa.txt> --outgroup EPI_ISL_18133029 --seed 39809473 --regional-context-max 250 --american-anchor-max 100
  ```

- **Maximum Likelihood Trees (IQ-TREE 3):** Infers phylogenies using ModelFinder (`-m TEST` by default) and 1000 ultrafast bootstraps (`-bb 1000` with `-bnni` for robust support). Full ML trees are computed for individual segments (e.g. `iqtree_fast_codon_segment` in `rules/01_ml_trees.smk`), and a pruned codon-partitioned tree is generated for HA (`iqtree_pruned_ha` in `rules/01_ml_trees.smk`).

  ```bash
  iqtree -s <aligned_segment.mafft> -pre <output_prefix> -seed 39809473 -nt 16 -bb 1000 -bnni
  ```

- **Root-to-Tip Regression (TreeTime):** Determines temporal signal and outlier sequences before Bayesian analysis using a clock filter of 4.0 IQRs (rule `run_root_to_tip_subsets_HA` in `rules/02_pre_beast.smk`).

  ```bash
  treetime clock --tree <ML.treefile> --aln <alignment.fasta> --dates <dates.tsv> --clock-filter 4.0 --outdir <output_dir>
  ```

### Bayesian Phylodynamics (BEAST 1.10.4)

Evaluates models and estimates time-scaled trees using BEAGLE (CPU only for maximum cross-compatibility). XML templates are stored in `template_beast/` (rule `run_final_beast` in `rules/03_beast.smk`).

```bash
beast -beagle_CPU -overwrite -threads 1 <model.xml>
```

**Bayesian Inference Priors (BEAST):** The final model employed an uncorrelated lognormal (UCLN) relaxed clock for the sequence alignment and a Bayesian skygrid coalescent prior. Geographic location and host category were modeled as discrete traits using strict clocks and continuous-time Markov-chain (CTMC) models with Bayesian stochastic search variable selection (BSSVS). Standard default priors generated by BEAUti v1.10.4 were used for all trait and clock parameters, including a Poisson prior on the BSSVS indicators to favor a sparse transition graph. A single MCMC chain was run for 200 million states, discarding the first 10% as burn-in to ensure convergence (ESS > 200). The complete XML specification is provided in this repository (`template_beast/final_run.xml`).

**Bayes Factor Calculations:** Bayes factors (BF) for the discrete trait transitions were calculated as the ratio of posterior odds to prior odds for each transition rate being non-zero. The prior odds were computed as $q / (1 - q)$, where $q = (K - 1) / E_{total}$ represents the expected prior probability of inclusion under the Poisson distribution, with $K$ being the number of discrete states and $E_{total}$ the total number of possible pairwise transitions in the graph. In cases where the MCMC sampling yielded a posterior probability of 1.0 for a transition (which mathematically leads to a division by zero in the posterior odds calculation), the Bayes factor is explicitly reported as 'Inf' (infinity) to denote absolute statistical support within the limits of the finite chain length.

- **Reassortment Typing (GenoFLU):** Used to classify the genomic constellations of all 8 segments and detect inter-lineage reassortments (rule `genoflu_multi` in `rules/build_inputs.smk`).

  ```bash
  python GenoFLU-multi.py <combined_segments.fasta> <output_tsv>
  ```
