# Historical analysis outputs

This directory contains preserved historical outputs; they have not been recalculated or independently validated for this upload.

- `01_oganov_distance_distribution_five_groups.csv` through `05_oganov_five_group_overall.csv` and `09_oganov_method_settings.txt`: five-group internal Oganov analysis.
- `07_oganov_distance_distribution_five_groups.*` and `08_oganov_similarity_distribution_five_groups.*`: distribution figures.
- `01_all_structures_to_C32_oganov.csv`: per-structure C32 comparison results.
- `C32_oganov_group_summary.csv`, `C32_oganov_high_similarity_bins.csv`, `Group_1_bar.*` through `Group_5_bar.*`, and `Fig6_*` / `Fig7_*`: tables and figures preserved from the historical C32 figure archive. Figure numbers are archive names, not confirmed final manuscript numbers.

The ten carbon structure input CSV files are now available in `../data/`. The C32 reference structure is not included, so the C32 comparison cannot yet be reproduced from this repository alone.

A zero cosine distance does not prove crystal equivalence. Oganov similarity results do not replace StructureMatcher checks against training data or Materials Project.

These carbon results do not constitute the complete paper-wide release of the 385 candidate novelty checks, seven-condition evaluations, or CrysGNN validation data. Those materials require their own verified source files and correspondence to the final manuscript.
