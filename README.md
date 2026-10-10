# Measuring the genetic content of diagnostic codes

[![reproduce](https://github.com/zxu69461-creator/genetic-content-of-diagnostic-codes/actions/workflows/reproduce.yml/badge.svg)](https://github.com/zxu69461-creator/genetic-content-of-diagnostic-codes/actions/workflows/reproduce.yml)

Aggregate data and analytical code for

> Xu Z†, Cheng J†, Chen X, Luo J, Guo X\*, Zhu X\*. Measuring the genetic content of diagnostic codes: an ankylosing spondylitis-related signal in seronegative rheumatoid arthritis. Manuscript submitted to *Human Genomics* (2026).
>
> † equal contribution; \* corresponding authors

## Contents

`additional_file_3/` is identical, file for file, to **Additional file 3** submitted with the manuscript (archive SHA-256 `f52f92935b8324e4ba20f7a5f6f0ac69976d3703ee24681eb261ab287c3bef1a`). Its [README](additional_file_3/README.md) describes the contents, data conventions and scope. Every file is listed with its SHA-256 in [`SHA256.json`](additional_file_3/SHA256.json), and the first check below verifies them.

| Folder | What it holds |
|---|---|
| `code/` | Core checks: the 1677-endpoint HLA-B\*27:05 atlas counts, the headline HLA estimate and the 22-locus / expanded-MHC 21-locus AS-equivalent estimates |
| `revision/` | The targeted analyses: fixed-panel regressions in FinnGen and UK Biobank, chromosome and locus influence, conditional null simulations, HLA covariance sensitivity, aggregate LD matrices and the recorded genome-wide model-check criteria |
| `tables/` | The 43 supplementary tables in machine-readable form (`table_index.json` maps numbers to titles) |
| `inputs/`, `provenance/` | Frozen aggregate inputs, public acquisition URLs and checksums |

## Running the checks

```bash
cd additional_file_3
python -m pip install -r requirements.txt   # numpy, pandas, scipy, pyarrow
python code/reproduce_core_results.py         # prints a JSON receipt; status PASS
python revision/code/analyze_revision.py      # recomputes the targeted analyses into revision/recomputed/
python revision/code/verify_revision.py       # independent check; status PASS
```

The checks need no network access and no project database, and finish in seconds on a laptop (validated with Python 3.11; versions in `provenance/validation_environment.json`). The `reproduce` workflow runs the same commands on every push; the badge above shows the latest result.

## Not included

Complete genome-wide summary statistics, individual-level reference genotypes and the external toolchains used for the genome-wide fits are not redistributed; obtain them from the original providers listed in the manuscript (Data availability; Table S1). Figure-drawing code depends on the full project database and is not part of this archive.

## Licence

The analytical code is released under the MIT License ([`additional_file_3/LICENSE`](additional_file_3/LICENSE)). Third-party summary data remain under their providers' terms; no new licence is asserted over them.
