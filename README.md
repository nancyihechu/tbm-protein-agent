# Agent 1: Template-Based Protein Structure Modeling

Agent 1 takes a protein sequence (FASTA), finds structural templates in the PDB,
builds 3D models with MODELLER, predicts the same protein with AlphaFold3, and
compares the two. It writes the models, scores, figures and a machine-readable
`metadata.json` for each target, which Agent 2 reads.

The targets are CASP15 proteins. All runs were done on the University of
Missouri Hellbender cluster with SLURM.

## Pipeline

```
 FASTA (targets/<T>/<T>.fasta)
   |
   v
 1. Template search          hhblits vs UniRef30 (3 iterations) -> <T>.a3m
    run_template_search_      hhsearch vs pdb70                 -> <T>.hhr
    target.slurm
   |
   v
 2. Screen templates         top_hits.py, screen_templates.py   -> screening_summary.csv
   |
   v
 3. Prepare templates        prepare_templates.py               -> templates/<pdb>_<chain>.pdb
   |
   v
 4. Alignment                build_pir.py (1 or more hits)      -> alignment/<T>.pir
   |
   +-------------------------------+
   v                               v
 5. MODELLER (5 models,          6. AlphaFold3 3.0.1 (GPU,
    DOPE + GA341)                   1 seed x 5 samples)
    run_modeller_target.slurm       run_af3_target.slurm
    -> tbm_models/                  -> af3_model/
   |                               |
   +---------------+---------------+
                   v
 7. Compare and record           finalize_target.py (calls compare_to_af3.py)
                                 -> comparison/, notes.txt, metadata.json
                   |
                   v
 8. Visuals (optional)           visualize_model.py, af3_quality_report.py
```

## Tools and versions

| Tool | Version | Used for | On Hellbender |
|---|---|---|---|
| HH-suite (hhblits, hhsearch) | 3.3.0 | Template search | `module load rosettafold/v1.1.0`, then `conda activate RoseTTAFold` (env in `/cluster/software/conda-envs-global/rosetta_1_2025/miniconda3`) |
| MODELLER | 10.8 | Homology models, DOPE, GA341, z-DOPE | conda env `tbm_tools` |
| PyMOL (open source) | 3.1.0 | Superposition images, file conversion | conda env `tbm_tools` |
| TM-align and TMscore | 20240303 (bioconda `tmalign`) | Structure comparison | conda env `tbm_tools` |
| Biopython | 1.88 | mmCIF parsing, template cleaning | conda env `tbm_tools` |
| Python, numpy, matplotlib | 3.11.16, 2.4.6, 3.11.2 | Scripts and plots | conda env `tbm_tools` |
| AlphaFold3 | 3.0.1 container | Comparison predictions | `module load alphafold3/alphafold3_v301_deepmind` |

Notes:

- The `tbm_tools` env is at `/cluster/pixstor/dkhf3-lab/users/nancy/envs/tbm_tools`
  and is activated with `module load miniconda3/v26.7.1-1_py314`,
  `eval "$(conda shell.bash hook)"`, `conda activate <path>`. `environment.yml`
  in this repo recreates it. MODELLER needs a free academic license key
  (set `KEY_MODELLER` before installing); the key is not stored in this repo.
- The AlphaFold3 3.0.1 module runs the Singularity image
  `/cluster/software/src/SINGULARITY_IMAGES/alphafold3/ALPHAFOLD_3.0.1/alphafold3_301.sif`,
  pulled from `docker://tuftsttsrt/alphafold3:3.0.1`. Inside the image, the
  Python package metadata reports version 3.0.0.
- The default `alphafold3` module (`v1_alphafold3`) is an unofficial
  reimplementation (kyegomez/AlphaFold3) and was not used.
- The demo run in `examples/` used the older `alphafold3/alphafold3_deepmind`
  module (package 3.0.0).

## Databases

| Database | Version or date | Path on Hellbender |
|---|---|---|
| UniRef30 (hhblits) | UniRef30_2020_06 (June 2020) | `/cluster/VAST/cietest/rosettafold_databases/UniRef30_2020_06/` |
| pdb70 (hhsearch) | Files dated 2020-03-31 and 2020-04-01 | `/cluster/VAST/cietest/alphafold_database/pdb70/` |
| PDB mmCIF (template coordinates) | Local copy | `/cluster/VAST/cietest/alphafold_database/pdb_mmcif/mmcif_files/` |
| AlphaFold3 databases | PDB mmCIF and seqres 2022-09-28, UniRef90 2022_05, MGnify 2022_05, UniProt 2021_04, BFD, Rfam 14.9, RNAcentral, NT 2023-02-23 | `/cluster/VAST/cietest/alphafold3_database/` |
| AlphaFold3 weights | `af3.bin` | `/cluster/VAST/cietest/alphafold3_model_parameters/` (not distributed) |

**Why this avoids the true CASP15 structures.** CASP15 ran in 2022. The pdb70
database used for the template search dates from April 2020, so experimental
structures of the CASP15 targets cannot be among its templates. UniRef30 is from
June 2020. AlphaFold3 was run with its default `max_template_date` of
2021-09-30 (the date used in the AlphaFold 3 paper), so its template search
ignores structures released after that date even though its PDB copy is from
2022-09-28.

## How to run

All commands are run from the repo root on Hellbender. Each target lives in
`targets/<T>/` and starts from `targets/<T>/<T>.fasta`.

```bash
# 1. Template search (one SLURM job per target, partition "general")
sbatch --job-name=ts_T1147 --output=targets/T1147/slurm-%j.out \
    SCRIPTS/run_template_search_target.slurm T1147

# 2. Screen: top hits, then PASS/FAIL
#    PASS = some hit has probability >= 95, identity >= 30 and coverage >= 70
python SCRIPTS/top_hits.py targets/T1147/hhsearch_output/T1147.hhr -n 3
python SCRIPTS/screen_templates.py T1194 T1183 T1159 T1147 T1120 T1133   # -> targets/screening_summary.csv

# Steps 3 to 8 use the tbm_tools env
module load miniconda3/v26.7.1-1_py314
eval "$(conda shell.bash hook)"
conda activate /cluster/pixstor/dkhf3-lab/users/nancy/envs/tbm_tools

# 3. Clean template chains (protein only, first model, no waters, ligands or hydrogens)
python SCRIPTS/prepare_templates.py 6N9A_B 6S84_C 2A6A_A -o targets/T1183/templates

# 4. PIR alignment from hhsearch hit numbers (several hits give a multi-template alignment)
python SCRIPTS/build_pir.py targets/T1183/hhsearch_output/T1183.hhr targets/T1183/T1183.fasta \
    --hit 1 3 5 -o targets/T1183/alignment/T1183.pir

# 5. MODELLER: target, first residue, last residue, template codes
sbatch --job-name=mod_T1183 --output=targets/T1183/slurm-%j.out \
    SCRIPTS/run_modeller_target.slurm T1183 1 200 6n9a_B 6s84_C 2a6a_A

# 6. AlphaFold3 (input JSON: targets/<T>/af3_model/<T>.json; partition dkhf3-lab-gpu,
#    account dkhf3-lab, 1 L40S GPU; retries with run_script_a3_unified_mem on GPU OOM)
sbatch --job-name=af3_T1183 --output=targets/T1183/af3_model/slurm-%j.out \
    SCRIPTS/run_af3_target.slurm T1183

# 7. Compare and write notes.txt and metadata.json
python SCRIPTS/finalize_target.py T1183 --start 1 --end 200 --hits 1 3 5 \
    --modeller-job <id> --af3-job <id> --template-choice "..." --note "..."

# Pairwise comparison on its own (any two PDB or mmCIF files)
python SCRIPTS/compare_to_af3.py <tbm_model.pdb> <af3_model.cif> --out-dir comparison_results

# 8. Visuals for one model (DOPE coloring) and the full report used for the demo
python SCRIPTS/visualize_model.py --model <model.pdb> --template <template.pdb> --out-dir <dir>
python SCRIPTS/af3_quality_report.py      # defaults to examples/test_target/

# Job status and the newest SLURM log in the repo
SCRIPTS/check_job.sh
```

The AlphaFold3 input JSON is:

```json
{"name": "T1183", "modelSeeds": [1], "dialect": "alphafold3", "version": 1,
 "sequences": [{"protein": {"id": "A", "sequence": "<sequence>"}}]}
```

## Folder layout

```
.
├── README.md
├── ALPHAFOLD3_OUTPUT_NOTICE.txt     Terms notice for the AlphaFold3 outputs
├── environment.yml                  tbm_tools conda env
├── SCRIPTS/
│   ├── fasta_reader.py              Read FASTA files
│   ├── run_template_search_target.slurm
│   ├── top_hits.py                  Parse .hhr files
│   ├── screen_templates.py          PASS/FAIL screen -> screening_summary.csv
│   ├── prepare_templates.py         Copy mmCIF, write cleaned chain PDB
│   ├── build_pir.py                 PIR alignment from 1 or more hhsearch hits
│   ├── run_modeller.py              MODELLER automodel, DOPE and GA341
│   ├── run_modeller_target.slurm
│   ├── run_af3_target.slurm         AlphaFold3 3.0.1
│   ├── compare_to_af3.py            TM-align comparison and superposition image
│   ├── finalize_target.py           Comparison, notes.txt, metadata.json
│   ├── visualize_model.py           Model on template, DOPE-colored image
│   ├── af3_quality_report.py        Presentation figures and tables
│   └── check_job.sh
├── targets/
│   ├── candidates.fasta             T1194, T1183, T1159, T1147, T1120
│   ├── candidates_backup.fasta      T1112, T1133
│   ├── screening_summary.csv        Screen of the 6 candidates
│   ├── top_hits_summary.csv         Top 3 hits for T1104, T1122, T1155
│   ├── T1147/ T1133/ T1183/         Modeled targets (layout below)
│   └── T1104/ T1122/ T1155/ T1194/ T1159/ T1120/   Screened only: FASTA and hhsearch_output/
└── examples/
    └── test_target/                 Demo run on T1112 (see its README.md)
```

Each modeled target:

```
targets/<T>/
├── <T>.fasta
├── hhsearch_output/   <T>.a3m, <T>.hhr
├── templates/         <pdb>.cif (copied), <pdb>_<chain>.pdb (cleaned)
├── alignment/         <T>.pir
├── tbm_models/        <T>.B9999000[1-5].pdb, trimmed PIR used by MODELLER
├── tbm_models.log     MODELLER log with the DOPE and GA341 summary
├── af3_model/         <T>.json (input), af3_run.log, <t>/ (AlphaFold3 output)
├── comparison/        TM-align output, superposition.png/.pse, per-residue distances
├── notes.txt          Templates, residue range, problems
└── metadata.json      Everything above in machine-readable form
```

## metadata.json

Paths inside `metadata.json` are relative to `targets/<T>/`.

| Field | Meaning |
|---|---|
| `target`, `description`, `length`, `fasta` | Target ID, FASTA header, sequence length, FASTA file |
| `created` | Date the file was written |
| `template_search` | Databases, `.hhr` file, and why the templates were chosen (`choice`) |
| `templates[]` | Per template: `pdb_chain`, `hhsearch_rank`, `probability`, `evalue`, `identity_pct`, `query_range`, `template_seqres_range`, `resolution_A` (null for NMR), `protein`, `file` |
| `modeled_range` | First and last residue in the TBM models (full-sequence numbering) |
| `alignment` | PIR file |
| `modeller` | `version`, `n_models`, `slurm_job`, `models[]` (`file`, `dope`, `ga341`), `best_model`, `best_dope`, `best_ga341` |
| `alphafold3` | `module`, `slurm_job`, `model` (top-ranked), `seeds`, `samples`, `ptm`, `ranking_score`, `fraction_disordered`, `mean_plddt` (all residues), `mean_plddt_modeled_range` |
| `comparison.tm_align` | Sequence-independent: `tm_score_norm_tbm`, `tm_score_norm_af3`, `rmsd`, `aligned_length`, `pairs_within_5A` |
| `comparison.residue_matched` | TMscore over the modeled range, each residue compared with the same residue: `tm_score`, `rmsd`, `gdt_ts`, `residues_in_common` |
| `problems[]` | Free-text issues, also listed in `notes.txt` |

Example (`targets/T1147/metadata.json`, `models[]` and `problems[]` shortened):

```json
{
  "target": "T1147",
  "description": "T1147 KEK_05752, MmpS5, Mycobacterium thermoresistibile, 103 residues",
  "length": 103,
  "fasta": "T1147.fasta",
  "created": "2026-10-05",
  "template_search": {
    "hhblits_db": "UniRef30_2020_06",
    "hhsearch_db": "pdb70 (2020-04)",
    "hhr": "hhsearch_output/T1147.hhr",
    "choice": "Only qualifying hit (screen PASS): prob 99.6, 55% identity."
  },
  "templates": [
    {
      "pdb_chain": "2LW3_A", "hhsearch_rank": 1, "probability": 99.63, "evalue": 9.3e-20,
      "identity_pct": 55, "query_range": [15, 103], "template_seqres_range": [2, 90],
      "resolution_A": null, "protein": "Putative membrane protein mmpS4",
      "file": "templates/2lw3_A.pdb"
    }
  ],
  "modeled_range": [1, 103],
  "alignment": "alignment/T1147.pir",
  "modeller": {
    "version": "10.8", "n_models": 5, "slurm_job": "18319194",
    "models": [{"file": "tbm_models/T1147.B99990001.pdb", "dope": -6472.051, "ga341": 0.999}],
    "best_model": "tbm_models/T1147.B99990001.pdb", "best_dope": -6472.051, "best_ga341": 0.999
  },
  "alphafold3": {
    "module": "alphafold3/alphafold3_v301_deepmind", "slurm_job": "18319106",
    "model": "af3_model/t1147/t1147_model.cif", "seeds": 1, "samples": 5,
    "ptm": 0.79, "ranking_score": 0.84, "fraction_disordered": 0.12,
    "mean_plddt": 91.21, "mean_plddt_modeled_range": 91.21
  },
  "comparison": {
    "dir": "comparison/T1147.B99990001_vs_t1147_model",
    "tm_align": {"tm_score_norm_tbm": 0.62892, "tm_score_norm_af3": 0.62892, "rmsd": 2.85,
                 "aligned_length": 90, "pairs_within_5A": 82},
    "residue_matched": {"tm_score": 0.5955, "rmsd": 8.448, "gdt_ts": 0.5534, "residues_in_common": 103}
  },
  "problems": ["Residues 1-14 have no template; MODELLER built them without template restraints."]
}
```

## Results

| Target | Template(s) | Identity | Query coverage | Residues modeled | Best DOPE (model) | AF3 mean pLDDT | AF3 pTM | TM-score vs AF3, residue-matched | TM-score vs AF3, TM-align |
|---|---|---|---|---|---|---|---|---|---|
| T1147 (103 aa) | 2LW3_A | 55% | 15-103 (86%) | 1-103 | -6472.1 (model 1) | 91.2 | 0.79 | 0.596 | 0.629 |
| T1133 (585 aa) | 3WJ9_B | 37% | 5-418 (71%) | 5-418 | -42136.5 (model 5) | 88.5 (95.6 over 5-418) | 0.78 | 0.893 | 0.917 |
| T1183 (200 aa) | 6N9A_B, 6S84_C, 2A6A_A | 24%, 25%, 25% | 1-199 (100%) | 1-200 | -21673.1 (model 5) | 93.2 | 0.90 | 0.744 | 0.751 |

- All 15 TBM models have GA341 of 0.999 or 1.000.
- TM-scores are normalized by the TBM model length. "Residue-matched" (TMscore)
  compares each residue only with the same residue in the AlphaFold3 model;
  "TM-align" finds the best structural match regardless of residue numbering.
- RMSD, residue-matched / TM-align: T1147 8.45 / 2.85 A (90 aligned),
  T1133 4.27 / 1.99 A (403 aligned), T1183 5.95 / 3.31 A (188 aligned).
- T1147: the core (residues 10-94) agrees with AlphaFold3. Residues 1-9 and
  95-103 account for most of the residue-matched RMSD; residues 95-103 are up to
  41 A from their AlphaFold3 positions.
- Demo run (T1112, `examples/test_target/`): 4EHX_A at 14% identity, residues
  103-460 modeled, TM-score vs AlphaFold3 0.554 residue-matched and 0.643 with
  TM-align.

## Target selection

Candidates were screened with the template search above. A target passes when
at least one hhsearch hit has probability >= 95, identity >= 30% and query
coverage >= 70%; the best such hit (highest identity, then coverage, then
probability) is its template (`best_template` in `screening_summary.csv`).

| Target | Length | Rank-1 hit | Best qualifying template | Result |
|---|---|---|---|---|
| T1194 | 168 | 1V7N_X (prob 16.1) | none | FAIL |
| T1183 | 200 | 6N9A_B (prob 100.0, 24%) | none | FAIL (identity only) |
| T1159 | 160 | 5HB6_B (prob 45.6) | none | FAIL |
| T1147 | 103 | 2LW3_A (prob 99.6, 55%) | 2LW3_A | PASS |
| T1120 | 235 | 6A6W_A (prob 77.8) | none | FAIL |
| T1133 | 585 | 6GSN_p (prob 99.9, 17%) | 3WJ9_B (rank 3, prob 99.9, 37%) | PASS |

An earlier batch (T1104, T1122, T1155, `targets/top_hits_summary.csv`) had no
hit above 89% probability. T1112 was the demo target.

T1147 and T1133 were modeled because they pass. **T1183 is an exception**: it
fails the 30% identity cutoff (its full-length templates have 24 to 25%
identity), but 22 hits reach 95% probability and several cover residues 1-199,
so it was modeled with a three-template alignment. For T1133 only 3WJ9_B
(*S. pombe* eIF2A) was used; the two higher-ranked hits (6GSN_p, 6GSM_p, chains
from cryo-EM ribosome complexes at 17% identity) were excluded.

## Known limitations

- **Low sequence identity.** T1183 templates are 24 to 25% identical and the
  T1112 demo template is 14%. At this level, parts of the alignment can be
  shifted along the sequence; in the demo, only 137 of 307 TM-align pairs had
  the same residue number.
- **Residues without a template.** T1133 residues 1-4 and 419-585 are not in the
  TBM model at all. T1147 residues 1-14 have no template and were built by
  MODELLER without template restraints. T1183 residue 200 has no template.
- **Template choices.** The three T1183 templates are the same protein
  (*Thermotoga maritima* TsaB) in different structures, so they add little
  structural variety. The T1147 template 2LW3 is an NMR ensemble; only its first
  model was used.
- **Old search databases.** pdb70 and UniRef30 date from 2020, so templates and
  homologs deposited later are missing.
- **AlphaFold3 was run locally on Hellbender** (module
  `alphafold3/alphafold3_v301_deepmind`, one seed, five samples), not on the
  AlphaFold Server.
- **The comparison is to AlphaFold3, not to the true structure.** TM-scores and
  RMSDs measure agreement with the AlphaFold3 prediction, and DOPE and GA341 are
  statistical estimates of model quality.

## For Agent 2

Start from `targets/<T>/metadata.json` for T1147, T1133 and T1183. All paths in
it are relative to `targets/<T>/`.

| What | Where | Format |
|---|---|---|
| Best TBM model | `modeller.best_model` | PDB, chain A, residues numbered as in the full FASTA sequence (T1133 starts at 5) |
| All TBM models and scores | `modeller.models[]` | PDB; `dope` (lower is better), `ga341` (0 to 1) |
| AlphaFold3 top model | `alphafold3.model` | mmCIF, chain A, residues 1 to length; per-atom pLDDT in the B-factor column |
| AlphaFold3 confidences | `af3_model/<t>/<t>_summary_confidences.json`, `<t>_confidences.json` | JSON (pTM, ranking score; full pLDDT and PAE) |
| Other AlphaFold3 samples | `af3_model/<t>/seed-1_sample-[0-4]/model.cif` | mmCIF |
| Comparison | `comparison.*` and `comparison/<dir>/summary.json` | JSON; `per_residue_distances.csv` for CA distances |
| Templates and alignment | `templates[]`, `alignment` | JSON, PIR |
| Caveats | `problems[]`, `notes.txt` | Text |
| Screening of all candidates | `targets/screening_summary.csv` | CSV: target, length, top_hit_pdb, best_template, probability, evalue, identity, coverage, combined_coverage, status |

`.a3m` files are hhblits alignments in A3M format and `.hhr` files are hhsearch
results. AlphaFold3 MSA features (`*_data.json`) and per-sample PAE matrices are
not in the repo because of their size.

## AlphaFold3 outputs

AlphaFold3 predictions in this repo are AlphaFold 3 Output and are provided
under the [AlphaFold 3 Output Terms of Use](https://github.com/google-deepmind/alphafold3/blob/main/OUTPUT_TERMS_OF_USE.md).
Each output folder keeps its `TERMS_OF_USE.md`, and `ALPHAFOLD3_OUTPUT_NOTICE.txt`
covers the derived files (trimmed structures, superpositions and figures).
If you use them, cite: Abramson, J. et al. Accurate structure prediction of
biomolecular interactions with AlphaFold 3. *Nature* (2024).
