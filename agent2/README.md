# Agent 2: experimental-reference evaluation

**Owner: Matthew. AI orchestration: ChatGPT (Astra).** This is the evaluation half
of the Project 2 workflow presented on October 2. It reads Nancy's Agent 1 model
bundle and calculates accuracy against independently verified experimental
structures. The original Agent 1 scripts, metadata, predictions, and comparisons
remain unchanged. All additions are contained in this `agent2/` folder.

The AI follows [AGENT_PROMPT.md](AGENT_PROMPT.md) to inspect inputs, run the tools,
interpret evidence, and prepare the report. The numerical work is deterministic
Python and the official TMscore program. This is a runnable tool workflow for the
ChatGPT/Astra agent; it does not make an unattended model API call or require an
API key. An LLM must not invent or estimate missing metric values.

## Assignment and presentation alignment

Source: `3_protein_TBM_2026_v3.pdf`, slides 80-85; revised presentation content
slides 8-10. Results presentation: **October 9, 2026**.

| Requirement | Agent 2 implementation / deliverable |
|---|---|
| Evaluate predictions against true structures | Score the DOPE-selected TBM model and ranked AF3 output against the same verified reference |
| Multiple quality metrics | Official TM-score and GDT-TS; Kabsch C-alpha RMSD; explicitly labeled C-alpha lDDT |
| Superimpose predicted and true structures | Fixed-correspondence superposed C-alpha PDBs, PNG overlays, transforms, per-residue deviations |
| Three CASP15 template-based single-chain targets | T1147, T1133, T1183; all officially TBM-easy, evaluated as the specified chain/domain |
| Compare with AlphaFold 3 web-server models | Three genuine October 6 server jobs, preserved and evaluated as `alphafold3_server` |
| Traceable results and limitations | CSV/JSON tables, residue mappings, checksums, versions, raw scoring logs, Markdown report |

Nancy's existing AF3 models were produced locally on Hellbender and are labeled
`alphafold3_local`. The new `alphafold3_server` rows fulfill the web-server
comparison requirement. Her TBM-versus-AF3 results describe prediction agreement;
Agent 2 adds each prediction's accuracy against experiment.

## Verified targets and references

| Target | Experimental chain | CASP evaluation domain | Domain length | Primary TBM model |
|---|---|---|---:|---|
| T1147 | 8EM5 A | 12-103 | 92 | T1147.B99990001.pdb |
| T1133 | 8DYS A | 4-427 | 424 | T1133.B99990005.pdb |
| T1183 | 8IFX B | 1-195 | 195 | T1183.B99990005.pdb |

The [official CASP domain summary](https://predictioncenter.org/casp15/domains_summary.cgi)
establishes the intervals, classification, and PDB accessions. See
[reference provenance](references/README.md) for sequence checks, chain choice,
missing coordinates, hashes, and the alternative 8IEY citation for T1183. These
are current RCSB experimental files, not a claim to reproduce the original CASP
assessors' complete native-file processing or ranking.

T1133 is particularly important: its full sequence has 585 residues, Nancy's TBM
contains 414 (5-418), and the chosen experimental chain has 409 observed C-alpha
atoms. The paired evaluated intersection contains 399 residues. The report shows
coverage and domain normalization alongside scores.

## Run locally

From the repository root, Python 3.11+:

```bash
python -m venv agent2/.venv
# Linux / macOS:
source agent2/.venv/bin/activate
# Windows PowerShell instead:
# .\agent2\.venv\Scripts\Activate.ps1
python -m pip install -r agent2/requirements.txt
python -m agent2 audit --config agent2/config_with_server.json
python -m agent2 evaluate --config agent2/config_with_server.json --tmscore /path/to/TMscore --require-complete
python -m unittest discover -s agent2/tests -v
```

`audit` checks the handoff, source hashes, reference assignment, and recorded
template accession exclusion. `evaluate` additionally verifies coordinate
correspondence and computes metrics/images. Outputs default to a fresh
`agent2/results/run-<UTC timestamp>/` folder. You may pass
`--out agent2/results/my_new_run`; an existing folder is always refused.
There is no overwrite flag. Source data are read-only inputs.

The official [pylelab/USalign TMscore](https://github.com/pylelab/USalign) executable
must be supplied with `--tmscore`, configured in `config.json`, or on PATH as
`TMscore`. Use **TMscore**, because the shared residue correspondence is fixed;
TMalign performs a different structural alignment. The runtime never downloads or
installs an executable. Without a valid tool, RMSD and lDDT-CA remain available,
TM-score/GDT-TS are null, and the report states the missing dependency.

On Hellbender, activate Nancy's existing `tbm_tools` environment and use its
`TMscore` executable if available. Agent 2 does not require MODELLER, AF3 GPU
inference, or PyMOL to evaluate the existing files.

The CLI exits 1 for input/scoring workflow errors. A report with clearly recorded
pending assignment items can be produced successfully; use `--require-complete`
to return exit code 2 for any pending item. This prevents a partial result from
being treated as a complete assignment by an automated submission check.

## Completed AlphaFold web-server baseline

All three jobs completed on October 6, 2026. The complete downloads retain five
predictions per target, confidence files, submitted requests, MSAs, templates,
and the server's output terms. See [server provenance and import instructions](baselines/server_submission_2026-10-06/README.md).

The prepared JSON files are planning artifacts. Actual jobs were submitted by
pasting the unchanged sequences into the signed-in Chrome server UI. Exported
job requests are the authoritative record of the inputs. Each confirms seed 1
and one protein copy; the import verifies exact sequence equality. Primary model
selection uses the highest exported server ranking score, with the lowest model
index breaking a tie, before any experimental-reference evaluation.

Run the completed configuration into a fresh output directory:

```bash
python -m agent2 evaluate --config agent2/config_with_server.json --tmscore /path/to/TMscore --require-complete
```

The [completed report](results/with_server_2026-10-06/REPORT.md) and
[metrics CSV](results/with_server_2026-10-06/metrics.csv) contain nine comparisons:
TBM, local AF3 and server AF3 for each target. The earlier six-comparison run
remains intact as a historical snapshot.

### Import future server runs

Submit each unchanged `targets/<T>/<T>.fasta` using the web server. Save the
downloaded outputs and the submitted FASTA under a new folder such as
`agent2/baselines/<T>/`. Do not put files in Nancy's `af3_model/` folders. Preserve
job identifiers, server ranking/selection evidence, and download provenance.
Choose the server-ranked primary model before inspecting its experimental score.

Replace that target's `af3_server: null` in a **copy** of `config.json` with:

```json
{
  "path": "agent2/baselines/T1147/model_0.cif",
  "submitted_fasta": "agent2/baselines/T1147/submitted.fasta",
  "chain": "A",
  "model_index": 0,
  "mapping_mode": "sequence",
  "job_id": "ACTUAL_SERVER_JOB_ID",
  "source_url": "ACTUAL_SERVER_JOB_OR_DOWNLOAD_URL"
}
```

Then run `python -m agent2 evaluate --config agent2/config_with_server.json
--tmscore /path/to/TMscore`. Values in this example are placeholders, not a job
submission. Sequence equality is checked against Nancy's FASTA. If automatic
sequence mapping is ambiguous, use `mapping_mode: "explicit"` and a `mapping`
object mapping target positions to author residue IDs; residue identities are
still checked. A new run freezes a common mask across every supplied method.

## Outputs to use for the October 9 presentation

```text
agent2/results/<run>/
  REPORT.md                    interpretation, coverage, pending items
  metrics.csv                  target x method accuracy table
  summary.json                 machine-readable report
  config_snapshot.json         exact configuration used
  provenance.json              code/input/tool provenance and dependency versions
  <target>/
    residue_mapping.csv        every target position, correspondences, exclusions
    frozen_mask.json           evaluated positions shared by all supplied methods
    reference_common_ca.pdb
    <method>/
      evaluation.json          values, definitions, raw-tool details
      model_common_ca.pdb
      model_superposed_ca.pdb
      superposition.png         mapped C-alpha overlay and error plot
      transform.json
      per_residue_deviation.csv
      tmscore_stdout.txt
      tmscore_stderr.txt
```

For presentation slide 8, use the measured accuracy table, one reference overlay,
and the coverage caveat. For slide 9, compare TBM, local AF3 and server AF3 accuracy
with consistent labels. For slide 10, report completed runs, reproducibility,
and limitations. Do not claim
one method generally outperforms another from three descriptive target examples.

Read [METRICS.md](METRICS.md) before interpreting or quoting any number.

The current [measured report](results/with_server_2026-10-06/REPORT.md)
contains nine completed prediction-versus-experiment comparisons. The
[initial report](results/reference_evaluation_2026-10-06/REPORT.md) preserves the
earlier state before web-server outputs arrived. [Validation](VALIDATION.md)
records the tested behavior.
The official scoring tool's source commit, build command, license, and binary
checksum are recorded in [tool provenance](tool_build_provenance.json).

For this Windows checkout, the validated executable is at
`agent2/.tools/TMscore.exe` and the prepared Python environment is at
`agent2/.venv/Scripts/python.exe`. Both are local, ignored tools; portable source
packages exclude them. Use `agent2/.venv/Scripts/python.exe -m agent2 evaluate
--config agent2/config_with_server.json --tmscore agent2/.tools/TMscore.exe
--require-complete` from this repository root for a new run.

To reproduce the tool on a machine with a C++ compiler, obtain the official
USalign source at commit `1fa25a958fe3095900eba2ef9b48562bdc462251` and compile
`TMscore.cpp` with the repository's documented compiler command. Our Windows
build used isolated Zig 0.15.2; the exact invocation is in the provenance file.
Retain the upstream license when redistributing a source or executable copy.

AlphaFold Server outputs and derivatives remain subject to the
[AlphaFold Server Output Terms of Use](https://alphafoldserver.com/output-terms),
which accompany each download. Our derivatives extract and renumber shared
C-alpha atoms, rigidly superpose them, and calculate metrics and figures; the
downloaded all-atom predictions are preserved unchanged.

## Preservation

Upstream snapshot: `d65f2b5ae2acb150704bf20347604c3de3b896c3` from
`https://github.com/nancyihechu/tbm-protein-agent`, copied October 6, 2026.
Work was prepared on local branch `agent2-evaluation`. Before changes, the parent
workspace received both an untouched ZIP and a verified full-history Git bundle
under `output/Project2_Repository_Backup/`. GitHub publication was authorized on
October 6; the final addition is confined to this `agent2/` folder.
The original source remains recoverable at the upstream commit above. This folder's
Git attributes preserve exact file bytes so recorded checksums survive checkout.

Every upstream tracked file, including the root README, is unchanged from the
original source commit. The earlier README addition was removed when the user
clarified that all additions should stay inside `agent2/`.
Keep reference structures under Agent 2 so they are not fed into Agent 1 template
search. Any future edits to Agent 1 need a new preserved copy or an explicit,
reviewable patch.
