# Validation on October 6, 2026

- **34 tests passed**, with the real official TMscore integration test enabled.
  Command: `python -m unittest discover -s agent2/tests -v`, with
  `AGENT2_TMSCORE_EXECUTABLE` pointing to the local official executable.
- Tests cover rigid-transform invariance, perturbation, reflection rejection,
  ambiguous residue mapping, insertion codes, missing C-alpha atoms, alternate
  occupancy, changed reference hashes, missing metrics, and output preservation.
- The official TMscore build passed identity and rigid-transform checks, including
  user-length normalization. Its source tree was unmodified; see
  [tool build provenance](tool_build_provenance.json).
- End-to-end evaluation succeeded for the preselected TBM, local AF3 and genuine
  web-server models of all three targets: **nine measured comparisons**, each with
  all four metrics. `--require-complete` returned exit code 0 and `COMPLETE`.
  Shared C-alpha counts were 92 (T1147), 399 (T1133), and 195 (T1183).
- The original six superposition/error-plot figures and three new server figures
  were visually inspected. Labels and legends are readable.
- Independent review reproduced all nine RMSD/lDDT values from original
  coordinates and reran official TMscore for all nine matching TM/GDT values.
  Frozen masks, derived PDB numbering, per-residue tables and checksums agree.
- Downloaded server ZIPs passed CRC checks. All three exported sequences, seed 1,
  one-copy inputs and selected CIF sequences were verified. All five samples per
  target, MSAs, templates, confidence files and terms are retained with hashes.
- The server primary models were selected from exported ranking scores before
  experimental scoring. All selected sample 0; T1147 used the recorded lowest-index
  rule to resolve a rounded tie. Default template dates omitted by the server are
  recorded as unverified rather than filled from intended submission settings.
- The main README remains append-only relative to Nancy's original. Its checkout
  bytes remain an exact prefix of the updated file. All 231 other upstream tracked
  files remain unchanged; see [preservation check](PRESERVATION_CHECK.json).
- The original full-history Git bundle passed `git bundle verify`. An untouched
  ZIP and original checkout README are also saved outside the working repository.
- GitHub publication was authorized after evaluation. The published addition is
  confined to `agent2/` and an appended root README section; the preservation
  record describes the pre-publication check. Three user-authorized AlphaFold
  Server jobs were submitted and completed.

The configured evaluation deliverables are complete. The
[current report](results/with_server_2026-10-06/REPORT.md) contains all nine rows;
the earlier [six-comparison report](results/reference_evaluation_2026-10-06/REPORT.md)
remains unchanged. Presentation assembly is separate. References and metric
conventions retain the limitations documented in [METRICS.md](METRICS.md).
