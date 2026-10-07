"""Read Agent 1 artifacts, freeze correspondence, score, and write separate outputs.

The LLM follows AGENT_PROMPT.md. Scientific values are computed here and in
metrics.py; this module does not call an LLM or a prediction service.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path, rows, fields):
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def repo_path(root, value):
    path = (root / value).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError(f"Input path leaves repository: {value}")
    return path


def read_fasta(path):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    headers = [line for line in lines if line.startswith(">")]
    sequence = "".join(line.strip() for line in lines if line.strip() and not line.startswith(">"))
    if len(headers) != 1 or not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
        raise ValueError(f"Expected one canonical protein FASTA: {path}")
    return sequence


def reserve_output(root, requested):
    base = (root / "agent2/results").resolve()
    name = requested or "agent2/results/run-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = repo_path(root, name)
    if path == base or not path.is_relative_to(base):
        raise ValueError("Outputs must be a new subdirectory of agent2/results")
    path.mkdir(parents=True, exist_ok=False)
    return path


def load_target(root, entry):
    metadata_path = repo_path(root, entry["metadata"])
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    target = entry["target"]
    if metadata["target"] != target:
        raise ValueError("Target ID disagrees with Agent 1 metadata")
    fasta_path = repo_path(root, str(metadata_path.parent.relative_to(root) / metadata["fasta"]))
    sequence = read_fasta(fasta_path)
    if len(sequence) != metadata["length"]:
        raise ValueError("FASTA length disagrees with Agent 1 metadata")
    start, end = entry["domain_range"]
    if not (1 <= start <= end <= len(sequence)) or entry["normalization_length"] != end - start + 1:
        raise ValueError("Invalid CASP domain interval/normalization length")
    models = []
    for method, relative, chain in [
        ("tbm", metadata["modeller"]["best_model"], entry["tbm_chain"]),
        ("alphafold3_local", metadata["alphafold3"]["model"], entry["af3_chain"]),
    ]:
        path = repo_path(root, str(metadata_path.parent.relative_to(root) / relative))
        if not path.is_file():
            raise FileNotFoundError(f"Missing {method}: {path}")
        models.append({"method": method, "path": path, "chain": chain, "model_index": 0,
                       "mapping_mode": "author_number_equals_target"})
    server = entry.get("af3_server")
    if server:
        # Server provenance is explicit; local cluster outputs never get relabeled.
        if not server.get("job_id") or not server.get("source_url"):
            raise ValueError("AF3 server output needs job_id and source_url")
        submitted = read_fasta(repo_path(root, server["submitted_fasta"]))
        if submitted != sequence:
            raise ValueError("AF3 server submitted sequence differs from Agent 1 FASTA")
        models.append({"method": "alphafold3_server", "path": repo_path(root, server["path"]),
                       "chain": server["chain"], "model_index": server.get("model_index", 0),
                       "mapping_mode": server.get("mapping_mode", "sequence"),
                       "mapping": server.get("mapping"), "job_id": server["job_id"],
                       "source_url": server["source_url"]})
    reference = entry["reference"]
    ref_path = repo_path(root, reference["path"])
    if reference.get("verified") is not True or not reference.get("source_url"):
        raise ValueError("Reference assignment must be verified with source evidence")
    if sha256(ref_path) != reference["sha256"]:
        raise ValueError("Reference checksum changed; verify assignment before updating the manifest")
    native_id = reference["pdb_id"].upper()
    templates = [item["pdb_chain"].split("_")[0].upper() for item in metadata["templates"]]
    if native_id in templates:
        raise ValueError("Experimental reference accession appears in Agent 1 template inputs")
    return metadata, sequence, fasta_path, models, ref_path


def audit_target(root, entry):
    metadata, sequence, fasta_path, models, ref_path = load_target(root, entry)
    return {
        "target": entry["target"], "status": "INPUTS_PRESENT", "target_length": len(sequence),
        "domain": entry["domain_id"], "domain_range": entry["domain_range"],
        "classification": entry["classification"], "reference": entry["reference"],
        "models": [{**{k: v for k, v in model.items() if k != "path"},
                    "path": str(model["path"].relative_to(root)), "sha256": sha256(model["path"])} for model in models],
        "fasta_sha256": sha256(fasta_path),
        "metadata_sha256": sha256(repo_path(root, entry["metadata"])),
        "tbm_selection": "Agent 1 best_model selected by DOPE before experimental evaluation",
        "af3_selection": "Agent 1 top-ranked model, fixed before experimental evaluation",
        "tbm_modeled_range": metadata["modeled_range"],
        "template_accession_check": "No exact reference accession in recorded templates; this is not a training-data audit",
        "af3_server_status": "PROVIDED" if entry.get("af3_server") else "NOT_PROVIDED",
        "agent1_problems": metadata.get("problems", []),
    }


def mapped_structure(path, chain, model_index, sequence, mode, explicit=None, altloc=None):
    from .metrics import read_ca_structure, map_to_target
    structure = read_ca_structure(path, chain, model_index, altloc=altloc)
    if mode == "author_number_equals_target":
        explicit = {str(residue["number"]): residue["residue_id"] for residue in structure["residues"]
                    if not residue["insertion_code"] and 1 <= residue["number"] <= len(sequence)}
        # Identity checks in map_to_target prevent trusting numbers alone.
        if len(explicit) != len(structure["residues"]):
            raise ValueError("Number-based mapping would silently discard residues; supply an explicit mapping")
    elif mode not in {"sequence", "explicit"}:
        raise ValueError(f"Unknown mapping mode: {mode}")
    if mode == "explicit" and not explicit:
        raise ValueError("Explicit mapping is required")
    return structure, map_to_target(sequence, structure, explicit)


def render_pair(path, target, method, positions, reference, fitted, errors):
    import numpy as np
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    fig = Figure(figsize=(12, 5.3), layout="constrained")
    FigureCanvasAgg(fig)
    ax = fig.add_subplot(121, projection="3d")
    # Break curves at missing residues instead of drawing imaginary connections.
    splits = np.where(np.diff(positions) != 1)[0] + 1
    for index, ids in enumerate(np.split(np.arange(len(positions)), splits)):
        ax.plot(*reference[ids].T, color="#E3AE4B", linewidth=2,
                label="Experimental reference" if index == 0 else None)
        ax.plot(*fitted[ids].T, color="#237CA8", linewidth=1.5,
                label=method if index == 0 else None)
    points = np.vstack([reference, fitted])
    center = (points.max(axis=0) + points.min(axis=0)) / 2
    radius = max(float(np.ptp(points, axis=0).max()) / 2, 1)
    ax.set(xlim=(center[0]-radius, center[0]+radius), ylim=(center[1]-radius, center[1]+radius),
           zlim=(center[2]-radius, center[2]+radius), title="C-alpha superposition; fixed shared residues")
    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()
    ax.legend(loc="lower left", fontsize=8)
    error_ax = fig.add_subplot(122)
    for ids in np.split(np.arange(len(positions)), splits):
        error_ax.plot(np.asarray(positions)[ids], errors[ids], color="#237CA8")
    error_ax.set(xlabel="Target residue", ylabel="C-alpha deviation after Kabsch fit (angstrom)",
                 title="Per-residue structural deviation", ylim=(0, max(float(max(errors))*1.1, 1)))
    error_ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"{target} | {method} versus experiment", fontsize=15, fontweight="bold")
    fig.savefig(path, dpi=180)


def evaluate_target(root, entry, output, tmscore):
    import numpy as np
    from .metrics import kabsch_rmsd, lddt_ca, run_tmscore, _write_ca_pdb
    metadata, sequence, fasta, models, ref_path = load_target(root, entry)
    destination = output / entry["target"]
    destination.mkdir()
    spec = entry["reference"]
    reference, reference_map = mapped_structure(ref_path, spec["chain"], spec.get("model_index", 0), sequence,
                                                spec["mapping_mode"], spec.get("mapping"), spec.get("altloc"))
    structures, maps = {}, {}
    for model in models:
        structures[model["method"]], maps[model["method"]] = mapped_structure(
            model["path"], model["chain"], model["model_index"], sequence,
            model["mapping_mode"], model.get("mapping"))
    start, end = entry["domain_range"]
    domain = set(range(start, end + 1))
    common = set(reference_map) & domain
    for mapping in maps.values():
        common &= set(mapping)
    positions = sorted(common)
    if len(positions) < 3:
        raise ValueError("Fewer than three experimentally observed positions shared by all supplied methods")
    # Freeze mask before fitting; every method is evaluated at these same positions.
    map_rows = []
    for position in range(1, len(sequence) + 1):
        row = {"target_position": position, "amino_acid": sequence[position-1],
               "in_casp_domain": position in domain, "evaluated": position in common,
               "reference_residue": reference_map.get(position, {}).get("residue_id", ""),
               "reference_altloc": reference_map.get(position, {}).get("altloc", ""),
               "reference_occupancy": reference_map.get(position, {}).get("occupancy", "")}
        for method, mapping in maps.items():
            row[method + "_residue"] = mapping.get(position, {}).get("residue_id", "")
            row[method + "_altloc"] = mapping.get(position, {}).get("altloc", "")
        reasons = []
        if position not in domain:
            reasons.append("outside_CASP_domain")
        if position not in reference_map:
            reasons.append("no_reference_CA")
        reasons.extend("no_" + method + "_CA" for method, mapping in maps.items() if position not in mapping)
        row["exclusion_reason"] = ";".join(reasons)
        map_rows.append(row)
    write_csv(destination / "residue_mapping.csv", map_rows, list(map_rows[0]))
    write_json(destination / "frozen_mask.json", {"target_positions": positions,
               "methods": list(maps), "normalization_length": entry["normalization_length"],
               "scope": "Intersection of experimental and all supplied model C-alpha positions within CASP domain"})
    ref_xyz = np.asarray([reference_map[p]["xyz"] for p in positions])
    ref_ca = destination / "reference_common_ca.pdb"
    _write_ca_pdb(ref_ca, positions, sequence, ref_xyz)
    rows = []
    for model in models:
        method = model["method"]
        method_dir = destination / method
        method_dir.mkdir()
        xyz = np.asarray([maps[method][p]["xyz"] for p in positions])
        fit = kabsch_rmsd(xyz, ref_xyz)
        lddt = lddt_ca(xyz, ref_xyz)
        ca_path = method_dir / "model_common_ca.pdb"
        _write_ca_pdb(ca_path, positions, sequence, xyz)
        _write_ca_pdb(method_dir / "model_superposed_ca.pdb", positions, sequence, fit["fitted"])
        score = run_tmscore(ca_path, ref_ca, entry["normalization_length"], len(positions), tmscore)
        (method_dir / "tmscore_stdout.txt").write_text(score.get("stdout", ""), encoding="utf-8")
        (method_dir / "tmscore_stderr.txt").write_text(score.get("stderr", score.get("diagnostic", "")), encoding="utf-8")
        write_json(method_dir / "transform.json", {"convention": "row-vector: fitted = xyz @ rotation + translation",
                   "rotation": fit["rotation_row_vector"].tolist(), "translation": fit["translation"].tolist()})
        write_csv(method_dir / "per_residue_deviation.csv", [
            {"target_position": p, "ca_deviation_angstrom": float(d)} for p, d in zip(positions, fit["distances_angstrom"])],
            ["target_position", "ca_deviation_angstrom"])
        render_pair(method_dir / "superposition.png", entry["target"], method, positions, ref_xyz,
                    fit["fitted"], fit["distances_angstrom"])
        complete_metrics = all(value is not None for value in [score.get("tm_score"), score.get("gdt_ts"), lddt["score"], fit["rmsd_angstrom"]])
        row = {"target": entry["target"], "method": method, "status": "SCORED" if complete_metrics else "PARTIAL_METRICS",
               "model": str(model["path"].relative_to(root)), "reference": spec["pdb_id"] + ":" + spec["chain"],
               "domain": entry["domain_id"], "target_length": len(sequence),
               "normalization_length": entry["normalization_length"], "matched_ca": len(positions),
               "observed_reference_ca_in_domain": len(set(reference_map) & domain),
               "model_ca_in_full_target": len(maps[method]),
               "evaluated_fraction_of_domain": len(positions) / entry["normalization_length"],
               "model_coverage_of_full_target": len(maps[method]) / len(sequence),
               "tm_score": score.get("tm_score"), "gdt_ts": score.get("gdt_ts"),
               "lddt_ca_common_mask": lddt["score"], "rmsd_ca_angstrom": fit["rmsd_angstrom"],
               "tmscore_status": score["status"]}
        write_json(method_dir / "evaluation.json", {"summary": row, "lddt_details": lddt,
                   "tmscore_details": score, "model_sha256": sha256(model["path"]),
                   "reference_sha256": sha256(ref_path), "reference_altloc": spec.get("altloc", "highest occupancy; ties require explicit selection"),
                   "missing_ca_records": structures[method]["missing_ca"],
                   "reference_missing_ca_records": reference["missing_ca"]})
        rows.append(row)
    if not entry.get("af3_server"):
        rows.append({"target": entry["target"], "method": "alphafold3_server", "status": "NOT_PROVIDED",
                     "reference": spec["pdb_id"] + ":" + spec["chain"], "domain": entry["domain_id"],
                     "tm_score": None, "gdt_ts": None, "lddt_ca_common_mask": None, "rmsd_ca_angstrom": None})
    return rows


def format_number(value):
    return "pending" if value is None else f"{value:.4f}"


def report(summary):
    lines = ["# Agent 2 experimental-reference evaluation", "", f"Status: **{summary['status']}**", "",
             "Owner: Matthew. Orchestration: ChatGPT (Astra), with deterministic scoring tools.", "",
             "## Input and assignment status", "",
             "| Target | Reference | CASP domain | Classification | AF3 server |",
             "|---|---|---|---|---|"]
    for item in summary["audit"]:
        if item.get("status") == "ERROR":
            lines.append(f"| {item['target']} | ERROR | {item['error']} | | |")
        else:
            lines.append(f"| {item['target']} | {item['reference']['pdb_id']}:{item['reference']['chain']} | "
                         f"{item['domain_range'][0]}-{item['domain_range'][1]} | {item['classification']} | {item['af3_server_status']} |")
    lines.extend(["", "## Accuracy against experiment", "",
                  "| Target | Method | Status | Paired CA/domain length | TM-score | GDT-TS (0-1) | lDDT-CA | RMSD (A) |",
                  "|---|---|---|---:|---:|---:|---:|---:|"])
    for row in summary["metrics"]:
        count = f"{row['matched_ca']}/{row['normalization_length']}" if "matched_ca" in row else "pending"
        scores = " | ".join(format_number(row.get(key)) for key in ["tm_score", "gdt_ts", "lddt_ca_common_mask", "rmsd_ca_angstrom"])
        lines.append(f"| {row['target']} | {row['method']} | {row['status']} | {count} | {scores} |")
    if not summary["metrics"]:
        lines.append("| No scores in audit mode | | | | | | | |")
    lines += ["", "## Interpretation and limits", "",
              "- TBM is Agent 1's DOPE-selected model; the local AF3 model is Agent 1's top-ranked output. Selection precedes experimental scoring.",
              "- Each target uses one fixed intersection of observed reference and supplied model C-alpha residues within its CASP domain. No score-based residue pruning.",
              "- TM-score and GDT-TS use the declared full CASP domain length. Missing positions remain in that denominator. These are our reproducible calculations, not official CASP assessor scores.",
              "- RMSD uses a Kabsch fit over the fixed paired C-alpha mask. lDDT-CA uses that mask and is not all-atom lDDT or pLDDT. See METRICS.md for definitions.",
              "- T1133's TBM covers 414/585 target residues; its experimental reference and CASP domain also cover only part of the full sequence. High regional agreement does not establish full-length accuracy.",
              "- Local AF3 models were generated on Hellbender. Only rows labeled alphafold3_server satisfy the web-server comparison requirement; inspect their input status above.",
              "- All methods use the same mask per run. Adding a new method can change the intersection; compare methods within one run.",
              "- Existing Agent 1 scores compare TBM with AF3; those remain separate from this experimental-reference evaluation.",
              "- These are retrospective predictions of CASP15 targets. Template-date restrictions do not establish training independence or reproduce a blind CASP submission.",
              "- Server and local AF3 runs may differ in model release, databases, sampling and other settings. Score differences cannot be attributed solely to server versus local execution."]
    if any(row.get("method") == "alphafold3_server" and row.get("status") == "SCORED" for row in summary["metrics"]):
        lines += ["", "## AlphaFold Server provenance and output notice", "",
                  "- Server primary models were selected by highest exported ranking_score before experimental scoring; ties use the lowest server model-file index. All five candidates and original download archives are retained.",
                  "- See the supplied server config and baseline import provenance for submitted sequences, seeds, job URLs, template settings and checksums. A missing template date in an exported request is a server default, not an explicitly verified date.",
                  "- AlphaFold Server outputs and derivatives are subject to the [AlphaFold Server Output Terms of Use](https://alphafoldserver.com/output-terms). The original terms files accompany the downloaded outputs. Derivatives here extract common C-alpha residues, renumber those pairs, rigidly superpose coordinates, and calculate metrics and figures; original predictions are unmodified.",
                  "- Cite Abramson et al., [Accurate structure prediction of biomolecular interactions with AlphaFold 3](https://doi.org/10.1038/s41586-024-07487-w), Nature (2024)."]
    lines += ["", "## Remaining deliverables", ""]
    lines.extend("- " + gap for gap in summary["pending"])
    if not summary["pending"]:
        lines.append("- None for the configured evaluation outputs; integrate the figures and tables into the October 9 presentation.")
    if summary["errors"]:
        lines += ["", "## Errors", ""] + ["- " + error for error in summary["errors"]]
    lines += ["", "## Evidence", "",
              "- [CASP15 target assignments](https://predictioncenter.org/casp15/targetlist.cgi)",
              "- [CASP15 domain definitions and classification](https://predictioncenter.org/casp15/domains_summary.cgi)",
              "- Source checksums, software versions and configuration are saved with this report.", ""]
    return "\n".join(lines)


def run(root, args):
    root = root.resolve()
    config_path = repo_path(root, args.config)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    targets = config["targets"]
    ids = [item["target"] for item in targets]
    if config.get("schema_version") != 1 or len(ids) != 3 or len(set(ids)) != 3:
        raise ValueError("Configuration must declare schema 1 and exactly three distinct targets")
    if any(not target.replace("_", "").isalnum() for target in ids):
        raise ValueError("Target IDs must be safe directory names")
    output = reserve_output(root, args.out)
    write_json(output / "config_snapshot.json", config)
    summary = {"schema_version": 1, "status": "AUDIT", "command": args.command,
               "created_utc": datetime.now(timezone.utc).isoformat(), "source_commit": config["source_commit"],
               "audit": [], "metrics": [], "pending": [], "errors": []}
    versions = {"python": platform.python_version(), "platform": platform.platform()}
    for name in ["biopython", "numpy", "matplotlib"]:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = "unavailable"
    source_hashes = {str(path.relative_to(root)): sha256(path) for path in (root / "agent2").glob("*.py")}
    write_json(output / "provenance.json", {"versions": versions, "config_sha256": sha256(config_path),
               "code_sha256": source_hashes, "tmscore_requested": args.tmscore or config.get("tmscore_executable")})
    for entry in targets:
        try:
            audit = audit_target(root, entry)
            summary["audit"].append(audit)
            if not entry.get("af3_server"):
                summary["pending"].append(f"{entry['target']}: obtain AlphaFold web-server output using the identical FASTA; record job provenance and evaluate it.")
            if args.command == "evaluate":
                rows = evaluate_target(root, entry, output, args.tmscore or config.get("tmscore_executable"))
                summary["metrics"].extend(rows)
                if any(row["status"] == "PARTIAL_METRICS" for row in rows):
                    summary["pending"].append(f"{entry['target']}: at least one required metric is unavailable. Inspect evaluation.json and raw logs; supply a compatible TMscore for TM/GDT or resolve missing lDDT contacts.")
        except (ValueError, OSError, KeyError, ImportError) as exc:
            error = f"{entry['target']}: {exc}"
            summary["errors"].append(error)
            if not any(item["target"] == entry["target"] for item in summary["audit"]):
                summary["audit"].append({"target": entry["target"], "status": "ERROR", "error": str(exc)})
    if args.command == "audit":
        summary["pending"].append("Run the experimental-reference evaluation to produce measured scores, mapping tables, and superposition images.")
    summary["status"] = "ERROR" if summary["errors"] else "PARTIAL_ASSIGNMENT" if summary["pending"] else "COMPLETE"
    write_json(output / "summary.json", summary)
    if summary["metrics"]:
        fields = list(dict.fromkeys(key for row in summary["metrics"] for key in row))
        write_csv(output / "metrics.csv", summary["metrics"], fields)
    (output / "REPORT.md").write_text(report(summary), encoding="utf-8")
    return summary, output
