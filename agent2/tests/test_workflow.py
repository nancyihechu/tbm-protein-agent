"""Integration checks for preservation and evaluation completeness."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from agent2.workflow import evaluate_target, load_target, reserve_output, run, sha256


REPO = Path(__file__).resolve().parents[1]


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "agent2").mkdir()
        self.config = json.loads((REPO / "config.json").read_text())

    def test_output_cannot_overwrite_or_escape_agent2_results(self):
        path = reserve_output(self.root, "agent2/results/first")
        (path / "keep.txt").write_text("preserve")
        with self.assertRaises(FileExistsError):
            reserve_output(self.root, "agent2/results/first")
        with self.assertRaises(ValueError):
            reserve_output(self.root, "targets/T1147")
        with self.assertRaises(ValueError):
            reserve_output(self.root, "../escape")
        self.assertEqual((path / "keep.txt").read_text(), "preserve")

    def test_changed_reference_blocks_scoring(self):
        config = deepcopy(self.config)
        config["targets"][0]["reference"]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "checksum"):
            load_target(REPO.parent, config["targets"][0])

    def test_local_af3_does_not_satisfy_web_server_requirement(self):
        (self.root / "agent2/config.json").write_text(json.dumps(self.config))
        args = SimpleNamespace(config="agent2/config.json", out="agent2/results/test", command="audit", tmscore=None)
        with patch("agent2.workflow.audit_target", side_effect=lambda _, e: {
            "target": e["target"], "status": "INPUTS_PRESENT", "reference": e["reference"],
            "domain_range": e["domain_range"], "classification": e["classification"], "af3_server_status": "NOT_PROVIDED"}):
            result, output = run(self.root, args)
        self.assertEqual(result["status"], "PARTIAL_ASSIGNMENT")
        self.assertEqual(sum("web-server" in item for item in result["pending"]), 3)
        self.assertEqual(result["metrics"], [])

    def test_reference_and_model_inputs_are_only_read(self):
        config = self.config["targets"][0]
        root = REPO.parent
        tracked = [root / config["metadata"], root / config["reference"]["path"]]
        before = [sha256(p) for p in tracked]
        metadata, sequence, _, models, _ = load_target(root, config)
        self.assertEqual(len(sequence), 103)
        self.assertEqual([item["method"] for item in models], ["tbm", "alphafold3_local"])
        self.assertEqual(metadata["target"], "T1147")
        self.assertEqual([sha256(p) for p in tracked], before)

    def test_three_distinct_targets_required(self):
        self.config["targets"][1] = self.config["targets"][0]
        (self.root / "agent2/config.json").write_text(json.dumps(self.config))
        args = SimpleNamespace(config="agent2/config.json", out=None, command="audit", tmscore=None)
        with self.assertRaisesRegex(ValueError, "distinct"):
            run(self.root, args)

    def test_missing_gdt_never_marks_evaluation_complete(self):
        valid_tm_missing_gdt = {"status": "AVAILABLE", "tm_score": 0.8, "gdt_ts": None}
        with patch("agent2.metrics.run_tmscore", return_value=valid_tm_missing_gdt), patch("agent2.workflow.render_pair"):
            rows = evaluate_target(REPO.parent, self.config["targets"][0], self.root, None)
        measured = [row for row in rows if row["method"] != "alphafold3_server"]
        self.assertEqual(len(measured), 2)
        self.assertTrue(all(row["status"] == "PARTIAL_METRICS" for row in measured))
        self.assertTrue(all(row["gdt_ts"] is None for row in measured))


if __name__ == "__main__":
    unittest.main()
