# Agent 1 demo app (T1147)

A Streamlit app that runs Agent 1 live for CASP15 target T1147 and shows the
results. It runs inside an interactive SLURM job on Hellbender and uses the
CPUs of that job (`$SLURM_CPUS_PER_TASK`). It does not submit SLURM jobs.

## What it does

**Live run** (the "Start live run" button), five steps:

1. Search for templates: hhblits against UniRef30, then hhsearch against pdb70.
2. Choose a template with the cutoffs in `SCRIPTS/screen_templates.py`
   (probability >= 95, identity >= 30%, coverage >= 70%).
3. Clean the template chain and build the PIR alignment.
4. Build 5 MODELLER models (DOPE and GA341).
5. Compare the best model with the AlphaFold3 model computed ahead of time
   (`targets/T1147/af3_model/`), using `SCRIPTS/compare_to_af3.py`.
   AlphaFold3 is not run live.

At startup the app tests `claude -p`. If it works, Claude Code runs the five
steps in non-interactive mode and its real output is shown as **Agent output**.
Claude is only allowed to run the five `app/pipeline.py step ...` commands. If
`claude -p` does not work, the app runs the steps directly and shows the
**Pipeline log**.

Every run is written to `runs/<timestamp>/` (status, logs, models, comparison).
Saved results in `targets/` are only read, never written.

**Replay saved run** replays a finished run (the path in `app/saved_run.txt`)
as a backup if a live run is not possible.

**Results tab**: metrics table, 3D viewer (TBM model colored by DOPE, AF3 model
colored by pLDDT, superposition), and the per-residue distance chart. 3Dmol.js
is shipped in `app/static/` (BSD-3-Clause, see `LICENSE-3Dmol.txt`), so the
page needs no internet.

## Start the app

1. Start an interactive job on a compute node, for example:

   ```bash
   srun --partition=general --cpus-per-task=8 --mem=32G --time=01:00:00 --pty bash
   ```

   The test run below used 8 CPUs. Ask for at least 30 minutes so the job
   does not end during the demo.

2. On the compute node, start the app and print the node name:

   ```bash
   cd /cluster/pixstor/dkhf3-lab/users/nancy/TBM_AGENT/app
   module load miniconda3/v26.7.1-1_py314
   eval "$(conda shell.bash hook)"
   conda activate /cluster/pixstor/dkhf3-lab/users/nancy/envs/tbm_tools
   echo "Node: $(hostname -s)   Port: 8501"
   streamlit run app.py
   ```

   The app listens on `localhost:8501` on the compute node only.

## Open it from a Windows laptop

In PowerShell on the laptop, run this (replace `c023` with the node name
printed above):

```powershell
ssh -L 8501:localhost:28501 nigby@hellbender-login.rnet.missouri.edu -t "ssh -L 28501:localhost:8501 c023"
```

Then open <http://localhost:8501> in a browser on the laptop. Keep the
PowerShell window open during the demo.

How it works: the first hop logs you in to the login node and forwards laptop
port 8501 to port 28501 on the login node. The second hop runs on the login
node, uses your cluster SSH key to reach the compute node, and forwards port
28501 to the app on port 8501. Port 28501 on the shared login node avoids
clashing with anyone else using 8501 there; if it is taken, pick another
number and change it in both places.

If you normally connect to Hellbender with a different address, use that
address instead of `hellbender-login.rnet.missouri.edu`.

## Requirements

- The `tbm_tools` conda env (see `environment.yml`) plus Streamlit, installed
  with `pip install streamlit` (version 1.65.0 was used).
- HH-suite from the RoseTTAFold conda env and the databases listed in the main
  README (paths are set in `app/pipeline.py`).
- Optional: Claude Code (`claude`) logged in on the cluster, for agent mode.

## Files

| File | Purpose |
|---|---|
| `app.py` | Streamlit app |
| `pipeline.py` | The five steps; each writes `status.json`, `pipeline.log` and `logs/<step>.log` in the run folder |
| `.streamlit/config.toml` | Light theme, one accent color, port 8501, localhost only |
| `static/3Dmol-min.js` | 3Dmol.js 2.4.2, bundled for offline use |
| `saved_run.txt` | Run folder used by "Replay saved run" |

The pipeline can also be run without the app:

```bash
python app/pipeline.py all --run-dir runs/test     # from the repo root
```
