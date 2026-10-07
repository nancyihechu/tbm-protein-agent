"""Agent 1 demo: live template-based modeling run for T1147, with results viewer.

Start from the app/ folder (see app/README.md):  streamlit run app.py
"""

import html
import re
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime

import altair as alt
import pandas as pd
import streamlit as st

APP = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(APP)
RUNS = os.path.join(REPO, "runs")
PIPELINE = os.path.join(APP, "pipeline.py")
SAVED_RUN_FILE = os.path.join(APP, "saved_run.txt")
JS_3DMOL = os.path.join(APP, "static", "3Dmol-min.js")
PY = sys.executable

ACCENT = "#2a78d6"
INK = "#14181c"
MUTED = "#5b6570"
LINE = "#e3e7eb"
AF3_GREY = "#b9c0c8"
PLDDT_BANDS = [(90, "#0053d6", "Very high (90+)"), (70, "#65cbf3", "Confident (70-90)"),
               (50, "#ffdb13", "Low (50-70)"), (0, "#ff7d45", "Very low (<50)")]

STEP_TITLES = ["Search for templates", "Choose a template", "Align the sequence to the template",
               "Build 3D models with MODELLER", "Compare with the AlphaFold3 model"]

AGENT_PROMPT = """You are Agent 1, a template-based protein structure modeling agent, giving a live demo
for CASP15 target T1147 (MmpS5, 103 residues).

Run these five commands with the Bash tool, one at a time and in this order. Wait for each one
to finish before starting the next. Use a Bash timeout of 900000 ms.

1. {py} {pipeline} step search --run-dir {run}
2. {py} {pipeline} step select --run-dir {run}
3. {py} {pipeline} step align --run-dir {run}
4. {py} {pipeline} step model --run-dir {run}
5. {py} {pipeline} step compare --run-dir {run}

After each command, write one short plain-language sentence for a non-expert audience about what
that step found, using only numbers from its output. If a command fails, stop and explain the
failure in one sentence. Do not run any other commands and do not edit files. After step 5,
finish with a two-sentence summary of the result."""

st.set_page_config(page_title="Agent 1 demo: T1147", layout="wide", initial_sidebar_state="collapsed")

st.markdown(f"""
<style>
html, body, [class*="css"] {{ font-size: 19px; }}
#MainMenu, footer, [data-testid="stToolbar"], [data-testid="stDecoration"] {{ visibility: hidden; height: 0; }}
.block-container {{ padding-top: 2.2rem; max-width: 1180px; }}
h1 {{ font-size: 2.3rem !important; font-weight: 700 !important; color: {INK}; margin-bottom: 0 !important; }}
h2, h3 {{ color: {INK}; }}
.sub {{ color: {MUTED}; font-size: 1.1rem; margin: .2rem 0 1.2rem; }}
.badge {{ display: inline-block; padding: .2rem .7rem; border-radius: 999px; font-size: .9rem;
          border: 1px solid {LINE}; color: {MUTED}; margin-right: .5rem; }}
.badge.on {{ border-color: {ACCENT}; color: {ACCENT}; }}
.tl {{ display: grid; gap: .2rem; margin: .4rem 0 1rem; }}
.st-row {{ display: grid; grid-template-columns: 2.4rem 1fr auto; gap: .8rem; align-items: start;
           padding: .75rem 0; border-bottom: 1px solid {LINE}; }}
.st-dot {{ width: 1.9rem; height: 1.9rem; border-radius: 50%; display: grid; place-items: center;
           font-weight: 700; font-size: .95rem; border: 2px solid {LINE}; color: {MUTED}; background: #fff; }}
.st-dot.running {{ border-color: {ACCENT}; color: {ACCENT}; animation: pulse 1.2s ease-in-out infinite; }}
.st-dot.done {{ background: {ACCENT}; border-color: {ACCENT}; color: #fff; }}
.st-dot.failed {{ background: #d03b3b; border-color: #d03b3b; color: #fff; }}
@keyframes pulse {{ 50% {{ box-shadow: 0 0 0 6px rgba(42,120,214,.15); }} }}
@media (prefers-reduced-motion: reduce) {{ .st-dot.running {{ animation: none; }} }}
.st-title {{ font-weight: 600; font-size: 1.08rem; color: {INK}; }}
.st-sum {{ color: {MUTED}; font-size: .98rem; margin-top: .15rem; line-height: 1.45; }}
.st-err {{ color: #b42323; font-size: .95rem; margin-top: .15rem; }}
.st-time {{ font-variant-numeric: tabular-nums; color: {MUTED}; font-size: .95rem; white-space: nowrap; }}
.total {{ font-size: 1.15rem; margin: .3rem 0 .8rem; color: {INK}; }}
.total b {{ font-variant-numeric: tabular-nums; }}
.logbox {{ border: 1px solid {LINE}; border-radius: 8px; background: #fafbfc; padding: .8rem 1rem;
           height: 430px; overflow-y: auto; font-family: ui-monospace, Menlo, Consolas, monospace;
           font-size: .82rem; line-height: 1.5; color: {INK}; white-space: pre-wrap; word-break: break-word; }}
.logbox .say {{ font-family: inherit; color: {INK}; font-size: .9rem; display: block; margin: .35rem 0; }}
.logbox .cmd {{ color: {ACCENT}; display: block; }}
.logbox .out {{ color: {MUTED}; display: block; }}
.panel-label {{ font-weight: 600; font-size: 1.05rem; margin-bottom: .4rem; }}
.metrics {{ width: 100%; border-collapse: collapse; font-size: 1rem; }}
.metrics td {{ padding: .5rem .6rem; border-bottom: 1px solid {LINE}; }}
.metrics td:last-child {{ text-align: right; font-variant-numeric: tabular-nums; font-weight: 600; }}
.legend {{ display: flex; flex-wrap: wrap; gap: .4rem 1.2rem; color: {MUTED}; font-size: .92rem; margin-top: .4rem; }}
.legend i {{ display: inline-block; width: .9rem; height: .9rem; border-radius: 3px; margin-right: .35rem;
             vertical-align: -2px; }}
.note {{ color: {MUTED}; font-size: .95rem; }}
</style>
""", unsafe_allow_html=True)


# ---------------------------------------------------------------- helpers

def clean_env():
    env = dict(os.environ)
    for key in ("CLAUDECODE", "CLAUDE_CODE_ENTRYPOINT", "CLAUDE_CODE_SSE_PORT"):
        env.pop(key, None)
    env["BASH_DEFAULT_TIMEOUT_MS"] = "900000"
    env["BASH_MAX_TIMEOUT_MS"] = "900000"
    return env


@st.cache_resource(show_spinner=False)
def claude_check():
    """Test whether 'claude -p' works on this node. Returns (ok, detail)."""
    exe = shutil.which("claude")
    if not exe:
        return False, "claude command not found"
    t0 = time.time()
    try:
        out = subprocess.run([exe, "-p", "Reply with exactly: OK", "--max-turns", "1"], capture_output=True,
                             text=True, timeout=120, env=clean_env(), cwd=REPO)
    except subprocess.TimeoutExpired:
        return False, "claude -p timed out"
    if out.returncode == 0 and "OK" in out.stdout:
        return True, f"claude -p answered in {time.time() - t0:.0f} s"
    return False, (out.stderr or out.stdout).strip()[:200] or f"exit {out.returncode}"


def read_json(path, default=None):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return default


def pid_alive(pid):
    try:
        os.kill(pid, 0)
        return True
    except (OSError, TypeError):
        return False


def run_active(run):
    proc = read_json(os.path.join(run, "process.json"), {})
    status = read_json(os.path.join(run, "status.json"), {})
    return pid_alive(proc.get("pid")) and status.get("state") not in ("done", "failed")


def latest_active_run():
    if not os.path.isdir(RUNS):
        return None
    for name in sorted(os.listdir(RUNS), reverse=True):
        run = os.path.join(RUNS, name)
        if os.path.exists(os.path.join(run, "process.json")) and run_active(run):
            return run
    return None


def saved_run():
    try:
        path = open(SAVED_RUN_FILE).read().strip()
    except OSError:
        return None
    path = path if os.path.isabs(path) else os.path.join(REPO, path)
    return path if read_json(os.path.join(path, "status.json"), {}).get("state") == "done" else None


def fmt_secs(s):
    if s is None:
        return ""
    s = int(round(s))
    return f"{s // 60} min {s % 60:02d} s" if s >= 60 else f"{s} s"


def start_run(agent_ok):
    os.makedirs(RUNS, exist_ok=True)
    run = os.path.join(RUNS, datetime.now().strftime("%Y%m%d-%H%M%S"))
    mode = "agent" if agent_ok else "direct"
    subprocess.run([PY, PIPELINE, "init", "--run-dir", run, "--mode", mode], check=True, cwd=REPO)
    if agent_ok:
        prompt = AGENT_PROMPT.format(py=PY, pipeline=PIPELINE, run=run)
        cmd = [shutil.which("claude"), "-p", prompt, "--output-format", "stream-json", "--verbose",
               "--max-turns", "30", "--allowedTools", f"Bash({PY} {PIPELINE} step:*)"]
        out = open(os.path.join(run, "agent_stream.jsonl"), "w")
        err = open(os.path.join(run, "agent_stderr.log"), "w")
    else:
        cmd = [PY, PIPELINE, "all", "--run-dir", run]
        out = open(os.path.join(run, "direct_stdout.log"), "w")
        err = subprocess.STDOUT
    proc = subprocess.Popen(cmd, stdout=out, stderr=err, cwd=REPO, env=clean_env(), start_new_session=True)
    with open(os.path.join(run, "process.json"), "w") as handle:
        json.dump({"pid": proc.pid, "mode": mode, "started": time.time()}, handle)
    return run


# ---------------------------------------------------------------- timeline and log panels

def timeline_html(status, reveal=None, now=None):
    steps = status["steps"]
    rows = []
    for i, step in enumerate(steps):
        state, summary, secs, err = step["state"], step["summary"], step["seconds"], step["error"]
        if reveal is not None:  # replay: show steps up to 'reveal' (float)
            if i < int(reveal):
                pass
            elif i == int(reveal) and reveal < len(steps):
                state, summary, secs, err = "running", "", None, ""
            else:
                state, summary, secs, err = "pending", "", None, ""
        elif state == "running" and step["start"]:
            secs = (now or time.time()) - step["start"]
        mark = {"done": "&#10003;", "failed": "!", "running": str(i + 1)}.get(state, str(i + 1))
        body = f'<div class="st-title">{html.escape(step["title"])}</div>'
        if summary:
            body += f'<div class="st-sum">{html.escape(summary)}</div>'
        elif state == "running":
            body += '<div class="st-sum">Working...</div>'
        if err:
            body += f'<div class="st-err">{html.escape(err)}</div>'
        rows.append(f'<div class="st-row"><div class="st-dot {state}">{mark}</div><div>{body}</div>'
                    f'<div class="st-time">{fmt_secs(secs)}</div></div>')
    return '<div class="tl">' + "".join(rows) + "</div>"


def agent_entries(run):
    """Turn Claude Code stream-json output into display entries (kind, text)."""
    entries = []
    path = os.path.join(run, "agent_stream.jsonl")
    if not os.path.exists(path):
        return entries
    short = lambda s: s.replace(f"{PY} ", "python ").replace(PIPELINE, "app/pipeline.py").replace(run, "<run>")
    for line in open(path):
        try:
            ev = json.loads(line)
        except ValueError:
            continue
        kind = ev.get("type")
        if kind == "system" and ev.get("subtype") == "init":
            entries.append(("out", f"Claude Code session started (model {ev.get('model', '?')})"))
        elif kind == "assistant":
            for block in ev.get("message", {}).get("content", []):
                if block.get("type") == "text" and block.get("text", "").strip():
                    entries.append(("say", block["text"].strip()))
                elif block.get("type") == "tool_use":
                    entries.append(("cmd", "$ " + short(block.get("input", {}).get("command", block.get("name", "")))))
        elif kind == "user":
            for block in ev.get("message", {}).get("content", []):
                if isinstance(block, dict) and block.get("type") == "tool_result":
                    content = block.get("content")
                    if isinstance(content, list):
                        content = " ".join(c.get("text", "") for c in content if isinstance(c, dict))
                    text = short(str(content or "")).strip()
                    entries.append(("out", text[:600] + (" ..." if len(text) > 600 else "")))
        elif kind == "result":
            secs = (ev.get("duration_ms") or 0) / 1000
            entries.append(("out", f"Agent finished ({ev.get('subtype', '')}) after {fmt_secs(secs)}, "
                                   f"{ev.get('num_turns', '?')} turns."))
    return entries


def log_html(run, mode, fraction=None):
    if mode == "agent":
        entries = agent_entries(run)
        if fraction is not None:
            entries = entries[:max(1, int(len(entries) * fraction))]
        bold = lambda t: re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", html.escape(t))
        parts = [f'<span class="{k}">{bold(t) if k == "say" else html.escape(t)}</span>' for k, t in entries]
        if not parts:
            parts = ['<span class="out">Starting Claude Code...</span>']
    else:
        try:
            lines = open(os.path.join(run, "pipeline.log")).read().splitlines()
        except OSError:
            lines = []
        if fraction is not None:
            lines = lines[:max(1, int(len(lines) * fraction))]
        parts = [f'<span class="out">{html.escape(l)}</span>' for l in lines[-200:]] or \
                ['<span class="out">Starting...</span>']
    return '<div class="logbox">' + "\n".join(parts) + "</div>"


def render_run_view(run, reveal=None, fraction=None):
    status = read_json(os.path.join(run, "status.json"))
    if not status:
        st.info("Preparing the run...")
        return status
    mode = status.get("mode", "direct")
    if reveal is not None:
        total = sum(s["seconds"] or 0 for s in status["steps"][:int(reveal)])
        label = f"Recorded time: <b>{fmt_secs(total)}</b>"
    else:
        end = status.get("finished") or time.time()
        label = f"Total time: <b>{fmt_secs(end - status['created'])}</b>"
        if status.get("state") == "failed":
            label += " (stopped on an error)"
    left, right = st.columns([1.15, 1], gap="large")
    with left:
        st.markdown(f'<div class="total">{label}</div>', unsafe_allow_html=True)
        st.markdown(timeline_html(status, reveal=reveal), unsafe_allow_html=True)
    with right:
        st.markdown(f'<div class="panel-label">{"Agent output" if mode == "agent" else "Pipeline log"}</div>',
                    unsafe_allow_html=True)
        st.markdown(log_html(run, mode, fraction), unsafe_allow_html=True)
    return status


@st.fragment(run_every=1.0)
def live_panel(run):
    status = render_run_view(run)
    if status and not run_active(run):
        if status.get("state") not in ("done", "failed"):
            # The process ended without finishing every step
            status["state"] = "failed"
            status["finished"] = status.get("finished") or time.time()
            with open(os.path.join(run, "status.json"), "w") as handle:
                json.dump(status, handle, indent=2)
        st.rerun(scope="app")


@st.fragment(run_every=0.6)
def replay_panel(run, t0):
    n = len(STEP_TITLES)
    reveal = min((time.time() - t0) / 1.6, n)
    render_run_view(run, reveal=reveal, fraction=reveal / n)
    if reveal >= n:
        st.session_state.replay_done = True
        st.rerun(scope="app")


# ---------------------------------------------------------------- results

def color_dope(values):
    lo, hi = min(values.values()), max(values.values())
    def mix(a, b, t):
        a, b = [int(a[i:i + 2], 16) for i in (1, 3, 5)], [int(b[i:i + 2], 16) for i in (1, 3, 5)]
        return "#" + "".join(f"{round(x + (y - x) * t):02x}" for x, y in zip(a, b))
    out = {}
    for r, v in values.items():
        t = (v - lo) / (hi - lo) if hi > lo else 0.5
        out[r] = mix(ACCENT, "#f1f1f1", t * 2) if t < 0.5 else mix("#f1f1f1", "#d03b3b", (t - 0.5) * 2)
    return out


def color_plddt(values):
    return {r: next(c for cut, c, _ in PLDDT_BANDS if v >= cut) for r, v in values.items()}


def ca_bfactors(pdb_text):
    vals = {}
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            vals[int(line[22:26])] = float(line[60:66])
    return vals


@st.cache_data(show_spinner=False)
def js_3dmol():
    return open(JS_3DMOL).read().replace("</script", "<\\/script")


def viewer(models, height=520):
    """models: list of (pdb_text, {resi: hex} or None, default_hex)."""
    specs = json.dumps([{"pdb": p, "colors": {str(k): v for k, v in (c or {}).items()}, "base": b}
                        for p, c, b in models])
    page = f"""
<div id="v" style="width:100%;height:{height}px;position:relative;border:1px solid {LINE};border-radius:8px;"></div>
<script>{js_3dmol()}</script>
<script>
const specs = {specs};
const v = $3Dmol.createViewer(document.getElementById("v"), {{backgroundColor: "white"}});
for (const s of specs) {{
  const m = v.addModel(s.pdb, "pdb");
  m.setStyle({{}}, {{cartoon: {{colorfunc: a => s.colors[a.resi] || s.base}}}});
}}
v.zoomTo(); v.render(); v.zoom(1.1, 600);
</script>"""
    st.iframe(page, height=height + 8)


def results_view(run, origin):
    res = read_json(os.path.join(run, "results.json"))
    models = read_json(os.path.join(run, "models.json"), [])
    if not res:
        st.markdown('<p class="note">No results yet. Start a live run or replay the saved run.</p>',
                    unsafe_allow_html=True)
        return
    vdir = os.path.join(run, res["target"], "viewer")
    st.markdown(f'<p class="note">Results from run {os.path.basename(run)} ({origin}). '
                f'The AlphaFold3 model was computed ahead of time and loaded from '
                f'<code>{html.escape(res["af3_model"])}</code>.</p>', unsafe_allow_html=True)

    left, right = st.columns([1, 1.5], gap="large")
    tm, rm = res["tm_align"], res["residue_matched"]
    rows = [
        ("Template", f'{res["template"]} ({res["identity"]}% identical)'),
        ("Template covers", f'residues {res["query_range"][0]}-{res["query_range"][1]} ({res["coverage"]:.0f}%)'),
        ("Residues modeled", f'{res["modeled_range"][0]}-{res["modeled_range"][1]}'),
        ("Best TBM model (DOPE)", f'{res["best_dope"]:.0f}'),
        ("Best TBM model (GA341)", f'{res["best_ga341"]:.3f}'),
        ("AF3 mean pLDDT, whole protein", f'{res["af3_mean_plddt"]:.1f}'),
        ("AF3 pTM", f'{res["af3_ptm"]:.2f}'),
        ("TM-score vs AF3 (TM-align)", f'{tm["tm_score_norm_tbm"]:.3f}'),
        ("TM-score vs AF3 (residue-matched)", f'{rm["tm_score"]:.3f}'),
        ("RMSD (TM-align)", f'{tm["rmsd"]:.2f} &#8491; over {tm["aligned_length"]}'),
        ("RMSD (residue-matched)", f'{rm["rmsd"]:.2f} &#8491;'),
    ]
    with left:
        st.markdown("### Metrics")
        st.markdown('<table class="metrics">' + "".join(f"<tr><td>{a}</td><td>{b}</td></tr>" for a, b in rows)
                    + "</table>", unsafe_allow_html=True)
        if models:
            st.markdown("<br>", unsafe_allow_html=True)
            st.markdown('<p class="note">All five models (DOPE, lower is better): ' +
                        ", ".join(f'{m["file"].split(".")[1][-1]}: {m["dope"]:.0f}' for m in models) + "</p>",
                        unsafe_allow_html=True)

    with right:
        st.markdown("### 3D structure")
        choice = st.radio("View", ["TBM model, colored by DOPE", "AF3 model computed ahead of time, colored by pLDDT",
                                   "Superposition"], horizontal=False, label_visibility="collapsed")
        tbm = open(os.path.join(vdir, "tbm_dope.pdb")).read()
        af3 = open(os.path.join(vdir, "af3.pdb")).read()
        if choice.startswith("TBM"):
            viewer([(tbm, color_dope(ca_bfactors(tbm)), "#cccccc")])
            legend = (f'<span><i style="background:{ACCENT}"></i>Favorable (low DOPE)</span>'
                      '<span><i style="background:#f1f1f1;border:1px solid #ddd"></i>Average</span>'
                      '<span><i style="background:#d03b3b"></i>Unfavorable (high DOPE)</span>')
        elif choice.startswith("AF3"):
            viewer([(af3, color_plddt(ca_bfactors(af3)), "#cccccc")])
            legend = "".join(f'<span><i style="background:{c}"></i>{l}</span>' for _, c, l in PLDDT_BANDS)
        else:
            sup = open(os.path.join(vdir, "tbm_superposed.pdb")).read()
            viewer([(af3, None, AF3_GREY), (sup, None, ACCENT)])
            legend = (f'<span><i style="background:{ACCENT}"></i>TBM model</span>'
                      f'<span><i style="background:{AF3_GREY}"></i>AF3 model computed ahead of time</span>')
        st.markdown(f'<div class="legend">{legend}</div>', unsafe_allow_html=True)

    st.markdown("### Distance between the two models at each residue")
    df = pd.read_csv(os.path.join(run, "per_residue.csv"))
    q0 = res["query_range"][0]
    base = alt.Chart(df).encode(x=alt.X("resi:Q", title="Residue", scale=alt.Scale(domain=[1, int(df.resi.max())])))
    layers = []
    if q0 > res["modeled_range"][0]:
        layers.append(alt.Chart(pd.DataFrame({"a": [res["modeled_range"][0]], "b": [q0 - 0.5],
                                              "label": ["No template"]}))
                      .mark_rect(color="#eef1f4").encode(x="a:Q", x2="b:Q"))
        layers.append(alt.Chart(pd.DataFrame({"x": [(res["modeled_range"][0] + q0) / 2], "y": [float(df.ca_deviation_A.max())],
                                              "t": ["no template"]}))
                      .mark_text(color=MUTED, fontSize=14, baseline="top").encode(x="x:Q", y="y:Q", text="t:N"))
    layers.append(alt.Chart(pd.DataFrame({"y": [5]})).mark_rule(color=MUTED, strokeDash=[5, 4]).encode(y="y:Q"))
    layers.append(base.mark_line(color=ACCENT, strokeWidth=2.5).encode(
        y=alt.Y("ca_deviation_A:Q", title="CA distance (Å)"),
        tooltip=[alt.Tooltip("resi:Q", title="Residue"), alt.Tooltip("ca_deviation_A:Q", title="Distance (Å)", format=".1f"),
                 alt.Tooltip("af3_plddt:Q", title="AF3 pLDDT", format=".0f")]))
    chart = alt.layer(*layers).properties(height=300).configure_axis(
        labelFontSize=15, titleFontSize=16, gridColor=LINE, labelColor=MUTED, titleColor=INK).configure_view(strokeWidth=0)
    st.altair_chart(chart, width="stretch")
    st.markdown('<p class="note">Distance between matching residues after superposing the TBM model on the '
                'AF3 model. Dashed line: 5 &#8491;.</p>', unsafe_allow_html=True)


# ---------------------------------------------------------------- page

agent_ok, agent_detail = claude_check()
if "run" not in st.session_state:
    st.session_state.run = latest_active_run()
st.session_state.setdefault("replay", None)
st.session_state.setdefault("replay_done", False)

st.markdown("<h1>Agent 1: Template-based modeling</h1>", unsafe_allow_html=True)
st.markdown('<div class="sub">Target T1147 &middot; MmpS5, <i>Mycobacterium thermoresistibile</i> &middot; '
            '103 residues</div>', unsafe_allow_html=True)
st.markdown(
    f'<span class="badge {"on" if agent_ok else ""}">'
    f'{"Runs through Claude Code (claude -p)" if agent_ok else "Runs the scripts directly"}</span>'
    f'<span class="badge">{os.uname().nodename.split(".")[0]} &middot; '
    f'{os.environ.get("SLURM_CPUS_PER_TASK", "?")} CPUs</span>', unsafe_allow_html=True)
if not agent_ok:
    st.caption(f"claude -p is not available here ({agent_detail}).")

tab_run, tab_results = st.tabs(["Run", "Results"])

with tab_run:
    active = st.session_state.run is not None and run_active(st.session_state.run)
    replaying = st.session_state.replay is not None and not st.session_state.replay_done
    c1, c2, _ = st.columns([1.3, 1, 3])
    if c1.button("Start live run", type="primary", disabled=active or replaying, width="stretch"):
        st.session_state.run = start_run(agent_ok)
        st.session_state.replay = None
        st.rerun()
    saved = saved_run()
    if c2.button("Replay saved run", disabled=active or replaying or not saved, width="stretch"):
        st.session_state.replay = {"run": saved, "t0": time.time()}
        st.session_state.replay_done = False
        st.rerun()

    if st.session_state.replay and not st.session_state.replay_done:
        st.markdown(f'<p class="note">Replay of saved run {os.path.basename(st.session_state.replay["run"])} '
                    "(not live; recorded times shown).</p>", unsafe_allow_html=True)
        replay_panel(st.session_state.replay["run"], st.session_state.replay["t0"])
    elif st.session_state.replay and st.session_state.replay_done:
        st.markdown(f'<p class="note">Replay of saved run {os.path.basename(st.session_state.replay["run"])} '
                    "(not live; recorded times shown).</p>", unsafe_allow_html=True)
        render_run_view(st.session_state.replay["run"], reveal=len(STEP_TITLES), fraction=1.0)
    elif st.session_state.run and active:
        live_panel(st.session_state.run)
    elif st.session_state.run:
        render_run_view(st.session_state.run)
    else:
        st.markdown('<p class="note">Press <b>Start live run</b> to search for templates, align, and build '
                    'five MODELLER models for T1147 on this node. AlphaFold3 is not run live; its model was '
                    'computed ahead of time.</p>', unsafe_allow_html=True)

with tab_results:
    if st.session_state.replay and st.session_state.replay_done:
        results_view(st.session_state.replay["run"], "saved run")
    elif st.session_state.run and read_json(os.path.join(st.session_state.run, "status.json"), {}).get("state") == "done":
        results_view(st.session_state.run, "live run")
    elif saved_run():
        results_view(saved_run(), "saved run")
    else:
        results_view(RUNS, "none")
