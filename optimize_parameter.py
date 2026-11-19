#!/usr/bin/env python3
"""
Optuna parameter refinement for one clustering algorithm on one sample.

Each Optuna trial writes the sampled parameters into a temporary copy of the
config, re-runs only the clustering + scoring rules of the Snakemake workflow
(reusing every upstream output), reads the chosen metric back and returns it to
Optuna.

Modes
-----
  --mode group   (default) refine the algorithm across many normalisation x
                 feature-selection workflows at once; the objective is the MEAN
                 of the chosen metric over those workflows. By default all
                 existing workflows are used; restrict with --norm / --feat.
  --mode single  refine exactly one workflow (--norm X --feat Y); the objective
                 is that workflow's metric value.

In both modes every trial records the mean AND the value of each individual
workflow.

Metrics
-------
  --metrics NMI | ARI | HOM | COM | CHAOS | PAS | ASW
  Higher is better for NMI, ARI, HOM, COM and ASW; lower is better for CHAOS
  and PAS. The optimisation direction and the --target test follow that.

  --source cluster    (default) read 09_<alg>_cluster_metrics.csv
                      (raw cluster IDs vs ground truth; all 7 metrics)
  --source annotated  read 09_<alg>_metrics.csv (after SingleR annotation;
                      NMI / ARI / HOM / COM only)

Early stop
----------
  --target V  stop the study as soon as any completed trial reaches V
              (>= V for higher-is-better metrics, <= V for CHAOS / PAS).
              Trials already running in parallel are allowed to finish.

Examples
--------
  python optimize_parameter.py --sample Mouse_brain_squidpy --algorithm SpaGCN \\
      --metrics PAS --mode group --trials 50 --cores 64 --config config.yaml --target 0.05

  python optimize_parameter.py --sample Mouse_brain_squidpy --algorithm BANKSY \\
      --metrics ARI --mode single --norm SCTransform --feat moransi --trials 40 --target 0.6

Outputs (prefix set by --output_prefix)
---------------------------------------
  <prefix>_trials.csv     one row per trial: parameters, state, mean metric,
                          number of workflows scored, one column per workflow
  <prefix>_workflows.csv  one row per trial x workflow with every metric found
  Both files are rewritten after every trial, so an interrupted run keeps its
  results.
"""
import argparse
import logging
import os
import shutil
import subprocess
import threading
from copy import deepcopy

import numpy as np
import optuna
import pandas as pd
import yaml

optuna.logging.set_verbosity(optuna.logging.WARNING)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
HIGHER_IS_BETTER = {"NMI": True, "ARI": True, "HOM": True, "COM": True, "ASW": True,
                    "CHAOS": False, "PAS": False}
SOURCE_METRICS = {"cluster": ["NMI", "ARI", "HOM", "COM", "CHAOS", "PAS", "ASW"],
                  "annotated": ["NMI", "ARI", "HOM", "COM"]}
SOURCE_FILE = {"cluster": "09_{alg}_cluster_metrics.csv",
               "annotated": "09_{alg}_metrics.csv"}
SOURCE_RULES = {"cluster": ["compute_cluster_metrics"],
                "annotated": ["annotate_clusters", "compute_metrics"]}

class TrialFailed(Exception):
    """Raised when a trial produced no usable score; Optuna marks it FAIL."""


# ---------------------------------------------------------------------------
# Search spaces
# ---------------------------------------------------------------------------
def suggest_params(trial, alg):
    params = {}
    if alg == "BANKSY":
        params["lambda"]     = trial.suggest_categorical("lambda",     [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
        params["k_geom"]     = trial.suggest_categorical("k_geom",     [10, 15, 25, 30, 40, 50, 60, 70, 80])
        params["resolution"] = trial.suggest_categorical("resolution", [0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0])

    elif alg == "BayesSpace":
        params["q"]       = trial.suggest_int("q", 5, 20)
        params["nrep"]    = trial.suggest_categorical("nrep",    [500, 1000, 2000, 3000, 4000, 5000])
        params["burn_in"] = trial.suggest_categorical("burn_in", [100, 300, 500, 700, 1000])
        params["gamma"]   = trial.suggest_categorical("gamma",   [1, 2, 3, 4, 5])

    elif alg == "SpatialLeiden":
        params["resolution"]  = trial.suggest_categorical("resolution",  [0.5, 1.0, 1.5, 2.0])
        params["layer_ratio"] = trial.suggest_categorical("layer_ratio", [0.5, 1.0, 1.8, 2.5, 3.0])
        params["n_neighbors"] = trial.suggest_categorical("n_neighbors", [6, 10, 15, 18, 24, 30])

    elif alg == "SpaGCN":
        params["n_clusters"] = trial.suggest_int("n_clusters", 5, 23)
        params["beta"]       = trial.suggest_int("beta", 10, 100)
        params["alpha"]      = trial.suggest_categorical("alpha",  [0.5, 1, 1.5, 2])
        params["refine"]     = trial.suggest_categorical("refine", [True, False])

    elif alg == "stlearn":
        params["n_neighbors"]     = trial.suggest_int("n_neighbors", 3, 30)
        params["n_pcs"]           = trial.suggest_categorical("n_pcs", [10, 15, 20, 25, 30, 40, 50])
        params["resolution"]      = trial.suggest_float("resolution", 0.1, 2.0, step=0.1)
        params["use_morphology"]  = trial.suggest_categorical("use_morphology", [True, False])

    return params


SUPPORTED_ALGORITHMS = ["BANKSY", "BayesSpace", "SpatialLeiden", "SpaGCN", "stlearn"]


def load_config(config_path="config.yaml"):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


# ---------------------------------------------------------------------------
# Workflow selection
# ---------------------------------------------------------------------------
def select_workflows(sample, mode, norms_cfg, feats_cfg, norm_arg, feat_arg):
    """Return the list of (norm, feat) workflows to refine.

    A workflow is usable only if its upstream directory data/<sample>/<norm>/<feat>
    exists (steps 0-3 already run).
    """
    if mode == "single":
        if not norm_arg or not feat_arg or len(norm_arg) != 1 or len(feat_arg) != 1:
            raise SystemExit("--mode single needs exactly one --norm and one --feat, "
                             "e.g. --norm SCTransform --feat moransi")
    norms = norm_arg if norm_arg else norms_cfg
    feats = feat_arg if feat_arg else feats_cfg

    workflows, missing = [], []
    for n in norms:
        for f in feats:
            if os.path.isdir(os.path.join("data", sample, n, f)):
                workflows.append((n, f))
            else:
                missing.append(f"{n}/{f}")
    if missing:
        logging.warning(f"Skipping {len(missing)} workflow(s) without upstream data: {', '.join(missing)}")
    if not workflows:
        raise SystemExit(f"No usable workflows for sample {sample} "
                         f"(checked data/{sample}/<norm>/<feat>).")
    return workflows


def wf_name(norm, feat):
    return f"{norm}/{feat}"


def read_metrics(path, metric_names):
    """Return {metric: value} from a metrics CSV (first row), or None if unusable."""
    if not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path)
    except Exception:
        return None
    if df.empty:
        return None
    col = {c.upper(): c for c in df.columns}
    row = df.iloc[0]
    out = {}
    for m in metric_names:
        if m in col:
            try:
                out[m] = float(row[col[m]])
            except (TypeError, ValueError):
                out[m] = np.nan
    if "TRUTH" in col:
        out["truth"] = row[col["TRUTH"]]
    return out


# ---------------------------------------------------------------------------
# Objective
# ---------------------------------------------------------------------------
class SpatialObjective:
    """
    Optuna objective that runs one Snakemake trial per call.

    Each trial gets its own isolated namespace (data/_opt_t{N}_{sample}/ etc.,
    backed by --scratch_dir) so parallel trials never write to the same paths
    and nothing is written into the shared analysis tree.
    """

    def __init__(self, base_config, args, workflows):
        self.base_config = base_config
        self.args = args
        self.workflows = workflows
        self.metric = args.metrics
        self.metric_file = SOURCE_FILE[args.source].format(alg=args.algorithm)
        self.trial_rows = []      # one dict per trial
        self.workflow_rows = []   # one dict per trial x workflow
        self._lock = threading.Lock()
        self.cores_per_trial = max(1, args.cores // max(1, args.n_parallel))

    # ------------------------------------------------------------------
    # Per-trial file-system isolation
    # ------------------------------------------------------------------
    def _setup_trial_namespace(self, trial_number):
        """
        Build an isolated namespace for one trial on scratch:

          {scratch}/_opt_t{N}_{sample}/data/     symlinks to data/{sample}/ files,
                                                 recursing into {norm}/{feat}/ but
                                                 not into algorithm folders
          {scratch}/_opt_t{N}_{sample}/results/  symlinks to results/{sample}/{norm}/{feat}/
                                                 files (e.g. 04_{feat}.csv gene lists,
                                                 needed by SpaGCN / stLearn / SpatialLeiden)
          {scratch}/_opt_t{N}_{sample}/logs/

        and expose it to Snakemake through working-directory symlinks
        data/, results/, logs/_opt_t{N}_{sample} and data/Ground_truth/_opt_t{N}_{sample}.
        """
        sample = self.args.sample
        trial_sample = f"_opt_t{trial_number}_{sample}"
        ws = os.path.join(self.args.scratch_dir, trial_sample)
        ws_data, ws_res, ws_logs = (os.path.join(ws, d) for d in ("data", "results", "logs"))
        for d in (ws_data, ws_res, ws_logs):
            os.makedirs(d, exist_ok=True)

        def _link_tree(src, dst, depth=0):
            # depth 0 = sample dir, 1 = {norm}/, 2 = {feat}/ (stop: no algorithm folders)
            try:
                entries = os.listdir(src)
            except OSError:
                return
            for item in entries:
                s, d = os.path.join(src, item), os.path.join(dst, item)
                if os.path.islink(d) or os.path.exists(d):
                    continue
                if os.path.isfile(s):
                    os.symlink(os.path.abspath(s), d)
                elif os.path.isdir(s) and depth < 2:
                    os.makedirs(d, exist_ok=True)
                    _link_tree(s, d, depth + 1)

        _link_tree(os.path.join("data", sample), ws_data)
        _link_tree(os.path.join("results", sample), ws_res)

        for wf_base, target in ((os.path.join("data", trial_sample), ws_data),
                                (os.path.join("results", trial_sample), ws_res),
                                (os.path.join("logs", trial_sample), ws_logs)):
            os.makedirs(os.path.dirname(wf_base), exist_ok=True)
            if not os.path.lexists(wf_base):
                os.symlink(os.path.abspath(target), wf_base)

        real_gt = os.path.join("data", "Ground_truth", sample)
        trial_gt = os.path.join("data", "Ground_truth", trial_sample)
        if os.path.isdir(real_gt) and not os.path.lexists(trial_gt):
            os.symlink(os.path.abspath(real_gt), trial_gt)
        return trial_sample

    def _teardown_trial_namespace(self, trial_sample):
        for base in ("data", "results", "logs"):
            link = os.path.join(base, trial_sample)
            if os.path.islink(link):
                os.unlink(link)
            elif os.path.exists(link):
                shutil.rmtree(link, ignore_errors=True)
        gt_link = os.path.join("data", "Ground_truth", trial_sample)
        if os.path.islink(gt_link):
            os.unlink(gt_link)
        ws = os.path.join(self.args.scratch_dir, trial_sample)
        if os.path.exists(ws):
            shutil.rmtree(ws, ignore_errors=True)

    # ------------------------------------------------------------------
    # Result bookkeeping
    # ------------------------------------------------------------------
    def _record(self, trial_number, params, state, mean_value, per_wf, all_metrics):
        row = {"trial_id": trial_number, "state": state,
               f"mean_{self.metric}": mean_value,
               "n_workflows_scored": sum(v is not None and not np.isnan(v) for v in per_wf.values()),
               "n_workflows": len(self.workflows)}
        row.update(params)
        for (n, f) in self.workflows:
            row[wf_name(n, f)] = per_wf.get(wf_name(n, f), np.nan)
        with self._lock:
            self.trial_rows.append(row)
            for (n, f) in self.workflows:
                r = {"trial_id": trial_number, "normalization": n, "feature_selection": f,
                     "algorithm": self.args.algorithm, "state": state}
                r.update(params)
                r.update(all_metrics.get(wf_name(n, f)) or {})
                self.workflow_rows.append(r)
            self.write_outputs()

    def write_outputs(self):
        """Rewrite both CSVs (called under the lock after every trial)."""
        if not self.trial_rows:
            return
        better_high = HIGHER_IS_BETTER[self.metric]
        t = pd.DataFrame(self.trial_rows).sort_values(
            f"mean_{self.metric}", ascending=not better_high, na_position="last")
        t.to_csv(f"{self.args.output_prefix}_trials.csv", index=False)
        pd.DataFrame(self.workflow_rows).sort_values(["trial_id", "normalization", "feature_selection"]) \
            .to_csv(f"{self.args.output_prefix}_workflows.csv", index=False)

    # ------------------------------------------------------------------
    # Optuna objective call
    # ------------------------------------------------------------------
    def __call__(self, trial):
        alg = self.args.algorithm
        params = suggest_params(trial, alg)
        run_config = deepcopy(self.base_config)
        section = alg   # each algorithm has its own block in config.yaml
        run_config.setdefault(section, {})
        run_config[section].update(params)

        temp_config = os.path.abspath(f"temp_opt_{self.args.sample}_{alg}_{trial.number}.yaml")
        with open(temp_config, "w") as fh:
            yaml.dump(run_config, fh)

        trial_sample = self._setup_trial_namespace(trial.number)
        per_wf, all_metrics = {}, {}
        try:
            targets = [f"results/{trial_sample}/{n}/{f}/{alg}/{self.metric_file}"
                       for (n, f) in self.workflows]
            cmd = ["snakemake", *targets,
                   "--cores", str(self.cores_per_trial),
                   "--configfile", temp_config,
                   "--rerun-incomplete", "--nolock", "--keep-going", "-p",
                   "--allowed-rules", f"cluster_{alg}", "reduce", "export_normalized", *SOURCE_RULES[self.args.source],
                   "-R", f"cluster_{alg}"]
            logging.info(f"Trial {trial.number}: {alg} | params={params} | "
                         f"cores={self.cores_per_trial} | workflows={len(targets)}")
            rc = subprocess.run(cmd).returncode
            if rc != 0:
                logging.warning(f"Trial {trial.number}: Snakemake exited with code {rc} "
                                f"(some workflows may have failed)")

            # Read every workflow's metrics BEFORE the namespace is deleted.
            for (n, f), path in zip(self.workflows, targets):
                m = read_metrics(path, SOURCE_METRICS[self.args.source])
                all_metrics[wf_name(n, f)] = m
                per_wf[wf_name(n, f)] = m.get(self.metric, np.nan) if m else np.nan
        finally:
            self._teardown_trial_namespace(trial_sample)
            if os.path.exists(temp_config):
                os.remove(temp_config)

        values = np.array([v for v in per_wf.values()], dtype=float)
        n_ok = int(np.sum(~np.isnan(values)))
        complete = n_ok == len(self.workflows) or (n_ok > 0 and self.args.allow_partial)
        mean_value = float(np.nanmean(values)) if n_ok else np.nan

        # Report the trial's numbers to the log in both modes.
        wf_txt = ", ".join(f"{k}={v:.4f}" for k, v in per_wf.items() if not np.isnan(v))
        logging.info(f"Trial {trial.number}: mean {self.metric} = {mean_value:.4f} "
                     f"({n_ok}/{len(self.workflows)} workflows scored) | {wf_txt}")

        trial.set_user_attr("per_workflow", {k: (None if np.isnan(v) else v) for k, v in per_wf.items()})
        trial.set_user_attr("n_workflows_scored", n_ok)
        self._record(trial.number, params, "COMPLETE" if complete else "FAIL", mean_value, per_wf, all_metrics)

        if not complete:
            raise TrialFailed(f"Trial {trial.number}: only {n_ok}/{len(self.workflows)} workflows "
                              f"produced {self.metric} (use --allow_partial to score partial trials)")
        return mean_value


# ---------------------------------------------------------------------------
# Early stop on target
# ---------------------------------------------------------------------------
def make_target_callback(target, higher_is_better, metric):
    def _cb(study, trial):
        if trial.state != optuna.trial.TrialState.COMPLETE or trial.value is None:
            return
        reached = trial.value >= target if higher_is_better else trial.value <= target
        if reached:
            logging.info(f"Target reached: trial {trial.number} {metric} = {trial.value:.4f} "
                         f"({'>=' if higher_is_better else '<='} {target}). Stopping the study.")
            study.set_user_attr("stopped_on_target", trial.number)
            study.stop()
    return _cb


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True)
    p.add_argument("--algorithm", required=True, choices=SUPPORTED_ALGORITHMS)
    p.add_argument("--metrics", "--metric", dest="metrics", default="NMI",
                   type=lambda s: s.upper(), choices=list(HIGHER_IS_BETTER),
                   help="Metric to refine (default NMI). CHAOS and PAS are minimised, the rest maximised.")
    p.add_argument("--source", choices=["cluster", "annotated"], default="cluster",
                   help="cluster: raw cluster IDs vs truth (09_<alg>_cluster_metrics.csv, all metrics); "
                        "annotated: after SingleR annotation (09_<alg>_metrics.csv, NMI/ARI/HOM/COM).")
    p.add_argument("--mode", choices=["group", "single"], default="group",
                   help="group: mean over many norm x feature workflows; single: one workflow.")
    p.add_argument("--norm", nargs="+", default=None,
                   help="Normalisation(s). single mode: exactly one. group mode: optional subset "
                        "(default: all in config).")
    p.add_argument("--feat", nargs="+", default=None,
                   help="Feature-selection method(s). single mode: exactly one. group mode: optional "
                        "subset (default: all in config).")
    p.add_argument("--target", type=float, default=None,
                   help="Stop as soon as a trial reaches this value (>= for NMI/ARI/HOM/COM/ASW, "
                        "<= for CHAOS/PAS).")
    p.add_argument("--allow_partial", action="store_true",
                   help="Score a trial on the workflows that succeeded even if some failed "
                        "(default: a trial with any failed workflow is marked FAIL).")
    p.add_argument("--trials", type=int, default=10)
    p.add_argument("--cores", type=int, default=16,
                   help="Total CPU cores; divided evenly among parallel trials.")
    p.add_argument("--n_parallel", type=int, default=1,
                   help="Number of Optuna trials to run at once (default 1).")
    p.add_argument("--seed", type=int, default=None, help="Seed for the TPE sampler (reproducible proposals).")
    p.add_argument("--config", default="config.yaml")
    p.add_argument("--output_prefix", default=None,
                   help="Prefix for the output CSVs (default optuna_<sample>_<alg>_<mode>_<metric>).")
    p.add_argument("--scratch_dir", default="/vast/scratch/users/xu.ya/opt_trials",
                   help="Base directory on scratch for trial workspaces.")
    args = p.parse_args()

    if args.metrics not in SOURCE_METRICS[args.source]:
        raise SystemExit(f"--metrics {args.metrics} is not available with --source {args.source} "
                         f"(available: {', '.join(SOURCE_METRICS[args.source])}).")
    if args.output_prefix is None:
        args.output_prefix = f"optuna_{args.sample}_{args.algorithm}_{args.mode}_{args.metrics}"

    higher = HIGHER_IS_BETTER[args.metrics]
    base_config = load_config(args.config)
    workflows = select_workflows(args.sample, args.mode,
                                 base_config.get("normalizations", []),
                                 base_config.get("features", {}).get("methods", []),
                                 args.norm, args.feat)

    objective = SpatialObjective(base_config, args, workflows)
    sampler = optuna.samplers.TPESampler(seed=args.seed)
    study = optuna.create_study(direction="maximize" if higher else "minimize", sampler=sampler)

    if args.algorithm == "BANKSY":
        study.enqueue_trial({"lambda": 0.2, "k_geom": 30, "resolution": 0.8})
        study.enqueue_trial({"lambda": 0.4, "k_geom": 30, "resolution": 0.8})
        study.enqueue_trial({"lambda": 0.1, "k_geom": 30, "resolution": 0.8})
        study.enqueue_trial({"lambda": 0.2, "k_geom": 70, "resolution": 0.8})

    callbacks = []
    if args.target is not None:
        callbacks.append(make_target_callback(args.target, higher, args.metrics))

    logging.info(f"Starting: sample={args.sample} algorithm={args.algorithm} mode={args.mode} "
                 f"metric={args.metrics} ({'maximise' if higher else 'minimise'}, source={args.source}) "
                 f"workflows={len(workflows)} trials={args.trials} target={args.target} "
                 f"n_parallel={args.n_parallel} cores_per_trial={objective.cores_per_trial}")
    logging.info(f"Workflows: {', '.join(wf_name(n, f) for n, f in workflows)}")

    try:
        study.optimize(objective, n_trials=args.trials, n_jobs=args.n_parallel,
                       callbacks=callbacks, catch=(TrialFailed,))
    except KeyboardInterrupt:
        logging.warning("Interrupted; writing results collected so far.")

    with objective._lock:
        objective.write_outputs()

    done = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
    print("\nOptimisation finished.")
    print(f"Trials run      : {len(study.trials)}  (completed {len(done)})")
    if "stopped_on_target" in study.user_attrs:
        print(f"Stopped early   : trial {study.user_attrs['stopped_on_target']} reached the target {args.target}")
    if not done:
        print("No trial completed; see the log above.")
        return
    best = study.best_trial
    label = "mean " + args.metrics if args.mode == "group" else args.metrics
    print(f"Best trial      : {best.number}")
    print(f"Best parameters : {best.params}")
    print(f"Best {label:<10}: {best.value:.4f}")
    print("Per workflow    :")
    for k, v in best.user_attrs.get("per_workflow", {}).items():
        print(f"    {k:<28} {('NA' if v is None else f'{v:.4f}')}")
    print(f"Results         : {args.output_prefix}_trials.csv, {args.output_prefix}_workflows.csv")


if __name__ == "__main__":
    main()
