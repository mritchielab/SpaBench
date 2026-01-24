#!/usr/bin/env python3
import argparse
import yaml
import subprocess
import pandas as pd
import numpy as np
import os
import sys
import optuna 
import logging
from copy import deepcopy

# --- Configuration & Setup ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def suggest_params(trial, alg):
    params = {}
    if alg == "BANKSY":
        params["lambda"] = trial.suggest_categorical("lambda", [0.1, 0.2,0.3, 0.4, 0.5, 0.6,0.7, 0.8, 0.9])
        params["k_geom"] = trial.suggest_categorical("k_geom", [10, 15, 25, 30, 40, 50, 60, 70, 80])
        params["resolution"] = trial.suggest_categorical("resolution", [0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0])
    elif alg == "BayesSpace":
        # q is the number of clusters. 
        # Usually, we set this to the ground truth number, but if you want to optimize:
        params["q"] = trial.suggest_int("q", 5, 20) 
        # nrep and burn_in affect convergence
        params["nrep"] = trial.suggest_categorical("nrep", [500, 1000, 2000, 3000, 4000, 5000])
        params["burn_in"] = trial.suggest_categorical("burn_in", [100, 300, 500, 700, 1000])
        params["gamma"] = trial.suggest_categorical("gamma", [1, 2, 3, 4, 5])
        # gamma is the spatial smoothing parameter (VERY important for BayesSpace)
        # params["gamma"] = trial.suggest_float("gamma", 1.0, 5.0)
    elif alg == "SpatialLeiden":
        params["resolution"] = trial.suggest_categorical("resolution", [0.5, 1.0, 1.5, 2.0])
        params["layer_ratio"] = trial.suggest_categorical("layer_ratio", [0.5, 1.0, 1.8, 2.5, 3.0])
        params["n_neighbors"] = trial.suggest_categorical("n_neighbors", [10, 15, 30])
    # Add other algorithms as needed
    return params

def load_config(config_path="config.yaml"):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def generate_target_list(sample, alg, norms, feats):
    targets = []
    for n in norms:
        for f in feats:
            t = f"results/{sample}/{n}/{f}/{alg}/09_{alg}_metrics.csv"
            targets.append(t)
    return targets

def clean_clustering_results(sample, alg, norms, feats):
    """
    Deletes clustering-related results only. 
    By deleting the metrics CSVs, we force Snakemake to rerun the clustering 
    pipeline for the new trial parameters.
    """
    for n in norms:
        for f in feats:
            alg_dir = f"results/{sample}/{n}/{f}/{alg}"
            if os.path.exists(alg_dir):
                logging.info(f"Cleaning clustering directory: {alg_dir}")
                for filename in os.listdir(alg_dir):
                    file_path = os.path.join(alg_dir, filename)
                    try:
                        if os.path.isfile(file_path) or os.path.islink(file_path):
                            os.remove(file_path)
                    except Exception as e:
                        logging.warning(f"Failed to delete {file_path}: {e}")

def get_mean_nmi(target_files):
    nmis = []
    for f in target_files:
        if os.path.exists(f):
            try:
                df = pd.read_csv(f)
                col_map = {c.lower(): c for c in df.columns}
                if 'nmi' in col_map:
                    val = float(df[col_map['nmi']].iloc[0])
                    nmis.append(val)
            except Exception:
                pass
    return np.mean(nmis) if nmis else 0.0

class SpatialObjective:
    def __init__(self, base_config, args, norms, feats):
        self.base_config = base_config
        self.args = args
        self.norms = norms
        self.feats = feats
        self.results_history = [] 

    def __call__(self, trial):
        params = suggest_params(trial, self.args.algorithm)
        run_config = deepcopy(self.base_config)
        
        # Ensure the algorithm section exists in config
        if "clustering" not in run_config:
            run_config["clustering"] = {}
        
        # Update config with Optuna suggestions
        for k, v in params.items():
            run_config["clustering"][k] = v
        
        temp_config = f"temp_opt_{self.args.sample}_{self.args.algorithm}_{trial.number}.yaml"
        with open(temp_config, 'w') as f:
            yaml.dump(run_config, f)

        targets = generate_target_list(self.args.sample, self.args.algorithm, self.norms, self.feats)
        
        # 1. Clean only the metrics to trigger the 'compute_metrics' rule
        clean_clustering_results(self.args.sample, self.args.algorithm, self.norms, self.feats)

        # 2. Refined Snakemake command
        cluster_rule = f"cluster_{self.args.algorithm}"
        
        cmd = [
            "snakemake",
            *targets,  # Use * to unpack the list of targets here
            "--cores", str(self.args.cores),
            "--configfile", temp_config,
            "--rerun-incomplete",
            "--nolock",
            "-p",
            "--allowed-rules", cluster_rule, "annotate_clusters", "compute_metrics", "reduce",
            "-R", cluster_rule
        ]

        logging.info(f"Trial {trial.number}: Running {self.args.algorithm} with {params}")
        
        success = True
        try:
            # We use check=True to catch execution errors
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            success = False
            logging.error(f"Trial {trial.number} FAILED! Error: {e}")

        mean_nmi = get_mean_nmi(targets) if success else 0.0
        
        if os.path.exists(temp_config): 
            os.remove(temp_config)

        record = {"trial_id": trial.number, "nmi": mean_nmi, "success": success}
        record.update(params)
        self.results_history.append(record)
        return mean_nmi

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True)
    p.add_argument("--algorithm", required=True)
    p.add_argument("--default_nmi", type=float, required=True)
    p.add_argument("--trials", type=int, default=10)
    p.add_argument("--cores", type=int, default=16)
    p.add_argument("--config", default="config.yaml")
    p.add_argument("--output_csv", default=None) 
    args = p.parse_args()

    if args.output_csv is None:
        args.output_csv = f"benchmark_results_{args.sample}_{args.algorithm}.csv"

    base_config = load_config(args.config)
    norms = base_config.get("normalizations", [])
    feats = base_config.get("features", {}).get("methods", [])

    objective = SpatialObjective(base_config, args, norms, feats)
    study = optuna.create_study(direction="maximize")

    if args.algorithm == "BANKSY":
        study.enqueue_trial({"lambda": 0.2, "k_geom": 30, "resolution": 0.8})
        study.enqueue_trial({"lambda": 0.4, "k_geom": 30, "resolution": 0.8})
        study.enqueue_trial({"lambda": 0.1, "k_geom": 30, "resolution": 0.8})
        study.enqueue_trial({"lambda": 0.2, "k_geom": 70, "resolution": 0.8})

    logging.info(f"Starting optimization. Results will save to {args.output_csv}")
    
    try:
        study.optimize(objective, n_trials=args.trials)
    except KeyboardInterrupt:
        pass

    if objective.results_history:
        df_results = pd.DataFrame(objective.results_history)
        df_results.sort_values(by="nmi", ascending=False, inplace=True)
        df_results.to_csv(args.output_csv, index=False)
        print(f"\nOptimization Complete. Best Parameters: {study.best_params}")
    else:
        print("\nOptimization finished but no results were recorded.")

if __name__ == "__main__":
    main()
