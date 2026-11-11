from pathlib import Path
import argparse, subprocess, sys, yaml, logging, json, time

ROOT = Path(__file__).resolve().parent
CONFIG = ROOT / "config.yaml"

def load_config():
    with CONFIG.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def main():
    parser = argparse.ArgumentParser(description="Automated January 2m air-temperature research workflow")
    parser.add_argument("--skip-gee", action="store_true", help="Use already downloaded temperature stacks")
    parser.add_argument("--skip-analysis", action="store_true", help="Run GEE extraction/download only")
    args = parser.parse_args()

    cfg = load_config()
    project_id = cfg["project"].get("gee_cloud_project", "")
    if not args.skip_gee and (not project_id or project_id.startswith("REPLACE_")):
        raise SystemExit(
            "Set project.gee_cloud_project in config.yaml to an Earth Engine-enabled "
            "Google Cloud project ID before running extraction."
        )
    log_dir = ROOT / cfg["paths"]["output_dir"]
    log_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(log_dir / "automation.log", encoding="utf-8"),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logging.info("Starting %s", cfg["project"]["name"])

    boundary = ROOT / cfg["project"]["boundary_file"]
    if not boundary.exists():
        logging.info("Study boundary missing; preparing the configured study boundary.")
        from runpy import run_path
        run_path(str(ROOT / "01_BOUNDARY" / "prepare_boundary.py"), run_name="__main__")
        if not boundary.exists():
            raise FileNotFoundError(f"Boundary preparation did not create {boundary}")

    if not args.skip_gee:
        logging.info("Stage 1: GEE extraction and download")
        from gee_extract import extract_download
        extract_download(ROOT, cfg)

    if not args.skip_analysis:
        logging.info("Stage 2: preprocessing, statistics, mapping and manuscript")
        from analysis import run_analysis
        run_analysis(ROOT, cfg)

    logging.info("Workflow finished. Check %s.", log_dir)
    print("\nCompleted. Main deliverables are in:", ROOT / cfg["paths"]["output_dir"])

if __name__ == "__main__":
    main()
