"""Compute autocorrelation-aware Mann–Kendall and FDR trend results."""
from pathlib import Path
import sys,yaml
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from analysis import run_analysis
if __name__=='__main__':run_analysis(ROOT,yaml.safe_load((ROOT/'config.yaml').read_text(encoding='utf-8')))
