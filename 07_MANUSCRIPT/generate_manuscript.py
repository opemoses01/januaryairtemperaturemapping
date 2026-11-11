"""Generate the data-driven Word manuscript section by rerunning analysis."""
from pathlib import Path
import sys,yaml
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from analysis import run_analysis
if __name__=='__main__':
    cfg=yaml.safe_load((ROOT/'config.yaml').read_text(encoding='utf-8'))
    run_analysis(ROOT,cfg)
    print('Generated results section and Word document in outputs/Manuscript/.')
