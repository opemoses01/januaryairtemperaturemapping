"""Run the Earth Engine extraction stage (also downloads compact raster stacks)."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from gee_extract import extract_download
import yaml
if __name__=='__main__':
    cfg=yaml.safe_load((ROOT/'config.yaml').read_text(encoding='utf-8'))
    extract_download(ROOT,cfg)
