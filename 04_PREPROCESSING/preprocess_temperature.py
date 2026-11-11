"""Validate downloaded stacks; analysis performs preprocessing consistently."""
from pathlib import Path
import rasterio,yaml
ROOT=Path(__file__).resolve().parents[1]; C=yaml.safe_load((ROOT/'config.yaml').read_text(encoding='utf-8'))
for label in ('main','baseline'):
    p=ROOT/C['paths']['data_dir']/f'january_temperature_{label}.tif'
    if not p.exists(): raise SystemExit(f'Missing {p}. Run 02_GEE_EXTRACTION/gee_temperature_export.py first.')
    with rasterio.open(p) as ds: print(p.name,ds.count,'years;',ds.width,'x',ds.height,ds.crs)
print('Validated. The analysis stage computes climatology, anomalies and derived rasters.')
