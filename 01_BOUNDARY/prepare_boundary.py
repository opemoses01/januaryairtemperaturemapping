"""Download Nigeria ADM1 from geoBoundaries and select the six study states."""
from pathlib import Path
import requests, geopandas as gpd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'01_BOUNDARY/southwest_nigeria.geojson'
STATES={'Oyo','Ogun','Osun','Lagos','Ekiti','Ondo'}
API='https://www.geoboundaries.org/api/current/gbOpen/NGA/ADM1/'

def main():
    meta=requests.get(API,timeout=60); meta.raise_for_status(); info=meta.json()
    url=info.get('gjDownloadURL')
    if not url: raise RuntimeError(f'geoBoundaries metadata did not provide a GeoJSON URL: {info}')
    response=requests.get(url,timeout=180); response.raise_for_status()
    src=ROOT/'work_nigeria_adm1.geojson'; src.write_bytes(response.content)
    gdf=gpd.read_file(src)
    field=next((f for f in ['shapeName','shapeName_en','name','NAME_1','admin1Name'] if f in gdf.columns),None)
    if field is None: raise ValueError(f'No recognized state-name field in {list(gdf.columns)}')
    name=gdf[field].astype(str).str.strip().str.replace(r'\s+State$','',regex=True)
    selected=gdf[name.isin(STATES)].copy(); selected['shapeName']=name[name.isin(STATES)].values
    missing=STATES-set(selected.shapeName)
    if missing: raise ValueError(f'States missing from boundary source: {sorted(missing)}; names found: {sorted(name.unique())}')
    selected=selected.to_crs('EPSG:4326'); selected.to_file(OUT,driver='GeoJSON')
    src.unlink(missing_ok=True)
    print(f'Wrote six state features to {OUT}; source: {info.get("boundaryID",API)}')

if __name__=='__main__': main()
