from pathlib import Path
from datetime import date, timedelta
from io import BytesIO
import getpass, json, logging, requests, zipfile
import geopandas as gpd
import ee
import rasterio

def _init_ee(project_id, cfg):
    try:
        ee.Initialize(project=project_id)
        return
    except Exception as init_error:
        client_id = cfg["project"].get("oauth_client_id")
        if not client_id:
            raise RuntimeError(
                "Earth Engine credentials are unavailable. Configure a supported OAuth client "
                "before running extraction."
            ) from init_error
        secret = getpass.getpass("Enter the rotated OAuth client secret (input is hidden): ")
        if not secret.strip():
            raise RuntimeError("No OAuth client secret was entered.") from init_error
        # Earth Engine's helper persists the refresh token in the user's standard
        # credentials directory, not in this project. Never put the secret in YAML.
        ee.oauth.CLIENT_ID = client_id
        ee.oauth.CLIENT_SECRET = secret.strip()
        ee.Authenticate(auth_mode="localhost:8085", force=True)
        ee.Initialize(project=project_id)

def _ee_geometry(shapefile):
    gdf = gpd.read_file(shapefile)
    if gdf.empty:
        raise ValueError("The boundary shapefile contains no features.")
    if gdf.crs is None:
        raise ValueError("The boundary shapefile has no CRS. Define its correct CRS before running.")
    gdf = gdf.to_crs("EPSG:4326")
    geojson = json.loads(gdf.to_json())
    fc = ee.FeatureCollection(geojson)
    return fc.geometry()

def _daily_mean_for_january(year, geometry, collection_id, band):
    """Average hourly observations within each day, then average daily means."""
    start = ee.Date.fromYMD(year, 1, 1)
    day_offsets = ee.List.sequence(0, 30)

    def daily_mean(day_offset):
        offset = ee.Number(day_offset)
        day_start = start.advance(offset, "day")
        day_end = day_start.advance(1, "day")
        return (ee.ImageCollection(collection_id)
                .filterDate(day_start, day_end)
                .select(band)
                .mean()
                .rename(f"jan_{year}"))

    # Mapping the 31 days server-side keeps the request graph small while
    # retaining equal weighting of daily hourly means.
    daily_images = day_offsets.map(daily_mean)
    monthly = ee.ImageCollection.fromImages(daily_images).mean()
    return monthly.subtract(273.15).rename(f"jan_{year}").clip(geometry).toFloat()

def _download_stack(images, years, geometry, out_path, cfg):
    stack = ee.Image.cat(images).rename([f"jan_{y}" for y in years]).toFloat()
    params = {
        "name": out_path.stem,
        "bands": [f"jan_{y}" for y in years],
        "region": geometry,
        "scale": cfg["gee"]["native_scale_m"],
        "crs": cfg["gee"]["crs"],
        "format": "GEO_TIFF"
    }
    url = stack.getDownloadURL(params)
    response = requests.get(url, timeout=600)
    response.raise_for_status()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = response.content
    if zipfile.is_zipfile(BytesIO(payload)):
        with zipfile.ZipFile(BytesIO(payload)) as archive:
            tif_names = [n for n in archive.namelist() if n.lower().endswith((".tif", ".tiff"))]
            if not tif_names:
                raise RuntimeError(f"Earth Engine returned a ZIP without a GeoTIFF for {out_path.name}.")
            payload = archive.read(tif_names[0])
    temp_path = out_path.with_suffix(out_path.suffix + ".partial")
    temp_path.write_bytes(payload)
    try:
        with rasterio.open(temp_path) as raster:
            if raster.count != len(years):
                raise RuntimeError(f"Expected {len(years)} raster bands, received {raster.count} for {out_path.name}.")
            if raster.width < 1 or raster.height < 1:
                raise RuntimeError(f"Downloaded raster has an empty grid: {out_path.name}.")
        temp_path.replace(out_path)
    except Exception:
        temp_path.unlink(missing_ok=True)
        raise
    if not out_path.exists() or out_path.stat().st_size < 1000:
        raise RuntimeError(f"Download appears incomplete: {out_path}")
    logging.info("Downloaded %s (%.2f MB)", out_path, out_path.stat().st_size / 1e6)

def extract_download(root, cfg):
    project_id = cfg["project"]["gee_cloud_project"]
    if not project_id or project_id.startswith("REPLACE_"):
        raise ValueError("Set project.gee_cloud_project in config.yaml to your Google Cloud Project ID.")
    _init_ee(project_id, cfg)

    boundary = root / cfg["project"]["boundary_file"]
    geometry = _ee_geometry(boundary)
    collection_id = cfg["gee"]["collection"]
    band = cfg["gee"]["band"]
    data_dir = root / cfg["paths"]["data_dir"]
    data_dir.mkdir(parents=True, exist_ok=True)

    main_years = list(range(cfg["study"]["main_start_year"], cfg["study"]["main_end_year"] + 1))
    baseline_years = list(range(cfg["study"]["baseline_start_year"], cfg["study"]["baseline_end_year"] + 1))

    # January is a 31-day month. The daily-mean approach is retained to match
    # the manuscript methodology rather than averaging an arbitrary time window.
    for label, years in [("main", main_years), ("baseline", baseline_years)]:
        target = data_dir / f"january_temperature_{label}.tif"
        if target.exists() and target.stat().st_size > 1000:
            try:
                with rasterio.open(target) as raster:
                    if raster.count == len(years) and raster.width > 0 and raster.height > 0:
                        logging.info("Found valid existing %s stack; skipping GEE download.", target.name)
                        continue
            except rasterio.errors.RasterioError:
                logging.warning("Existing %s is not a readable GeoTIFF; replacing it.", target.name)
        logging.info("Building %s stack: %s-%s", label, years[0], years[-1])
        images = [
            _daily_mean_for_january(y, geometry, collection_id, band)
            for y in years
        ]
        _download_stack(images, years, geometry, target, cfg)

    metadata = {
        "collection": collection_id,
        "band": band,
        "unit": "degrees Celsius",
        "main_years": main_years,
        "baseline_years": baseline_years,
        "aggregation": "hourly observations -> daily mean -> January mean of daily means",
        "source_resolution_m": cfg["gee"]["native_scale_m"],
        "crs": cfg["gee"]["crs"]
    }
    (data_dir / "extraction_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
