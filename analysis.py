from pathlib import Path
import json, logging, math, re, warnings
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.transform import xy
from rasterio.features import geometry_mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
from scipy.stats import false_discovery_control, theilslopes, linregress
from statsmodels.stats.diagnostic import acorr_ljungbox
from pyproj import CRS, Transformer
from shapely.geometry import Polygon
from shapely.ops import transform as transform_geometry
import pymannkendall as mk
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from docx import Document
from docx.shared import Inches, Pt

EXPECTED_STATES = ["Oyo", "Ogun", "Osun", "Lagos", "Ekiti", "Ondo"]

def read_stack(path):
    with rasterio.open(path) as src:
        arr = src.read().astype("float32")
        profile = src.profile.copy()
        transform = src.transform
        crs = src.crs
        nodata = src.nodata
    if nodata is not None:
        arr[arr == nodata] = np.nan
    arr[~np.isfinite(arr)] = np.nan
    return arr, profile, transform, crs

def save_raster(path, array, profile):
    p = profile.copy()
    p.update(count=1, dtype="float32", nodata=-9999, compress="deflate")
    out = np.where(np.isfinite(array), array, -9999).astype("float32")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(path, "w", **p) as dst:
        dst.write(out, 1)

def resample_for_display(array, profile, target_crs, resolution_m, resampling=Resampling.bilinear):
    """Make an explicitly interpolated display copy; this does not add information."""
    bounds = rasterio.transform.array_bounds(array.shape[0], array.shape[1], profile["transform"])
    dst_transform, width, height = calculate_default_transform(
        profile["crs"], target_crs, array.shape[1], array.shape[0], *bounds,
        resolution=resolution_m
    )
    result = np.full((height, width), np.nan, dtype="float32")
    reproject(
        source=array.astype("float32"), destination=result,
        src_transform=profile["transform"], src_crs=profile["crs"], src_nodata=np.nan,
        dst_transform=dst_transform, dst_crs=target_crs, dst_nodata=np.nan,
        resampling=resampling
    )
    out_profile = profile.copy()
    out_profile.update(
        crs=target_crs, transform=dst_transform, width=width, height=height,
        count=1, dtype="float32", nodata=-9999, compress="deflate"
    )
    return result, out_profile

def detect_state_field(gdf, requested_field, mapping):
    if requested_field and requested_field in gdf.columns:
        return requested_field
    aliases = {s.lower() for s in EXPECTED_STATES}
    for col in gdf.columns:
        if col == gdf.geometry.name:
            continue
        vals = set(gdf[col].dropna().astype(str).str.strip().str.lower())
        if len(vals.intersection(aliases)) >= 4:
            return col
    # Search likely state-name columns with common suffixes.
    for col in gdf.columns:
        if col == gdf.geometry.name:
            continue
        vals = set(gdf[col].dropna().astype(str).str.strip().str.lower())
        hits = sum(any(v == s or v == s + " state" for v in vals) for s in aliases)
        if hits >= 4:
            return col
    raise ValueError(
        "Could not automatically identify the state-name field. "
        "Set study.state_name_field in config.yaml to the correct shapefile attribute."
    )

def normalise_states(gdf, field, mapping):
    def clean(v):
        raw = str(v).strip()
        if raw in mapping:
            return mapping[raw]
        low = raw.lower()
        for s in EXPECTED_STATES:
            if low in (s.lower(), (s + " state").lower()):
                return s
        return mapping.get(raw.upper(), raw)
    gdf = gdf.copy()
    gdf["_state_norm"] = gdf[field].apply(clean)
    return gdf[gdf["_state_norm"].isin(EXPECTED_STATES)].copy()

def weighted_mean(values, weights):
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not good.any():
        return np.nan
    w = weights[good]
    return float(np.sum(values[good] * w) / np.sum(w))

def weighted_sd(values, weights):
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if good.sum() < 2:
        return np.nan
    v, w = values[good], weights[good]
    mean = np.sum(v*w) / np.sum(w)
    return float(np.sqrt(np.sum(w*(v-mean)**2) / np.sum(w)))

def weighted_median(values, weights):
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not good.any():
        return np.nan
    v, w = values[good], weights[good]
    order = np.argsort(v)
    v, w = v[order], w[order]
    return float(v[np.searchsorted(np.cumsum(w), w.sum()/2, side="left")])

def equal_area_coverage_weights(shape, transform, source_crs, region_geom, states, equal_area_crs):
    """Return m² overlap weights for the region and each state, including partial cells."""
    rows, cols = shape
    src = CRS.from_user_input(source_crs)
    dst = CRS.from_user_input(equal_area_crs)
    projector = Transformer.from_crs(src, dst, always_xy=True).transform
    region_equal = transform_geometry(projector, region_geom)
    state_equal = {
        row["_state_norm"]: transform_geometry(projector, row.geometry)
        for _, row in states.iterrows()
    }
    region_weights = np.zeros(shape, dtype="float64")
    state_weights = {name: np.zeros(shape, dtype="float64") for name in state_equal}
    bx0, by0, bx1, by1 = region_geom.bounds
    for r in range(rows):
        for c in range(cols):
            corners = [
                xy(transform, r, c, offset="ul"),
                xy(transform, r, c, offset="ur"),
                xy(transform, r, c, offset="lr"),
                xy(transform, r, c, offset="ll"),
            ]
            xs, ys = zip(*corners)
            if max(xs) <= bx0 or min(xs) >= bx1 or max(ys) <= by0 or min(ys) >= by1:
                continue
            cell = transform_geometry(projector, Polygon(corners))
            if cell.is_empty or cell.area <= 0:
                continue
            region_area = cell.intersection(region_equal).area
            if region_area <= 0:
                continue
            region_weights[r, c] = region_area
            for name, state_geom in state_equal.items():
                state_area = cell.intersection(state_geom).area
                if state_area > 0:
                    state_weights[name][r, c] = state_area
    return region_weights, state_weights

MAP_CRS = CRS.from_epsg(32631)  # WGS 84 / UTM zone 31N; grid north is shown.

def warp_map_array_to_utm(array, profile, categorical=False):
    """Reproject a display layer to UTM 31N for map graticules and metric axes."""
    source_crs = CRS.from_user_input(profile["crs"])
    height, width = array.shape
    bounds = rasterio.transform.array_bounds(height, width, profile["transform"])
    dst_transform, dst_width, dst_height = calculate_default_transform(
        source_crs, MAP_CRS, width, height, *bounds
    )
    output = np.full((dst_height, dst_width), np.nan, dtype="float32")
    reproject(
        source=np.asarray(array, dtype="float32"), destination=output,
        src_transform=profile["transform"], src_crs=source_crs, src_nodata=np.nan,
        dst_transform=dst_transform, dst_crs=MAP_CRS, dst_nodata=np.nan,
        resampling=Resampling.nearest if categorical else Resampling.bilinear,
    )
    return output, dst_transform

def mask_to_study_area(array, transform, states):
    """Mask reprojected display pixels outside the six-state study boundary."""
    projected = states.to_crs(MAP_CRS)
    inside = geometry_mask(
        [geom for geom in projected.geometry if geom is not None and not geom.is_empty],
        out_shape=array.shape, transform=transform, invert=True, all_touched=False,
    )
    result = array.copy()
    result[~inside] = np.nan
    return result

def decorate_utm_axes(ax, states, title):
    """Apply consistent UTM grid, state labels, grid-north arrow and scale bar."""
    projected = states.to_crs(MAP_CRS)
    projected.boundary.plot(ax=ax, color="#263746", linewidth=0.75, zorder=4)
    name_col = "_state_norm" if "_state_norm" in projected.columns else "NAME_1"
    for _, row in projected.iterrows():
        point = row.geometry.representative_point()
        ax.annotate(
            str(row[name_col]), (point.x, point.y), ha="center", va="center",
            fontsize=8, fontweight="bold", color="#172b3a", zorder=5,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.82),
        )
    minx, miny, maxx, maxy = projected.total_bounds
    padx, pady = (maxx-minx)*0.055, (maxy-miny)*0.055
    ax.set_xlim(minx-padx, maxx+padx)
    ax.set_ylim(miny-pady, maxy+pady)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("UTM Zone 31N Easting (m)")
    ax.set_ylabel("UTM Zone 31N Northing (m)")
    ax.xaxis.set_major_locator(MultipleLocator(50000))
    ax.yaxis.set_major_locator(MultipleLocator(50000))
    for axis in (ax.xaxis, ax.yaxis):
        formatter = ScalarFormatter(useOffset=False)
        formatter.set_scientific(False)
        axis.set_major_formatter(formatter)
    ax.grid(True, which="major", color="#526575", alpha=0.30, linewidth=0.55, linestyle=":", zorder=3)
    ax.tick_params(labelsize=7.5, colors="#34495e")
    ax.set_title(title, fontsize=11, weight="bold", pad=10)
    # Arrow points to projected grid north, consistent with the UTM coordinate grid.
    ax.annotate("", xy=(0.95, 0.94), xytext=(0.95, 0.83), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="#172b3a", lw=1.4), zorder=7)
    ax.text(0.95, 0.955, "N", transform=ax.transAxes, ha="center", va="bottom",
            fontsize=9, fontweight="bold", color="#172b3a", zorder=7)
    # A metric scale bar for the projected map frame.
    x0 = minx + padx*0.35
    y0 = miny + pady*0.65
    bar = 50000.0
    ax.plot([x0, x0+bar], [y0, y0], color="white", linewidth=4, zorder=6)
    ax.plot([x0, x0+bar], [y0, y0], color="#172b3a", linewidth=2, zorder=7)
    ax.text(x0+bar/2, y0+pady*0.15, "50 km", ha="center", va="bottom", fontsize=7.5,
            color="#172b3a", zorder=7, bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.0))

def make_map(array, profile, gdf, title, out, cmap="viridis", vmin=None, vmax=None,
             diverging=False, label="Temperature (°C)", dpi=300):
    map_array, transform = warp_map_array_to_utm(array, profile)
    map_array = mask_to_study_area(map_array, transform, gdf)
    vals = map_array[np.isfinite(map_array)]
    if vals.size and vmin is None: vmin = float(np.nanpercentile(vals, 2))
    if vals.size and vmax is None: vmax = float(np.nanpercentile(vals, 98))
    if diverging and vmin is not None and vmax is not None:
        lim = max(abs(vmin), abs(vmax))
        norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
        image_kwargs = {"norm": norm}
    else:
        image_kwargs = {"vmin": vmin, "vmax": vmax}
    h, w = map_array.shape
    extent = [transform.c, transform.c+w*transform.a, transform.f+h*transform.e, transform.f]
    fig, ax = plt.subplots(figsize=(9, 8))
    im = ax.imshow(map_array, extent=extent, origin="upper", cmap=cmap, **image_kwargs)
    decorate_utm_axes(ax, gdf, title)
    cb = fig.colorbar(im, ax=ax, shrink=0.78, pad=0.03)
    cb.set_label(label)
    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

def make_categorical_map(array, profile, gdf, title, out, cmap, norm, ticks, ticklabels, dpi=300):
    map_array, transform = warp_map_array_to_utm(array, profile, categorical=True)
    map_array = mask_to_study_area(map_array, transform, gdf)
    h, w = map_array.shape
    extent = [transform.c, transform.c+w*transform.a, transform.f+h*transform.e, transform.f]
    fig, ax = plt.subplots(figsize=(9, 8))
    im = ax.imshow(map_array, extent=extent, origin="upper", cmap=cmap, norm=norm,
                   interpolation="nearest")
    decorate_utm_axes(ax, gdf, title)
    cb = fig.colorbar(im, ax=ax, ticks=ticks, shrink=0.78, pad=0.03)
    cb.ax.set_yticklabels(ticklabels)
    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

def run_analysis(root, cfg):
    root = Path(root)
    outroot = root / cfg["paths"]["output_dir"]
    for folder in ["Rasters", "Rasters/1km_resampled_display", "Maps/Annual", "Maps/Annual_1km_resampled",
                   "Maps/Anomalies", "Maps/Anomalies_1km_resampled", "Maps/Trends", "Maps/Trends_1km_resampled",
                   "Graphs", "Tables", "Manuscript"]:
        (outroot / folder).mkdir(parents=True, exist_ok=True)

    main_path = root / cfg["paths"]["data_dir"] / "january_temperature_main.tif"
    base_path = root / cfg["paths"]["data_dir"] / "january_temperature_baseline.tif"
    if not main_path.exists() or not base_path.exists():
        raise FileNotFoundError("Missing GEE temperature stack(s). Run without --skip-gee first.")

    main, profile, transform, crs = read_stack(main_path)
    baseline, base_profile, _, base_crs = read_stack(base_path)
    if main.shape[1:] != baseline.shape[1:] or crs != base_crs or transform != base_profile["transform"]:
        raise ValueError("Main and baseline raster grids do not match. Re-export both stacks using the same boundary/scale/CRS.")

    main_years = list(range(cfg["study"]["main_start_year"], cfg["study"]["main_end_year"] + 1))
    base_years = list(range(cfg["study"]["baseline_start_year"], cfg["study"]["baseline_end_year"] + 1))
    baseline_label = f"{base_years[0]}_{base_years[-1]}"
    baseline_label_display = f"{base_years[0]}–{base_years[-1]}"
    if main.shape[0] != len(main_years) or baseline.shape[0] != len(base_years):
        raise ValueError("Raster band count does not match configured year ranges.")

    boundary = gpd.read_file(root / cfg["project"]["boundary_file"])
    if boundary.crs is None:
        raise ValueError("Boundary CRS is undefined.")
    boundary = boundary.to_crs(crs)
    field = detect_state_field(boundary, cfg["study"].get("state_name_field", ""), cfg["study"].get("state_name_mapping", {}))
    states = normalise_states(boundary, field, cfg["study"].get("state_name_mapping", {}))
    missing_states = sorted(set(EXPECTED_STATES) - set(states["_state_norm"].unique()))
    if missing_states:
        raise ValueError(f"Shapefile is missing separate state features for: {missing_states}. A state attribute and six state geometries are required.")
    # Dissolve multipart/duplicate records into one geometry per state and
    # restrict all subsequent processing to the six requested states.
    states = states.dissolve(by="_state_norm", as_index=False)
    boundary = states.copy()
    region_geom = boundary.geometry.union_all()
    region_weights, state_weights = equal_area_coverage_weights(
        main.shape[1:], transform, crs, region_geom, states,
        cfg["analysis"]["equal_area_crs"]
    )
    inside = region_weights > 0
    main[:, ~inside] = np.nan
    baseline[:, ~inside] = np.nan

    # Main descriptive rasters.
    climatology = np.nanmean(baseline, axis=0)
    annual_mean = np.nanmean(main, axis=0)
    annual_min = np.nanmin(main, axis=0)
    annual_max = np.nanmax(main, axis=0)
    annual_sd = np.nanstd(main, axis=0, ddof=1)
    # CV uses absolute temperature in kelvin so the denominator has a physical zero.
    annual_cv = annual_sd / (annual_mean + 273.15) * 100.0
    for name, arr in [
        (f"january_climatology_{baseline_label}.tif", climatology),
        ("january_mean_2000_2026.tif", annual_mean),
        ("january_min_2000_2026.tif", annual_min),
        ("january_max_2000_2026.tif", annual_max),
        ("january_standard_deviation_2000_2026.tif", annual_sd),
        ("january_coefficient_of_variation_absolute_pct_2000_2026.tif", annual_cv),
    ]:
        save_raster(outroot / "Rasters" / name, arr, profile)

    # Per-year anomalies relative to each pixel's configured January baseline.
    anomalies = main - climatology[None, :, :]
    for idx, year in enumerate(main_years):
        # Standalone annual January temperature GeoTIFF, preserving the source grid.
        save_raster(outroot / "Rasters" / f"january_temperature_{year}.tif", main[idx], profile)
        save_raster(outroot / "Rasters" / f"anomaly_january_{year}.tif", anomalies[idx], profile)

    # Regional and state annual summaries.
    regional_rows, state_rows = [], []
    for idx, year in enumerate(main_years):
        regional_rows.append({
            "year": year,
            "regional_mean_C": weighted_mean(main[idx], region_weights),
            "regional_min_C": float(np.nanmin(main[idx][inside])),
            "regional_max_C": float(np.nanmax(main[idx][inside])),
            "regional_sd_C": weighted_sd(main[idx], region_weights),
            "regional_anomaly_C": weighted_mean(anomalies[idx], region_weights),
            "valid_grid_cells": int(np.isfinite(main[idx][inside]).sum())
        })
        for _, row in states.iterrows():
            state_w = state_weights[row["_state_norm"]]
            mask = state_w > 0
            state_rows.append({
                "state": row["_state_norm"], "year": year,
                "mean_temperature_C": weighted_mean(main[idx], state_w),
                "mean_anomaly_C": weighted_mean(anomalies[idx], state_w),
                "sd_temperature_C": weighted_sd(main[idx], state_w)
            })
    regional = pd.DataFrame(regional_rows)
    annual_min_max_avg = regional[["year", "regional_min_C", "regional_max_C", "regional_mean_C"]].rename(
        columns={
            "regional_min_C": "january_minimum_C",
            "regional_max_C": "january_maximum_C",
            "regional_mean_C": "january_area_weighted_average_C"
        }
    )
    state_ts = pd.DataFrame(state_rows)
    regional_sd = float(regional["regional_mean_C"].std(ddof=1))
    regional_mean = float(regional["regional_mean_C"].mean())
    regional_cv_abs = regional_sd / (regional_mean + 273.15) * 100
    state_descriptive = []
    for state, frame in state_ts.groupby("state"):
        values = frame["mean_temperature_C"].dropna()
        state_mean = float(values.mean()) if not values.empty else np.nan
        state_sd = float(values.std(ddof=1)) if len(values) > 1 else np.nan
        state_descriptive.append({
            "state": state,
            "mean_annual_january_temperature_C": state_mean,
            "minimum_annual_state_mean_C": float(values.min()) if not values.empty else np.nan,
            "maximum_annual_state_mean_C": float(values.max()) if not values.empty else np.nan,
            "sd_of_annual_state_means_C": state_sd,
            "cv_absolute_temperature_pct": state_sd / (state_mean + 273.15) * 100 if np.isfinite(state_sd) else np.nan,
            "mean_annual_anomaly_C": float(frame["mean_anomaly_C"].mean())
        })
    state_descriptive_df = pd.DataFrame(state_descriptive)
    regional.to_csv(outroot / "Tables" / "regional_annual_temperature.csv", index=False)
    annual_min_max_avg.to_csv(outroot / "Tables" / "annual_january_min_max_average.csv", index=False)
    state_ts.to_csv(outroot / "Tables" / "state_annual_temperature.csv", index=False)
    state_descriptive_df.to_csv(outroot / "Tables" / "state_descriptive_statistics.csv", index=False)

    # Grid-cell Mann-Kendall and Sen slope; use Hamed-Rao modified MK when
    # statistically significant lag-1 autocorrelation is detected.
    nrows, ncols = main.shape[1:]
    slope = np.full((nrows, ncols), np.nan, dtype="float32")
    mk_z = np.full_like(slope, np.nan)
    mk_p = np.full_like(slope, np.nan)
    mk_q = np.full_like(slope, np.nan)
    mk_tau = np.full_like(slope, np.nan)
    lag1 = np.full_like(slope, np.nan)
    autocorr_adjusted = np.zeros((nrows, ncols), dtype="uint8")
    alpha = float(cfg["analysis"]["alpha"])
    ac_alpha = float(cfg["analysis"]["autocorrelation_alpha"])
    min_years = int(cfg["analysis"]["minimum_valid_years"])
    p_values, locations = [], []
    logging.info("Running pixel-wise Mann-Kendall and Sen slope analysis...")
    for r in range(nrows):
        for c in range(ncols):
            x = main[:, r, c]
            valid = np.isfinite(x)
            if valid.sum() < min_years:
                continue
            x = x[valid].astype(float)
            years = np.array(main_years)[valid].astype(float)
            # Lag-one autocorrelation diagnostic on consecutive valid years.
            if len(x) > 3:
                try:
                    ac = float(np.corrcoef(x[:-1], x[1:])[0, 1])
                    acp = float(acorr_ljungbox(x, lags=[1], return_df=True)["lb_pvalue"].iloc[0])
                except Exception:
                    ac, acp = np.nan, 1.0
            else:
                ac, acp = np.nan, 1.0
            lag1[r, c] = ac
            try:
                if np.isfinite(acp) and acp < ac_alpha:
                    result = mk.hamed_rao_modification_test(x)
                    autocorr_adjusted[r, c] = 1
                else:
                    result = mk.original_test(x)
                sen = theilslopes(x, years)
                slope[r, c] = float(sen.slope)
                mk_z[r, c] = float(result.z)
                mk_p[r, c] = float(result.p)
                mk_tau[r, c] = float(result.Tau)
                p_values.append(float(result.p))
                locations.append((r, c))
            except Exception as exc:
                logging.warning("Trend test failed at pixel (%s,%s): %s", r, c, exc)

    if p_values:
        q_values = false_discovery_control(np.asarray(p_values), method="bh")
        for (r, c), q in zip(locations, q_values):
            mk_q[r, c] = q

    use_fdr = bool(cfg["analysis"]["use_fdr_for_final_classification"])
    sig = (mk_q < float(cfg["analysis"]["fdr_alpha"])) if use_fdr else (mk_p < alpha)
    trend_class = np.zeros((nrows, ncols), dtype="float32")
    # Distinguish cells with no valid trend estimate from tested, non-significant cells.
    trend_class[~inside | ~np.isfinite(slope)] = np.nan
    trend_class[sig & (slope > 0)] = 1
    trend_class[sig & (slope < 0)] = -1
    # No significant trend is 0; masked outside study area remains nodata.
    slope_valid = inside & np.isfinite(slope)
    intensity_cutoffs = np.nanquantile(slope[slope_valid], [1/3, 2/3])
    slope_intensity = np.full((nrows, ncols), np.nan, dtype="float32")
    slope_intensity[slope_valid & (slope < intensity_cutoffs[0])] = 1
    slope_intensity[slope_valid & (slope >= intensity_cutoffs[0]) & (slope < intensity_cutoffs[1])] = 2
    slope_intensity[slope_valid & (slope >= intensity_cutoffs[1])] = 3

    trend_arrays = {
        "sen_slope_C_per_year.tif": slope,
        "sen_slope_C_per_decade.tif": slope * 10.0,
        "mann_kendall_Z.tif": mk_z,
        "mann_kendall_p.tif": mk_p,
        "mann_kendall_FDR_q.tif": mk_q,
        "mann_kendall_tau.tif": mk_tau,
        "lag1_autocorrelation.tif": lag1,
        "trend_classification.tif": trend_class,
        "sen_slope_intensity_tertiles.tif": slope_intensity
    }
    for name, arr in trend_arrays.items():
        save_raster(outroot / "Rasters" / name, arr, profile)
    intensity_rows = []
    for code, label in [(1, "Lower tertile"), (2, "Middle tertile"), (3, "Upper tertile")]:
        values = slope[slope_intensity == code]
        intensity_rows.append({
            "intensity_class": label,
            "sen_slope_min_C_per_decade": float(np.min(values)),
            "sen_slope_max_C_per_decade": float(np.max(values)),
            "grid_cells": int(values.size),
            "percent_of_analysed_cells": float(100*values.size/slope_valid.sum())
        })
    pd.DataFrame(intensity_rows).to_csv(outroot / "Tables" / "sen_slope_intensity_tertiles.csv", index=False)

    # State-level trend summaries and weighted area proportions.
    state_trends = []
    for _, row in states.iterrows():
        state_w = state_weights[row["_state_norm"]]
        valid = (state_w > 0) & np.isfinite(slope) & np.isfinite(mk_p)
        warm = valid & sig & (slope > 0)
        cool = valid & sig & (slope < 0)
        nonsig = valid & ~sig
        total_w = float(np.sum(state_w[valid]))
        state_trends.append({
            "state": row["_state_norm"],
            "analysed_area_km2": total_w / 1e6,
            "mean_sen_slope_C_per_decade": weighted_mean(slope*10, state_w),
            "median_sen_slope_C_per_decade": weighted_median(slope*10, state_w),
            "significant_warming_area_km2": float(np.sum(state_w[warm])/1e6),
            "significant_cooling_area_km2": float(np.sum(state_w[cool])/1e6),
            "non_significant_area_km2": float(np.sum(state_w[nonsig])/1e6),
            "significant_warming_percent": 100*float(np.sum(state_w[warm]))/total_w if total_w else np.nan,
            "significant_cooling_percent": 100*float(np.sum(state_w[cool]))/total_w if total_w else np.nan,
            "non_significant_percent": 100*float(np.sum(state_w[nonsig]))/total_w if total_w else np.nan,
            "mean_MK_p": weighted_mean(mk_p, state_w),
            "mean_FDR_q": weighted_mean(mk_q, state_w)
        })
    state_trends_df = pd.DataFrame(state_trends)
    state_trends_df.to_csv(outroot / "Tables" / "state_trend_summary.csv", index=False)
    trend_area_rows = []
    trend_valid = inside & np.isfinite(slope) & np.isfinite(mk_p)
    trend_total_area = float(np.sum(region_weights[trend_valid]))
    for code, label in [(-1, "significant_cooling"), (0, "not_significant"), (1, "significant_warming")]:
        selected = trend_valid & (trend_class == code)
        area = float(np.sum(region_weights[selected]))
        trend_area_rows.append({
            "trend_class": label,
            "area_km2": area / 1e6,
            "percent_of_analysed_area": 100*area/trend_total_area if trend_total_area else np.nan
        })
    trend_area_df = pd.DataFrame(trend_area_rows)
    trend_area_df.to_csv(outroot / "Tables" / "regional_trend_area.csv", index=False)

    # Annual temperature and anomaly maps.
    annual_vmin = float(np.nanpercentile(main, 2))
    annual_vmax = float(np.nanpercentile(main, 98))
    anomaly_lim = float(np.nanpercentile(np.abs(anomalies), 98))
    display_resolution_m = int(cfg["analysis"].get("display_resolution_m", 1000))
    display_crs = cfg["analysis"]["equal_area_crs"]
    boundary_display = boundary.to_crs(display_crs)
    for idx, year in enumerate(main_years):
        make_map(main[idx], profile, boundary, f"January 2 m air temperature — {year}",
                 outroot / "Maps" / "Annual" / f"temperature_{year}.png",
                 cmap="turbo", vmin=annual_vmin, vmax=annual_vmax,
                 label="Air temperature (°C)", dpi=cfg["analysis"]["map_dpi"])
        make_map(anomalies[idx], profile, boundary, f"January temperature anomaly — {year}",
                 outroot / "Maps" / "Anomalies" / f"anomaly_{year}.png",
                 cmap="RdBu_r", vmin=-anomaly_lim, vmax=anomaly_lim, diverging=True,
                 label="Temperature anomaly (°C)", dpi=cfg["analysis"]["map_dpi"])
        display_arr, display_profile = resample_for_display(
            main[idx], profile, display_crs, display_resolution_m
        )
        display_name = f"january_temperature_{year}_{display_resolution_m}m_resampled.tif"
        save_raster(
            outroot / "Rasters" / "1km_resampled_display" / display_name,
            display_arr, display_profile
        )
        with rasterio.open(outroot / "Rasters" / "1km_resampled_display" / display_name, "r+") as dst:
            dst.update_tags(
                source_resolution="ERA5-Land nominal 11132 m",
                output_pixel_size_m=str(display_resolution_m),
                resampling="bilinear",
                interpretation="resampled display grid; no additional observational detail"
            )
        make_map(
            display_arr, display_profile, boundary_display,
            f"January 2 m air temperature — {year} ({display_resolution_m/1000:g} km resampled display)",
            outroot / "Maps" / "Annual_1km_resampled" / f"temperature_{year}_{display_resolution_m}m_resampled.png",
            cmap="turbo", vmin=annual_vmin, vmax=annual_vmax,
            label="Air temperature (°C); bilinear resampling", dpi=cfg["analysis"]["map_dpi"]
        )
        anomaly_display, anomaly_profile = resample_for_display(
            anomalies[idx], profile, display_crs, display_resolution_m
        )
        anomaly_name = f"anomaly_january_{year}_{display_resolution_m}m_resampled.tif"
        save_raster(outroot / "Rasters" / "1km_resampled_display" / anomaly_name,
                    anomaly_display, anomaly_profile)
        make_map(
            anomaly_display, anomaly_profile, boundary_display,
            f"January temperature anomaly — {year} ({display_resolution_m/1000:g} km resampled display)",
            outroot / "Maps" / "Anomalies_1km_resampled" / f"anomaly_{year}_{display_resolution_m}m_resampled.png",
            cmap="RdBu_r", vmin=-anomaly_lim, vmax=anomaly_lim, diverging=True,
            label="Temperature anomaly (°C); bilinear resampling", dpi=cfg["analysis"]["map_dpi"]
        )

    make_map(climatology, profile, boundary, f"January mean 2 m air temperature climatology ({baseline_label_display})",
             outroot / "Maps" / f"january_climatology_{baseline_label}.png", cmap="turbo",
             label="Air temperature (°C)", dpi=cfg["analysis"]["map_dpi"])
    make_map(annual_sd, profile, boundary, "Interannual standard deviation of January temperature (2000–2026)",
             outroot / "Maps" / "january_temperature_standard_deviation.png", cmap="magma",
             label="Standard deviation (°C)", dpi=cfg["analysis"]["map_dpi"])
    make_map(annual_cv, profile, boundary, "Interannual coefficient of variation of January temperature (2000–2026)",
             outroot / "Maps" / "january_temperature_coefficient_of_variation.png", cmap="magma",
             label="CV (%) using absolute temperature (K)", dpi=cfg["analysis"]["map_dpi"])
    make_map(slope*10, profile, boundary, "Sen's slope of January temperature (°C/decade)",
             outroot / "Maps" / "Trends" / "sen_slope_C_per_decade.png", cmap="RdYlBu_r",
             diverging=True, label="Sen's slope (°C/decade)", dpi=cfg["analysis"]["map_dpi"])
    make_map(mk_z, profile, boundary, "Mann–Kendall Z statistic",
             outroot / "Maps" / "Trends" / "mann_kendall_Z.png", cmap="RdBu_r",
             diverging=True, label="MK Z statistic", dpi=cfg["analysis"]["map_dpi"])
    make_map(mk_q if use_fdr else mk_p, profile, boundary,
             "False discovery rate q-values" if use_fdr else "Mann–Kendall p-values",
             outroot / "Maps" / "Trends" / "trend_significance.png", cmap="viridis_r",
             vmin=0, vmax=0.1, label="Adjusted q-value" if use_fdr else "p-value",
             dpi=cfg["analysis"]["map_dpi"])

    def save_resampled_layer(array, stem, title, cmap="viridis", label="", vmin=None, vmax=None,
                             diverging=False, resampling=Resampling.bilinear):
        display, display_prof = resample_for_display(
            array, profile, display_crs, display_resolution_m, resampling=resampling
        )
        tif = outroot / "Rasters" / "1km_resampled_display" / f"{stem}_{display_resolution_m}m_resampled.tif"
        save_raster(tif, display, display_prof)
        with rasterio.open(tif, "r+") as dst:
            dst.update_tags(
                source_resolution="ERA5-Land nominal 11132 m",
                output_pixel_size_m=str(display_resolution_m),
                resampling="nearest" if resampling == Resampling.nearest else "bilinear",
                interpretation="resampled display grid; no additional observational detail"
            )
        make_map(
            display, display_prof, boundary_display,
            f"{title} ({display_resolution_m/1000:g} km resampled display)",
            outroot / "Maps" / "Trends_1km_resampled" / f"{stem}_{display_resolution_m}m_resampled.png",
            cmap=cmap, vmin=vmin, vmax=vmax, diverging=diverging,
            label=label + "; resampled display", dpi=cfg["analysis"]["map_dpi"]
        )
        return display, display_prof

    for stem, array, title, cmap_name, label, lo, hi, diverging in [
        (f"january_climatology_{baseline_label}", climatology,
         f"January temperature climatology ({baseline_label_display})", "turbo", "Temperature (°C)", None, None, False),
        ("january_temperature_standard_deviation", annual_sd,
         "Interannual standard deviation of January temperature (2000–2026)", "magma", "Standard deviation (°C)", None, None, False),
        ("january_temperature_coefficient_of_variation", annual_cv,
         "Interannual coefficient of variation of January temperature (2000–2026)", "magma", "CV (%)", None, None, False),
        ("sen_slope_C_per_decade", slope*10,
         "Sen slope of January temperature", "RdYlBu_r", "Sen slope (°C/decade)", None, None, True),
        ("mann_kendall_Z", mk_z,
         "Mann–Kendall Z statistic", "RdBu_r", "MK Z statistic", None, None, True),
        ("mann_kendall_p_values", mk_p,
         "Mann–Kendall raw p-values", "viridis_r", "Raw p-value", 0, 0.1, False),
        ("mann_kendall_FDR_q_values", mk_q,
         "Benjamini–Hochberg FDR q-values", "viridis_r", "FDR q-value", 0, 0.1, False),
    ]:
        save_resampled_layer(array, stem, title, cmap_name, label, lo, hi, diverging)

    trend_display, trend_display_profile = resample_for_display(
        trend_class, profile, display_crs, display_resolution_m, resampling=Resampling.nearest
    )
    trend_tif = outroot / "Rasters" / "1km_resampled_display" / f"trend_classification_{display_resolution_m}m_resampled.tif"
    save_raster(trend_tif, trend_display, trend_display_profile)
    with rasterio.open(trend_tif, "r+") as dst:
        dst.update_tags(
            source_resolution="ERA5-Land nominal 11132 m", output_pixel_size_m=str(display_resolution_m),
            resampling="nearest", interpretation="resampled display grid; categorical classes; no additional observational detail"
        )
    from matplotlib.colors import ListedColormap, BoundaryNorm
    cmap = ListedColormap(["#2166ac", "#eeeeee", "#b2182b"])
    norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], cmap.N)
    make_categorical_map(
        trend_display, trend_display_profile, boundary_display,
        f"January trend classes ({display_resolution_m/1000:g} km resampled display)",
        outroot / "Maps" / "Trends_1km_resampled" / f"trend_classification_{display_resolution_m}m_resampled.png",
        cmap, norm, [-1, 0, 1], ["Significant cooling", "Not significant", "Significant warming"],
        dpi=cfg["analysis"]["map_dpi"],
    )

    intensity_display, intensity_profile = resample_for_display(
        slope_intensity, profile, display_crs, display_resolution_m, resampling=Resampling.nearest
    )
    intensity_tif = outroot / "Rasters" / "1km_resampled_display" / f"sen_slope_intensity_tertiles_{display_resolution_m}m_resampled.tif"
    save_raster(intensity_tif, intensity_display, intensity_profile)
    with rasterio.open(intensity_tif, "r+") as dst:
        dst.update_tags(
            source_resolution="ERA5-Land nominal 11132 m", output_pixel_size_m=str(display_resolution_m),
            resampling="nearest", interpretation="relative Sen slope tertiles; display-only classes, not significance classes"
        )
    intensity_cmap = ListedColormap(["#fee08b", "#f46d43", "#a50026"])
    intensity_norm = BoundaryNorm([0.5, 1.5, 2.5, 3.5], intensity_cmap.N)
    make_categorical_map(
        intensity_display, intensity_profile, boundary_display,
        f"Relative Sen slope intensity ({display_resolution_m/1000:g} km resampled display)",
        outroot / "Maps" / "Trends_1km_resampled" / f"sen_slope_intensity_tertiles_{display_resolution_m}m_resampled.png",
        intensity_cmap, intensity_norm, [1, 2, 3], ["Lower tertile", "Middle tertile", "Upper tertile"],
        dpi=cfg["analysis"]["map_dpi"],
    )

    # Discrete trend classes.
    cmap = ListedColormap(["#2166ac", "#eeeeee", "#b2182b"])
    norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], cmap.N)
    make_categorical_map(
        trend_class, profile, boundary, "Spatial classification of January temperature trends",
        outroot / "Maps" / "Trends" / "trend_classification.png", cmap, norm, [-1, 0, 1],
        ["Significant cooling", "Not significant", "Significant warming"],
        dpi=cfg["analysis"]["map_dpi"],
    )

    # Time-series plots.
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(regional["year"], regional["regional_mean_C"], marker="o", linewidth=1.5)
    ax.set(title="Regional January mean 2 m air temperature (2000–2026)",
           xlabel="Year", ylabel="Temperature (°C)")
    ax.grid(alpha=0.25); fig.tight_layout()
    regional_temperature_plot = outroot / "Graphs" / "regional_temperature_timeseries_current.png"
    fig.savefig(regional_temperature_plot, dpi=300); plt.close(fig)

    # Descriptive OLS line; Sen's slope and Mann–Kendall remain the inferential trend procedures.
    regional_fit = linregress(regional["year"].to_numpy(dtype=float), regional["regional_mean_C"].to_numpy(dtype=float))
    fitted = regional_fit.intercept + regional_fit.slope * regional["year"].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(regional["year"], regional["regional_mean_C"], marker="o", linewidth=1.5, label="Regional January mean")
    ax.plot(regional["year"], fitted, color="#b2182b", linewidth=2,
            label=f"OLS fit: {regional_fit.slope*10:.3f} °C/decade")
    ax.set(title="Regional January mean temperature and linear trend (2000–2026)",
           xlabel="Year", ylabel="Temperature (°C)")
    ax.legend(); ax.grid(alpha=0.25); fig.tight_layout()
    fig.savefig(outroot / "Graphs" / "regional_temperature_decadal_trend.png", dpi=300); plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(regional["year"], regional["regional_anomaly_C"], marker="o", linewidth=1.5)
    ax.axhline(0, linewidth=1)
    ax.set(title="Regional January temperature anomalies (2000–2026)",
           xlabel="Year", ylabel="Anomaly (°C)")
    ax.grid(alpha=0.25); fig.tight_layout()
    regional_anomaly_plot = outroot / "Graphs" / "regional_anomaly_timeseries_current.png"
    fig.savefig(regional_anomaly_plot, dpi=300); plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 6))
    for state in EXPECTED_STATES:
        subset = state_ts[state_ts["state"] == state]
        ax.plot(subset["year"], subset["mean_temperature_C"], marker=".", label=state)
    ax.set(title="State-level January mean 2 m air temperature (2000–2026)",
           xlabel="Year", ylabel="Temperature (°C)")
    ax.legend(ncol=3); ax.grid(alpha=0.25); fig.tight_layout()
    state_temperature_plot = outroot / "Graphs" / "state_temperature_timeseries_current.png"
    fig.savefig(state_temperature_plot, dpi=300); plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 6))
    for state in EXPECTED_STATES:
        subset = state_ts[state_ts["state"] == state]
        ax.plot(subset["year"], subset["mean_anomaly_C"], marker=".", label=state)
    ax.axhline(0, linewidth=1)
    ax.set(title="State-level January temperature anomalies (2000–2026)",
           xlabel="Year", ylabel="Anomaly (°C)")
    ax.legend(ncol=3); ax.grid(alpha=0.25); fig.tight_layout()
    state_anomaly_plot = outroot / "Graphs" / "state_anomaly_timeseries_current.png"
    fig.savefig(state_anomaly_plot, dpi=300); plt.close(fig)

    # Additional CSV tables and Excel workbook.
    summary = pd.DataFrame([{
        "period": f"{main_years[0]}-{main_years[-1]}",
        "regional_mean_of_annual_means_C": float(regional["regional_mean_C"].mean()),
        "minimum_annual_regional_mean_C": float(regional["regional_mean_C"].min()),
        "maximum_annual_regional_mean_C": float(regional["regional_mean_C"].max()),
        "sd_of_annual_regional_means_C": float(regional["regional_mean_C"].std(ddof=1)),
        "cv_absolute_temperature_pct": regional_cv_abs,
        "mean_annual_anomaly_C": float(regional["regional_anomaly_C"].mean()),
        "trend_test_cells": len(p_values),
        "autocorrelation_adjusted_cells": int(autocorr_adjusted.sum()),
        "significance_basis": "FDR q-value" if use_fdr else "raw MK p-value"
    }])
    summary.to_csv(outroot / "Tables" / "regional_summary.csv", index=False)
    with pd.ExcelWriter(outroot / "Tables" / "January_Temperature_Results.xlsx", engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Regional summary", index=False)
        regional.to_excel(writer, sheet_name="Regional annual", index=False)
        annual_min_max_avg.to_excel(writer, sheet_name="Annual min max average", index=False)
        state_ts.to_excel(writer, sheet_name="State annual", index=False)
        state_descriptive_df.to_excel(writer, sheet_name="State descriptive", index=False)
        state_trends_df.to_excel(writer, sheet_name="State trends", index=False)
        trend_area_df.to_excel(writer, sheet_name="Regional trend area", index=False)

    # Results text and Word manuscript. Text is computed only from actual outputs.
    warmest = regional.loc[regional["regional_mean_C"].idxmax()]
    coolest = regional.loc[regional["regional_mean_C"].idxmin()]
    mean_trend = float(np.nanmedian(slope[inside]) * 10)
    results_text = (
        f"Across the regional annual series, the highest January mean temperature occurred in "
        f"{int(warmest['year'])} ({warmest['regional_mean_C']:.2f} °C), while the lowest occurred "
        f"in {int(coolest['year'])} ({coolest['regional_mean_C']:.2f} °C). The mean of the annual "
        f"regional January temperatures was {regional['regional_mean_C'].mean():.2f} °C, with a "
        f"standard deviation of {regional['regional_mean_C'].std(ddof=1):.2f} °C; the coefficient "
        f"of variation calculated with absolute temperature in kelvin was {regional_cv_abs:.2f}%. "
        f"The median "
        f"pixel-level Sen's slope was {mean_trend:.3f} °C per decade. These summaries describe "
        f"the processed ERA5-Land estimates and should be interpreted alongside the spatial maps, "
        f"state-level results and significance outputs."
    )
    (outroot / "Manuscript" / "auto_results_text.txt").write_text(results_text, encoding="utf-8")

    template = root / cfg["project"].get("manuscript_template", "")
    if template.exists():
        doc = Document(str(template))
    else:
        doc = Document()
        doc.add_heading("Spatiotemporal Analysis of January 2-metre Air Temperature Variability and Long-Term Trends in Southwestern Nigeria (2000–2026)", 0)
    # Correct the one inconsistent baseline-period statement in the source text.
    for paragraph in doc.paragraphs:
        text = paragraph.text
        fixed = text.replace("1991–2026 January climatological reference period", f"{baseline_label_display} January climatological reference period")
        fixed = fixed.replace("1991-2026 January climatological reference period", f"{base_years[0]}-{base_years[-1]} January climatological reference period")
        if fixed != text:
            paragraph.text = fixed
    doc.add_heading("2.10. Implemented statistical procedures", level=2)
    doc.add_paragraph(
        f"The regional and state time series were summarized from the annual January rasters. "
        f"Grid-cell trends were assessed using the Mann–Kendall test; lag-one serial dependence "
        f"was screened with the Ljung–Box test at α = {ac_alpha:.2f}, and the Hamed–Rao "
        f"modified Mann–Kendall test was used where that screen was significant. Sen slopes "
        f"were estimated from pairwise annual values and reported in °C per decade. "
        f"For spatial trend classification, the Benjamini–Hochberg false-discovery-rate "
        f"adjustment was {'applied' if use_fdr else 'not applied'} at α = "
        f"{float(cfg['analysis']['fdr_alpha'] if use_fdr else alpha):.2f}; positive and negative "
        f"significant slopes were classified as warming and cooling, respectively. Regional "
        f"Regional and state zonal summaries and trend-area proportions use cell-polygon "
        f"intersection areas calculated in the configured equal-area projection, including "
        f"fractional coverage at state and regional boundaries. The coefficient of variation "
        f"is calculated as the standard deviation of annual January means divided by mean "
        f"absolute temperature in kelvin, multiplied by 100, avoiding a Celsius-zero denominator. "
        f"The analysis boundary is the supplied Southwest.shp shapefile, using its configured "
        f"state-name field; no boundary source is inferred."
    )
    doc.add_page_break()
    doc.add_heading("3. Results", level=1)
    doc.add_heading("3.1. Regional January Temperature Variation", level=2)
    doc.add_paragraph(results_text)
    doc.add_picture(str(regional_temperature_plot), width=Inches(6.0))
    doc.add_heading("3.2. January Temperature Anomalies", level=2)
    doc.add_paragraph(
        f"Relative to the {baseline_label_display} January climatological baseline, the mean regional anomaly "
        f"across {main_years[0]}–{main_years[-1]} was "
        f"{regional['regional_anomaly_C'].mean():.2f} °C. Positive values represent warmer-than-"
        f"baseline January conditions, whereas negative values represent cooler-than-baseline conditions."
    )
    doc.add_picture(str(regional_anomaly_plot), width=Inches(6.0))
    doc.add_heading("3.3. Spatial Temperature Trends", level=2)
    doc.add_paragraph(
        f"The pixel-wise trend analysis evaluated {len(p_values)} grid cells. "
        f"{int(autocorr_adjusted.sum())} cells were assigned the Hamed–Rao modified Mann–Kendall "
        f"test following the lag-one autocorrelation diagnostic. The final trend classification "
        f"used {'FDR-adjusted q-values' if use_fdr else 'unadjusted Mann–Kendall p-values'} "
        f"at the configured significance threshold."
    )
    doc.add_picture(str(outroot / "Maps" / "Trends" / "sen_slope_C_per_decade.png"), width=Inches(6.0))
    doc.add_picture(str(outroot / "Maps" / "Trends" / "trend_classification.png"), width=Inches(6.0))
    doc.add_heading("3.4. State-Level Trend Comparison", level=2)
    table = doc.add_table(rows=1, cols=5)
    try:
        table.style = "Light Shading Accent 1"
    except KeyError:
        # Some supplied manuscript templates do not include the built-in table style.
        pass
    headers = ["State", "Mean Sen slope (°C/decade)", "Warming (%)", "Cooling (%)", "Not significant (%)"]
    for cell, text in zip(table.rows[0].cells, headers):
        cell.text = text
    for _, r in state_trends_df.iterrows():
        cells = table.add_row().cells
        vals = [r["state"], f"{r['mean_sen_slope_C_per_decade']:.3f}",
                f"{r['significant_warming_percent']:.1f}",
                f"{r['significant_cooling_percent']:.1f}",
                f"{r['non_significant_percent']:.1f}"]
        for cell, val in zip(cells, vals): cell.text = str(val)
    doc.add_paragraph(
        "The accompanying Excel workbook and CSV tables contain the complete regional and "
        "state-level numerical summaries. All numerical findings in this generated section are "
        "calculated from the downloaded raster data; they should be reviewed before journal submission."
    )
    doc.save(outroot / "Manuscript" / "AirTemperature_Manuscript_with_Results.docx")
    logging.info("Analysis completed. Results and manuscript saved to %s", outroot)
