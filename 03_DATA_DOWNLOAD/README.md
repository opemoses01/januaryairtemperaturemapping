# Downloaded input temperature data

The GEE extraction stage writes its initial analysis rasters directly into this folder:

- `january_temperature_main.tif`: one band per year, January 2000–2026.
- `january_temperature_baseline.tif`: one band per year, January 1991–2000 as set in `config.yaml`.
- `extraction_metadata.json`: dataset, band, years, units, aggregation, scale and CRS.

Each band is an annual January mean in degrees Celsius, computed in Earth Engine from the ERA5-Land hourly 2 m temperature collection. These are the input datasets used by the Python analysis stage. They are downloaded here automatically when `RUN_ALL.py` completes the GEE extraction stage; authentication must succeed first.
