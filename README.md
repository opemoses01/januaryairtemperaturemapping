# January Air Temperature Analysis for Southwestern Nigeria

Python workflow for extracting January 2 m air temperature from ERA5-Land through Google Earth Engine, downloading annual raster stacks, and generating trend analyses, maps, tables, and an optional manuscript draft.

## Setup

1. Install Python 3.11 or later and create a virtual environment.
2. Install dependencies with `python -m pip install -r 00_CONFIGURATION/requirements.txt`.
3. Edit `config.yaml`: set `project.gee_cloud_project` and, if needed for authentication, `project.oauth_client_id`. The repository configuration uses placeholders. Never commit OAuth client secrets, refresh tokens, or credential files.
4. Run `python RUN_ALL.py`.

Earth Engine access requires an enabled Google Cloud project and valid Earth Engine authorization. The workflow prompts for an OAuth client secret only when needed and saves credentials outside this repository.

## Analysis settings

Defaults use ERA5-Land hourly `temperature_2m`, January 2000–2026, with 1991–2000 as the baseline. The six-state study boundary is in `01_BOUNDARY/studyarea/`. Update the study period and boundary settings in `config.yaml` as needed.

Use `python RUN_ALL.py --skip-gee` to run analysis with previously downloaded raster stacks, or `python RUN_ALL.py --skip-analysis` to perform extraction only.

## Outputs and manuscript

Input raster stacks and generated results are written to `03_DATA_DOWNLOAD/` and `08_OUTPUTS/`. These generated files are excluded from version control by `.gitignore`. The manuscript template is not included; set `project.manuscript_template` to a local Word document path if you want the optional manuscript output.

The 1 km maps are resampled display products and add no observational detail beyond the native ERA5-Land grid. Review generated statistics, maps, and manuscript text before publication, and cite ERA5-Land appropriately.
