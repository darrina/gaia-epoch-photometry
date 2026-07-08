# Gaia DR3 Epoch Photometry Variability Analyzer

Identifies astronomical objects whose **BP or RP flux changed by more than 100 %** over the Gaia observation period, using the first 20 files from the [Gaia DR3 epoch photometry archive](https://cdn.gea.esac.esa.int/?prefix=Gaia/gdr3/Photometry/epoch_photometry/).

---

## How it works

For every `source_id` in the input files the script:

1. Extracts the `bp_flux` and `rp_flux` observation arrays.
2. Discards missing, null, NaN, and infinite values.
3. Computes **min** and **max** of the remaining valid fluxes for each band (BP and RP).
4. Calculates the percentage change for each band:

   ```
   percentage_change = ((max_flux − min_flux) / min_flux) × 100
   ```

   (computed only when `min_flux > 0`, to avoid division-by-zero / sign ambiguity)

5. Takes the **larger** of the two band percentage changes as the final `percentage_change`.
6. Emits the source when `percentage_change > 100 %`.

### Output columns

| Column | Description |
|---|---|
| `source_id` | Gaia DR3 source identifier |
| `bp_min_flux` | Minimum valid BP flux across all observations |
| `bp_max_flux` | Maximum valid BP flux across all observations |
| `rp_min_flux` | Minimum valid RP flux across all observations |
| `rp_max_flux` | Maximum valid RP flux across all observations |
| `percentage_change` | Largest variability (%) across BP and RP bands |

---

## Installation

Python **3.8 +** required (uses `f-strings`, `typing` generics, and `argparse` features present since 3.8).

```bash
pip install -r requirements.txt
```

Dependencies are minimal — only the standard `requests` library is needed (plus `pytest` for tests).

---

## Usage

### Download and process from the Gaia CDN (default)

```bash
python process_gaia.py
```

The script fetches the XML file listing from the Gaia CDN, selects the first 20 `EpochPhotometry_*.csv.gz` files (sorted alphabetically), downloads them, and writes qualifying sources to `results.csv`.

### Process local files (e.g. the challenge template `.data/in` directory)

```bash
python process_gaia.py --data-dir .data/in
```

If the `--data-dir` directory contains `EpochPhotometry_*.csv.gz` files the CDN download is skipped entirely.

### All options

```
usage: process_gaia.py [-h] [--output FILE] [--data-dir DIR] [--verbose]

options:
  -h, --help            show this help message and exit
  --output FILE, -o FILE
                        Path to the output CSV file (default: results.csv)
  --data-dir DIR, -d DIR
                        Directory with local EpochPhotometry_*.csv.gz files.
                        Falls back to CDN download when empty or absent.
                        (default: .data/in)
  --verbose, -v         Enable verbose (DEBUG-level) logging
```

---

## Running the tests

```bash
pip install pytest
python -m pytest tests/ -v
```

---

## Data source

- **Archive**: [Gaia DR3 Epoch Photometry CDN](https://cdn.gea.esac.esa.int/?prefix=Gaia/gdr3/Photometry/epoch_photometry/)
- **Column descriptions**: [Gaia DR3 data model — epoch_photometry](https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_photometry/ssec_dm_epoch_photometry.html)
