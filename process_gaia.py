#!/usr/bin/env python3
"""
Gaia DR3 Epoch Photometry Analyzer
====================================
Identifies astronomical objects whose BP or RP flux changed by more than 100%
over the observation period, using the first 20 files from the Gaia DR3 epoch
photometry archive.

Algorithm (per source_id):
  1. Extract bp_flux and rp_flux arrays from each file row.
  2. Discard missing, null, NaN, or infinite flux values.
  3. Compute min and max of the remaining valid values for each band.
  4. Percentage change = ((max_flux - min_flux) / min_flux) * 100
     (only when min_flux > 0, to avoid division-by-zero / sign ambiguity)
  5. Take the larger percentage_change across BP and RP.
  6. Emit the source if percentage_change > 100 %.

Output CSV columns:
  source_id, bp_min_flux, bp_max_flux, rp_min_flux, rp_max_flux,
  percentage_change

Usage:
  python process_gaia.py [--output results.csv] [--data-dir .data/in]

Comments / feedback on the challenge:
  - The Gaia epoch photometry files are large gzip-compressed CSVs where each
    row holds a single source and the flux measurements are stored as
    space-separated arrays inside a single CSV cell (e.g., "1.23 4.56 7.89").
  - Python (without heavy scientific dependencies) makes the min/max computation over
    variable-length arrays very clean and efficient.
  - A streaming approach (row-by-row parsing) is used so that even very large
    files can be processed without loading everything into RAM at once.
  - Writing code in Python as recommended by the bonus criteria was
    straightforward for this scientific data-processing task.
"""

import argparse
import csv
import gzip
import io
import logging
import math
import os
import re
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
from typing import BinaryIO, Dict, Generator, List, Optional, Tuple, Union

import requests

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Gaia DR3 epoch photometry CDN base URL
CDN_BASE = "https://cdn.gea.esac.esa.int"
CDN_LIST_URL = f"{CDN_BASE}/?prefix=Gaia/gdr3/Photometry/epoch_photometry/"

# Number of files to process (as specified in the challenge benchmark)
NUM_FILES = 20

# Percentage-change threshold
THRESHOLD = 100.0

# CSV output columns
OUTPUT_FIELDS = [
    "source_id",
    "bp_min_flux",
    "bp_max_flux",
    "rp_min_flux",
    "rp_max_flux",
    "percentage_change",
]

# HTTP request settings
# Gaia epoch photometry files can be several hundred MB each.
# 300 s gives ~1 MB/s on a slow connection, which is conservative but safe.
REQUEST_TIMEOUT = 300  # seconds per download
MAX_RETRIES = 3
RETRY_DELAY = 5  # seconds between retries
DOWNLOAD_CHUNK_SIZE = 1024 * 1024  # 1 MiB

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# File listing
# ---------------------------------------------------------------------------


def fetch_file_listing() -> List[str]:
    """
    Query the Gaia CDN (Google Cloud Storage bucket) for the list of epoch
    photometry CSV files and return them sorted alphabetically.

    Returns a list of full download URLs for the first NUM_FILES files.
    """
    logger.info("Fetching file listing from %s", CDN_LIST_URL)

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = requests.get(CDN_LIST_URL, timeout=60)
            resp.raise_for_status()
            break
        except requests.RequestException as exc:
            logger.warning("Attempt %d/%d failed: %s", attempt, MAX_RETRIES, exc)
            if attempt == MAX_RETRIES:
                raise
            time.sleep(RETRY_DELAY)

    # Parse S3-compatible XML listing
    try:
        root = ET.fromstring(resp.text)
    except ET.ParseError as exc:
        raise ValueError(f"Failed to parse CDN XML listing: {exc}") from exc

    # Namespace can vary; strip it for robustness
    ns_match = re.match(r"\{(.+?)\}", root.tag)
    ns = f"{{{ns_match.group(1)}}}" if ns_match else ""

    files: List[str] = []
    for contents in root.findall(f"{ns}Contents"):
        key_el = contents.find(f"{ns}Key")
        if key_el is None or not key_el.text:
            continue
        key = key_el.text
        # Accept only epoch photometry CSV gzip files
        filename = key.rsplit("/", 1)[-1]
        if filename.startswith("EpochPhotometry_") and filename.endswith(".csv.gz"):
            files.append(key)

    files.sort()
    logger.info("Found %d CSV files in the archive.", len(files))

    selected = files[:NUM_FILES]
    urls = [f"{CDN_BASE}/{key}" for key in selected]
    for url in urls:
        logger.info("  Selected: %s", url)

    return urls


# ---------------------------------------------------------------------------
# File download
# ---------------------------------------------------------------------------


def download_file(url: str) -> str:
    """
    Download a file with retries, streaming it to a temporary file.

    Returns the temporary file path. The caller is responsible for removing it,
    including if downstream processing later fails after this function returns.
    """
    for attempt in range(1, MAX_RETRIES + 1):
        temp_path: Optional[str] = None
        try:
            logger.info("[%d/%d] Downloading %s", attempt, MAX_RETRIES, url)
            resp = requests.get(url, timeout=REQUEST_TIMEOUT, stream=True)
            resp.raise_for_status()
            total_bytes = 0
            with tempfile.NamedTemporaryFile(delete=False, suffix=".csv.gz") as fh:
                temp_path = fh.name
                for chunk in resp.iter_content(chunk_size=DOWNLOAD_CHUNK_SIZE):
                    if chunk:
                        fh.write(chunk)
                        total_bytes += len(chunk)
            logger.info("  Downloaded %.1f MB", total_bytes / 1_048_576)
            return temp_path
        except requests.RequestException as exc:
            if temp_path and os.path.exists(temp_path):
                os.remove(temp_path)
            logger.warning("Download attempt %d failed: %s", attempt, exc)
            if attempt == MAX_RETRIES:
                raise
            time.sleep(RETRY_DELAY)
    raise RuntimeError(f"All {MAX_RETRIES} download attempts failed for {url}")


# ---------------------------------------------------------------------------
# Flux array parsing
# ---------------------------------------------------------------------------


def parse_flux_array(value: str) -> List[float]:
    """
    Parse a CSV cell that contains a variable-length array of flux values.

    Gaia epoch photometry CSV files encode arrays as space-separated floats
    within a single cell, optionally wrapped in square brackets.  Example:
      "1.23456e+03 -4.56e+02 7.89e+03"
      "[1.23456e+03 -4.56e+02 7.89e+03]"

    Returns a list of finite, non-NaN floats.  Invalid tokens are silently
    skipped as per the problem specification.
    """
    if not value:
        return []

    # Remove brackets if present
    stripped = value.strip().lstrip("[").rstrip("]")

    valid: List[float] = []
    for token in stripped.split():
        token = token.strip().rstrip(",")
        if not token or token.lower() in ("nan", "null", "none", "inf", "-inf"):
            continue
        try:
            f = float(token)
            if math.isfinite(f):
                valid.append(f)
        except ValueError:
            pass

    return valid


# ---------------------------------------------------------------------------
# Per-band statistics
# ---------------------------------------------------------------------------


def band_stats(fluxes: List[float]) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """
    Given a list of valid flux values, return (min_flux, max_flux, pct_change).

    percentage_change = ((max - min) / min) * 100
    Only computed when min_flux > 0 (avoids division-by-zero and sign issues).
    Requires at least 2 values; returns (None, None, None) otherwise.
    """
    if len(fluxes) < 2:
        return None, None, None

    mn = min(fluxes)
    mx = max(fluxes)

    if mn <= 0:
        # Cannot compute a meaningful positive ratio when min ≤ 0
        return float(mn), float(mx), None

    pct = ((mx - mn) / mn) * 100.0
    return float(mn), float(mx), pct


# ---------------------------------------------------------------------------
# Row processing
# ---------------------------------------------------------------------------


def process_row(row: Dict[str, str]) -> Optional[Dict]:
    """
    Process one CSV row (one source_id) and return a result dict if the source
    meets the > THRESHOLD % variability criterion, otherwise return None.
    """
    source_id = row.get("source_id", "").strip()
    if not source_id:
        return None

    bp_fluxes = parse_flux_array(row.get("bp_flux", ""))
    rp_fluxes = parse_flux_array(row.get("rp_flux", ""))

    bp_min, bp_max, bp_pct = band_stats(bp_fluxes)
    rp_min, rp_max, rp_pct = band_stats(rp_fluxes)

    # Collect valid percentage changes from both bands
    pcts = [p for p in [bp_pct, rp_pct] if p is not None]
    if not pcts:
        return None

    percentage_change = max(pcts)
    if percentage_change <= THRESHOLD:
        return None

    return {
        "source_id": source_id,
        "bp_min_flux": bp_min,
        "bp_max_flux": bp_max,
        "rp_min_flux": rp_min,
        "rp_max_flux": rp_max,
        "percentage_change": round(percentage_change, 6),
    }


# ---------------------------------------------------------------------------
# File processing
# ---------------------------------------------------------------------------


def iter_csv_rows(data: Union[bytes, BinaryIO]) -> Generator[Dict[str, str], None, None]:
    """
    Decompress a gzip-compressed CSV bytes object or binary file and yield rows
    as dicts.
    """
    if isinstance(data, bytes):
        # Wrap in-memory gzip bytes in a file-like object for gzip.open().
        source = io.BytesIO(data)
    else:
        source = data

    with gzip.open(source, "rt", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            yield row


def process_file_data(data: Union[bytes, BinaryIO], filename: str) -> List[Dict]:
    """
    Process the contents of a single epoch photometry CSV.gz file from raw
    bytes or a binary file object.

    When passing a file object, it must be opened in binary mode and positioned
    at the start of the gzip-compressed data.
    Returns a list of result dicts for qualifying sources.
    """
    results: List[Dict] = []
    row_count = 0
    match_count = 0

    for row in iter_csv_rows(data):
        row_count += 1
        result = process_row(row)
        if result is not None:
            results.append(result)
            match_count += 1

    logger.info(
        "  %s: %d sources processed, %d with >%.0f%% variability",
        filename, row_count, match_count, THRESHOLD,
    )
    return results


# ---------------------------------------------------------------------------
# Local file mode
# ---------------------------------------------------------------------------


def find_local_files(data_dir: str) -> List[str]:
    """
    Return sorted paths to CSV.gz epoch photometry files in data_dir.
    Returns the first NUM_FILES found (alphabetical order).
    """
    if not os.path.isdir(data_dir):
        return []

    files = sorted(
        os.path.join(data_dir, f)
        for f in os.listdir(data_dir)
        if f.endswith(".csv.gz") and "EpochPhotometry" in f
    )
    selected = files[:NUM_FILES]
    logger.info("Found %d local files in '%s'.", len(selected), data_dir)
    for p in selected:
        logger.info("  Local file: %s", p)
    return selected


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def main(args: argparse.Namespace) -> None:
    output_path = args.output
    data_dir = args.data_dir

    all_results: List[Dict] = []
    t_start = time.time()

    # ---- Determine source of files ----------------------------------------
    local_files = find_local_files(data_dir) if data_dir else []

    if local_files:
        logger.info("Using %d local file(s) from '%s'.", len(local_files), data_dir)
        for path in local_files:
            fname = os.path.basename(path)
            logger.info("Processing local file: %s", fname)
            with open(path, "rb") as fh:
                data = fh.read()
            results = process_file_data(data, fname)
            all_results.extend(results)

    else:
        # Download from Gaia CDN
        logger.info("No local files found. Downloading from Gaia CDN.")
        urls = fetch_file_listing()

        for i, url in enumerate(urls, start=1):
            fname = url.rsplit("/", 1)[-1]
            logger.info("--- File %d/%d: %s ---", i, len(urls), fname)
            temp_path = download_file(url)
            try:
                with open(temp_path, "rb") as fh:
                    results = process_file_data(fh, fname)
                all_results.extend(results)
            finally:
                if os.path.exists(temp_path):
                    os.remove(temp_path)

    # ---- Write output -------------------------------------------------------
    logger.info("Writing %d results to '%s'.", len(all_results), output_path)
    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(all_results)

    elapsed = time.time() - t_start
    logger.info(
        "Done. %d qualifying sources written to '%s' in %.1f s.",
        len(all_results), output_path, elapsed,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Gaia DR3 Epoch Photometry Analyzer — finds sources whose "
            "BP or RP flux changed by >100%% across observations."
        )
    )
    parser.add_argument(
        "--output", "-o",
        default="results.csv",
        metavar="FILE",
        help="Path to the output CSV file (default: results.csv)",
    )
    parser.add_argument(
        "--data-dir", "-d",
        default=".data/in",
        metavar="DIR",
        help=(
            "Directory containing local EpochPhotometry_*.csv.gz files. "
            "If files are found here, the CDN download is skipped. "
            "(default: .data/in)"
        ),
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose (DEBUG-level) logging",
    )
    parsed = parser.parse_args()

    if parsed.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    main(parsed)
