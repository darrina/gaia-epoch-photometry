"""
Tests for process_gaia.py

Tests cover:
  - parse_flux_array  : tokenising flux strings with various encodings
  - band_stats        : per-band min/max/percentage_change computation
  - process_row       : end-to-end per-source filtering logic
"""

import math
import gzip
import io
import csv
import sys
import os

# Make sure the parent directory is on the path so we can import process_gaia
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
import process_gaia

from process_gaia import (
    DOWNLOAD_CHUNK_SIZE,
    parse_flux_array,
    band_stats,
    download_file,
    process_row,
    iter_csv_rows,
    process_file_data,
)


# ---------------------------------------------------------------------------
# parse_flux_array
# ---------------------------------------------------------------------------


class TestParseFluxArray:
    def test_empty_string(self):
        assert parse_flux_array("") == []

    def test_none_like_string(self):
        assert parse_flux_array("   ") == []

    def test_space_separated(self):
        result = parse_flux_array("1.0 2.0 3.0")
        assert result == [1.0, 2.0, 3.0]

    def test_bracket_wrapped(self):
        result = parse_flux_array("[1.0 2.0 3.0]")
        assert result == [1.0, 2.0, 3.0]

    def test_scientific_notation(self):
        result = parse_flux_array("1.23e3 4.56e-2")
        assert math.isclose(result[0], 1230.0)
        assert math.isclose(result[1], 0.0456)

    def test_nan_skipped(self):
        result = parse_flux_array("1.0 nan 3.0")
        assert result == [1.0, 3.0]

    def test_null_skipped(self):
        result = parse_flux_array("1.0 null 3.0")
        assert result == [1.0, 3.0]

    def test_inf_skipped(self):
        result = parse_flux_array("1.0 inf -inf 3.0")
        assert result == [1.0, 3.0]

    def test_negative_values_kept(self):
        result = parse_flux_array("-5.0 2.0 3.0")
        assert -5.0 in result
        assert len(result) == 3

    def test_single_value(self):
        result = parse_flux_array("42.0")
        assert result == [42.0]

    def test_comma_stripped(self):
        # Some encodings may include trailing commas after each token
        result = parse_flux_array("1.0, 2.0, 3.0")
        assert result == [1.0, 2.0, 3.0]


# ---------------------------------------------------------------------------
# band_stats
# ---------------------------------------------------------------------------


class TestBandStats:
    def test_empty_list_returns_nones(self):
        assert band_stats([]) == (None, None, None)

    def test_single_value_returns_nones(self):
        assert band_stats([5.0]) == (None, None, None)

    def test_basic_percentage_change(self):
        # min=1, max=3 → pct = ((3-1)/1)*100 = 200%
        mn, mx, pct = band_stats([1.0, 3.0])
        assert mn == 1.0
        assert mx == 3.0
        assert math.isclose(pct, 200.0)

    def test_no_change(self):
        mn, mx, pct = band_stats([5.0, 5.0])
        assert math.isclose(pct, 0.0)

    def test_exactly_100_percent(self):
        # min=2, max=4 → 100%
        _, _, pct = band_stats([2.0, 4.0])
        assert math.isclose(pct, 100.0)

    def test_min_zero_pct_is_none(self):
        mn, mx, pct = band_stats([0.0, 5.0])
        assert mn == 0.0
        assert mx == 5.0
        assert pct is None

    def test_min_negative_pct_is_none(self):
        mn, mx, pct = band_stats([-3.0, 5.0])
        assert mn == -3.0
        assert mx == 5.0
        assert pct is None

    def test_all_positive_many_values(self):
        fluxes = [10.0, 20.0, 5.0, 15.0]  # min=5, max=20
        mn, mx, pct = band_stats(fluxes)
        assert mn == 5.0
        assert mx == 20.0
        assert math.isclose(pct, ((20 - 5) / 5) * 100)


# ---------------------------------------------------------------------------
# process_row
# ---------------------------------------------------------------------------


def _make_row(source_id="1234567890", bp="", rp=""):
    return {"source_id": source_id, "bp_flux": bp, "rp_flux": rp}


class TestProcessRow:
    def test_no_source_id_returns_none(self):
        assert process_row(_make_row(source_id="")) is None

    def test_no_flux_data_returns_none(self):
        assert process_row(_make_row(bp="", rp="")) is None

    def test_single_bp_value_no_result(self):
        # Only 1 value → cannot compute variability
        assert process_row(_make_row(bp="5.0", rp="")) is None

    def test_below_threshold_returns_none(self):
        # 50% change (min=2, max=3)
        row = _make_row(bp="2.0 3.0", rp="")
        assert process_row(row) is None

    def test_exactly_100_percent_returns_none(self):
        # Exactly 100% is NOT > 100%, so should be excluded
        row = _make_row(bp="1.0 2.0", rp="")
        assert process_row(row) is None

    def test_above_threshold_bp_only(self):
        # 200% change (min=1, max=3)
        row = _make_row(bp="1.0 3.0", rp="")
        result = process_row(row)
        assert result is not None
        assert result["source_id"] == "1234567890"
        assert math.isclose(result["bp_min_flux"], 1.0)
        assert math.isclose(result["bp_max_flux"], 3.0)
        assert result["rp_min_flux"] is None
        assert result["rp_max_flux"] is None
        assert math.isclose(result["percentage_change"], 200.0)

    def test_above_threshold_rp_only(self):
        # RP changes 300% (min=1, max=4), BP single value → no bp stats
        row = _make_row(bp="5.0", rp="1.0 4.0")
        result = process_row(row)
        assert result is not None
        assert math.isclose(result["rp_min_flux"], 1.0)
        assert math.isclose(result["rp_max_flux"], 4.0)
        assert math.isclose(result["percentage_change"], 300.0)

    def test_max_of_bp_and_rp_chosen(self):
        # BP: min=1, max=3 → 200%; RP: min=1, max=6 → 500%; result=500%
        row = _make_row(bp="1.0 3.0", rp="1.0 6.0")
        result = process_row(row)
        assert result is not None
        assert math.isclose(result["percentage_change"], 500.0)

    def test_nan_values_ignored(self):
        # Effective fluxes after NaN removal: [1.0, 3.0] → 200%
        row = _make_row(bp="1.0 nan 3.0", rp="")
        result = process_row(row)
        assert result is not None
        assert math.isclose(result["percentage_change"], 200.0)

    def test_all_nan_bp_rp_usable(self):
        # All BP values are NaN, RP has 500% change
        row = _make_row(bp="nan nan nan", rp="1.0 6.0")
        result = process_row(row)
        assert result is not None
        assert math.isclose(result["percentage_change"], 500.0)

    def test_min_zero_falls_back_to_rp(self):
        # BP min=0 → no pct; RP min=1, max=4 → 300% → result=300%
        row = _make_row(bp="0.0 5.0", rp="1.0 4.0")
        result = process_row(row)
        assert result is not None
        assert math.isclose(result["percentage_change"], 300.0)

    def test_both_bands_min_zero_returns_none(self):
        row = _make_row(bp="0.0 5.0", rp="0.0 8.0")
        assert process_row(row) is None


# ---------------------------------------------------------------------------
# iter_csv_rows / process_file_data (integration)
# ---------------------------------------------------------------------------


def _make_csv_gz(rows, fieldnames):
    """Helper: build a gzip-compressed CSV bytes object from a list of row dicts."""
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
    data = buf.getvalue().encode("utf-8")
    gz_buf = io.BytesIO()
    with gzip.open(gz_buf, "wb") as fh:
        fh.write(data)
    return gz_buf.getvalue()


class TestIterCsvRows:
    def test_basic_iteration(self):
        gz = _make_csv_gz(
            [{"source_id": "1", "bp_flux": "1.0 2.0", "rp_flux": ""}],
            ["source_id", "bp_flux", "rp_flux"],
        )
        rows = list(iter_csv_rows(gz))
        assert len(rows) == 1
        assert rows[0]["source_id"] == "1"

    def test_accepts_binary_file_object(self, tmp_path):
        gz = _make_csv_gz(
            [{"source_id": "1", "bp_flux": "1.0 2.0", "rp_flux": ""}],
            ["source_id", "bp_flux", "rp_flux"],
        )
        path = tmp_path / "test.csv.gz"
        path.write_bytes(gz)

        with path.open("rb") as fh:
            rows = list(iter_csv_rows(fh))

        assert len(rows) == 1
        assert rows[0]["source_id"] == "1"


class TestDownloadFile:
    def test_streams_response_to_temp_file(self, monkeypatch):
        chunks = [b"abc", b"", b"def", b"ghi"]

        class MockResponse:
            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size):
                assert chunk_size == DOWNLOAD_CHUNK_SIZE
                yield from chunks

            @property
            def content(self):
                raise AssertionError("download_file should not access resp.content")

        monkeypatch.setattr(
            "process_gaia.requests.get",
            lambda *args, **kwargs: MockResponse(),
        )

        temp_path = download_file("https://example.com/test.csv.gz")
        try:
            with open(temp_path, "rb") as fh:
                assert fh.read() == b"abcdefghi"
        finally:
            os.remove(temp_path)

    def test_skips_empty_chunks_when_streaming(self, monkeypatch):
        chunks = [b"abc", b"", b"def"]

        class MockResponse:
            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size):
                assert chunk_size == DOWNLOAD_CHUNK_SIZE
                yield from chunks

        monkeypatch.setattr(
            "process_gaia.requests.get",
            lambda *args, **kwargs: MockResponse(),
        )

        temp_path = download_file("https://example.com/test.csv.gz")
        try:
            with open(temp_path, "rb") as fh:
                assert fh.read() == b"abcdef"
        finally:
            os.remove(temp_path)

    def test_removes_temp_file_on_streaming_error(self, monkeypatch, tmp_path):
        temp_path = tmp_path / "download.csv.gz"

        class FakeTempFile:
            def __init__(self, path):
                self.name = str(path)
                self._fh = None

            def __enter__(self):
                self._fh = open(self.name, "wb")
                return self

            def write(self, data):
                return self._fh.write(data)

            def __exit__(self, exc_type, exc, tb):
                self._fh.close()

        class MockResponse:
            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size):
                assert chunk_size == DOWNLOAD_CHUNK_SIZE
                yield b"abc"
                raise process_gaia.requests.RequestException("stream interrupted")

        monkeypatch.setattr(
            "process_gaia.requests.get",
            lambda *args, **kwargs: MockResponse(),
        )
        monkeypatch.setattr(
            "process_gaia.tempfile.NamedTemporaryFile",
            lambda *args, **kwargs: FakeTempFile(temp_path),
        )

        with pytest.raises(process_gaia.requests.RequestException, match="stream interrupted"):
            download_file("https://example.com/test.csv.gz")

        assert not temp_path.exists()

    def test_removes_temp_file_on_http_error(self, monkeypatch):
        class MockResponse:
            def raise_for_status(self):
                raise process_gaia.requests.HTTPError("bad response")

        monkeypatch.setattr(
            "process_gaia.requests.get",
            lambda *args, **kwargs: MockResponse(),
        )

        def fail_named_tempfile(*args, **kwargs):
            raise AssertionError("temp file should not be created on HTTP error")

        monkeypatch.setattr("process_gaia.tempfile.NamedTemporaryFile", fail_named_tempfile)

        with pytest.raises(process_gaia.requests.HTTPError):
            download_file("https://example.com/test.csv.gz")


class TestProcessFileData:
    def test_no_qualifying_sources(self):
        gz = _make_csv_gz(
            [{"source_id": "1", "bp_flux": "1.0 1.5", "rp_flux": ""}],
            ["source_id", "bp_flux", "rp_flux"],
        )
        results = process_file_data(gz, "test.csv.gz")
        assert results == []

    def test_one_qualifying_source(self):
        gz = _make_csv_gz(
            [{"source_id": "42", "bp_flux": "1.0 5.0", "rp_flux": ""}],
            ["source_id", "bp_flux", "rp_flux"],
        )
        results = process_file_data(gz, "test.csv.gz")
        assert len(results) == 1
        assert results[0]["source_id"] == "42"
        assert math.isclose(results[0]["percentage_change"], 400.0)

    def test_mixed_sources(self):
        rows = [
            {"source_id": "1", "bp_flux": "10.0 10.5", "rp_flux": ""},  # <100%
            {"source_id": "2", "bp_flux": "1.0 10.0", "rp_flux": ""},   # 900%
            {"source_id": "3", "bp_flux": "", "rp_flux": "2.0 2.5"},    # <100%
        ]
        gz = _make_csv_gz(rows, ["source_id", "bp_flux", "rp_flux"])
        results = process_file_data(gz, "test.csv.gz")
        assert len(results) == 1
        assert results[0]["source_id"] == "2"
