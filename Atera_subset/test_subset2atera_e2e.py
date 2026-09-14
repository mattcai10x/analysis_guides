"""End-to-end test for subset2atera.py against a small, synthetic Atera bundle.

This does NOT use spatialdata_io's own atera test fixture (`spatialdata_io/tests
/test_atera.py`'s `_build_atera_bundle`): that fixture deliberately follows the
"best-effort guess" schema for transcripts.zarr.zip (flat x/y/z/feature_name
arrays) and a generic anndata.write_zarr() for cell_feature_matrix.zarr.zip --
which is exactly the schema subset2atera.py's design notes explain is NOT what a
real Atera bundle looks like (see the module docstring of subset2atera.py). Using
it here would validate the wrong thing.

Instead, this builds a bundle using `atera_dataset_tools`'s own real-schema
read/write functions directly (the same functions `atera_dataset_tools.crop` uses
against production data), with 12 cells laid out along a line and cell/transcript/
cell-feature-matrix data kept consistent across files by construction.

Run directly:
    PYTHONPATH=<atera-dataset-tools>:<spatialdata-io>/src:<sopa> \
        python3 test_subset2atera_e2e.py
or under pytest (same PYTHONPATH requirement) -- see the README section in this
directory for the exact command used in development.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import uuid

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from atera_dataset_tools import validate  # noqa: E402
from atera_dataset_tools.formats import cell_feature_matrix as cfm  # noqa: E402
from atera_dataset_tools.formats import cells as cells_mod  # noqa: E402
from atera_dataset_tools.formats import transcripts as tr_mod  # noqa: E402
from atera_dataset_tools.formats.transcripts import CELL_ID_SENTINEL, UNKNOWN_CODEWORD_INDEX  # noqa: E402

PIXEL_SIZE = 0.2125
N_CELLS = 12
GENE_NAMES = ["GeneA", "GeneB", "GeneC"]
TRANSCRIPT_GRID_SIZE = 20.0
TRANSCRIPTS_PER_CELL = 8


def _cell_id_str(prefix: np.ndarray, suffix: np.ndarray) -> list[str]:
    """The real Xenium/Atera hex-nibble-to-letter cell-ID string convention
    (spatialdata_io.readers.xenium.cell_id_str_from_prefix_suffix_uint32),
    reimplemented inline here so this test has no import-time dependency on
    spatialdata_io just to build the fixture (subset2atera.py itself does
    depend on it, at run time, for the real pipeline)."""
    out = []
    for p, s in zip(prefix, suffix):
        nibbles = [(int(p) >> shift) & 0xF for shift in (28, 24, 20, 16, 12, 8, 4, 0)]
        out.append("".join(chr(ord("a") + n) for n in nibbles) + f"-{int(s)}")
    return out


def build_synthetic_bundle(root: str) -> None:
    os.makedirs(root, exist_ok=True)

    # -- cells.zarr.zip: a real CellSegmentationDataset fixture (has the
    #    no-nucleus / multi-nucleus edge cases baked in by make_fixture) --
    cells_path = os.path.join(root, "cells.zarr.zip")
    cells_mod.make_fixture(cells_path, number_cells=N_CELLS, max_vertices=6, grid_size=50.0, seed=0)
    cell_seg = cells_mod.read(cells_path)
    prefix, suffix = cell_seg.cell_id[:, 0], cell_seg.cell_id[:, 1]
    cell_id_strs = _cell_id_str(prefix, suffix)
    cell_cx = cell_seg.cell_summary[:, 0]
    cell_cy = cell_seg.cell_summary[:, 1]

    # -- transcripts.zarr.zip: real RnaDatasetCsV3 schema, ~TRANSCRIPTS_PER_CELL
    #    transcripts scattered near each cell's centroid (+ a few unassigned,
    #    scattered across the whole domain) --
    rng = np.random.default_rng(0)
    codeword_count = 4  # last one left unassigned, matching atera_dataset_tools convention
    codeword_gene_mapping = np.array([0, 1, 2, UNKNOWN_CODEWORD_INDEX], dtype=np.uint32)
    codeword_gene_names = GENE_NAMES + [""]
    codeword_category = np.zeros((codeword_count, len(tr_mod.CATEGORY_COLUMNS)), dtype=bool)
    codeword_category[:-1, tr_mod.CATEGORY_COLUMNS.index("is_gene")] = True
    codeword_category[-1, tr_mod.CATEGORY_COLUMNS.index("is_unassigned_codeword")] = True
    gene_category = np.zeros((len(GENE_NAMES), len(tr_mod.CATEGORY_COLUMNS)), dtype=bool)
    gene_category[:, tr_mod.CATEGORY_COLUMNS.index("is_gene")] = True

    xs, ys, zs, qvs, cws, cids_prefix, cids_suffix, overlaps = [], [], [], [], [], [], [], []
    for c in range(N_CELLS):
        for _ in range(TRANSCRIPTS_PER_CELL):
            xs.append(cell_cx[c] + rng.uniform(-1.5, 1.5))
            ys.append(cell_cy[c] + rng.uniform(-1.5, 1.5))
            zs.append(0.0)
            qvs.append(rng.uniform(0.0, 40.0))
            cws.append(rng.integers(0, 3))
            cids_prefix.append(prefix[c])
            cids_suffix.append(suffix[c])
            overlaps.append(rng.integers(0, 2))
    # A handful of unassigned-to-any-cell transcripts, spread across the domain.
    n_unassigned = N_CELLS * 2
    for _ in range(n_unassigned):
        xs.append(rng.uniform(cell_cx.min() - 5, cell_cx.max() + 5))
        ys.append(rng.uniform(cell_cy.min() - 5, cell_cy.max() + 5))
        zs.append(0.0)
        qvs.append(rng.uniform(0.0, 40.0))
        cws.append(rng.integers(0, 3))
        cids_prefix.append(CELL_ID_SENTINEL)
        cids_suffix.append(CELL_ID_SENTINEL)
        overlaps.append(0)

    transcripts = tr_mod.Transcripts(
        number_genes=len(GENE_NAMES),
        gene_names=GENE_NAMES,
        codeword_count=codeword_count,
        codeword_gene_mapping=codeword_gene_mapping,
        codeword_gene_names=codeword_gene_names,
        codeword_category=codeword_category,
        gene_category=gene_category,
        x_position=np.asarray(xs, dtype=np.float32),
        y_position=np.asarray(ys, dtype=np.float32),
        z_position=np.asarray(zs, dtype=np.float32),
        quality_score=np.asarray(qvs, dtype=np.float16),
        codeword_index=np.asarray(cws, dtype=np.uint32),
        cell_id=np.stack([cids_prefix, cids_suffix], axis=1).astype(np.uint32),
        overlaps_nucleus=np.asarray(overlaps, dtype=np.uint8),
        kit_info=[{"kit_name": "synthetic_kit", "kit_type": "test", "codeword_count": codeword_count}],
        dataset_uuid=str(uuid.uuid4()),
    )
    tr_mod.write(os.path.join(root, "transcripts.zarr.zip"), transcripts, grid_size=TRANSCRIPT_GRID_SIZE)
    total_input_transcripts = len(xs)

    # -- cell_feature_matrix.zarr.zip (CSR) + csc_cell_feature_matrix.zarr.zip (CSC) --
    counts = rng.integers(0, 6, size=(N_CELLS, len(GENE_NAMES))).astype(np.int32)
    obs = pd.DataFrame(index=pd.Index(cell_id_strs, name="barcode"))
    var = pd.DataFrame(index=pd.Index(GENE_NAMES, name="feature_name"))
    import anndata

    adata = anndata.AnnData(X=counts, obs=obs, var=var)
    cfm.write(
        os.path.join(root, "cell_feature_matrix.zarr.zip"), adata,
        obs_index_col="barcode", var_index_col="feature_name", major_axis="cell",
    )
    cfm.write(
        os.path.join(root, "csc_cell_feature_matrix.zarr.zip"), adata,
        obs_index_col="barcode", var_index_col="feature_name", major_axis="feature",
    )

    # -- tiny morphology images (single-resolution; no real pyramid) --
    import tifffile

    os.makedirs(os.path.join(root, "morphology_2d"), exist_ok=True)
    os.makedirs(os.path.join(root, "morphology_3d"), exist_ok=True)
    max_x_px = int(np.ceil((cell_cx.max() + 10) / PIXEL_SIZE))
    max_y_px = int(np.ceil((cell_cy.max() + 10) / PIXEL_SIZE))
    img2d = rng.integers(0, 255, size=(max_y_px, max_x_px), dtype=np.uint8)
    tifffile.imwrite(
        os.path.join(root, "morphology_2d", "ch0000_dapi.ome.tif"), img2d, photometric="minisblack"
    )
    img3d = rng.integers(0, 255, size=(2, max_y_px, max_x_px), dtype=np.uint8)
    tifffile.imwrite(
        os.path.join(root, "morphology_3d", "ch0000_dapi_3d.ome.tif"), img3d, photometric="minisblack"
    )

    # -- experiment.spatial manifest --
    manifest = {
        "major_version": 1,
        "minor_version": 2,
        "chemistry_version": "Atera v1 (synthetic)",
        "run_name": "synthetic_test_run",
        "run_start_time": "2026-01-01T00:00:00Z",
        "slide_name": "synthetic_slide",
        "region_name": "synthetic_region",
        "sample_type": "fresh_frozen",
        "pixel_size": PIXEL_SIZE,
        "z_step_size": 1.5,
        "total_num_targets": len(GENE_NAMES),
        "num_cells": N_CELLS,
        "panel_config": "synthetic_panel",
        "explorer_files": {
            "transcripts_zarr_filepath": "transcripts.zarr.zip",
            "cells_zarr_filepath": "cells.zarr.zip",
            "cell_feature_zarr_filepath": "cell_feature_matrix.zarr.zip",
            "cell_feature_viz_zarr_filepath": "csc_cell_feature_matrix.zarr.zip",
        },
        "images": {
            "morphology_2d_filepath": "morphology_2d/ch0000_dapi.ome.tif",
            "morphology_3d_filepath": "morphology_3d/ch0000_dapi_3d.ome.tif",
        },
    }
    with open(os.path.join(root, "experiment.spatial"), "w") as f:
        json.dump(manifest, f, indent=4)

    return {
        "n_cells": N_CELLS,
        "n_transcripts": total_input_transcripts,
        "cell_cx": cell_cx,
        "cell_cy": cell_cy,
    }


def build_polygon_geojson(path: str, max_x: float) -> None:
    """A rectangle covering the left half of the cell layout (x < max_x/2)."""
    import geopandas as gpd
    from shapely.geometry import box

    half = max_x / 2.0
    poly = box(-5.0, -20.0, half, 80.0)
    gpd.GeoDataFrame({"geometry": [poly]}, crs="EPSG:4326").to_file(path, driver="GeoJSON")


def run_subset2atera(input_dir: str, polygon_path: str, output_dir: str, extra_args: list[str] | None = None) -> str:
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "subset2atera.py")
    cmd = [sys.executable, script, "-i", input_dir, "-p", polygon_path, "-o", output_dir]
    if extra_args:
        cmd += extra_args
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    print(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(f"subset2atera.py failed (exit {result.returncode}); see stdout/stderr above")
    return result.stdout + result.stderr


def test_dry_run(tmp_path):
    root = str(tmp_path / "bundle")
    info = build_synthetic_bundle(root)
    polygon_path = str(tmp_path / "polygon.geojson")
    build_polygon_geojson(polygon_path, max_x=float(info["cell_cx"].max()))

    output_dir = str(tmp_path / "out_dryrun")
    out = run_subset2atera(root, polygon_path, output_dir, extra_args=["--dry-run"])
    assert "[dry-run]" in out
    assert not os.path.exists(os.path.join(output_dir, "experiment.spatial"))


def test_end_to_end(tmp_path):
    root = str(tmp_path / "bundle")
    info = build_synthetic_bundle(root)
    polygon_path = str(tmp_path / "polygon.geojson")
    build_polygon_geojson(polygon_path, max_x=float(info["cell_cx"].max()))

    output_dir = str(tmp_path / "out")
    run_subset2atera(root, polygon_path, output_dir, extra_args=["--keep-tmp"])

    # -- structural validity (real atera_dataset_tools validators) --
    assert validate.validate_cells(os.path.join(output_dir, "cells.zarr.zip")) == []
    assert validate.validate_transcripts(os.path.join(output_dir, "transcripts.zarr.zip")) == []
    assert validate.validate_cell_feature_matrix(os.path.join(output_dir, "cell_feature_matrix.zarr.zip")) == []
    assert validate.validate_cell_feature_matrix(os.path.join(output_dir, "csc_cell_feature_matrix.zarr.zip")) == []

    # -- experiment.spatial is valid JSON with the expected shape --
    with open(os.path.join(output_dir, "experiment.spatial")) as f:
        out_manifest = json.load(f)
    assert out_manifest["pixel_size"] == PIXEL_SIZE
    assert out_manifest["run_name"] == "synthetic_test_run"  # copied forward

    # -- it actually subsetted something, not a no-op --
    out_cells = cells_mod.read(os.path.join(output_dir, "cells.zarr.zip"))
    out_transcripts = tr_mod.read(os.path.join(output_dir, "transcripts.zarr.zip"))
    assert 0 < out_cells.number_cells < info["n_cells"]
    assert 0 < out_transcripts.number_rnas < info["n_transcripts"]
    assert out_manifest["num_cells"] == out_cells.number_cells
    assert out_manifest["num_transcripts"] == out_transcripts.number_rnas

    out_adata = cfm.read(os.path.join(output_dir, "cell_feature_matrix.zarr.zip"))
    assert out_adata.n_obs == out_cells.number_cells

    # -- morphology images were cropped (smaller than source), not just copied --
    import tifffile

    src_img = tifffile.imread(os.path.join(root, "morphology_2d", "ch0000_dapi.ome.tif"))
    out_img = tifffile.imread(os.path.join(output_dir, "morphology_2d", "ch0000_dapi.ome.tif"))
    assert out_img.shape[-1] < src_img.shape[-1]  # narrower in X (half the domain was cut)

    # -- Stage-2 tmp dir preserved with --keep-tmp --
    assert os.path.isdir(os.path.join(output_dir, ".stage2_tmp"))

    print(f"input cells={info['n_cells']} transcripts={info['n_transcripts']}")
    print(f"output cells={out_cells.number_cells} transcripts={out_transcripts.number_rnas}")


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        class _P:
            def __init__(self, base):
                self._base = base

            def __truediv__(self, other):
                return os.path.join(self._base, other)

        tmp_path = _P(td)
        test_dry_run(tmp_path)
        test_end_to_end(tmp_path)
    print("ALL OK")
