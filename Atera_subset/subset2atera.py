# ---
# title: "atera subset tool"
# date: "13-September-2026"
# author: "Matthew Cai @10x Genomics"
# ---
#
# Spatially subset an Atera ("Gen2") output bundle using a polygon GeoJSON, producing
# a smaller output bundle in the same format, for visualization in ziggy.
#
# This is the Atera-format sibling of `Xenium_subset/subset2zarr.py`. It is NOT a
# mechanical swap of `xenium()`/`sopa.io.explorer.write()` for
# `atera()`/`sopa.io.atera.write()`: a real Atera bundle is far too large to load
# fully into memory/dask before filtering (see the module docstring in the repo's
# design notes / the task write-up this was built from -- a real bundle is
# `transcripts.zarr.zip` ~43GB, `cell_feature_matrix.zarr.zip` ~6.9GB,
# `csc_cell_feature_matrix.zarr.zip` ~6.9GB (redundant, viz-only),
# `binned_transcripts.zarr.zip` ~5.3GB (viz-only, derived), `cells.zarr.zip` ~1.3GB;
# ~64GB total). Instead, this script filters at the raw-zarr level *before*
# constructing any `SpatialData` object, using `atera_dataset_tools` primitives that
# touch each large file at most once, and only builds an in-memory `SpatialData`
# object out of the already-small, already-cropped intermediate data.
#
# --------------------------------------------------------------------------------
# IMPORTANT DEVIATION FROM A NAIVE spatialdata_io.atera()/sopa.io.atera.write() PORT
# --------------------------------------------------------------------------------
# While building this pipeline, three concrete schema incompatibilities were found
# between the "real" Atera schema (implemented in `atera_dataset_tools`, confirmed
# against a real production bundle per that package's own module docstrings) and
# the schema assumed by `spatialdata_io.readers.atera.atera()` / `sopa.io.atera.write()`
# (which the latter two openly flag as a "best-effort guess ... not confirmed
# against a real production bundle" for the transcript file, and read/write the
# cell-feature matrix as a generic `anndata.write_zarr()`/`anndata.read_zarr()`
# store rather than the real custom X/obs/var-descriptor layout):
#
#   1. transcripts.zarr.zip: the real schema (`atera_dataset_tools.formats.transcripts`,
#      `RnaDatasetCsV3`) is a flat `grid/` group keyed by `(grid_x, grid_y)` tiles,
#      each holding `location` (x,y,z), `quality_score`, `codeword_identity` (an index
#      into a root-level `codeword_gene_mapping` table) and `cell_id` arrays.
#      `spatialdata_io.readers.atera._read_transcripts_dataframe` instead expects
#      flat `x`/`y`/`z`/`feature_name` arrays directly under the root or under each
#      tile group -- it will not find real transcript data in a real bundle (it
#      would either silently misparse or raise `ValueError`/`KeyError`).
#   2. cell_feature_matrix.zarr.zip / csc_cell_feature_matrix.zarr.zip: the real
#      schema is a custom `X`/`obs`/`var` zarr layout (`atera_dataset_tools.formats
#      .cell_feature_matrix`, a pure-Python reimplementation of turing's compiled
#      writer). `spatialdata_io.readers.atera._get_table_and_circles` reads this
#      file with `anndata.read_zarr()`, which expects AnnData's own on-disk zarr
#      layout -- a real bundle's file will not open as a valid AnnData store.
#   3. cells.zarr.zip: `sopa.io.atera.write_polygons` (the writer side) asserts a
#      strict 1-to-1 cell<->nucleus pairing and cannot write cells with zero or
#      multiple nuclei -- a case `atera_dataset_tools.formats.cells.make_fixture`
#      (and the real spec) explicitly documents as a case Atera data has.
#
# Given this, this script uses `spatialdata_io.readers.atera.atera()` only for the
# parts of the schema that really are compatible (images, labels, cell/nucleus
# boundary shapes, and the manifest-derived `SpatialData.attrs`), and uses
# `atera_dataset_tools`'s own (real-schema, production-confirmed) read/write/crop
# functions directly for transcripts, the cell-feature matrices, and cells.zarr.zip.
# `sopa.io.atera.write()` is not called at all in the main pipeline -- see the
# report this script shipped with for the full rationale. `spatialdata.SpatialData
# .write()` (the *internal* SpatialData zarr representation used by `--zarr_out`,
# unrelated to the Atera bundle format) has no such issue and is used normally.
#
# --------------------------------------------------------------------------------
# STAGED PIPELINE (memory/runtime budget is the primary design constraint)
# --------------------------------------------------------------------------------
# Stage 0: parse the polygon GeoJSON, compute its bounding box.
#
# Stage 1 (cheap, in-memory-trivial): read only cells.zarr.zip's /cell_summary +
#   /cell_id from the ORIGINAL bundle (via the new, lightweight
#   `atera_dataset_tools.formats.cells.read_cell_summary`, added alongside this
#   script -- see that module for why it never touches /masks or /polygon_sets).
#   Test each cell's centroid for polygon containment directly: this is already
#   exact, not a bbox approximation, so it gives the final `keep_cell_ids` set
#   immediately.
#
# Stage 2 (tile-level / targeted crops on the ORIGINAL large files -> small
#   intermediates in a deterministic, resumable tmp directory; each large file is
#   touched at most once):
#   a. crop_transcripts_to_bbox (tile-level prefilter, no decompression of kept
#      tiles' point data).
#   b. crop_cells_by_ids using the exact keep_cell_ids from Stage 1.
#   c. subset_cell_feature_matrix on cell_feature_matrix.zarr.zip, then (strictly
#      sequentially, never held in memory simultaneously) on
#      csc_cell_feature_matrix.zarr.zip.
#   d. Morphology OME-TIFF(s): windowed reads via tifffile's zarr-store interface
#      (`TiffFile(...).series[...].levels[i].aszarr()`), one pyramid level at a
#      time, so the full slide image is never materialized in memory. See
#      `crop_morphology_image()`'s docstring for exactly what is shipped vs. ideal.
#   e. binned_transcripts.zarr.zip and csc_cell_feature_matrix.zarr.zip are both
#      viz-only/derived (not needed for *correctness*), but ARE included by
#      default: the whole point of this script is a bundle that's actually
#      visualizable in ziggy, not just structurally valid. csc_* was never
#      skipped (see 2c above). binned_transcripts.zarr.zip is cropped via
#      `crop_binned_transcripts_to_bbox`, which streams one gene at a time
#      (see `atera_dataset_tools.crop`'s module docstring) rather than loading
#      the whole ~5.3GB file at once -- pass `--skip-binned-transcripts` to
#      omit it anyway (e.g. if ziggy's density view isn't needed for a given
#      use case and you'd rather skip the extra I/O).
#
# Stage 3 (small data now -- safe to fully materialize):
#   a. `spatialdata_io.atera()` on the Stage-2 bundle for images/labels/shapes only.
#   b. Exact polygon (not bbox) filtering of the cropped transcripts (cheap now).
#   c. Manually build the `table`/`points["transcripts"]` elements from
#      `atera_dataset_tools` reads (see deviation notes above) and sync the table
#      against the kept cell/nucleus boundaries, mirroring `subset2zarr.py` step 5.
#   d. Write the final bundle: transcripts.zarr.zip and the cell-feature matrices
#      via `atera_dataset_tools` (real schema); cells.zarr.zip and morphology
#      images are the already-correct Stage-2 intermediates (re-cropped only if
#      the sync step in 3c dropped anything further, which should not normally
#      happen since Stage 1 was already exact).
#   e. Patch `experiment.spatial` with copied-forward manifest fields + recomputed
#      `num_cells`/`transcripts_per_cell`/`num_transcripts`/`num_transcripts_high_quality`.
#   f. Optional full zarr dump of the final small `SpatialData` object (`-z/--zarr_out`).

from __future__ import annotations

import argparse
import json
import os
import resource
import shutil
import sys
import time
import warnings

import geopandas as gpd
import numpy as np
import pandas as pd
import zarr

warnings.filterwarnings("ignore")

# `zarr.config` is a zarr-v3-only knob; `atera_dataset_tools` (and this script's
# own zarr usage) targets zarr v2 (zarr.Blosc, ZipStore, open_group), so this is
# best-effort only -- skip it quietly on a v2-only zarr install.
if hasattr(zarr, "config"):
    zarr.config.set({"array.rectilinear_chunks": True})

from atera_dataset_tools import crop as adt_crop  # noqa: E402
from atera_dataset_tools.formats import cell_feature_matrix as adt_cfm  # noqa: E402
from atera_dataset_tools.formats import cells as adt_cells  # noqa: E402
from atera_dataset_tools.formats import transcripts as adt_transcripts  # noqa: E402
from atera_dataset_tools.formats.transcripts import UNKNOWN_CODEWORD_INDEX  # noqa: E402

# --------------------------------------------------------------------------------
# A note on cell-ID string conventions (a second, smaller schema wrinkle found
# while building this pipeline, distinct from the three noted above): the true
# on-disk Xenium/Atera cell ID string is a hex-nibble-to-letter encoding of the
# (prefix, suffix) uint32 pair (`spatialdata_io.readers.xenium
# .cell_id_str_from_prefix_suffix_uint32`, e.g. "aaaaaaaa-1" -- see that
# function's docstring/the 10x xoa-output-zarr#cellID spec it cites), which is
# what `spatialdata_io.atera()` uses to build `cell_boundaries.index`. But
# `atera_dataset_tools.crop.crop_cells_by_ids`'s own `keep_cell_ids` token format
# is a *plain-decimal* "<prefix>-<suffix>" string (see its docstring) -- a
# different, simpler convention that cannot even parse a real hex-letter ID
# (`int("aaaaaaaa")` raises). Rather than gamble on which convention a real
# `cell_feature_matrix.zarr.zip`'s `obs` barcode column actually uses on disk,
# this script sidesteps the question entirely: it selects cells by *row
# position* (unambiguous for `crop_cells_by_ids`, which explicitly also accepts
# a plain integer row index), and -- assuming `cells.zarr.zip` and
# `cell_feature_matrix.zarr.zip` are row-aligned per cell, which a single
# turing pipeline run emitting both files together should guarantee -- reads
# back `cell_feature_matrix.zarr.zip`'s own real obs-index strings at those row
# positions (via the new `atera_dataset_tools.formats.cell_feature_matrix
# .read_obs_index`, added alongside this script) rather than reconstructing or
# guessing them. This makes the pipeline correct regardless of which string
# convention the real barcodes turn out to use.
# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
# Manifest field names (mirrors spatialdata_io.readers.atera.AteraKeys / the real
# experiment.spatial layout -- kept as plain strings here so this script has no
# hard dependency on spatialdata_io for anything other than the atera() reader call).
# --------------------------------------------------------------------------------
MANIFEST_FILENAME = "experiment.spatial"
EXPLORER_FILES_KEY = "explorer_files"
IMAGES_KEY = "images"
TRANSCRIPTS_ZARR_KEY = "transcripts_zarr_filepath"
TRANSCRIPTS_VIZ_ZARR_KEY = "transcripts_viz_zarr_filepath"
CELLS_ZARR_KEY = "cells_zarr_filepath"
CELL_FEATURE_ZARR_KEY = "cell_feature_zarr_filepath"
CELL_FEATURE_VIZ_ZARR_KEY = "cell_feature_viz_zarr_filepath"
MORPHOLOGY_2D_KEY = "morphology_2d_filepath"
MORPHOLOGY_3D_KEY = "morphology_3d_filepath"

# Manifest fields copied forward verbatim from the source experiment.spatial (the
# rest -- num_cells, transcripts_per_cell, num_transcripts, num_transcripts_high_quality
# -- are always recomputed from the subsetted data; see patch_manifest()).
MANIFEST_KEYS_TO_COPY = [
    "major_version",
    "minor_version",
    "chemistry_version",
    "run_name",
    "run_start_time",
    "slide_name",
    "region_name",
    "sample_type",
    "pixel_size",
    "z_step_size",
    "total_num_targets",
    "panel_config",
]

_MEM_T0 = time.monotonic()


def _peak_rss_gb() -> float:
    """Peak RSS of this process so far, in GiB (ru_maxrss is KiB on Linux, bytes on macOS)."""
    ru = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    divisor = (1024 * 1024) if sys.platform.startswith("linux") else (1024 * 1024 * 1024)
    return ru / divisor


def log_checkpoint(label: str) -> None:
    elapsed = time.monotonic() - _MEM_T0
    print(f"[{elapsed:8.1f}s elapsed | peak RSS {_peak_rss_gb():6.2f} GiB] {label}")


# --------------------------------------------------------------------------------
# Stage 0
# --------------------------------------------------------------------------------


def load_polygon(polygon_path: str, *, units: str, pixel_size: float):
    """Load the selection polygon and convert it to microns if needed.

    Every downstream spatial comparison in this script -- cell centroids
    (`/cell_summary`), transcript positions, tile bboxes -- is in microns.
    Polygons exported from an image-space viewer (drawing against the
    morphology image, which is natively pixel-indexed) are in *pixels*
    instead, and nothing about a GeoJSON file's contents distinguishes the
    two: the coordinates are just floats, and a pixel-space polygon over a
    real tissue image can easily have a bounding box that superficially
    looks plausible in micron space too (both are "a few thousand units"),
    so this mismatch does not fail loudly -- it just silently selects zero
    (or the wrong) cells. That's exactly what happened building this
    pipeline: a real bundle's dry run reported tile-level bbox overlap but
    zero exact cell matches, which turned out to be a pixel/micron unit
    mismatch, not a location or shapely bug. Hence `--polygon-units` is a
    required flag rather than a silently-assumed default.
    """
    from shapely import affinity

    gdf = gpd.read_file(polygon_path)
    geom = gdf.geometry[0]
    if units == "pixels":
        # pixel_size is microns-per-pixel (same manifest field used elsewhere
        # in this script, e.g. the morphology-image crop and the transcripts
        # Scale transform) -- so microns = pixels * pixel_size. Scale about
        # the true origin (0, 0), not shapely's default "center" origin,
        # since this is an absolute unit conversion, not a resize.
        geom = affinity.scale(geom, xfact=pixel_size, yfact=pixel_size, origin=(0, 0))
    return geom


# --------------------------------------------------------------------------------
# Stage 1: exact cell selection from /cell_summary + /cell_id only
# --------------------------------------------------------------------------------


def select_cells(cells_zarr_path: str, target_polygon) -> tuple[np.ndarray, int, int]:
    """Return (keep_row_indices, n_cells_total, n_cells_kept) using only the cheap
    cell_summary/cell_id arrays -- never touches /masks or /polygon_sets.

    Row indices (not cell-ID strings) are the canonical selection here -- see the
    module-level note on cell-ID string conventions above for why.
    """
    from shapely import vectorized

    summary = adt_cells.read_cell_summary(cells_zarr_path)
    cx, cy = summary.cell_centroid_x, summary.cell_centroid_y
    keep_mask = vectorized.contains(target_polygon, cx, cy)
    keep_rows = np.where(keep_mask)[0]
    return keep_rows, summary.number_cells, len(keep_rows)


def cfm_barcodes_for_rows(cell_feature_matrix_src: str, keep_rows: np.ndarray) -> list[str]:
    """Read cell_feature_matrix.zarr.zip's real obs-index strings at `keep_rows`
    (cheap: only the small obs-index array is read, never `X`). Assumes
    row-alignment with cells.zarr.zip -- see the module-level note above."""
    obs_index = adt_cfm.read_obs_index(cell_feature_matrix_src)
    return list(np.asarray(obs_index)[keep_rows].astype(str))


# --------------------------------------------------------------------------------
# Stage 2: targeted crops of the large source files -> small tmp intermediates
# --------------------------------------------------------------------------------


def _skip_if_present(dst: str, force: bool) -> bool:
    if force:
        return False
    return os.path.exists(dst)


def crop_morphology_image(
    src_path: str,
    dst_path: str,
    bbox_px: tuple[float, float, float, float],
    margin_px: float = 32.0,
) -> None:
    """Windowed crop of one (potentially pyramidal) OME-TIFF, without ever loading
    a full resolution level into memory.

    Shipped implementation: for each pyramid level present in the source series
    (``TiffFile(...).series[0].levels``), open that level's zarr-store view via
    ``.aszarr()``, compute the proportionally-scaled crop window for that level, and
    slice only that window out (a tiled TIFF decodes only the touched tiles; a
    stripped/non-tiled TIFF decodes only the touched scanlines -- either way, far
    less than the whole level). All cropped levels are written back as a single
    new *pyramidal* OME-TIFF (subimages appended as literal separate pages via
    `tifffile.imwrite(..., data=[level0, level1, ...])`), preserving the same
    number of resolution levels as the source.

    Known gap vs. an ideal implementation: this does not re-derive OME-XML pyramid
    metadata (`SubIFDs`, per-level `PhysicalSize*`, etc.) -- the levels are written
    as plain sequential pages/series, which round-trips the pixel data and level
    count but may not be recognized as a "true" OME pyramid by every downstream
    tool. Memory cost of what is shipped: bounded by the crop window's pixel count
    at each level (at most a few tens of MB for a typical ROI), never the full
    slide (which is exactly the ~10s-of-GB-per-channel object this function exists
    to avoid loading).
    """
    import tifffile

    min_x, min_y, max_x, max_y = bbox_px
    min_x -= margin_px
    min_y -= margin_px
    max_x += margin_px
    max_y += margin_px

    with tifffile.TiffFile(src_path) as tf:
        series = tf.series[0]
        levels = list(getattr(series, "levels", [series]))
        base_shape = levels[0].shape
        base_axes = levels[0].axes  # e.g. "YX", "CYX", "CZYX"
        y_ax = base_axes.index("Y")
        x_ax = base_axes.index("X")
        base_h, base_w = base_shape[y_ax], base_shape[x_ax]

        cropped_levels = []
        for i, level in enumerate(levels):
            # NOTE(tifffile-zarr3): tifffile >=2026.5.2 rewrote ZarrTiffStore
            # for zarr format 3 / NGFF 0.5. Calling .aszarr() on an individual
            # per-level TiffPageSeries object with no level= kwarg does NOT
            # give you just that level -- confirmed empirically: ZarrTiffStore
            # reads self._data from `arg.levels`, which apparently returns the
            # *same* full pyramid-level list regardless of which per-level
            # object `arg` is, so every call produced a >1-level, NGFF
            # multiscales *group* (9 numeric keys "0".."8" for this file's
            # 9-level pyramid) rather than a bare per-level Array -- and
            # `za[tuple(slice, ...)]` on a Group raises TypeError (expects a
            # string key). ZarrTiffStore's own documented `level` parameter
            # ("Pyramidal level to wrap") is exactly for this: passing it
            # explicitly makes self._data a single-item list, which takes the
            # non-multiscales path and returns a plain, directly-sliceable
            # Array store -- matching the old (pre-2026.5.2) behavior exactly.
            store = series.aszarr(level=i)
            try:
                za = zarr.open(store, mode="r")
                lvl_axes = level.axes
                ly, lx = lvl_axes.index("Y"), lvl_axes.index("X")
                lvl_h, lvl_w = level.shape[ly], level.shape[lx]
                scale_y = lvl_h / base_h
                scale_x = lvl_w / base_w

                y0 = max(0, int(np.floor(min_y * scale_y)))
                y1 = min(lvl_h, int(np.ceil(max_y * scale_y)))
                x0 = max(0, int(np.floor(min_x * scale_x)))
                x1 = min(lvl_w, int(np.ceil(max_x * scale_x)))
                y0, y1 = min(y0, y1), max(y0, y1)
                x0, x1 = min(x0, x1), max(x0, x1)
                y1, x1 = max(y1, y0 + 1), max(x1, x0 + 1)

                sl = [slice(None)] * len(lvl_axes)
                sl[ly] = slice(y0, y1)
                sl[lx] = slice(x0, x1)
                cropped_levels.append(np.asarray(za[tuple(sl)]))
            finally:
                if hasattr(store, "close"):
                    store.close()

    os.makedirs(os.path.dirname(os.path.abspath(dst_path)) or ".", exist_ok=True)
    if len(cropped_levels) == 1:
        tifffile.imwrite(dst_path, cropped_levels[0], photometric="minisblack")
    else:
        # Sequential pages, one per (cropped) pyramid level. See docstring: this
        # preserves level count/pixels, not full OME pyramid SubIFD metadata.
        #
        # NOTE(tifffile-zarr3): newer tifffile (confirmed against a version
        # after 2026.5.2) refuses TiffWriter(..., append=True) with
        # "cannot append to file containing metadata" whenever the base file
        # was written with tifffile's default embedded "shaped" metadata --
        # a real safety check (that metadata encodes an array shape/page
        # count that a blind append would leave stale), not a bug. Since this
        # code was never trying to maintain synced OME/shaped pyramid
        # metadata across the appended pages anyway (see comment above),
        # metadata=None on the first write keeps the base file a plain,
        # metadata-less multi-page TIFF that's genuinely appendable, rather
        # than forcing past the check with append='force'.
        tifffile.imwrite(
            dst_path, cropped_levels[0], photometric="minisblack", metadata=None
        )
        with tifffile.TiffWriter(dst_path, append=True) as tw:
            for lvl in cropped_levels[1:]:
                tw.write(lvl, photometric="minisblack")


def maybe_crop_binned_transcripts(src: str, dst: str, bbox, skip: bool) -> bool:
    """Returns True if a binned_transcripts.zarr.zip was written to dst.

    Included by default (see `skip`'s default in `main()`): binned_transcripts.zarr.zip
    is viz-only/derived, but this script's whole point is a bundle that's actually
    visualizable in ziggy, so it isn't optional in practice. Uses
    `crop_binned_transcripts_to_bbox`, which streams one gene at a time (reads,
    filters, and writes each gene's nested-zip archive independently) rather than
    loading the whole file into memory -- see `atera_dataset_tools.crop`'s module
    docstring. Possible future optimization, not yet implemented: regenerating the
    density raster directly from the already-small cropped transcripts, rather
    than touching the original file at all -- would save the one remaining full
    pass over binned_transcripts.zarr.zip, but isn't needed for correctness or to
    keep memory bounded (the streaming crop already does that).
    """
    if skip:
        return False
    adt_crop.crop_binned_transcripts_to_bbox(src, dst, bbox)
    return True


def stage2_crop(input_dir: str, manifest: dict, bbox, keep_rows: np.ndarray, tmp_dir: str, args) -> dict:
    """Crop every large source file into `tmp_dir`, skipping files that already
    exist there unless `args.force_redo` is set. Returns a dict of the manifest
    keys we care about -> absolute path of the tmp intermediate.
    """
    os.makedirs(tmp_dir, exist_ok=True)
    explorer_files = manifest[EXPLORER_FILES_KEY]
    images_manifest = manifest.get(IMAGES_KEY, {})
    out = {}

    # -- transcripts.zarr.zip: tile-level bbox prefilter --
    src = os.path.join(input_dir, explorer_files[TRANSCRIPTS_ZARR_KEY])
    dst = os.path.join(tmp_dir, "transcripts.zarr.zip")
    if not _skip_if_present(dst, args.force_redo):
        adt_crop.crop_transcripts_to_bbox(src, dst, bbox)
    else:
        print(f"  [resume] {dst} already exists, skipping re-crop")
    out["transcripts"] = dst
    log_checkpoint("Stage 2a: transcripts.zarr.zip cropped to bbox tiles")

    # -- cells.zarr.zip: exact cell crop, selected by row index (unambiguous;
    #    see the module-level note on cell-ID string conventions) --
    keep_row_tokens = [str(int(r)) for r in keep_rows]
    src = os.path.join(input_dir, explorer_files[CELLS_ZARR_KEY])
    dst = os.path.join(tmp_dir, "cells.zarr.zip")
    if not _skip_if_present(dst, args.force_redo):
        adt_crop.crop_cells_by_ids(src, dst, keep_row_tokens)
    else:
        print(f"  [resume] {dst} already exists, skipping re-crop")
    out["cells"] = dst
    log_checkpoint("Stage 2b: cells.zarr.zip cropped to exact keep_rows")

    # -- cell_feature_matrix.zarr.zip (CSR), then csc_* (CSC) -- sequentially --
    # `subset_cell_feature_matrix` matches by obs-index *string*, not row index,
    # so first read back this file's own real barcode strings at `keep_rows`
    # (cheap: obs-index only, never touches X) -- see module-level note above.
    cfm_rel = explorer_files[CELL_FEATURE_ZARR_KEY]
    src = os.path.join(input_dir, cfm_rel)
    keep_cells_cfm = cfm_barcodes_for_rows(src, keep_rows)
    dst = os.path.join(tmp_dir, "cell_feature_matrix.zarr.zip")
    if not _skip_if_present(dst, args.force_redo):
        adt_crop.subset_cell_feature_matrix(src, dst, keep_cells=keep_cells_cfm)
    else:
        print(f"  [resume] {dst} already exists, skipping re-crop")
    out["cell_feature_matrix"] = dst
    log_checkpoint("Stage 2c: cell_feature_matrix.zarr.zip (CSR) subset")

    if CELL_FEATURE_VIZ_ZARR_KEY in explorer_files:
        # csc_cell_feature_matrix.zarr.zip is presumed row-aligned with the CSR
        # file above (same per-cell ordering, both derived from the same table);
        # reuse the same barcode strings rather than re-deriving them.
        src = os.path.join(input_dir, explorer_files[CELL_FEATURE_VIZ_ZARR_KEY])
        dst = os.path.join(tmp_dir, "csc_cell_feature_matrix.zarr.zip")
        if not _skip_if_present(dst, args.force_redo):
            adt_crop.subset_cell_feature_matrix(src, dst, keep_cells=keep_cells_cfm)
        else:
            print(f"  [resume] {dst} already exists, skipping re-crop")
        out["csc_cell_feature_matrix"] = dst
        log_checkpoint("Stage 2c: csc_cell_feature_matrix.zarr.zip (CSC) subset")

    # -- morphology images: windowed pyramid crop --
    pixel_size = float(manifest["pixel_size"])
    bbox_px = tuple(v / pixel_size for v in bbox)
    for key, out_key in ((MORPHOLOGY_2D_KEY, "morphology_2d"), (MORPHOLOGY_3D_KEY, "morphology_3d")):
        if key not in images_manifest:
            continue
        rel = images_manifest[key]
        src = os.path.join(input_dir, rel)
        dst = os.path.join(tmp_dir, rel)
        if not _skip_if_present(dst, args.force_redo):
            crop_morphology_image(src, dst, bbox_px)
        else:
            print(f"  [resume] {dst} already exists, skipping re-crop")
        out[out_key] = dst
        out[out_key + "_rel"] = rel
    log_checkpoint("Stage 2d: morphology image(s) windowed-cropped")

    # -- binned_transcripts.zarr.zip: skipped by default --
    if TRANSCRIPTS_VIZ_ZARR_KEY in explorer_files:
        src = os.path.join(input_dir, explorer_files[TRANSCRIPTS_VIZ_ZARR_KEY])
        dst = os.path.join(tmp_dir, "binned_transcripts.zarr.zip")
        if not _skip_if_present(dst, args.force_redo):
            wrote = maybe_crop_binned_transcripts(src, dst, bbox, args.skip_binned_transcripts)
        else:
            wrote = True
            print(f"  [resume] {dst} already exists, skipping re-crop")
        if wrote:
            out["binned_transcripts"] = dst
        log_checkpoint(
            "Stage 2e: binned_transcripts.zarr.zip "
            + ("skipped (--skip-binned-transcripts)" if not wrote else "cropped")
        )

    return out


# --------------------------------------------------------------------------------
# Stage 3: everything left is small -- build SpatialData, exact-filter transcripts,
# sync table/boundaries, write final bundle.
# --------------------------------------------------------------------------------


def _resolve_gene_names(t: "adt_transcripts.Transcripts") -> np.ndarray:
    gene_ix = np.full(t.codeword_index.shape, -1, dtype=np.int64)
    valid_cw = t.codeword_index < len(t.codeword_gene_mapping)
    mapped = t.codeword_gene_mapping[np.where(valid_cw, t.codeword_index, 0).astype(np.int64)]
    has_gene = valid_cw & (mapped != UNKNOWN_CODEWORD_INDEX) & (mapped < t.number_genes)
    gene_ix[has_gene] = mapped[has_gene].astype(np.int64)
    gene_names = np.asarray(t.gene_names + ["Unassigned"])
    gene_ix[~has_gene] = len(t.gene_names)
    return gene_names[gene_ix]


def filter_transcripts_to_polygon(cropped_transcripts_path: str, target_polygon):
    """Read the (already tile-bbox-cropped, small) transcripts.zarr.zip and keep
    only the transcripts whose exact (x, y) fall inside `target_polygon`.
    Returns the filtered `Transcripts` dataclass and the resolved gene-name array
    (aligned to the *filtered* rows) for downstream reuse (points df + metrics).
    """
    from shapely import vectorized

    t = adt_transcripts.read(cropped_transcripts_path)
    mask = vectorized.contains(target_polygon, t.x_position.astype(np.float64), t.y_position.astype(np.float64))

    filtered = adt_transcripts.Transcripts(
        number_genes=t.number_genes,
        gene_names=t.gene_names,
        codeword_count=t.codeword_count,
        codeword_gene_mapping=t.codeword_gene_mapping,
        codeword_gene_names=t.codeword_gene_names,
        codeword_category=t.codeword_category,
        gene_category=t.gene_category,
        x_position=t.x_position[mask],
        y_position=t.y_position[mask],
        z_position=t.z_position[mask],
        quality_score=t.quality_score[mask],
        codeword_index=t.codeword_index[mask],
        cell_id=t.cell_id[mask],
        overlaps_nucleus=t.overlaps_nucleus[mask],
        kit_info=t.kit_info,
        dataset_uuid=t.dataset_uuid,
        spatial_units=t.spatial_units,
        coordinate_space=t.coordinate_space,
        data_format=t.data_format,
    )
    gene_names_per_row = _resolve_gene_names(t)[mask]
    return filtered, gene_names_per_row


def build_table(cropped_cfm_path: str, cropped_cells_path: str, region: str = "cell_boundaries"):
    """Build a TableModel-ready AnnData from the (already cropped, small)
    cell_feature_matrix.zarr.zip + cells.zarr.zip, mirroring the obs enrichment
    that `spatialdata_io.readers.atera._get_table_and_circles` does for the
    schema-compatible parts (cell_summary), without going through that function's
    `anndata.read_zarr()` call (see module-level deviation notes: the real
    cell_feature_matrix schema is not a generic AnnData zarr store).
    """
    from spatialdata.models import TableModel

    adata = adt_cfm.read(cropped_cfm_path)
    summary = adt_cells.read_cell_summary(cropped_cells_path)
    # Positional (row-index) alignment, not a label/string join -- see the
    # module-level note on cell-ID string conventions: cropped_cfm_path and
    # cropped_cells_path are two independently-cropped-but-row-aligned files
    # (both produced from the same keep_rows selection in Stage 2), and there is
    # no guarantee their respective cell-ID *strings* use the same convention.
    if summary.number_cells != adata.n_obs:
        raise ValueError(
            f"cells.zarr.zip has {summary.number_cells} cells but "
            f"cell_feature_matrix.zarr.zip has {adata.n_obs}; the two are expected "
            "to be row-aligned per cell (see module docstring) and are not."
        )
    cell_summary_df = pd.DataFrame(summary.cell_summary, columns=adt_cells.CELL_SUMMARY_COLUMNS)
    cell_summary_df.index = adata.obs.index

    adata.obs["cell_id"] = adata.obs.index.astype(str)
    for col in ("z_level", "nucleus_count", "cell_area"):
        adata.obs[col] = cell_summary_df[col].to_numpy()
    adata.obsm["spatial"] = cell_summary_df[["cell_centroid_x", "cell_centroid_y"]].to_numpy()
    adata.obs["region"] = region
    adata.obs["region"] = adata.obs["region"].astype("category")

    return TableModel.parse(adata, region=region, region_key="region", instance_key="cell_id")


def build_points(filtered_transcripts, gene_names_per_row, pixel_size: float):
    from spatialdata.models import PointsModel
    from spatialdata.transformations.transformations import Scale
    from spatialdata_io.readers.xenium import cell_id_str_from_prefix_suffix_uint32

    from atera_dataset_tools.formats.transcripts import CELL_ID_SENTINEL

    t = filtered_transcripts
    has_cell = ~((t.cell_id[:, 0] == CELL_ID_SENTINEL) & (t.cell_id[:, 1] == CELL_ID_SENTINEL))
    cell_id_str = np.full(t.number_rnas, "", dtype=object)
    if has_cell.any():
        cell_id_str[has_cell] = cell_id_str_from_prefix_suffix_uint32(
            t.cell_id[has_cell, 0], t.cell_id[has_cell, 1]
        )

    df = pd.DataFrame(
        {
            "x": t.x_position,
            "y": t.y_position,
            "z": t.z_position,
            "feature_name": pd.Categorical(gene_names_per_row),
            "quality_score": t.quality_score.astype(np.float32),
            "cell_id": cell_id_str,
        }
    )
    transform = Scale([1.0 / pixel_size, 1.0 / pixel_size], axes=("x", "y"))
    return PointsModel.parse(
        df,
        coordinates={"x": "x", "y": "y", "z": "z"},
        feature_key="feature_name",
        instance_key="cell_id",
        transformations={"global": transform},
    )


def sync_table_and_boundaries(sdata) -> None:
    """Mirrors subset2zarr.py step 5: keep only cells present in both the table
    and the cell_boundaries GeoDataFrame. Expected to be a no-op here (Stage 1
    already picked one exact, canonical `keep_cell_ids` set applied consistently
    to every cropped file), but kept as cheap defensive parity with the original
    script -- it protects against any edge case where the manual `build_table`/
    `build_points` construction above disagrees with `cells.zarr.zip`.
    """
    sdata["table"].obs["region"] = "cell_boundaries"
    sdata.set_table_annotates_spatialelement("table", region="cell_boundaries", instance_key="cell_id")

    boundary_ids = set(sdata["cell_boundaries"].index.astype(str))
    obs_cell_ids = sdata["table"].obs["cell_id"].astype(str)
    mask = obs_cell_ids.isin(boundary_ids)
    sdata["table"] = sdata["table"][mask].copy()

    common_ids = list(sdata["table"].obs["cell_id"].astype(str))
    sdata["cell_boundaries"] = sdata["cell_boundaries"].loc[
        sdata["cell_boundaries"].index.astype(str).isin(common_ids)
    ]
    if "nucleus_boundaries" in sdata.shapes:
        sdata["nucleus_boundaries"] = sdata["nucleus_boundaries"].loc[
            sdata["nucleus_boundaries"]["cell_id"].astype(str).isin(common_ids)
        ]
    return common_ids


def patch_manifest(input_dir: str, output_dir: str, num_cells: int, transcripts_per_cell: int,
                    num_transcripts: int, num_transcripts_high_quality: int, written_files: dict) -> None:
    source_path = os.path.join(input_dir, MANIFEST_FILENAME)
    target_path = os.path.join(output_dir, MANIFEST_FILENAME)

    source_data = {}
    if os.path.exists(source_path):
        with open(source_path) as f:
            source_data = json.load(f)

    extracted = {k: source_data[k] for k in MANIFEST_KEYS_TO_COPY if k in source_data}

    target_data = dict(extracted)
    target_data.update(
        {
            "num_cells": num_cells,
            "transcripts_per_cell": transcripts_per_cell,
            "num_transcripts": num_transcripts,
            "num_transcripts_high_quality": num_transcripts_high_quality,
        }
    )

    explorer_files = {TRANSCRIPTS_ZARR_KEY: "transcripts.zarr.zip", CELLS_ZARR_KEY: "cells.zarr.zip",
                      CELL_FEATURE_ZARR_KEY: "cell_feature_matrix.zarr.zip"}
    if "csc_cell_feature_matrix" in written_files:
        explorer_files[CELL_FEATURE_VIZ_ZARR_KEY] = "csc_cell_feature_matrix.zarr.zip"
    if "binned_transcripts" in written_files:
        explorer_files[TRANSCRIPTS_VIZ_ZARR_KEY] = "binned_transcripts.zarr.zip"
    target_data[EXPLORER_FILES_KEY] = explorer_files

    images = {}
    if "morphology_2d_rel" in written_files:
        images[MORPHOLOGY_2D_KEY] = written_files["morphology_2d_rel"]
    if "morphology_3d_rel" in written_files:
        images[MORPHOLOGY_3D_KEY] = written_files["morphology_3d_rel"]
    target_data[IMAGES_KEY] = images

    with open(target_path, "w") as f:
        json.dump(target_data, f, indent=4)


# --------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description="Subset an Atera output bundle using a polygon.")
    parser.add_argument("-i", "--input", required=True, help="Path to input Atera bundle directory")
    parser.add_argument("-p", "--polygon", required=True, help="Path to selection GeoJSON")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("-z", "--zarr_out", help="Optional: path to save the subsetted SpatialData .zarr")
    parser.add_argument(
        "--polygon-units", required=True, choices=["microns", "pixels"],
        help=(
            "Units the polygon GeoJSON's coordinates are in. Required (no default) because "
            "this is silently wrong in a way that doesn't error -- a pixel-space polygon's "
            "bbox can look like a perfectly plausible micron-space region, and you just get "
            "zero (or the wrong) cells selected instead of a loud failure. Everything else in "
            "this bundle (cell centroids, transcript positions) is in microns; pass 'pixels' "
            "if your polygon was drawn/exported against the morphology image instead."
        ),
    )
    parser.add_argument("--dry-run", action="store_true",
                         help="Run Stage 0+1 only: report how many cells/tiles would be kept, then exit.")
    parser.add_argument("--skip-binned-transcripts", dest="skip_binned_transcripts", action="store_true",
                         default=False,
                         help="Skip binned_transcripts.zarr.zip entirely (default: False -- it is cropped "
                              "and included by default so the output bundle is visualizable in ziggy's "
                              "density view; the crop streams one gene at a time, so this is cheap).")
    parser.add_argument("--include-binned-transcripts", dest="skip_binned_transcripts", action="store_false",
                         help="Crop binned_transcripts.zarr.zip from the original file (this is the default; "
                              "this flag exists to override an earlier --skip-binned-transcripts).")
    parser.add_argument("--keep-tmp", action="store_true",
                         help="Do not delete the Stage-2 intermediate directory afterward.")
    parser.add_argument("--tmp-dir", default=None,
                         help="Deterministic Stage-2 intermediate directory (default: <output>/.stage2_tmp).")
    parser.add_argument("--force-redo", dest="force_redo", action="store_true",
                         help="Re-run Stage 2 crops even if their tmp-dir output already exists.")
    args = parser.parse_args()

    log_checkpoint("Stage 0: starting")

    with open(os.path.join(args.input, MANIFEST_FILENAME)) as f:
        manifest = json.load(f)
    pixel_size = float(manifest["pixel_size"])

    target_polygon = load_polygon(args.polygon, units=args.polygon_units, pixel_size=pixel_size)
    bbox = target_polygon.bounds  # (min_x, min_y, max_x, max_y) -- always microns from here on

    cells_zarr_path = os.path.join(args.input, manifest[EXPLORER_FILES_KEY][CELLS_ZARR_KEY])
    keep_rows, n_total, n_kept = select_cells(cells_zarr_path, target_polygon)
    log_checkpoint(f"Stage 1: {n_kept}/{n_total} cells kept (exact centroid-in-polygon test)")

    if args.dry_run:
        grid_attrs = adt_transcripts.read_grid_attrs(
            os.path.join(args.input, manifest[EXPLORER_FILES_KEY][TRANSCRIPTS_ZARR_KEY])
        )
        grid_size = float(grid_attrs["grid_size"])
        keep_tiles = adt_crop._tiles_in_bbox(grid_attrs["grid_keys"], grid_size, bbox)
        print(f"[dry-run] polygon bbox: {bbox}")
        print(f"[dry-run] cells: {n_kept}/{n_total} would be kept")
        print(f"[dry-run] transcript tiles: {len(keep_tiles)}/{len(grid_attrs['grid_keys'])} would be kept "
              f"(tile-level bbox prefilter, grid_size={grid_size})")
        print("[dry-run] exiting before Stage 2 (no large files touched).")
        return

    if n_kept == 0:
        raise SystemExit("No cells fall inside the given polygon; refusing to produce an empty bundle.")

    tmp_dir = args.tmp_dir or os.path.join(args.output, ".stage2_tmp")
    os.makedirs(args.output, exist_ok=True)
    stage2_files = stage2_crop(args.input, manifest, bbox, keep_rows, tmp_dir, args)

    # spatialdata_io.atera() needs its own experiment.spatial manifest inside
    # tmp_dir (it only reads `pixel_size`, `explorer_files.cells_zarr_filepath`,
    # and `images.*` for the elements requested below -- transcripts/cell_feature
    # are disabled and their manifest keys are never consulted).
    tmp_manifest = {
        "pixel_size": pixel_size,
        "explorer_files": {CELLS_ZARR_KEY: "cells.zarr.zip"},
        "images": {},
    }
    if "morphology_2d_rel" in stage2_files:
        tmp_manifest["images"][MORPHOLOGY_2D_KEY] = stage2_files["morphology_2d_rel"]
    if "morphology_3d_rel" in stage2_files:
        tmp_manifest["images"][MORPHOLOGY_3D_KEY] = stage2_files["morphology_3d_rel"]
    with open(os.path.join(tmp_dir, MANIFEST_FILENAME), "w") as f:
        json.dump(tmp_manifest, f)

    # -- Stage 3a: SpatialData for the schema-compatible parts only --
    from spatialdata_io.readers.atera import atera as read_atera

    # `scale_factors=None` skips building a multiscale image/labels pyramid: the
    # whole point of this pipeline is that the cropped ROI is already small, so a
    # multi-resolution pyramid is unlikely to be worth it (and a tiny crop can be
    # too small to downsample at all, which would otherwise raise inside
    # `multiscale_spatial_image`).
    sdata = read_atera(
        tmp_dir,
        cells_boundaries=True,
        nucleus_boundaries=True,
        cells_as_circles=False,
        cells_labels=True,
        nucleus_labels=True,
        cell_feature_matrix=False,
        transcripts=False,
        morphology_2d="morphology_2d" in stage2_files,
        morphology_3d="morphology_3d" in stage2_files,
        image_models_kwargs={"scale_factors": None},
        labels_models_kwargs={"scale_factors": None},
    )
    log_checkpoint("Stage 3a: SpatialData built (images/labels/shapes only)")

    # -- Stage 3b: exact polygon filter on transcripts (cheap: already tile-cropped) --
    filtered_transcripts, gene_names_per_row = filter_transcripts_to_polygon(
        stage2_files["transcripts"], target_polygon
    )
    sdata.points["transcripts"] = build_points(filtered_transcripts, gene_names_per_row, pixel_size)
    log_checkpoint(
        f"Stage 3b: transcripts exact-filtered to polygon "
        f"({filtered_transcripts.number_rnas} kept)"
    )

    # -- Stage 3c: table + sync against boundaries --
    sdata.tables["table"] = build_table(stage2_files["cell_feature_matrix"], stage2_files["cells"])
    common_ids = sync_table_and_boundaries(sdata)
    log_checkpoint(f"Stage 3c: table synced against boundaries ({len(common_ids)} cells)")

    # -- Rechunk images to 1024x1024 (parity with subset2zarr.py step 4) --
    # With `scale_factors=None` (used above, since our crops are already small
    # and may be too small to build any pyramid level from at all), a parsed
    # image can come back as a single-scale plain DataArray instead of a
    # multiscale DataTree/Dataset -- handle both shapes.
    if sdata.images:
        for img_key, data_tree in sdata.images.items():
            if hasattr(data_tree, "keys"):
                for scale_key in data_tree.keys():
                    da = data_tree[scale_key]["image"]
                    data_tree[scale_key]["image"] = da.chunk({"c": -1, "y": 1024, "x": 1024})
            else:
                sdata.images[img_key] = data_tree.chunk({"c": -1, "y": 1024, "x": 1024})

    # -- Stage 3d: write the final bundle --
    os.makedirs(args.output, exist_ok=True)

    # cells.zarr.zip / cell_feature_matrix(+csc): re-crop only if sync dropped
    # cells further than Stage 1 already did; otherwise a plain copy.
    if len(common_ids) < n_kept:
        print(f"  [sync] table sync dropped {n_kept - len(common_ids)} more cell(s); re-cropping to match")
        # `common_ids` are cell_feature_matrix barcode strings (see build_table);
        # translate to row positions within the (already small) Stage-2-cropped
        # cells.zarr.zip via the same row-alignment assumption used throughout,
        # rather than assuming `common_ids` parse as crop_cells_by_ids tokens.
        stage2_cfm_obs = adt_cfm.read_obs_index(stage2_files["cell_feature_matrix"])
        common_ids_set = set(common_ids)
        keep_rows_stage2 = [str(i) for i, bc in enumerate(stage2_cfm_obs) if str(bc) in common_ids_set]
        adt_crop.crop_cells_by_ids(
            stage2_files["cells"], os.path.join(args.output, "cells.zarr.zip"), keep_rows_stage2
        )
        adt_crop.subset_cell_feature_matrix(
            stage2_files["cell_feature_matrix"], os.path.join(args.output, "cell_feature_matrix.zarr.zip"),
            keep_cells=common_ids,
        )
        if "csc_cell_feature_matrix" in stage2_files:
            adt_crop.subset_cell_feature_matrix(
                stage2_files["csc_cell_feature_matrix"],
                os.path.join(args.output, "csc_cell_feature_matrix.zarr.zip"),
                keep_cells=common_ids,
            )
    else:
        shutil.copy2(stage2_files["cells"], os.path.join(args.output, "cells.zarr.zip"))
        shutil.copy2(stage2_files["cell_feature_matrix"], os.path.join(args.output, "cell_feature_matrix.zarr.zip"))
        if "csc_cell_feature_matrix" in stage2_files:
            shutil.copy2(
                stage2_files["csc_cell_feature_matrix"],
                os.path.join(args.output, "csc_cell_feature_matrix.zarr.zip"),
            )

    # transcripts.zarr.zip: write the exact-polygon-filtered Transcripts directly.
    grid_attrs = adt_transcripts.read_grid_attrs(stage2_files["transcripts"])
    adt_transcripts.write(
        os.path.join(args.output, "transcripts.zarr.zip"), filtered_transcripts,
        grid_size=float(grid_attrs["grid_size"]),
    )

    written_files = dict(stage2_files)
    for key in ("morphology_2d", "morphology_3d"):
        if key in stage2_files:
            rel = stage2_files[key + "_rel"]
            dst = os.path.join(args.output, rel)
            os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
            shutil.copy2(stage2_files[key], dst)

    if "binned_transcripts" in stage2_files:
        shutil.copy2(stage2_files["binned_transcripts"], os.path.join(args.output, "binned_transcripts.zarr.zip"))

    log_checkpoint("Stage 3d: final bundle files written")

    # -- Stage 3e: patch experiment.spatial --
    table = sdata["table"]
    num_cells = table.shape[0]
    transcripts_per_cell = int(np.median(np.asarray(table.X.sum(axis=1)).ravel())) if num_cells else 0
    num_transcripts = filtered_transcripts.number_rnas
    num_transcripts_high_quality = int((filtered_transcripts.quality_score >= 20).sum())

    patch_manifest(
        args.input, args.output, num_cells, transcripts_per_cell, num_transcripts,
        num_transcripts_high_quality, written_files,
    )
    log_checkpoint("Stage 3e: experiment.spatial patched")

    # -- Stage 3f: optional full SpatialData zarr dump --
    if args.zarr_out:
        sdata.write(args.zarr_out, overwrite=True)
        log_checkpoint(f"Stage 3f: full SpatialData zarr written to {args.zarr_out}")

    if not args.keep_tmp:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        print(f"  [cleanup] removed {tmp_dir}")

    print("Done!")


if __name__ == "__main__":
    main()
