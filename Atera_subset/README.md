# Atera subset tool

`subset2atera.py` is the Atera ("Gen2") equivalent of `../Xenium_subset/subset2zarr.py`:
given an Atera output bundle directory + a polygon GeoJSON, it produces a
spatially-subsetted Atera bundle for visualization in ziggy.

Unlike a naive port of the Xenium script, this one is designed around the fact
that a real Atera bundle is large (`transcripts.zarr.zip` can be tens of GB) and
cannot be loaded fully into memory before filtering. See the module docstring at
the top of `subset2atera.py` for the full staged-pipeline design and the schema
caveats discovered while building it (also covered in this repo's PR/task notes).

## Dependencies

This script depends on three sibling packages that (as of writing) are only
available as local branch checkouts, not yet released:

- `atera-dataset-tools` (branch `cells-zarr-support`)
- `spatialdata-io` (branch `atera-support`)
- `sopa` (branch `atera-support`) -- imported only incidentally by `spatialdata_io`;
  `subset2atera.py` itself does not call into `sopa.io.atera` (see the module
  docstring's "IMPORTANT DEVIATION" section for why).

Point `PYTHONPATH` at all three checkouts, e.g.:

```bash
export PYTHONPATH="/path/to/atera-dataset-tools:/path/to/spatialdata-io/src:/path/to/sopa"
```

## Usage

`--polygon-units` is required -- see "Polygon units: microns vs. pixels" below
before running this for the first time against a real bundle.

```bash
python3 subset2atera.py -i <input_bundle_dir> -p <polygon.geojson> -o <output_dir> --polygon-units microns

# Cheap sanity check first (Stage 0+1 only, touches no large files):
python3 subset2atera.py -i <input_bundle_dir> -p <polygon.geojson> -o <output_dir> --polygon-units microns --dry-run

# Keep the Stage-2 intermediates around (for debugging, or to resume Stage 3
# after a crash without re-running the expensive crops):
python3 subset2atera.py -i <input_bundle_dir> -p <polygon.geojson> -o <output_dir> --polygon-units microns --keep-tmp

# Also dump the full assembled SpatialData object as a .zarr:
python3 subset2atera.py -i <input_bundle_dir> -p <polygon.geojson> -o <output_dir> --polygon-units microns -z <out.zarr>
```

Run again with the same `-o` (and without `--force-redo`) to resume: Stage 2
crops whose tmp-dir output already exists are skipped.

### Polygon units: microns vs. pixels

Every spatial comparison in this script (cell centroids, transcript positions,
tile bboxes) is in microns -- that's the unit `cells.zarr.zip`/
`transcripts.zarr.zip` store natively. A polygon drawn or exported against the
*morphology image* (e.g. from an image-space annotation tool) is typically in
**pixels** instead, and a GeoJSON file gives no indication either way -- the
coordinates are just floats.

This is a real trap, not a hypothetical one: the first real-bundle dry run
against this script reported plausible tile-level bbox overlap but exactly
zero cells kept, which turned out to be exactly this -- a pixel-space polygon
compared against micron-space cell centroids, not a location or code bug.
Because a pixel-space bbox over real tissue can easily *look* like a
plausible micron-space region (both are "a few thousand units"), this doesn't
fail loudly on its own -- it just silently selects the wrong (often zero)
cells. That's why `--polygon-units {microns,pixels}` is a required flag
rather than a silently-assumed default: get it wrong and you find out
immediately (zero cells, refused), rather than getting a subtly-wrong output
bundle.

If your polygon is in pixels, pass `--polygon-units pixels`; the script
converts it to microns using the bundle's own `pixel_size` (microns/pixel)
from `experiment.spatial`, before any cell/transcript comparison happens.

## Testing

There is no real ~64GB Atera bundle available in this environment, so
`test_subset2atera_e2e.py` builds a small, structurally-real synthetic bundle
(using `atera_dataset_tools`'s own read/write functions -- the real schema, not
`spatialdata_io`'s own best-effort test fixture, which intentionally follows a
*different*, unconfirmed schema for transcripts/cell-feature-matrix; see that
test's module docstring) and runs the full pipeline against it end-to-end.

```bash
export PYTHONPATH="/path/to/atera-dataset-tools:/path/to/spatialdata-io/src:/path/to/sopa"
python3 -m pytest test_subset2atera_e2e.py -q -s
# or, standalone (no pytest):
python3 test_subset2atera_e2e.py
```

It asserts:
- the output bundle passes `atera_dataset_tools.validate`'s structural checks
  for `cells.zarr.zip`, `transcripts.zarr.zip`, both cell-feature matrices, and
  `binned_transcripts.zarr.zip`;
- `experiment.spatial` is valid JSON with copied-forward + recomputed fields;
- the output actually has fewer cells and fewer transcripts than the input
  (not a no-op), and `binned_transcripts.zarr.zip` has fewer-or-equal tiles
  than the source (not a passthrough copy);
- the cropped morphology image is narrower than the source image.

## Included by default: `binned_transcripts.zarr.zip` and `csc_cell_feature_matrix.zarr.zip`

Both are viz-only/derived -- not needed for structural correctness -- but are
included in the output by default anyway, since the point of this script is a
bundle that's actually visualizable in ziggy, not just a valid one.
`csc_cell_feature_matrix.zarr.zip` was never skippable in the first place (it's
cropped the same way as the CSR matrix, just sequentially after it, to keep
peak memory to one ~6.9GB matrix at a time rather than two). `binned_transcripts
.zarr.zip` (~5.3GB) is cropped via `crop_binned_transcripts_to_bbox`, which as
of this session streams one gene at a time (reads, filters, and writes each
gene's nested-zip archive independently, then discards it) rather than loading
every gene's every tile into memory at once -- the latter was this format's
original crop implementation and would have meant a ~5.3GB peak just to crop
it. Pass `--skip-binned-transcripts` to omit it if the density view isn't
needed for a given use case.

### What this validates vs. what it can't

This sandbox can exercise structural correctness and the full control flow
(dry-run, Stage 1 exact cell selection, Stage 2 crops, Stage 3 SpatialData
assembly + sync + write, resume/`--force-redo`, `--zarr_out`) end-to-end, at a
scale of a dozen cells / ~100 transcripts / tiny images. It cannot validate:
- actual memory behavior at the real ~64GB bundle scale (the memory-safety
  claims for `crop_transcripts_to_bbox`/`subset_cell_feature_matrix`/the
  windowed morphology crop rest on reading those functions' implementations,
  not on having run them against anything close to that size);
- whether a *real* production bundle's `cell_feature_matrix.zarr.zip` obs-index
  strings really do use the same convention as `cells.zarr.zip`'s derived
  hex-letter cell-ID string once matched through `spatialdata_io`'s boundary
  reader (the pipeline is designed to be robust to this either way -- see the
  "cell-ID string conventions" note in `subset2atera.py` -- but this could only
  be truly confirmed against a real bundle);
- whether a real, multi-resolution pyramidal OME-TIFF crops/writes back out
  correctly (the synthetic fixture here only has single-resolution images).
