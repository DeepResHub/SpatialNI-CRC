# Minimal WSI-to-SpatialNI Example

This example follows the repository order:

1. WSI preprocessing
2. patch classification
3. SpatialNI feature construction
4. WSI patch class map export

By default, `classify_wsi_with_tiatoolbox()` performs Reinhard stain
normalization. You must provide a reference H&E patch image through
`reference_image`.

```python
from pathlib import Path

from patch_classification import classify_wsi_with_tiatoolbox
from spatialni_feature_engine import (
    build_patient_spatialni_from_predictions,
    plot_wsi_patch_class_map,
)

project_dir = Path(__file__).resolve().parent

wsi_path = project_dir / "wsi" / "TCGA-A6-2676-01Z-00-DX1.c465f6e0-b47c-48e9-bdb1-67077bb16c67.svs"
reference_image = project_dir / "outputs" / "TCGA-A6-2676-01" / "reference_patch.png"
output_dir = project_dir / "outputs" / "TCGA-A6-2676-01"
output_dir.mkdir(parents=True, exist_ok=True)

predictions = classify_wsi_with_tiatoolbox(
    wsi_path=wsi_path,
    output_csv=output_dir / "patch_predictions.csv",
    resolution=0.5,
    units="mpp",
    batch_size=16,
    device="cuda",
    reference_image=reference_image,
)

spatialni = build_patient_spatialni_from_predictions(
    predictions,
    patient_id="TCGA-A6-2676-01",
)
spatialni.to_csv(output_dir / "patient_spatialni.csv", index=False)

plot_wsi_patch_class_map(
    predictions,
    output_png=output_dir / "wsi_patch_class_map.png",
    wsi_path=wsi_path,
    patch_resolution=0.5,
    overview_resolution=8.0,
    units="mpp",
    title="TCGA-A6-2676-01 patch class map",
)
```

If no stain-normalization reference image is available and you only want to test
that the example slide can run end-to-end, explicitly disable normalization:

```python
predictions = classify_wsi_with_tiatoolbox(
    wsi_path=wsi_path,
    output_csv=output_dir / "patch_predictions.csv",
    resolution=0.5,
    units="mpp",
    batch_size=16,
    device="cuda",
    stain_normalization=False,
)
```

## Physical resolution and scope

All new patch extraction uses the native WSI MPP to resample to 0.5 micrometres/pixel.
A 512-pixel patch spans 256 micrometres; the full 3 x 3 window spans 768 micrometres per side.
Coordinates are expressed in the resampled patch grid. Supply a stain reference image for the standard pipeline.
The updated example outputs are generated at 0.5 MPP. Rerunning the example generates run_manifest.json, which records the software versions and file checksums for that run.
This public example generates features and a tissue map from user-supplied WSI data.
Private patient-level clinical data are not distributed, so the repository does not directly regenerate the manuscript clinical figures.

## Reproduce the bundled full-slide example

From the repository root, run:

```powershell
$env:NUMBA_DISABLE_INTEL_SVML = "1"
python run_public_example.py --wsi "wsi/TCGA-A6-2676-01Z-00-DX1.c465f6e0-b47c-48e9-bdb1-67077bb16c67.svs"
```

Use `--device cpu` for CPU inference. Numerical differences across devices and
software versions are possible. The runner derives a deterministic stain-reference
patch from the middle row of the Otsu-filtered coordinate table and saves it as
`reference_patch.png`. This is a public software demonstration reference, not the
manuscript stain reference. The runner processes all selected tissue patches,
exports predictions, patient SpatialNI/composite scores and the tissue map, and
records input, output and code checksums in `run_manifest.json`.
The saved CSV uses a portable reference filename, without local absolute paths.
