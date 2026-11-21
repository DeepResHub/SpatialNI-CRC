"""Run the complete public WSI example at 0.5 micrometres/pixel."""
import argparse
import hashlib
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import platform
import time

import numpy as np
from PIL import Image


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wsi', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path(__file__).parent/'outputs/TCGA-A6-2676-01')
    parser.add_argument('--device', default='cuda')
    args = parser.parse_args()
    from wsi_preprocessing import open_wsi, notebook_style_patch_locations, read_patches_from_locations
    from patch_classification import classify_wsi_with_tiatoolbox
    from spatialni_feature_engine import build_patient_spatialni_from_predictions, plot_wsi_patch_class_map

    start = time.time()
    args.output.mkdir(parents=True, exist_ok=True)
    reader = open_wsi(str(args.wsi))
    locations = notebook_style_patch_locations(args.wsi, resolution=.5, units='mpp')
    # Deterministic, public-demo-only stain reference; not the manuscript reference.
    ref_index = len(locations)//2
    ref_row = locations.iloc[[ref_index]]
    image = read_patches_from_locations(str(args.wsi), ref_row, resolution=.5, units='mpp')[0]
    reference = args.output/'reference_patch.png'
    Image.fromarray(image).save(reference)
    predictions = classify_wsi_with_tiatoolbox(
        args.wsi, reference_image=reference, resolution=.5, units='mpp',
        batch_size=16, device=args.device)
    assert len(predictions)==len(locations)
    assert np.array_equal(predictions[['x','y']].values, locations[['x','y']].values)
    assert not predictions.duplicated(['x','y']).any()
    predictions['reference_image'] = 'reference_patch.png'
    predictions.to_csv(args.output/'patch_predictions.csv', index=False)
    features = build_patient_spatialni_from_predictions(predictions, patient_id='TCGA-A6-2676-01')
    features.to_csv(args.output/'patient_spatialni.csv', index=False)
    plot_wsi_patch_class_map(predictions, args.output/'wsi_patch_class_map.png', wsi_path=args.wsi,
                            patch_resolution=.5, overview_resolution=8., units='mpp')
    weights=Path.home()/'.tiatoolbox/models/wide_resnet101_2-kather100k.pth'
    manifest = {
        'wsi_filename':args.wsi.name, 'wsi_sha256':sha256(args.wsi),
        'native_mpp':reader.info.mpp.tolist(), 'target_mpp':.5, 'units':'mpp',
        'patch_size':[512,512], 'stride':[512,512], 'patch_count':len(predictions),
        'mask':'otsu', 'stain_normalization':'Reinhard', 'reference_file':'reference_patch.png',
        'reference_origin':{'grid_row_zero_based':ref_index,'x':int(ref_row.x.iloc[0]),'y':int(ref_row.y.iloc[0]),
                            'scope':'Public software demonstration only; not the manuscript stain reference.'},
        'model':'wide_resnet101_2-kather100k', 'model_sha256':sha256(weights),
        'device':args.device, 'python':platform.python_version(),
        'versions':{n:metadata.version(n) for n in ['torch','torchvision','tiatoolbox','numpy','pandas','zarr','numcodecs']},
        'NUMBA_DISABLE_INTEL_SVML':os.environ.get('NUMBA_DISABLE_INTEL_SVML'),
        'elapsed_seconds':round(time.time()-start,2),
        'output_sha256':{n:sha256(args.output/n) for n in ['reference_patch.png','patch_predictions.csv','patient_spatialni.csv','wsi_patch_class_map.png']},
        'code_sha256':{n:sha256(Path(__file__).parent/n) for n in ['run_public_example.py','patch_classification.py','wsi_preprocessing.py','spatialni_feature_engine.py']},
        'scope':'Full example WSI inference in the recorded development environment; not validation of a fresh environment or reproduction of private clinical figures.'}
    (args.output/'run_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(manifest,indent=2),flush=True)


if __name__=='__main__':
    main()
