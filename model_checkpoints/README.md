# Model Checkpoints

This directory contains the Mega/CNN checkpoint archive supplied with the project snapshot.

These files are included so trained weights are not lost, but the preferred reproducibility path for the background-matched audit is still to use locked isolate-level prediction CSVs when available. Prediction CSVs allow reviewers to rerun the audit without reloading the neural network.

## Checkpoint Set

Path: `model_checkpoints/mega_cnn_archive_2026-04-22/`

| File | Approx size | SHA256 |
| --- | ---: | --- |
| `maldi_amr_seed0.pt` | 3.7 MB | `7e72071ceeface0563be5f1df02e66a48b401735a9464f4b9871abe95c903599` |
| `maldi_amr_seed1.pt` | 3.7 MB | `e1400061d2c6f31f2dc94fb6020ecd1e5eb092fb69ff287b145ec4c1ed7508a7` |
| `maldi_amr_seed2.pt` | 3.7 MB | `2830d8367c28b673fe39a97c5acb9a491638080d0b3aef93f5d553f39521cc1e` |
| `maldi_amr_seed3.pt` | 3.7 MB | `e31a04ae280455e931a35b2ed9482a5288ea7b45c18c8719fbfb9b78a0c23dc6` |
| `maldi_amr_seed4.pt` | 3.7 MB | `a0a4fc4bca3b00428a64e49226cdf6bf021b1f9323889c2ebb9ad9658c2dc565` |
| `maldi_mae_pretrained.pt` | 16 MB | `f2c18c8b35841a372fdfc171474ade72f47b24bb7ba81cc80a58ee9092c906b9` |

## Notes

- The five `maldi_amr_seed*.pt` files are the AMR ensemble seed checkpoints.
- `maldi_mae_pretrained.pt` is the masked autoencoder pretraining checkpoint.
- These are binary artifacts; for a final public paper release, consider mirroring them to Zenodo or OSF and keeping GitHub as the code/documentation repository.
