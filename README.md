# netPCF: Geometry-aware Pair Correlation Functions for Spatial Biology

This repository contains the scripts and datasets required to reproduce the analyses presented in the paper:

**netPCF: Geometry-aware Pair Correlation Functions for Spatial Biology**

The implementation of **netPCF** is available in the **SpaceNet** Python package under the `point_patterns` submodule.

## Installation

To install the latest version of SpaceNet, run:

```bash
pip install spacenet
```

For more information about SpaceNet, visit: https://www.spacenet-python.com

## Repository Contents

* **analysis_scripts**: 
  Python notebooks containing all analyses presented in the manuscript.

* **data**
  Synthetic and processed point cloud datasets used throughout the study.

  * **imc_breast_carcinoma_kuett_2022**:
    A 3D point cloud of a breast carcinoma sample acquired using imaging mass cytometry (IMC), with cell type annotations from Kütt et al. (2022).
    https://doi.org/10.1038/s43018-021-00301-w
  * **murine_embryo_data_murphy_2022**: 
    3D point clouds representing the surfaces of murine embryos at embryonic days E9.5 and E11.5, from Murphy et al. (2022).
    https://doi.org/10.1242/dev.200312
  * **synthetic_data**:
    Collections of synthetic 2D and 3D point cloud datasets generated for this study.

* **misc**: 
  Helper functions used throughout the analysis notebooks.

## Reference

If you use **netPCF**, the datasets, or any part of the analysis in your research, please cite:

Joshua W. Moore, Joshua A. Bull, and Helen M. Byrne. Netpcf: geometry-aware pair correlation functions for spatial biology. bioRxiv, 2026. doi:10.64898/2026.07.02.736020.
