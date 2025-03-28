# MassLite

MassLite v1.0

Author: Zhu Zou
Email: xyqzouzhu@gmail.com
Created for Single Cell Mass Spectrometry data pretreatment
Specifically efficient for Single-probe SCMS data

Description:
-------------
MassLite is a GUI application for processing mass spectrometry data files 
(imzML/mzML). It enables peak alignment, cell segmentation, background removal, 
and export of structured intensity tables.

Main Functionalities:
- Load and display TICs from imzML/mzML files
- Group scans by marker m/z ranges
- Perform peak alignment with ppm-level tolerance
- Normalize and filter spectra
- Export scan-wise or cell-wise intensity tables, aligned or unaligned

Dependencies:
- numpy, pandas, matplotlib, seaborn, tkinter, pyimzML, pymzml, umap-learn, sklearn, scipy
