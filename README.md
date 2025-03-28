# MassLite

MassLite v1.0

Author: Zhu Zou 

Email: xyqzouzhu@gmail.com

Description:
-------------
MassLite is a GUI application for processing mass spectrometry data files, specifically compatible with SCMS data without chromatographic separation. It enables peak alignment, cell segmentation, background removal, 
and export of structured intensity tables.

Main Functionalities:
- Read SCMS data from imzML/mzML files
- Group cell scans by biomarkers
- Perform peak alignment with ppm-level tolerance without binning
- Normalize and filter spectra
- Export scan-wise or cell-wise intensity tables, aligned or unaligned

Dependencies:
- numpy, pandas, matplotlib, seaborn, tkinter, pyimzML, pymzml, umap-learn, sklearn, scipy
