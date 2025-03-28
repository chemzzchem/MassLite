# MassLite

MassLite v1.0

Author: Zhu Zou 

Email: xyqzouzhu@gmail.com

Description:
-------------
MassLite is a GUI application for processing mass spectrometry data files, specifically compatible with SCMS data without chromatographic separation. It enables peak alignment, cell segmentation, background removal, 
and export of structured intensity tables.

Please cite the following publication if you are using MassLite: 
> Zou, Z.; Peng, Z.; Bhusal, D.; Wije Munige, S.; Yang, Z. **MassLite: An Integrated Python Platform for Single Cell Mass Spectrometry Metabolomics Data Pretreatment with Graphical User Interface and Advanced Peak Alignment Method**, _Analytica Chimica Acta_ (2024). [https://doi.org/10.1016/j.aca.2024.343124]


Main Functionalities:
- Read SCMS data from imzML/mzML files
- Group cell scans by biomarkers
- Perform peak alignment with ppm-level tolerance without binning
- Normalize and filter spectra
- Export scan-wise or cell-wise intensity tables, aligned or unaligned

Dependencies:
- numpy, pandas, matplotlib, seaborn, tkinter, pyimzML, pymzml, umap-learn, sklearn, scipy

Usage:
-------------
Download the public Python file and run. Follow the instructions on the GUI. The workflow starts from top to bottom, left to right. The code was generated with Python 3.11 in VSCode environment.

Detail:
-------------
**File Input** from mzML and imzML formats of data using pymzml and pyimzml package, both in profile and centroid settings. Profile data will be converted into centroid data during the analysis.

******

**Void Scan Filter** is designed for improvised SCMS data acquisition process with void scans when data is collected when no meaningful sample is injected. The unsupervised ML method can quickly target the void scans and exclude them in the follow-up analysis. 

*****

**Cell Scan Grouper** can automatically trace the source of scans into individual cells by the biomarker specified from user, by default an abundant lipid PC(34:1) commonly observed in cells. The result can be visualized through EIC.

*****

**Peak Alignment** is done with two advancement. First, the method adopted is binning-free so that peaks are truly comparing among each other regarding their similarity, which is m/z values in this case. The method is independent of order of files and not biased by the position or size of the bins. Second, an algebraic transformation is included to convert absolute mass difference in Da into relative mass difference in ppm. This transformation enabled quick update of the distance matrix for hierarchical clustering when performing peak alignment tasks.

*****

**Background Detection** is a byproduct of cell scan grouper after alignment. By setting up the threshold, all peaks with their maximum detected among one of the background scans will be redeemed as background signal and can be substantially removed from further analysis. The method is highly effective agains data with a large portion of background scans where background peaks could hardly be all picked out easily, especially useful for SCMS experiment associated with improvised data acquisition process.

*****

**Result Filter** and **Exportation** can be done with options to drop specific peaks, with or without peak alignment, in terms of cells or scans. The generated CSV can be used as input for further analysis such as t-test, ANOVA, random forest, etc.

License:
-------------
MassLite is licensed under Apache-2.0 license. The code is provided as is with no warranty.
