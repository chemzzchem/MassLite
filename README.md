# MassLite

**MassLite v1.0**  

Author: Zhu Zou  

Email: xyqzouzhu@gmail.com

MassLite is a GUI-based Python application for processing **Single Cell Mass Spectrometry (SCMS)** data, designed specifically for datasets acquired without chromatographic separation. It enables advanced peak alignment, biomarker-driven scan grouping, background signal removal, and export of structured intensity tables — all in an intuitive graphical environment.

---

## Citation

If you use MassLite in your research, please cite:

> Zou, Z.; Peng, Z.; Bhusal, D.; Wije Munige, S.; Yang, Z.  
> **MassLite: An Integrated Python Platform for Single Cell Mass Spectrometry Metabolomics Data Pretreatment with Graphical User Interface and Advanced Peak Alignment Method**.  
> _Analytica Chimica Acta_ (2024). [https://doi.org/10.1016/j.aca.2024.343124](https://doi.org/10.1016/j.aca.2024.343124)

---

## Key Features

- Load `.imzML` and `.mzML` files in both **profile** and **centroid** mode (profile data is automatically centroided)
- Group scans by cell-specific biomarkers (e.g., PC(34:1))
- Perform **binning-free** peak alignment with ppm-level resolution
- Detect and remove background peaks using scan-wise TIC dynamics
- Filter and normalize peaks using user-defined thresholds
- Export results as aligned/unfiltered CSV files (cell-wise or scan-wise)

---

## Dependencies

MassLite is written in **Python 3.11** and requires the following packages:

- `numpy`, `pandas`, `matplotlib`, `seaborn`
- `tkinter` 
- `pyimzML`, `pymzml`
- `scipy`, `sklearn`, `umap-learn` 

Usage:
-------------
Download the public Python file and run. Follow the instructions on the GUI. The workflow starts from top to bottom, left to right. The code was generated with Python 3.11 in VSCode environment.

Detail:
-------------
**File Input** from mzML and imzML formats of data using pymzml and pyimzml package, both in profile and centroid settings. Profile data will be converted into centroid data during the analysis.

******

**Void Scan Filter** addresses "void scans" in SCMS experiments when data was acquired without sample injection. MassLite uses unsupervised machine learning (K-means, UMAP) to effectively identify these scans, thus exclude them to improve downstream accuracy.

*****

**Cell Scan Grouper** can automatically detect and segment individual single cells by user-defined biomarker, which can be an abundant lipid PC(34:1) commonly observed in cells by default without specification. The result are visualized through Extracted Ion Chromatograms(EICs).

*****

**Peak Alignment** is implemented with a hierarchical clustering-based approach and two key innovations. First, the method adopted is binning-free to avoid intrinsic drawbacks from traditional binning method. It is independent of order of files and not biased by the position or size of the bins. Second, an algebraic transformation is included to convert absolute mass difference in Dalton into relative mass difference in ppm. This transformation enabled quick update of the distance matrix for hierarchical clustering when performing the peak alignment task.

*****

**Background Detection** can be performed along with cell scan grouper after alignment. Peaks with their maximum detected among one of the background scans will be redeemed as background signal and can be substantially removed from further analysis. This is especially useful in noisy or improvised acquisition scenarios.

*****

**Result Filter** and **Exportation** can be done with user-defined parameters. The generated CSV can be used as input for downstream analysis such as t-test, ANOVA, random forest, etc.

License:
-------------
MassLite is licensed under Apache-2.0 license. 

The code is provided as-is, with no warranty.
