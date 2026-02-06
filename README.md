## Maintained By

**[Sebastiano Montante, PhD](https://github.com/semontante)**  
📧 Email: [s.montante@exeter.ac.uk](mailto:s.montante@exeter.ac.uk)  

# Introduction

Welcome to the AlgaeScan github repository! 

![Project Logo](https://github.com/exeter-cytomics/AlgaeScan/raw/main/intro_img/logo_tool.png) 
![License](https://img.shields.io/github/license/exeter-cytomics/AlgaeScan) 
![Issues](https://img.shields.io/github/issues/exeter-cytomics/AlgaeScan)
[![PDF](https://img.shields.io/badge/AlgaeScan_manual-PDF-red)](https://github.com/semontante/flowMagic/blob/main/pk_manual/AlgaeScan_manual.pdf)
[![PDF](https://img.shields.io/badge/AlgaeScan_vignette-PDF-red)](https://github.com/semontante/flowMagic/blob/main/pk_manual/AlgaeScan_vignette.pdf)

Algae species represent crucial indicators of ecosystem health and integrity, as they form the 
foundation of most aquatic food chains. Nearly all aquatic animals rely on these primary producers. 
In addition, the composition and abundance of algae species provide insights into water quality, 
as different species thrive under specific chemical and physical conditions.
Traditional manual identification of algae species is time-consuming and requires specialized expertise. 

This repository stores our fast, customizable and user-friendly pipeline (i.e., with no complex parameters tuning) designed to automatically identify 
algae species in biological samples analyzed through the flow cytometry technology.

# Overview

The AlgaeScan pipeline is composed of two machine learning models trained on spectral flow cytometry data to perform distinct tasks:

    1) Model A: ML model trained to separate algae vs non-algae data.
    2) Model B: ML model trained to identify algae species present in the training set.

# Installation

The pipeline can be installed within R using the devtools library:

```R
library(devtools)
install_github("exeter-cytomics/AlgaeScan")
```

Alternatively, the user can download the AlgaeScan pipeline and install it from their local directory:

```R
install.packages("path/to/AlgaeScan.tar.gz",repos=NULL,type="source")
```

# Test script

The script below can be used to test the correct installation of the AlgaeScan package.





