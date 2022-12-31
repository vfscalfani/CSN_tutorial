# CSN_tutorial

This repository contains the Supporting Information code for:

Scalfani, V.F., Patel, V.D. & Fernandez, A.M. Visualizing chemical space networks with RDKit and NetworkX.
*J Cheminform*, **2022**, *14*, 87. https://doi.org/10.1186/s13321-022-00664-x

The original Jupyter Notebooks associated with the manuscript are in the *CSN_Jupyter_Notebooks/* folder. The `glucocorticoid_recepter_2034_2.csv` ChEMBL dataset (Additional File 1 in manuscript) is also provided in the *Dataset/* folder. Please read the dataset_license file for the dataset reuse terms.

The *Less_Memory_Calculations/* folder contains an alternative script for the CSN calculations that uses less memory.
This script/method was not part of the original article; we added it here as it was useful to us for running the
calculations on a Raspberry Pi 400 with only 4 GB RAM.

The `CSN_calculations_lessMem.py` script takes about 5 hours to run on the Raspberry Pi 400 (using 3 cores)
and about 25 minutes on a 12th generation Intel Core i9 desktop processor (using 22 cores).

