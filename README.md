# CSN_tutorial

This repository contains the Supporting Information code for:

Scalfani, V.F., Patel, V.D. & Fernandez, A.M. Visualizing chemical space networks with RDKit and NetworkX.
*J Cheminform*, **2022**, *14*, 87. https://doi.org/10.1186/s13321-022-00664-x

```bibtex
@article{scalfani2022visualizing,
  title={Visualizing chemical space networks with RDKit and NetworkX},
  author={Scalfani, Vincent F and Patel, Vishank D and Fernandez, Avery M},
  journal={Journal of Cheminformatics},
  volume={14},
  number={1},
  pages={87},
  year={2022},
  publisher={Springer}
}
```

The original Jupyter Notebooks associated with the manuscript are in the *CSN_Jupyter_Notebooks/* folder. The `glucocorticoid_recepter_2034_2.csv` ChEMBL dataset (Additional File 1 in manuscript) is also provided in the *Dataset/* folder. Please read the dataset_license file for the dataset reuse terms.

The *Less_Memory_Calculations/* folder contains an alternative script for the CSN calculations that uses less memory.
This script/method was not part of the original article; we added it here as it was useful to us for running the
calculations on a Raspberry Pi 400 with only 4 GB RAM.

Approximate run times for the `CSN_calculations_lessMem.py` script:

| Hardware                                    | Number of Cores used | Rounded Run Time |
|---------------------------------------------|:--------------------:|-----------------:|
| 12th generation Intel Core i9, 64 GB RAM    | 22                   |    25 min        |
| Raspberry Pi 5, 8 GB RAM                    | 3                    |    3 hours       |
| Raspberry Pi 400, 4 GB RAM                  | 3                    |    5 hours       |

