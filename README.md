# ptv-binning-codes
Object-oriented Matlab framework for the ensemble averaging of Particle Tracking Velocimetry (PTV) data

## Motivation
This repository contains the codes that I used to obtain ensemble averages of Particle Tracking Velocimetry (PTV) data for my [MSc Thesis](https://repository.tudelft.nl/islandora/object/uuid%3A51c6c715-43e3-4cc7-9808-be6d62927d6f?collection=education), which was later published as a [journal paper](https://iopscience.iop.org/article/10.1088/1361-6501/ac41dd/meta).

The original codes were developed for the following papers, please cite them if you use the codes of this repository to time-average your data:
- Jan F G Schneiders, Fulvio Scarano, Constantin Jux and Andrea Sciacchitano. "Coaxial volumetric velocimetry", 2018, [https://doi.org/10.1088/1361-6501/aab07d](https://doi.org/10.1088/1361-6501/aab07d)
- Constantin Jux, Andrea Sciacchitano, Jan F. G. Schneiders and Fulvio Scarano. "Robotic volumetric PIV of a full-scale cyclist", 2018, [https://doi.org/10.1088/1361-6501/aab07d](https://doi.org/10.1007/s00348-018-2524-1)

During my MSc Thesis I rewrote the codes to implement phase-averaging capabilities, an object oriented structure, and a consistent naming convention for the variables. My MSc Thesis work ended up in the following paper, please cite it if you use the codes to phase-average your data:
- Francesco M A Mitrotta, Jurij Sodja and Andrea Sciacchitano. "On the combined flow and structural measurements via robotic volumetric PTV", 2022, [https://doi.org/10.1088/1361-6501/ac41dd](https://doi.org/10.1088/1361-6501/ac41dd)

## Installation
You need to have [Matlab](https://www.mathworks.com/products/matlab.html) already installed in your system to use the codes of this repository.

1. Download the repository to a local folder (e.g. ~/ptv-binning-codes/) by running: 
```console
git clone https://github.com/fmamitrotta/ptv-binning-codes.git
```
2. Run Matlab and add the folder (~/ptv-binning-codes/) to your Matlab path.

## Usage
The folder (~/ptv-binning-codes/Binning/) contains all the essential classes, functions and script for the ensemble averaging of your PTV data. In such folder, the script `StbBinnerMain.m` provides a skeleton that you can modify to your needs and run to perform your ensemble average.

The folder (~/ptv-binning-codes/MitrottaThesis/) includes some additional codes and examples of how I used the binning codes to perform an ensemble average of my data in a phase-average fashion.

## Contributing
Please don't hesistate to throw feedback and suggestions. Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[GPL-3.0](https://choosealicense.com/licenses/gpl-3.0/)
