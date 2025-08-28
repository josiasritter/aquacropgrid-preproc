# aquacropgrid_dev
This repository is work-in-progress for developing code to run the Python library AquaCrop-OSPy in gridded format. Follow these steps to set up the project locally on your machine.

**Prerequisites**

Conda: You'll need either Miniconda or Anaconda installed. Miniconda is generally recommended for a lighter installation.


**Steps to download and preprocess climate data**
1. Go to Copernicus Climate Data Store website (https://cds.climate.copernicus.eu/)
  - Create a CDS account (register)
  - Access your API Token (Your profile -> API Token). You will need the token as input argument for running the script (step 5 below)
  - Accept “Terms of use” at the bottom of the following two pages:
      - ERA5: https://cds.climate.copernicus.eu/datasets/derived-era5-single-levels-daily-statistics?tab=download
      - ERA5-Land https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land?tab=download
2. Clone this code Repository
  - From command line (recommended), clone this repository to your local machine using Git (make sure you have Git installed on your system; you can download it from git-scm.com):
      - git clone https://github.com/josiasritter/aquacropgrid_dev.git
  - Alternatively, manually download this repository to the desired location on your computer
3. Navigate your current working directory to the newly created project directory:
      - cd aquacropgrid_dev
4. Create a new conda environment and install the required packages:
      - conda env create -f environment.yml
      - conda activate aquacropgrid_dev
5. Run script “climate.py”
  - First, manually adjust the input arguments at the top of the script and save the changes (the script will later be refactored into a function “climate” that accepts these input arguments)
  - python preprocessing/climate.py
