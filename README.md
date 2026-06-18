# aquacropgrid-preproc
This work-in-progress repository creates an automatic data download, preprocessing, and harmonisation pipeline for global gridded datasets to serve as inputs to agricultural crop and irrigation models. More specifically, this covers all necessary data to run FAO's AquaCrop model over large regions (e.g. countries or river basins) in gridded format: 
- Climate inputs (precipitation, temperature, reference evapotranspiration) in daily time step, both for past (from AgERA5) and future climate (from NASA NEX)
- Soil properties (from ISRIC Soilgrids)
- Planting and harvesting calendars (from GGCMI)
- Crop areas (from MAPSPAM)

Follow the steps below to set up and run the project locally on your machine.

**Prerequisites**

Conda: You'll need either Miniconda or Anaconda installed. Miniconda is generally recommended for a lighter installation.


**Steps to download and preprocess input datasets from global data sources (climate, soil, crop areas, planting and harvesting calendars)**
1. In case you would like to prepare past (1979-recent) climate data, go to Copernicus Climate Data Store website (https://cds.climate.copernicus.eu/) and do the following:
  - Create a CDS account (register)
  - Access your API Token (Your profile -> API Token). You will need the token as input argument for running the script (see step 5 below)
  - When logged into your CDS account, accept “Terms of use” at the bottom of the following page:
      - AgERA5: https://cds.climate.copernicus.eu/datasets/sis-agrometeorological-indicators?tab=download
2. Clone this repository:
  - From command line (recommended), clone this repository to your local machine using Git (make sure you have Git installed on your system; you can download it from git-scm.com):
      - git clone https://github.com/josiasritter/aquacropgrid-preproc 
  - Alternatively, manually download this repository to the desired location on your computer
3. Navigate working directory to the newly created project directory (via command line):
      - cd aquacropgrid-preproc
4. Create a new conda environment and install the required packages (via command line):
      - conda env create -f environment.yml
      - conda activate aquacropgrid-preproc
5. Run script “preproc_main.py”:
  - Manually adjust the input arguments at the top of the script and save the changes
  - Run in command line: python preproc_main.py
