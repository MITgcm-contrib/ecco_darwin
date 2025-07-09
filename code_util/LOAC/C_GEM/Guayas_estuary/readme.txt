# Config file and boundary conditions for C-GEM config of the Guayas estuary (Ecuador).
# Base version from code_python_FCO2, modfied with following references:

# Nguyen et al., 2021 https://doi.org/10.1016/j.scitotenv.2021.147261
# Twilley et al., 2001 https://doi.org/10.1007/978-3-662-04482-7_18
# Cifuentes et al., 1996 https://doi.org/10.1006/ecss.1996.0103
# Twilley et al., 1998 https://doi.org/10.1016/S1462-9011(98)00012-4
# Belliard et al., 2022 https://doi.org/10.1016/j.ecss.2022.107766

# and Guayas River centerlines from https://landscape.jpl.nasa.gov/

# First version is averaging conditions between dry and rainy seasons.

####################
#   Instructions   #
####################

git clone https://github.com/MITgcm-contrib/ecco_darwin.git
mkdir C_GEM
cd C_GEM
cp ../ecco_darwin/code_util/LOAC/C_GEM/code_python_FCO2/* .
cp -f ../ecco_darwin/code_util/LOAC/C_GEM/Guayas_estuary/* .

# Work in python virtual environment
# Create python virtual environment (for first time)
python3 -m venv CGEM_env
# Get in python virtual environment
source CGEM_env/bin/activate
# Install libraries (for first time)
pip install numpy numba
# run model
python3 main.py
# Leave virtual environment
deactivate

############################################################
Please, be mindful and do not push model outputs to GitHub repo.
############################################################