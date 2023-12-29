# ADDZYME
**ADDZYME** is a machine learning-based algorithm for prediction of the effect of an additive on the activity of an enzyme. It uses 30 ERT-Baseline models to report relative activity. For more information refer to our published paper.
![addzyme](https://github.com/miladrayka/addzyme/blob/main/addzyme.PNG)

# Contact

Milad Rayka, milad.rayka@yahoo.com

# Citation
Under publishing.

# Installation

Below packages should be installed for using ADDZYME. Dependencies:

- python = 3.8.16

- numpy = 1.23.5

- streamlit = 1.25.0

- molvs = 0.1.1

- joblib = 1.2.0

- rdkit-pypi  (2022/9/5 release)

- scikit-learn = 1.2.2

For installing first make a virtual environment and activate it.

On windows:

>    python py -m venv env
>    .\env\Scripts\activate

On macOS and Linux:

>    python3 -m venv env
>    source env/bin/activate

Which env is the location to create the virtual environment. Now you can install packages:

>    pip install *package_name*==*version"

# Usage

For using ADDZYME, after activating your environment, e.g., env, change your directory to ./codes and type the following command:

>   streamlit run addzyme_gui.py
