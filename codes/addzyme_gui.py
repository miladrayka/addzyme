
import glob
import joblib
import numpy as np
import streamlit as st
from molvs import standardize_smiles
from rdkit.Chem import Descriptors, MolFromSmiles
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator

st.set_page_config(
    page_title="ADDZYME GUI",
    page_icon="./title_logo.png",
    layout="centered",
    initial_sidebar_state="collapsed",
)

#Title

st.title("ADDZYME")
st.header("Additive Effects on Enzyme Activity Predictor")
st.write("**ADDZYME** is software for predicting effect of an additive on an enzyme activity.")

#Sidebar

st.sidebar.header("About")
st.sidebar.image("./github_logo.png", width=30)
st.sidebar.write("[GitHub](https://github.com/miladrayka), Developed by *[Milad Rayka](https://scholar.google.com/citations?user=NxF2f0cAAAAJ&hl=en)*.")
st.sidebar.divider()
st.sidebar.write("""**ADDZYME** is a software for predicting effect of an additive on an enzyme activity. It is developed at 
*Baqiyatallah University of Medical Sciences*.""")
st.sidebar.divider()
st.sidebar.write("""**Reference**:\n
Paper is *under production.*""")

#Body

with st.expander("See Information"):
    st.write("""**ADDZYME** is a machine learning-based algorithm for prediction of the effect of an additive on the activity of 
                an enzyme. It uses 30 ERT-Baseline models to report relative activity 
                For more information refer to our published paper.""")


substrate_smiles = st.text_input("Enter SMILES of the substrate:", 
                                placeholder="C1=CC=C(C=C1)C=O",
                                help="Use [PubChem](https://pubchem.ncbi.nlm.nih.gov/) for finding SMILES ")
if substrate_smiles:
    
    try:
        standard_substrate_smiles = standardize_smiles(substrate_smiles)
        substrate_molwt = Descriptors.MolWt(MolFromSmiles(standard_substrate_smiles))

        substrate_min_weight, substrate_max_weight = 120.063, 399.092

        if (substrate_molwt >= substrate_min_weight) & (substrate_molwt <= substrate_max_weight):
            st.success("Substrate is indomain.")
        else:
            st.warning(f"Substrate is outdomain. Molecular weight of substrate should be between {substrate_min_weight} and {substrate_max_weight}")
        
    except:
        assert substrate_smiles == None, "Enter valid SMILES. See the help."



additive_smiles = st.text_input("Enter SMILES of the additive:", 
                                placeholder="C1=CC=C(C=C1)C=O",
                                help="Use [PubChem](https://pubchem.ncbi.nlm.nih.gov/) for finding SMILES ")

if additive_smiles:
    
    try:
        standard_additive_smiles = standardize_smiles(additive_smiles)
        additive_molwt = Descriptors.MolWt(MolFromSmiles(standard_additive_smiles))

        additive_min_weight, additive_max_weight = 28.010, 276.116

        if (additive_molwt >= additive_min_weight) & (additive_molwt <= additive_max_weight):
            st.success("Additive is indomain.")
        else:
            st.warning(f"Additive is outdomain. Molecular weight of additive should be between {additive_min_weight} and {additive_max_weight}")
        
    except:
        assert additive_smiles == None, "Enter valid SMILES. See the help."
    

ec_number = ['3.1.1.1',
             '3.1.1.13',
             '3.1.1.2',
             '3.1.1.20',
             '3.1.1.25',
             '3.1.1.3',
             '3.1.1.43',
             '3.1.1.5',
             '3.1.1.59',
             '3.1.1.6',
             '3.1.1.60',
             '3.1.1.74',
             '3.1.1.79',
             '3.1.1.81',
             '3.1.3.12',
             '3.1.3.2',
             '3.1.4.46',
             '3.1.8.1',]

ec_number_index = dict(zip(ec_number, range(0, 18)))
one_hot = [0] * 18

selected_ec_number = st.selectbox("Select EC Number:", ec_number, help="Select one of the following EC number.")

if selected_ec_number:
    
    one_hot[ec_number_index[selected_ec_number]] = 1
    
temp = st.slider("Enter temperature:", min_value=20, max_value=75, step=5, help="Temperature applicability of domain is between 20 and 75.")
pH = st.slider("Enter pH:", min_value=5.0, max_value=9.0, step=0.5, help="pH applicability of domain is between 5 and 9.")

TEMP_MEAN, TEMP_STD = 37.769, 14.533
PH_MEAN, PH_STD =  7.519, 0.903

if temp:
    standard_temp = temp * TEMP_STD + TEMP_MEAN
    
if pH:
    standard_pH = pH * PH_STD + PH_MEAN

start = st.button("Predict Activity Status")

def load_and_predict(models, fv):

    predicted_values = []

    for model in models:
    
        reg = joblib.load(model)
        predicted_value = reg.predict(fv)
        predicted_values.append(predicted_value)
        
    return predicted_values

if start:
    
    with st.spinner("Please Wait..."):

        with open("./rdkit_descriptors.txt", "r") as file:
          f = file.read()

        rdkit_descriptors = f.replace("\n", " ").split(", ")
        DescCalc = MolecularDescriptorCalculator(rdkit_descriptors)

        substrate_mol = MolFromSmiles(standard_substrate_smiles)
        additive_mol = MolFromSmiles(standard_additive_smiles)

        substrate_desc = list(DescCalc.CalcDescriptors(substrate_mol))
        additive_desc = list(DescCalc.CalcDescriptors(additive_mol))

        baseline_fv = np.concatenate([substrate_desc, additive_desc, [standard_pH, standard_temp, *one_hot]])

        baseline_models = glob.glob(r".\models\ERT_Baseline_*.joblib")

        baseline_predicted_values = load_and_predict(baseline_models, baseline_fv.reshape(1, -1))

        l1 = np.array(baseline_predicted_values).reshape(30, 2)
        
        RA_MEAN = -17.321
        RA_STD =  71.558
        
        final_prediction = l1[:, 0] * RA_STD + RA_MEAN
        
        positive_cases = sum(final_prediction > 5)
        negative_cases = sum(final_prediction < -5)
        neutral_cases = sum((final_prediction >= -5) & (final_prediction <= 5))
        
        st.success(f"{positive_cases} cases of 30 model predict increase in activity.")
        st.success(f"{negative_cases} cases of 30 model predict decrase in activity.")
        st.success(f"{neutral_cases} cases of 30 model predict no change in activity.")
        
        
