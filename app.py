
import streamlit as st
from rdkit import Chem
import pandas as pd
from io import BytesIO
import base64

# Priority-based SMARTS definitions
SMARTS_PATTERNS = [
    ("Carboxylic Acid", Chem.MolFromSmarts("C(=O)[OH]")),
    ("Amide", Chem.MolFromSmarts("C(=O)N")),
    ("Ketone", Chem.MolFromSmarts("C(=O)[C;!$(C(=O)N);!$(C(=O)O)]")),
    ("Primary Amine", Chem.MolFromSmarts("[NX3;H2][#6]")),
    ("Secondary Amine", Chem.MolFromSmarts("[NX3;H1][#6]")),
    ("Tertiary Amine", Chem.MolFromSmarts("[NX3;H0][#6]")),
    ("Hydroxyl", Chem.MolFromSmarts("[OX2H]")),
    ("Ether", Chem.MolFromSmarts("C-O-C")),
    ("Aromatic Ring", Chem.MolFromSmarts("a")),
    ("Fluorine", Chem.MolFromSmarts("[F]")),
    ("Chlorine", Chem.MolFromSmarts("[Cl]")),
    ("Alkane", Chem.MolFromSmarts("[CH3,CH2]")),
    ("Alkene", Chem.MolFromSmarts("C=C")),
    ("Alkyne", Chem.MolFromSmarts("C#C")),
]

def detect_groups(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {name: "Invalid" for name, _ in SMARTS_PATTERNS}
    return {name: len(mol.GetSubstructMatches(pattern)) for name, pattern in SMARTS_PATTERNS}

def to_excel(df):
    output = BytesIO()
    writer = pd.ExcelWriter(output, engine='xlsxwriter')
    df.to_excel(writer, index=False, sheet_name='Results')
    writer.close()
    return output.getvalue()

st.title("Functional Group Detector")

smiles_input = st.text_area("Paste SMILES (one per line):")

if st.button("Detect Functional Groups"):
    smiles_list = [s.strip() for s in smiles_input.split("\n") if s.strip()]
    results = []
    for smi in smiles_list:
        row = {"SMILES": smi}
        row.update(detect_groups(smi))
        results.append(row)
    df = pd.DataFrame(results)
    st.dataframe(df)

    excel_data = to_excel(df)
    b64 = base64.b64encode(excel_data).decode()
    href = f'<a href="data:application/octet-stream;base64,{b64}" download="fg_results.xlsx">📥 Download Excel</a>'
    st.markdown(href, unsafe_allow_html=True)
