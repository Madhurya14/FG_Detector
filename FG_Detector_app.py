import streamlit as st
import pandas as pd
from rdkit import Chem

# Define the SMARTS patterns
SMARTS_PATTERNS = [
    ("Amide", Chem.MolFromSmarts("C(=O)N")),
    ("Carboxylic Acid", Chem.MolFromSmarts("C(=O)[OH]")),
    ("Ketone", Chem.MolFromSmarts("C(=O)[C;!$(C(=O)N);!$(C(=O)O)]")),
    ("Primary Amine", Chem.MolFromSmarts("[NX3;H2][#6]")),
    ("Secondary Amine", Chem.MolFromSmarts("[NX3;H1][#6]")),
    ("Tertiary Amine", Chem.MolFromSmarts("[NX3;H0][#6]")),
    ("Alcohol", Chem.MolFromSmarts("[OX2H][#6]")),
    ("Phenol", Chem.MolFromSmarts("c[OH]")),
    ("Ether", Chem.MolFromSmarts("[OD2]([#6])[#6]")),
    ("Alkene", Chem.MolFromSmarts("C=C")),
    ("Alkyne", Chem.MolFromSmarts("C#C")),
    ("Nitrile", Chem.MolFromSmarts("C#N")),
    ("Nitro", Chem.MolFromSmarts("[$([NX3](=O)=O)]")),
    ("Halide", Chem.MolFromSmarts("[F,Cl,Br,I]"))
]

def detect_groups(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ["Invalid SMILES"] * len(SMARTS_PATTERNS)
    return [len(mol.GetSubstructMatches(patt)) for _, patt in SMARTS_PATTERNS]

st.title("🔬 Functional Group Detector from SMILES")

uploaded_file = st.file_uploader("Upload a CSV file with SMILES column", type=["csv"])

if uploaded_file:
    df = pd.read_csv(uploaded_file)
    if "SMILES" not in df.columns:
        st.error("The uploaded file must have a column named 'SMILES'")
    else:
        st.success("File uploaded successfully. Processing...")

        fg_data = []
        for smi in df["SMILES"]:
            fg_data.append(detect_groups(smi))

        fg_df = pd.DataFrame(fg_data, columns=[name for name, _ in SMARTS_PATTERNS])
        output_df = pd.concat([df, fg_df], axis=1)

        st.write("### Detected Functional Groups")
        st.dataframe(output_df)

        st.download_button("📥 Download Results as CSV", output_df.to_csv(index=False), file_name="fg_results.csv", mime="text/csv")