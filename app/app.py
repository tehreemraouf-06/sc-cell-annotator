import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import streamlit as st, scanpy as sc, joblib, numpy as np, pandas as pd, tempfile
import matplotlib.pyplot as plt
from src.utils import to_matrix, get_log_adata

st.set_page_config(page_title="Cell Type Annotator", layout="wide")
st.title("🧬 Single-Cell Cell Type Annotator (PBMC)")
st.write("Upload a log-normalized `.h5ad` file with **gene symbols** as var names, "
         "or try the demo data. Low-confidence cells are labeled *Unassigned*.")

@st.cache_resource
def load_model():
    return joblib.load(os.path.join(os.path.dirname(__file__), '..', 'models', 'annotator.joblib'))
bundle = load_model(); model, genes = bundle['model'], bundle['genes']

threshold = st.slider("Confidence threshold", 0.0, 0.95, 0.6, 0.05)
file = st.file_uploader("Upload .h5ad", type=["h5ad"])
use_demo = st.button("Use demo data (pbmc68k_reduced)")

ad = None
if file is not None:
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        tmp.write(file.read()); path = tmp.name
    ad = sc.read_h5ad(path)
elif use_demo:
    ad = sc.datasets.pbmc68k_reduced()

if ad is not None:
    overlap = len(set(get_log_adata(ad).var_names) & set(genes))
    st.info(f"{ad.n_obs} cells | {overlap}/{len(genes)} model genes found in your data")
    if overlap < 0.3 * len(genes):
        st.warning("Low gene overlap: predictions may be unreliable. Check that var_names are gene symbols.")

    X = to_matrix(ad, genes)
    proba = model.predict_proba(X)
    conf = proba.max(axis=1)
    pred = model.classes_[proba.argmax(axis=1)].astype(object)
    pred[conf < threshold] = "Unassigned"
    ad.obs['predicted'] = pd.Categorical(pred); ad.obs['confidence'] = conf

    if 'X_umap' not in ad.obsm:
        sc.pp.pca(ad); sc.pp.neighbors(ad); sc.tl.umap(ad)

    c1, c2 = st.columns(2)
    with c1:
        fig, ax = plt.subplots(figsize=(5, 4))
        sc.pl.umap(ad, color='predicted', ax=ax, show=False)
        st.pyplot(fig)
    with c2:
        st.bar_chart(ad.obs['predicted'].value_counts())
    st.download_button("Download predictions (CSV)",
                       ad.obs[['predicted', 'confidence']].to_csv().encode(), "predictions.csv")
