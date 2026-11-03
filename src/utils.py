import numpy as np
import scipy.sparse as sp
from sklearn.preprocessing import StandardScaler

def get_log_adata(ad):
    """Use .raw (all genes, log-normalized) if present."""
    return ad.raw.to_adata() if ad.raw is not None else ad

def to_matrix(ad, genes):
    """Cells x genes matrix on a fixed gene list, z-scored per dataset.
    Genes missing from `ad` are filled with 0."""
    ad = get_log_adata(ad)
    idx = {g: i for i, g in enumerate(ad.var_names)}
    X = ad.X.toarray() if sp.issparse(ad.X) else np.asarray(ad.X)
    out = np.zeros((ad.n_obs, len(genes)), dtype=np.float32)
    for j, g in enumerate(genes):
        if g in idx:
            out[:, j] = X[:, idx[g]]
    return StandardScaler().fit_transform(out)
