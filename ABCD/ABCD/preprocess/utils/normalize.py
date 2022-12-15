import numpy as np

def normalize(X, add_genomics = False):
    X_ts = X[:, 0:7, :]
    mean = np.mean(X_ts, axis = (0, 2), keepdims =True)
    var = np.var(X_ts, axis = (0, 2), keepdims = True)
    X_ts = (X_ts - mean) / np.sqrt(var)
    X_ = X[:, 7:, :]
    if add_genomics:
        X_33 = X_[:,[32, 34, 35, 36, 37],:]
    else:
        X_33 = X_[:,[32, 34, 35],:]
    mean_33 = np.mean(X_33, axis = 0, keepdims =True)
    var_33 = np.var(X_33, axis = 0, keepdims = True)
    X_33 = (X_33 - mean_33) / np.sqrt(var_33)

    if add_genomics:
        X_[:,[32, 34, 35, 36, 37],:] = X_33
    else:
        X_[:,[32, 34, 35],:] = X_33
    X_out = np.concatenate([X_ts, X_], axis = 1)
    return X_out