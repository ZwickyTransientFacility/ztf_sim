import sklearn
from sklearn import cross_validation, ensemble, preprocessing, pipeline
from sklearn_pandas import DataFrameMapper, cross_val_score
from sklearn.externals import joblib
import pandas as pd
import numpy as np


def train_sky_model():

    df = pd.read_csv('../data/ptf-iptf_diq.csv.gz')
    # note that this is by pid, so there are multiple entries per image...

    # IPAC stores negative moonillf, but astroplan.moon_illumination does not
    df['moonillf'] = np.abs(df['moonillf'])

    # returns dataframes!
    X_train, X_test, y_train, y_test = cross_validation.train_test_split(
        df, df['sky_brightness'], test_size=0.2)

    # don't really need to standardize for RF, but preprocessing is nice
    mapper = DataFrameMapper([
        ('moonillf', preprocessing.StandardScaler()),
        ('moonalt',  preprocessing.StandardScaler()),
        ('moon_dist',  preprocessing.StandardScaler()),
        ('azimuth',  preprocessing.StandardScaler()),
        ('altitude',  preprocessing.StandardScaler()),
        ('sunalt',  preprocessing.StandardScaler()),
        ('filterkey',  None)])

    clf = pipeline.Pipeline([
        ('featurize', mapper),
        ('rf', ensemble.RandomForestRegressor(n_jobs=-1))])

    clf.fit(X_train, y_train)
    print clf.score(X_test, y_test)

    joblib.dump('../data/sky_model.pkl')

    return clf
