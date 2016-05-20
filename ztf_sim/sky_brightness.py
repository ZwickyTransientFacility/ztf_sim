import sklearn
from sklearn import cross_validation, ensemble, preprocessing, pipeline
from sklearn_pandas import DataFrameMapper, cross_val_score
from sklearn.externals import joblib
import pandas as pd
import numpy as np


class SkyBrightness(object):

    def __init__(self):
        self.clf_r = joblib.load('../data/sky_model/sky_model_r.pkl')
        self.clf_g = joblib.load('../data/sky_model/sky_model_g.pkl')

    def predict(self, df):
        """df is a dataframe with columns:
        mooonillf: 0-1
        moonalt: degrees
        moon_dist: degrees
        azimuth: degrees
        altitude: degrees
        sunalt: degrees
        filterkey: 1, 2"""
        # TODO: make this return a dataframe?

        # for now, only can run on one filterkey
        assert(len(set(df['filterkey'])) == 1)
        if df['filterkey'][0] == 1:
            return self.clf_g.predict(df)
        if df['filterkey'][0] == 2:
            return self.clf_r.predict(df)


def train_sky_model(filter_name='r'):

    filterid_map = {'r': 2, 'g': 1}

    df = pd.read_csv('../data/ptf-iptf_diq.csv.gz')
    # note that this is by pid, so there are multiple entries per image...

    df = df[df['filterkey'] == filterid_map[filter_name]]

    # IPAC stores negative moonillf, but astroplan.moon_illumination does not
    df['moonillf'] = np.abs(df['moonillf'])

    # returns dataframes!
    X_train, X_test, y_train, y_test = cross_validation.train_test_split(
        df, df['sky_brightness'], test_size=0.2)

    # don't really need to standardize for RF, but preprocessing is nice
    # preprocessing through sklearn_pandas raises a deprecation warning
    # from sklearn, so skip it.
    mapper = DataFrameMapper([
        ('moonillf', None),  # preprocessing.StandardScaler()),
        ('moonalt',  None),
        ('moon_dist',  None),
        ('azimuth',  None),
        ('altitude',  None),
        ('sunalt',  None)])
    #('filterkey',  None)])

    clf = pipeline.Pipeline([
        ('featurize', mapper),
        ('rf', ensemble.RandomForestRegressor(n_jobs=-1))])

    clf.fit(X_train, y_train)
    print clf.score(X_test, y_test)

    joblib.dump(clf, '../data/sky_model/sky_model_{}.pkl'.format(filter_name))

    return clf
