import sklearn
from sklearn import cross_validation, ensemble, preprocessing, pipeline, neighbors
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

        # for now, only can run on one filterkey
        assert(len(set(df['filterkey'])) == 1)
        if df['filterkey'][0] == 1:
            y = self.clf_g.predict(df)
        if df['filterkey'][0] == 2:
            y = self.clf_r.predict(df)
        return pd.Series(y, index=df.index, name='sky_brightness')


class FakeSkyBrightness(object):

    def __init__(self):
        pass

    def predict(self, df):
        y = np.ones(len(df)) * 20.
        return pd.Series(y, index=df.index, name='sky_brightness')


def train_sky_model(filter_name='r', df=None):

    filterid_map = {'r': 2, 'g': 1}

    if df is None:
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
        ('moonillf', preprocessing.StandardScaler()),
        ('moonalt',   preprocessing.StandardScaler()),
        ('moon_dist', preprocessing.StandardScaler()),
        ('azimuth',  preprocessing.StandardScaler()),
        ('altitude', preprocessing.StandardScaler()),
        ('sunalt',   preprocessing.StandardScaler())])
    #('filterkey',  None)])

    clf = pipeline.Pipeline([
        ('featurize', mapper),
        ('knr', neighbors.KNeighborsRegressor(n_neighbors=15, weights='distance', algorithm='auto'))])
    #('lm', linear_model.BayesianRidge())])
    #('rf', ensemble.RandomForestRegressor(n_jobs=-1))])

    clf.fit(X_train, y_train)
    print clf.score(X_test, y_test)

    joblib.dump(clf, '../data/sky_model/sky_model_{}.pkl'.format(filter_name))

    return clf
