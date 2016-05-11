import sklearn
from sklearn import cross_validation, ensemble
import pandas as pd


def train_sky_model():

    df = pd.read_csv('../data/ptf-iptf_diq.csv.gz')
    # note that this is by pid, so there are multiple entries per image...

    # delete uninformative columns:
    for col in [u'Unnamed: 0', u'pid', u'mjd', u'fwhmsex',
                u'airmass', u'limiting_mag']:
        del df[col]

    Y = df['sky_brightness'].as_matrix()
    del df['sky_brightness']
    X = df.as_matrix()

    X_train, X_test, y_train, y_test = cross_validation.train_test_split(
        X, Y, test_size=0.2)

    clf = ensemble.RandomForestRegressor(n_estimators=500)
    clf.fit(X_train, y_train)
    print clf.score(X_test, y_test)

    return clf
