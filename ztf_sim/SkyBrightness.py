"""Sky brightness model."""

import sklearn
from sklearn import model_selection, ensemble, preprocessing, pipeline
from sklearn import neighbors, svm, linear_model
from sklearn_pandas import DataFrameMapper
import joblib
import xgboost as xgb
import pandas as pd
import numpy as np
from .constants import FILTER_NAME_TO_ID, BASE_DIR


class SkyBrightness(object):
    """XGBoost-based sky brightness predictor trained on PTF/iPTF data.

    Loads one pre-trained pipeline (DataFrameMapper + XGBRegressor) per filter
    from ``data/sky_model/``. Predictions are returned as sky surface brightness
    in AB mag arcsec⁻².

    Attributes
    ----------
    clf_g : sklearn.pipeline.Pipeline
        Trained g-band sky brightness model.
    clf_r : sklearn.pipeline.Pipeline
        Trained r-band sky brightness model.
    clf_i : sklearn.pipeline.Pipeline
        Trained i-band sky brightness model.
    """

    def __init__(self):
        """Load pre-trained sky brightness models from disk.

        Models are read from ``data/sky_model/sky_model_{g,r,i}.pkl``.
        """
        self.clf_r = joblib.load(BASE_DIR + '../data/sky_model/sky_model_r.pkl')
        self.clf_g = joblib.load(BASE_DIR + '../data/sky_model/sky_model_g.pkl')
        self.clf_i = joblib.load(BASE_DIR + '../data/sky_model/sky_model_i.pkl')

    def predict(self, df):
        """Predict sky surface brightness for a set of pointings.

        Dispatches to the appropriate per-filter model based on
        ``filter_id`` in *df*.

        Parameters
        ----------
        df : pandas.DataFrame
            Pointings to evaluate. Required columns:

            * ``filter_id`` — int, 1 = g, 2 = r, 3 = i
            * ``moonillf`` — float, Moon illuminated fraction (0–1)
            * ``moonalt`` — float, Moon altitude in degrees
            * ``moon_dist`` — float, angular distance to Moon in degrees
            * ``azimuth`` — float, pointing azimuth in degrees
            * ``altitude`` — float, pointing altitude in degrees
            * ``sunalt`` — float, Sun altitude in degrees

        Returns
        -------
        pandas.Series of float
            Predicted sky brightness in AB mag arcsec⁻², indexed like *df*.
            Values are ``NaN`` for rows where the filter ID is not recognised.
        """

        filter_ids = df['filter_id'].unique()
        assert(np.sum(filter_ids > 3) == 0)

        sky = pd.Series(np.nan, index=df.index, name='sky_brightness')
        wg = (df['filter_id'] == FILTER_NAME_TO_ID['g'])
        if np.sum(wg):
            sky[wg] = self.clf_g.predict(df[wg])
        wr = (df['filter_id'] == FILTER_NAME_TO_ID['r'])
        if np.sum(wr):
            sky[wr] = self.clf_r.predict(df[wr])
        wi = (df['filter_id'] == FILTER_NAME_TO_ID['i'])
        if np.sum(wi):
            sky[wi] = self.clf_i.predict(df[wi])

        return sky


class FakeSkyBrightness(object):
    """Constant sky brightness stub for testing.

    Returns 20.0 AB mag arcsec⁻² for all inputs regardless of observing
    conditions.
    """

    def __init__(self):
        """Initialise the fake sky brightness model (no-op)."""
        pass

    def predict(self, df):
        """Return a constant sky brightness for every pointing.

        Parameters
        ----------
        df : pandas.DataFrame
            Pointings (any columns; only the index is used).

        Returns
        -------
        pandas.Series of float
            Series of 20.0 (AB mag arcsec⁻²) with the same index as *df*.
        """
        y = np.ones(len(df)) * 20.
        return pd.Series(y, index=df.index, name='sky_brightness')


def train_sky_model(filter_name='r', df=None):
    """Train and save an XGBoost sky brightness model for one filter.

    Fits a ``sklearn`` pipeline (``DataFrameMapper`` standardiser followed by
    ``XGBRegressor``) on PTF/iPTF DIQ data and serialises the result with
    ``joblib`` to ``data/sky_model/sky_model_{filter_name}.pkl``.

    Parameters
    ----------
    filter_name : str, optional
        Filter to train. One of ``'g'``, ``'r'``, ``'i'``. Default is
        ``'r'``.
    df : pandas.DataFrame or None, optional
        Training data. If ``None``, the function reads
        ``data/ptf-iptf_diq.csv.gz``. Required columns: ``filterkey``,
        ``sky_brightness``, ``moonillf``, ``moonalt``, ``moon_dist``,
        ``azimuth``, ``altitude``, ``sunalt``.

    Returns
    -------
    sklearn.pipeline.Pipeline
        The fitted pipeline (also saved to disk).

    Notes
    -----
    A 20 % hold-out test set is used to print the R² score. The iPTF i-band
    uses filter key 4 (mapped internally).
    """

    # PTF used 4 for i-band
    filterid_map = {'r': 2, 'g': 1, 'i': 4}

    if df is None:
        df = pd.read_csv(BASE_DIR + '../data/ptf-iptf_diq.csv.gz')
    # note that this is by pid, so there are multiple entries per image...

    df = df[df['filterkey'] == filterid_map[filter_name]].copy()

    # IPAC stores negative moonillf, but astroplan.moon_illumination does not
    df.loc[:, 'moonillf'] = np.abs(df['moonillf'])

    # returns dataframes!
    X_train, X_test, y_train, y_test = model_selection.train_test_split(
        df, df['sky_brightness'], test_size=0.2)

    # don't really need to standardize for RF, but preprocessing is nice
    # preprocessing through sklearn_pandas raises a deprecation warning
    # from sklearn, so skip it.
    mapper = DataFrameMapper([
        (['moonillf'], preprocessing.StandardScaler()),
        (['moonalt'],   preprocessing.StandardScaler()),
        (['moon_dist'], preprocessing.StandardScaler()),
        (['azimuth'],  preprocessing.StandardScaler()),
        (['altitude'], preprocessing.StandardScaler()),
        (['sunalt'],   preprocessing.StandardScaler())])
    #('filterkey',  None)])

    clf = pipeline.Pipeline([
        ('featurize', mapper),
        ('xgb', xgb.XGBRegressor())])
    #('svr', svm.SVR(kernel='poly',degree=2))])
    #('knr', neighbors.KNeighborsRegressor(n_neighbors=15, weights='distance', algorithm='auto'))])
    #('lm', linear_model.BayesianRidge())])
    #('rf', ensemble.RandomForestRegressor(n_jobs=-1))])

    clf.fit(X_train, y_train.values.reshape(-1, 1))
    print(clf.score(X_test, y_test.values.reshape(-1, 1)))

    joblib.dump(clf, BASE_DIR + '../data/sky_model/sky_model_{}.pkl'.format(filter_name))

    return clf
