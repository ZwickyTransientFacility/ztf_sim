# Installing ztf_sim

## Prerequisites

While you can install `ztf_sim` on your laptop, it can take a day or more to run a one-year simulation.  Accordingly you may prefer to install it on a remote server.

Since `ztf_sim` has several dependencies, we recommend installing in a [conda environment](http://conda.pydata.org/docs/using/envs.html).  

`ztf_sim` requires python 3.6 or later.

Packages needed:

* numpy
* scipy
* pandas
* [sqlalchemy](http://www.sqlalchemy.org/)
* [astropy](http://www.astropy.org/)
* [astroplan](http://www.astropy.org/)
* [scikit-learn](http://scikit-learn.org/)
* [sklearn_pandas](https://github.com/paulgb/sklearn-pandas)
* [xgboost](https://xgboost.readthedocs.io/)
* [transitions](https://github.com/tyarkoni/transitions)
* optional, for profiling: [pyinstrument](https://github.com/joerick/pyinstrument)

In the future we plan to provide appropriate recipes for installing these libraries in one go.

## Installing

`ztf_sim` is not yet packaged, so you'll need to download the source code to a convenient location:

    git clone https://github.com/ZwickyTransientFacility/ztf_sim.git
