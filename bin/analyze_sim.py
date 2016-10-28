#!/usr/bin/env python

"""Calculate basic efficiency statistics for a ztf_sim run."""

import sys
# hack to get the path right
sys.path.append('..')
import numpy as np
from collections import OrderedDict
from sqlalchemy import create_engine
import pandas as pd
from astropy.time import Time
import astropy.coordinates as coord
import astropy.units as u
from ztf_sim.constants import *
from ztf_sim.utils import *


def calc_stats(sim_name):

    df = df_read_from_sqlite(sim_name, tablename='Summary', directory='sims')

    stats = OrderedDict()

    stats['Simulation Name'] = sim_name

    stats['Number of Nights'] = len(df.night.unique())

    mjds = df.expMJD.apply(np.floor).unique()
    t_mjds = Time(mjds, format='mjd')
    hours_of_darkness = approx_hours_of_darkness(t_mjds).to(u.hour).value

    stats['Average Hours of Darkness'] = np.mean(hours_of_darkness)

    stats['Average Nightly Number of Exposures'] = \
        len(df) / stats['Number of Nights']

    # time used for science: exposing or slewing/filter changing
    stats['Fraction of time usable'] = \
        (df.visitExpTime.sum() + df.slewTime.sum()) / 3600. / \
        np.sum(hours_of_darkness)

    stats['Average Open Shutter Time (h)'] = \
        df.visitExpTime.sum() / 3600. / stats['Number of Nights']

    # fraction of usable science time
    stats['Open Shutter Fraction'] = \
        df.visitExpTime.sum() / (df.visitExpTime.sum() + df.slewTime.sum())

    stats['Mean Slew Time (s)'] = df.slewTime.mean()
    stats['Mean Slew Distance (deg)'] = np.degrees(df.slewDist.mean())

    stats['Median Airmass'] = df.airmass.median()

    # program breakdown
    pgrp = df.groupby('propID')
    stats['Program Fraction'] = (pgrp['fieldID'].agg(len) / len(df)).to_dict()

    # filter breakdown
    fgrp = df.groupby('filter')
    stats['Filter Fraction'] = (fgrp['fieldID'].agg(len) / len(df)).to_dict()

    # fraction of completed sequences

    # summed figure of merit

    for k, v in stats.iteritems():
        print('{}\t{}'.format(k, v))
    return stats

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: ./analyze_sim.py sim_name')
    else:
        calc_stats(sys.argv[1])
