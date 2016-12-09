
import json
import numpy as np
from ObservingProgram import ObservingProgram
from fields import Fields
from constants import *


class ZTFConfiguration(object):

    def __init__(self, config_file):

        self.load_configuration(config_file)
        self.check_configuration()

    def load_configuration(self, config_file):
        with open(config_file, 'r') as f:
            config = json.load(f)
        # TODO: construct and validate a schema
        self.config = config

    def check_configuration(self):

        if (np.sum([prog['observing_fraction'] for prog in
                    self.config['observing_programs']]) != 1.0):
            raise ValueError('Observing fractions must sum to 1')

        # could do this via schema validation
        for prog in self.config['observing_programs']:
            if prog['program_name'] not in PROGRAM_NAMES:
                raise ValueError('{} not in known programs'.format(
                    prog['program_name']))

    def build_observing_programs(self):

        OPs = []
        f = Fields()
        for prog in self.config['observing_programs']:
            field_ids = f.select_field_ids(**prog['field_selections'])
            OP = ObservingProgram(PROGRAM_NAME_TO_ID[prog['program_name']],
                                  prog['subprogram_name'], prog[
                                      'observing_fraction'],
                                  field_ids, prog['filter_ids'],
                                  prog['internight_gap_days'] * u.day,
                                  prog['n_visits_per_night'],
                                  prog['intranight_gap_min'] * u.min,
                                  prog['intranight_half_width_min'] * u.min,
                                  nightly_priority=prog['nightly_priority'],
                                  filter_choice=prog['filter_choice'])
            OPs.append(OP)

        return OPs
