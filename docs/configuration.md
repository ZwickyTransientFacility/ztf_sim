
# Configuring ztf_sim

*NOTE THAT THE DOCUMENTATION BELOW IS OUTDATED AND SOME OPTIONS NO LONGER WORK AS DESCRIBED*

Here is an example configuration file, which is in JSON format:

```
{
	"run_name": "test",
	"start_time": "2016-03-20 02:30:00",
	"weather_year": "None",
	"survey_duration_days": 1.0,
	"block_programs": false,
	"observing_programs": [
		{"program_name": "MSIP",
		 "subprogram_name": "all_sky",
		 "program_observing_fraction": 0.4,
		 "subprogram_fraction": 1.0,
		 "field_selections":{"dec_range":[-30,90],
				     "grid_id":0},
		 "filter_choice": "rotate",
		 "filter_ids": [1, 2],
		 "internight_gap_days": 3,
		 "n_visits_per_night": 2,
		 "intranight_gap_min": 60,
		 "intranight_half_width_min": 20,
		 "nightly_priority": "oldest"},
		{"program_name": "collaboration",
		 "subprogram_name": "all_sky",
		 "program_observing_fraction": 0.4,
		 "subprogram_fraction": 1.0,
		 "field_selections":{"dec_range":[-30,90],
				     "grid_id":0},
		 "filter_choice": "sequence",
		 "filter_ids": [1, 2, 2, 1],
		 "internight_gap_days": 1,
		 "n_visits_per_night": 4,
		 "intranight_gap_min": 60,
		 "intranight_half_width_min": 20,
		 "nightly_priority": "mean_observable_airmass"},
		{"program_name": "Caltech",
		 "subprogram_name": "all_sky",
		 "program_observing_fraction": 0.2,
		 "subprogram_fraction": 1.0,
		 "field_selections":{"dec_range":[-30,90],
				     "grid_id":0},
		 "filter_choice": "rotate",
		 "filter_ids": [1, 2],
		 "internight_gap_days": 3,
		 "n_visits_per_night": 3,
		 "intranight_gap_min": 90,
		 "intranight_half_width_min": 20,
		 "nightly_priority": "random"}
	]
}
```




## Run configuration

Basic parameters are set at the top level.

`run_name`: A string to identify this run.  The output database will have this name.

`start_time`: An astropy.Time-compliant string defining when the simulation should begin.

`weather_year`: Which year of historical PTF data to use to determine weather losses during the simulation.  May be any of `2009-2015`.  Set to `"None"` for no weather losses.

`survey_duration_days`: A float specifying how many days to simulate.

`block_programs`: A Boolean specifying whether one and only one program should take observations during one observing block.

## Observing Programs

`ztf_sim` can (in theory) support arbitary sub-surveys, each of which is specified as a component under `observing_programs`.  Each sub-survey is defined by the following elements:

`program_name`: There are three top-level observing programs: `"MSIP"`, `"collaboration"`, and `"Caltech"`, which split the available observing time 40/40/20.  

`subprogram_name`: A string identifying the sub-survey, such as `"Collaboration Galactic Plane"`.

`program_observing_fraction`: A float < 1.0 that specifies how much of the total time to be allocated to this program.  All MSIP and collaboration sub-programs should have this set to 0.4, and Caltech sub-programs should have it set to 0.2.

`subprogram_fraction`: A float < 1.0 that specifies how much of the _program_ time to be allocated to this sub-survey.  The sum of all `subprogram_fractions` for a given program (e.g., MSIP) should equal one.

`field_selections`: A nested set of criteria for assigning fields in the pre-defined field grid to the survey footprint.  All selections are optional. These arguments are passed to `Fields.select_fields`.  Options include:

* `ra_range`: Range of Right Ascension values (degrees)
* `dec_range`:  Range of Declination values (degrees)
* `l_range`: Range of Galactic longitude values (degrees)
* `b_range`: Range of Galactic latitude values (degrees)
* `abs_b_range`: Range of absolute Galactic latitude values (e.g., [20,90] to select extragalactic sky) (degrees)
* `ecliptic_lon_range`: Range of ecliptic longitude values (degrees)
* `ecliptic_lat_range`: Range of ecliptic latitude values (degrees)
* `grid_id`: 0 or 1.  If specified, only use fields from the primary or offset pointing grid.

It is also possible to select fields by the total number of observations (possibly subdivided by filter or programs) and or the last observation date (possibly subdivided by filter or programs).  These are less likely to be used for a new survey configuration.

`filter_choice`: Option for how to determine what filter(s) to use on a given night for this subprogram.  Currently implemented options include:

* `rotate`: rotate through the filter ids specified in `filter_ids`, with one and only one filter used for all observations per night
* `sequence`: every night, use the `filter_ids` in the order specified.  `filter_ids` should then have length `n_visits_per_night`.

`filter_ids`: list of filter ids (1=*g*, 2=*r*) to be used, according to the scheme defined by `filter_choice`.

`internight_gap_days`: Number of days before an observed field should be revisited.  Integer >= 1.

`n_visits_per_night`: Number of images to take per field per night.  Integer >= 1.

`intranight_gap_min`: Center of allowed time between exposures in minutes (e.g., 60 minutes)

`intranight_half_width_min`: Tolerance on window between exposures during the night.  If the previous exposure 
of the field occured at *T0*, the next exposure can occur at any time between *T0 + intranight_gap_min - intranight_half_width_min* and *T0 + intranight_gap_min + intranight_half_width_min*.

`nightly_priority`: Scheme for selecting a subset of the fields in the total sub-program footprint on any given night.  The number of fields selected per night is set by the `observing_fraction` and the number of hours in the night.  Currently implemented options include:

* `oldest`: Every night, select the set of fields that has been observed least recently by this sub-program.
*  `mean_observable_airmass`: Select the subset of fields with the lowest mean observable airmass (`> MAX_AIRMASS`) during the night.
*  `rotate`: rotate between stripes in right ascension.  The number of strips is set by `internight_gap_days`.
*  `random`: Select a random subset of fields.
