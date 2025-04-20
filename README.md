# ztf_sim
Fork of the :telescope: Scheduling library for the Zwicky Transient Facility adapted for the La Silla Schmidt Southern Survey (LS4).

Implements the Integer Linear Programming scheduling algorithm described in 
[Bellm et al. 2019](https://dx.doi.org/10.1088/1538-3873/ab0c2a) (PASP 131, 1000).

## Installation

Follow these [instructions](https://github.com/steveschulze/ls4_scheduler_general/README.md).


<!-- You will need a license for the [Gurobi](http://www.gurobi.com/) optimizer on the machine you want to run simulations on.  [Free academic licenses](http://www.gurobi.com/academia/for-universities) are readily available.

You are strongly encouraged to use `conda` and a conda environment for the installation.


```
conda create -n ztf_sim_test

conda activate ztf_sim_test

conda install python=3.11


conda config --add channels conda-forge 
conda config --add channels http://conda.anaconda.org/gurobi
conda install pip numpy scipy astropy astroplan pandas scikit-learn xgboost sqlalchemy gurobi


pip install sklearn_pandas 
pip install transitions

pip install -e git+https://github.com/ZwickyTransientFacility/ztf_sim.git#egg=ztf_sim
```

(To remove the environment, use `conda remove --name ztf_sim_test --all`.) -->


## Configuration

### Scheduler Configuration

The scheduler configuration determines which observing programs will run (fields, filters, cadences, etc.)  

An example set of scheduler configuration file is provided in `sims/`.

You can copy them locally with 
```
wget https://raw.githubusercontent.com/ZwickyTransientFacility/ztf_sim/master/sims/example_scheduler_config.json
wget https://raw.githubusercontent.com/ZwickyTransientFacility/ztf_sim/master/sims/survey_180501.json
wget https://raw.githubusercontent.com/ZwickyTransientFacility/ztf_sim/master/sims/reference_building.json
```

### Simulation Configuration

The simulation configuration determines which nights to simulate, which historical weather to use, and whether to overwrite the existing simulated database.

An example configuration file is provided in `config/default.cfg`.

You can copy it locally with `wget https://raw.githubusercontent.com/ZwickyTransientFacility/ztf_sim/master/config/default.cfg` 

## Running

`run_ztf_sim --help` summarizes the argument of the command line driver for the simulations.  Assuming you've copied the configuration files to your current directory, you should now be able to run

```
run_ztf_sim example_scheduler_config.json default.cfg
```

which will write an sqlite database file named `example_ztf_schedule.db` to the current directory.

## Output

The simulated schedule is written to a SQLite database in the LSST [Operations Simulator format](https://www.lsst.org/scientists/simulations/opsim/summary-table-column-descriptions-v335) with a few additional columns.  Example code for reading and summarizing the simulated schedule is located in the `bin/` directory, and can be run as

```
analyze_ztf_sim example_ztf_schedule.db
```

The simulator also writes a logfile, which may be monitored while the simulation is running with `tail -f example_ztf_schedule_log.txt`. The filename of database and logfile are generated from the keyword `run_name` in `example_scheduler_config.json`.

## How does this version differ from the original ZTF Scheduler?

- Docstrings were added with CoPilot. `README.md` was modified.
- The directories `notebooks` and `tests` were removed.
- `constants.py`: All telescope specific properties were re-configured to the La Silla Schmidt telescope. This also includes the systemwide variables P48_loc, P48_Observer, P48_slew_pars, and P48_slew_pars_goal.
- `Fields.py`: The declination cut (Line 83) was changed from dec >-32 deg (Palomar) to <+20 deg (La Silla).
<!-- - `Fields.py`: Lines 87-90 define the primary and secondary grids. -->
- `Fields.py`: Line 53 points to the file containing the pointing information. This was changed to `data/LS4_Fields_cleaned.txt`.
- The declination cut (Line 83) was changed from dec >-32 deg (Palomar) to <+20 deg (La Silla).
- `QueueManager.py`: ListQueueManager (L1363-1385) defines coordinate ranges to avoid. Not defined for the other optimisers. Toggled it off.
- `utils.py`: The TCS constraints of the P48 telescope (L419-L446) are toggled off.
- Two warnings are suppressed for convenience. They can be activated by toggling off Lines 10-13 in `/PATH/TO/LS4/SCHEDULER/bin/run_ztf.py`.<br>

  1. **Future warning**
   ```
   /PATH/TO/LS4/SCHEDULER/utils.py:60: FutureWarning: Setting the location attribute post initialization will be disallowed in a future version of Astropy. Instead you should set the location when creating the Time object. In the future, this will raise an AttributeError.
   time.location = P48_loc
   ```

  2. **Geolocation warning**
   ```
   WARNING: NonRotationTransformationWarning: transforming other coordinates from <ICRS Frame> to <GCRS Frame (obstime=2018-05-04 22:45:00, obsgeoloc=(-3808151.42273061, 4071448.43160014, -3093186.1485597) m, obsgeovel=(-296.88686399, -277.29877698, 0.51159639) m / s)>. Angular separation can depend on the direction of the transformation. [astropy.coordinates.baseframe]
   /PATH/TO/LS4/SCHEDULER/utils.py:60: FutureWarning: Setting the location attribute post initialization will be disallowed in a future version of Astropy. Instead you should set the location when creating the Time object. In the future, this will raise an AttributeError.
   ```


<!-- - are defined here. To avoid larger code changes, we overwrite the specs of P48 and of ZTF camera with the parameters of the LS4 camera. 2) Also contains the mappings between programme names and IDs and filters and fids. -->


<!-- IMPORTANT: observatory specific  -->

<!-- - Fields.py: L85-L89 define primary, secondary, etc grid -->

<!-- - Coordinates in the output database are reported in u.radian not u.deg. -->


<!-- - field_selection_function.py: in the case of cutting on ra, dec, Galactic latitude, and field ID is not adequate for your science case, you can define much more complex field selection functions. -->

<!-- - The simulator can throw many depreciation warnings. Activate L18-L21 in simulate.py to toggle most of them off. -->


<!-- - Don't know why but running `run_ztf_sim' can easily crash even when we do a simulation with the same setup. Unclear why this happens. It does not dependent on the computer architecture. -->