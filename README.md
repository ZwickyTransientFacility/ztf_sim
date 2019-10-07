# ztf_sim
:telescope: Scheduling library for the Zwicky Transient Facility.

Implements the Integer Linear Programming scheduling algorithm described in 
[Bellm et al. 2019](https://dx.doi.org/10.1088/1538-3873/ab0c2a) (PASP 131, 1000).

## Installation

You will need a license for the [Gurobi](http://www.gurobi.com/) optimizer on the machine you want to run simulations on.  [Free academic licenses](http://www.gurobi.com/academia/for-universities) are readily available.

You are strongly encouraged to use `conda` and a conda environment for the installation.


```
conda create -n ztf_sim_test

conda activate ztf_sim_test

conda install python=3.7


conda config --add channels conda-forge 
conda config --add channels http://conda.anaconda.org/gurobi
conda install pip numpy scipy astropy astroplan pandas scikit-learn xgboost sqlalchemy gurobi


pip install sklearn_pandas 
pip install transitions

pip install -e git+https://github.com/ZwickyTransientFacility/ztf_sim.git#egg=ztf_sim
```

(To remove the environment, use `conda remove --name ztf_sim_test --all`.)


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

The simulator also writes a logfile, which may be monitored while the simulation is running with `tail -f example_ztf_schedule_log.txt`
