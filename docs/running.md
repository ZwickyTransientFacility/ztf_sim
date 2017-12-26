# Running ztf_sim

Once you have a run configuration file, save it in the `ztf_sim/sims/` directory.  For clarity, we will use the packaged `test_config.json` as an example.

A run configuration file is also required; see `ztf_sim/config/default.cfg`.

To run the simulation, change to `ztf_sim/ztf_sim`.  If you are planning a long run (10 days or more), you may want to start a `screen`, `tmux`, or vnc session if you are running a remote server. 

Next, start `ipython`, then import the main observing script:

        In [1]: from ztf_sim.simulate import simulate
    
Optionally, start %pdb for interactive debugging if there's a crash: 

	In [2]: %pdb
	
Now start the observing code by passing it the name of the configuration files:

	In [3]: simulate("test_config.json", run_config_file='default.cfg')

If all goes well, the program will start running.  You can watch the progress by in the log file, which is `ztf_sim/sims/<run_name>_log.txt`.  `tail -f <filename>` will print the lines to your terminal as they are written to the file.
