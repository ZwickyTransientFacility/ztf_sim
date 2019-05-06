from ztf_sim.Scheduler import Scheduler
from astropy.time import Time

def test_list_queue():
    op_config = '/Users/ebellm/Eric/ZTF/survey_simulator/ztf_sim/sims/test_config.json'
    run_config_fullpath='/Users/ebellm/Eric/ZTF/survey_simulator/ztf_sim/config/default.cfg'
    scheduler = Scheduler(op_config,run_config_fullpath)
    scheduler.set_queue_manager(queue_manager='list')
    lod = [{'field_id':3881,'program_id':1,'subprogram_name':'test','filter_id':2},{'field_id':311,'program_id':1,'subprogram_name':'test','filter_id':1},{'field_id':312,'program_id':1,'subprogram_name':'test','filter_id':2}]
    scheduler.Q.load_list_queue(lod)
    assert( len(scheduler.Q.queue) == 3)
    scheduler.Q.load_list_queue(lod,append=True)
    assert( len(scheduler.Q.queue) == 6)
    #scheduler.Q.queue
    #scheduler.Q.next_obs({'current_time':Time.now()})

