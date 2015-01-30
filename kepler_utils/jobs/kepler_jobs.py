from celery import Celery
import math
import subprocess

from .generate import Generator, Simulation
from . import celeryconfig

app = Celery('kepler_jobs', backend = 'amqp', broker='amqp://guest@localhost//')
app.config_from_object (celeryconfig)


# This module is set up to use Celery to automate job running 
# To purge all queued jobs, use celery amqp queue.purge <QUEUE NAME>

@app.task
def run (name, original, run_location = '.', command = './kepler', force = False, query = None, tags = None, goal = "presn", **kwargs):
    generator = Generator (original)
    sim = Simulation (name, generator, run_location, command, force)
    try:
        p = sim.run (query = query, **kwargs)
        output = p.communicate () [0]
        ret = p.wait ()
    except TypeError as e:
        print (e)
        pass
    except NameError as e:
        print (e)
        pass
    except Exception as e:
        print (e)
        print (type (e))
        p.terminate ()
    sim.rebase (tags)
    latest = sim.getLatest ()
    if latest != goal:
        raise RuntimeError ("Did not run to goal")
    return latest