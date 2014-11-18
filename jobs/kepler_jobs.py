from celery import Celery
import math
import subprocess

from . import generate
from . import celeryconfig

app = Celery('kepler_jobs', backend = 'amqp', broker='amqp://guest@localhost//')
app.config_from_object (celeryconfig)

@app.task
def run (name, original, run_location = '.', command = './kepler', force = False, query = None, tags = None, **kwargs):
    generator = generate.Generator (original)
    sim = generate.Simulation (name, generator, run_location, command, force)
    print ("Running with arguments:", kwargs)
    try:
        print ("Trying to run...")
        p = sim.run (query = query, **kwargs)
        print ("Waiting...")
        output = p.communicate()[0]
        ret = p.wait()
    except TypeError as e:
        print (e)
        pass
    except NameError as e:
        print (e)
        pass
    sim.rebase (tags)
    return None
