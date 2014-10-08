from celery import Celery
import math
import subprocess

from . import generate
from . import celeryconfig

app = Celery('kepler_jobs', backend = 'amqp', broker='amqp://guest@localhost//')
app.config_from_object (celeryconfig)

@app.task
def run (name, original, run_location = '.', command = './kepler', force = False, **kwargs):
    print ("RUN")
    generator = generate.Generator (original)
    sim = generate.Simulation (name, generator, run_location, command, force)
    print (kwargs)
    p = sim.run (**kwargs)
    output = p.communicate()[0] 
    ret = p.wait()
    sim.rebase ()
    return ret
