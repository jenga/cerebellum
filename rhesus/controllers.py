from numpy import arange
from rhesus.utils import SimWriter
from collections import OrderedDict

class Simulation(object):

    def __init__(self, dt=0.01, fn=None):
        self.dt = dt
        self.models = OrderedDict()
        self.states = OrderedDict()
        self.spikes = OrderedDict()
        self.fn = fn
        self.writer = SimWriter()

    def addModel(self, model):
        self.models[model.name] = model
        self.states[model.name] = []
        self.spikes[model.name] = []

    def integrate(self,t):
        for name, model in self.models.items():
            # need to pull state out of this somehow
            model.integrate(t, self.dt)
            self.states[name].append(model.vrs)
            if 'spike' in model.outputs:
                if model.outputs['spike'] == 1:
                    self.spikes[model.name].append(t)

    def simulate(self, length=1000, buffer=10000):
        index = 0
        if self.fn:
            self.writer.create(self.fn, self.models, self.spikes.keys(), length/self.dt)
        for i in arange(0, length, self.dt):
            index += 1
            self.integrate(i)
            if index % buffer == 0 and self.fn:
                self.save()
                for key, value in self.states.items():
                    self.states[key] = []
        if self.fn:
            if index % buffer != 0 and self.fn:
                self.save()
            for name, data in self.spikes.items():
                self.writer.writeOutputs(name, data)
            self.writer.finish()

    def save(self):
        for name, data in self.states.items():
            self.writer.writeModel(name, data)
        #self.writer.finish()
    
class Simulator(object):

    def __init__(self, sim):
        self.simulation = sim

