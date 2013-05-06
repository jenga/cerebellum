import h5py

class SimWriter(object):

    def __init__(self):
        self.h5file = None;
        self.models = None;
        self.indexes = {} 

    def create(self, fn, models, outputs, size):
        self.h5file = h5py.File(fn)
        self.models = self.h5file.create_group('Models')
        for key, model in models.items():
            self.indexes[key] = 0
            grp = self.models.create_group(model.name)
            for st in model.vrs.keys():
                grp.create_dataset(st, (size,), maxshape=(None,))
        self.outputs = self.h5file.create_group('Outputs')
    
    def writeModel(self, name, data):
        for key in data[0].keys():
            d = [st[key] for st in data]
            dset = self.models[name][key] 
            dset[self.indexes[name]:self.indexes[name]+len(d)] = d
        self.indexes[name] += len(d)

    def writeOutputs(self, name, data): 
        if (len(data) > 0):
            self.outputs.create_group(name)
            dset = self.outputs[name].create_dataset('Spikes', (len(data),))
            dset[...] = data
    
    def finish(self):
        self.h5file.close()
