from rhesus.controllers import Simulation
from models.models import PC, DCN
from models.links import InhSynapse

pc1 = PC('PC1')
dcn1 = DCN('DCN1')
inh1 = InhSynapse('INH', pc1, dcn1)
sim = Simulation()
sim.addModel(pc1)
sim.addModel(dcn1)
sim.simulate()
vm = [st['Vs'] for st in sim.states['PC1']]
plot(vm)
