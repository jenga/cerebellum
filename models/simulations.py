from rhesus.controllers import Simulation
from models.models import  Noise, DCN, IClamp, PC, ComboNoise, PCST, PCSpont
from models.links import InhSynapse, NoiseLink, INoise, ISource

class PCDCNSimulation(Simulation):

    def __init__(self, dt=0.01, fn=None, gca=4.25, gih=0.8, corr=0.5):
        super().__init__(dt=dt, fn=fn)
        nc = Noise('common')
        dcn = DCN('dcn', gca=gca, gih=gih)
        icl = IClamp('dcn_i', on=0, dur=0, amp=4.0)
        ISource('dcn_icl', icl, dcn)
        noises = []
        combos = []
        self.addModel(nc)
        self.addModel(icl)
        self.addModel(dcn)
        for i in range(0,10):
            noise = Noise('n' + str(i))
            noises.append(noise)
            temp = ComboNoise('cn' + str(i), mu=0.6, c=corr)
            NoiseLink('ln' + str(i) + 'c', nc, temp)
            NoiseLink('ln' + str(i), noise, temp)
            combos.append(temp)
            self.addModel(noise)
            self.addModel(temp)
            pc = PC('pc' + str(i))
            INoise('in' + str(i), temp, pc)
            InhSynapse('inh' + str(i), pc, dcn)
            self.addModel(pc)

class PCSTSimulation(Simulation):

    def __init__(self, dt=0.01, fn=None, gca=4.25, gih=0.8):
        super().__init__(dt=dt, fn=fn)
        dcn = DCN('dcn', gca=gca, gih=gih)
        icl = IClamp('dcn_i', on=0, dur=0, amp=4.0)
        ISource('dcn_icl', icl, dcn)
        pc = PCST('pc')
        InhSynapse('inh', pc, dcn)
        self.addModel(icl)
        self.addModel(dcn)
        self.addModel(pc)

class PCST10Simulation(Simulation):

   def __init__(self, dt=0.01, fn=None, gca=4.25, gih=0.8):
       super().__init__(dt=dt, fn=fn)
       dcn = DCN('dcn', gca=gca, gih=gih)
       icl = IClamp('dcn_i', on=0, dur=0, amp=4.0)
       ISource('dcn_icl', icl, dcn)
       self.addModel(icl)
       self.addModel(dcn)
       for i in range(0,10):
           pc = PCST('pc'+str(i))
           InhSynapse('inh'+str(i), pc, dcn)
           self.addModel(pc)

class PCSTIM10Simulation(Simulation):

   def __init__(self, dt=0.01, fn=None, gca=4.25, gih=0.8):
       super().__init__(dt=dt, fn=fn)
       dcn = DCN('dcn', gca=gca, gih=gih)
       icl = IClamp('dcn_i', on=0, dur=0, amp=4.0)
       ISource('dcn_icl', icl, dcn)
       self.addModel(icl)
       self.addModel(dcn)
       for i in range(0,10):
           pc = PCSpont('pc'+str(i))
           InhSynapse('inh'+str(i), pc, dcn)
           self.addModel(pc)		  
