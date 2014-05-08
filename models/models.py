from rhesus.models.node import Node
from numpy import exp, pi

class PC(Node):

    '''Define model constants'''
    C = 1.5
    Cd = 1.5
    R = 0.75
    g_Na = 40.0
    g_K = 8.75
    g_Kd = 12.0
    g_IH = 0.03
    g_leak = 0.032
    n_tau = 15.0
    ih_tau = 100.0
    g_inh = 1.0
    g_exc = 0.0
    
    E_Na = 45.0
    E_K = -95.0
    E_leak = -77.0
    E_IH = -20.0
    E_exc = 0.0
    E_inh = -77.0

    def __init__(self, name):
        super().__init__(name)
        self.outputs = {'spike': 0}
        self.inputs = {'isource': [], 'exc': [], 'inh': []}
        self.vrs = {'Vs': -65, 'Vd': -65, 'h': 0, 'ih': 0, 'nd': 0}
        self.spiked = False

    def integrate(self, t, dt):
        v = self.vrs
        isources = self.inputs['isource']
        curr = 0.0
        for i in isources:
            curr = curr + i.transfer(dt)
        exc = 0.0
        for e in self.inputs['exc']:
            exc = exc + (e.transfer(dt) * (v['Vd'] - self.E_exc))
        inh = 0.0
        for i in self.inputs['inh']:
            inh = inh + (self.g_inh * i.transfer(dt) * (v['Vd'] - self.E_inh))
        dv = (((v['Vd'] - v['Vs'])) / self.R + curr -
                self.g_Na * self.m_inf(v['Vs']) * v['h'] * (v['Vs'] - self.E_Na) - 
                self.g_K * (1-v['h']) * (v['Vs'] - self.E_K) - 
                self.g_leak * (v['Vs'] - self.E_leak) - 
                self.g_IH * v['ih'] * (v['Vs'] - self.E_IH)) / self.C
        dvd = (((v['Vs'] - v['Vd']))/self.R - exc - inh -
                self.g_leak*(v['Vd'] - self.E_leak) - 
                self.g_Kd*v['nd']*(v['Vd'] - self.E_K))/self.Cd
        dh = (self.h_inf(v['Vs']) - v['h'])/self.h_tau(v['Vs'])
        dih = (self.ih_inf(v['Vs']) - v['ih'])/self.ih_tau
        dnd = (self.n_inf(v['Vd']) - v['nd'])/self.n_tau
        self.vrs = {'Vs': self.vrs['Vs'] + dv*dt,
                    'Vd': self.vrs['Vd'] + dvd*dt,
                    'h': self.vrs['h'] + dh*dt,
                    'ih': self.vrs['ih'] + dih*dt,
                    'nd': self.vrs['nd'] + dnd*dt}
        if not self.spiked and self.vrs['Vs'] > 0:
            self.spiked = True
            self.outputs['spike'] = 1
        elif self.spiked and self.vrs['Vs'] > 0:
            self.outputs['spike'] = 0
        elif self.spiked and self.vrs['Vs'] <= 0:
            self.spiked = False

    def m_inf(self, vm):
        return 1/(1+exp((vm - (-40))/-3))

    def h_inf(self, vm):
        return 1/(1+exp((vm - (-40))/3))

    def ih_inf(self, vm):
        return 1/(1+exp((vm - (-80))/3))

    def n_inf(self, vm):
        return 1/(1+exp((vm - (-35))/-3))

    def h_tau(self, vm):
        return 295.4/(4.0*pow(vm + 50, 2) + 400.0) + 0.012


class DCN(Node):

    E_K = -95.0
    E_leak = -77.0
    E_Ca = 128.0
    E_IH = -20.0
    E_Na = 45.0
    E_NaCa = -40.0
    E_exc = 0.0
    E_inh = -77.0
    
    g_K = 17.5
    g_Kd = 25.0
    g_Na = 58.0
    g_Ca = 4.25
    g_leak = 0.1
    g_IH = 0.8
    g_NaCa = 0.05

    tau_m = 0.09
    tau_n = 0.6
    tau_ihs = 75.0
    tau_ihf = 20.0
    tau_mt = 7.0
    tau_ht = 37.0
    tau_nd = 25.0

    g = 1.5
    kap = 0.5
    C = 3.0
    Cd = 3.0

    def __init__(self, name, gca=4.25, gih=0.8):
        super().__init__(name)
        self.outputs = {'spike': 0}
        self.inputs = {'isource': [], 'exc': [], 'inh': []}
        self.vrs = {'mt': 0.024517160789596117, 'Vs': -46.433814101677193, 'Vd': -49.858712712574359, 'ihf': 0.032103453134627931, 'ihs': 0.034306183881818941, 'h': 0.95868834690418592, 'nd': 0.00089539845217886688, 'm': 0.0067449017530451256, 'ht': 0.16175366140431852, 'n': 2.7053891749215154e-05, 'inh': 0.0, 'inh_g': 0.0}
        self.spiked = False
        self.g_Ca = gca
        self.g_IH = gih

    def integrate(self, t, dt):
        v = self.vrs
        curr = 0
        for i in self.inputs['isource']:
            curr = curr + i.transfer(dt)
        exc = 0.0
        for e in self.inputs['exc']:
            exc = exc + (e.transfer(dt) * (v['Vd'] - self.E_exc))
        inh = 0.0
        inh_g = 0.0
        for i in self.inputs['inh']:
            temp = i.transfer(dt)
            inh = inh + (temp * (v['Vd'] - self.E_inh))
            inh_g += temp
        dvs = (((v['Vd'] - v['Vs']) * self.g) / self.kap + curr -
                self.g_Na * v['m'] * v['h'] * (v['Vs'] - self.E_Na) -
                self.g_K * v['n'] * (v['Vs'] - self.E_K) -
                self.g_IH * v['ihs'] * v['ihf'] * (v['Vs'] - self.E_IH) -
                self.g_NaCa * (v['Vs'] - self.E_NaCa) -
                self.g_leak * (v['Vs'] - self.E_leak)) / self.C
        dvd = (((v['Vs'] - v['Vd']) * self.g) / (1 - self.kap) - exc - inh -
             self.g_Kd * v['nd'] * (v['Vd'] - self.E_K) - 
             self.g_Ca * v['mt'] * v['ht'] * (v['Vd'] - self.E_Ca) -
             self.g_leak * (v['Vd'] - self.E_leak)) / self.Cd
        dm = ((self.m_inf(v['Vs']) - v['m']) / self.tau_m)
        dh = ((self.h_inf(v['Vs']) - v['h']) / self.tau_h(v['Vs']))
        dn = ((self.n_inf(v['Vs']) - v['n']) / self.tau_n)
        dihs = ((self.ih_inf(v['Vs']) - v['ihs']) / self.tau_ihs)
        dihf = ((self.ih_inf(v['Vs']) - v['ihf']) / self.tau_ihf)
        dmt = ((self.mt_inf(v['Vd']) - v['mt']) / self.tau_mt)
        dht = ((self.ht_inf(v['Vd']) - v['ht']) / self.tau_ht)
        dnd = ((self.nd_inf(v['Vd']) - v['nd']) / self.tau_nd)
        self.vrs = {'Vs': self.vrs['Vs'] + dvs * dt,
                    'Vd': self.vrs['Vd'] + dvd * dt,
                    'm': self.vrs['m'] + dm * dt,
                    'h': self.vrs['h'] + dh * dt,
                    'n': self.vrs['n'] + dn * dt,
                    'ihs': self.vrs['ihs'] + dihs * dt,
                    'ihf': self.vrs['ihf'] + dihf * dt,
                    'mt': self.vrs['mt'] + dmt * dt,
                    'ht': self.vrs['ht'] + dht * dt,
                    'nd': self.vrs['nd'] + dnd * dt,
                    'inh': inh,
                    'inh_g': inh_g}
        if not self.spiked and self.vrs['Vs'] > 0:
            self.spiked = True
            self.outputs['spike'] = 1
        elif self.spiked and self.vrs['Vs'] > 0:
            self.outputs['spike'] = 0
        elif self.spiked and self.vrs['Vs'] <= 0:
            self.spiked = False

    def m_inf(self, v):
        return 1/(1 + exp(-(v + 33) / 2.8))

    def h_inf(self, v):
        return 1/(1 + exp(-(v + 40) / -2.9))

    def tau_h(self, v):
        return (((2 * 232)/ pi) * (20 / (2 * (v + 39) ** 2 + 20 ** 2))) * 0.2

    def n_inf(self, v):
        return 1 / (1 + exp(-(v + 20) / 2.7))

    def ih_inf(self, v):
        return 1 / (1 + exp((v + 92.1) / 11.4))

    def mt_inf(self, v):
        return 1 / (1 + exp((v + 31.3) / -4.5))

    def ht_inf(self, v):
        return 1 / (1 + exp((v + 63.8) / 6.9))

    def nd_inf(self, v):
        return 1 / (1 + exp(-(v + 30) / 3.5))

class IClamp(Node):

    def __init__(self, name, on=100.0, dur=200.0, amp=1.0, delta=0.0, freq=0.0):
        super().__init__(name)
        self.outputs = {'I': 0}
        self.vrs = {'I': 0, 'other': 0}
        self.on = on
        self.dur = dur
        self.amp = amp
        self.delta = delta
        self.freq = freq

    def integrate(self, t, dt):
        if (self.on == 0 and self.dur == 0):
            self.outputs = {'I': self.amp}
            self.vrs = {'I': self.amp}
        else: 
            if (t > self.on and t <= self.on + self.dur):
                self.outputs = {'I': self.amp}
                self.vrs = {'I': self.amp}
            else:
                self.outputs = {'I': 0}
                self.vrs = {'I': 0}

from math import sqrt
import random

class Noise(Node):

    def __init__(self, name, mu=0, sd=1, tau=5, corr_mu=30, corr_sd=100):
        super().__init__(name)
        self.outputs = {'s': 0}
        self.vrs = {'s': 0}
        self.mu = mu
        self.sd = sd
        self.tau = tau
        self.corr_mu = corr_mu
        self.corr_sd = corr_sd

    def integrate(self, t, dt):
        s = self.vrs['s']
        ds = -(s / self.tau) * dt + sqrt(2.0 / self.tau) * random.gauss(self.mu, self.sd) * sqrt(dt)
        self.vrs = {'s': s + ds}
        self.outputs = {'s': s + ds}

class ComboNoise(Node):

    def __init__(self, name, mu=0, sd=1, c=0.75):
        super().__init__(name)
        self.outputs = {'c': 0}
        self.inputs = {'s': []}
        self.vrs = {'s': 0, 'c': 0}
        self.mu = mu
        self.sd = sd
        self.c = c

    def integrate(self, t, dt):
        c = self.mu + self.sd * (sqrt(self.c) * self.inputs['s'][0].transfer(dt) + sqrt(1 - self.c) * self.inputs['s'][1].transfer(dt))
        self.vrs = {'c': c}
        self.outputs = {'c': c}

class PCST(Node):

    def __init__(self,name):
        super().__init__(name)
        self.outputs = {'spike': 0}
        self.inputs = {'isource': [], 'exc': [], 'inh': []}
        self.vrs = {'Vs':0}
        with open('Spont1.txt','r') as f:
             spikes=f.readline()
        spikes.strip('/n')
        spikes=spikes.split('\t')
        spikes=[float(x) for x in spikes]
        spikes=[x*1000 for x in spikes]
        self.spiketimes=[int(x) for x in spikes]
     
    def integrate(self,t,dt):
        if t in self.spiketimes:
           self.outputs['spike'] = 1
        else:
           self.outputs['spike'] = 0

class PCSpont1(Node):

    def __init__(self,name):
        super().__init__(name)
        self.outputs = {'spike': 0}
        self.inputs = {'isource': [], 'exc': [], 'inh': []}
        self.vrs = {'Vs':0}
        with open('Spont1.txt','r') as f:
             spikes=f.readline()
        spikes.strip('/n')
        spikes=spikes.split('\t')
        spikes=[float(x) for x in spikes]
        spikes=[x*1000 for x in spikes]
        self.spiketimes=[int(x) for x in spikes]

    def integrate(self,t,dt):
        if t in self.spiketimes:
           self.outputs['spike'] = 1
        else:
           self.outputs['spike'] = 0

class PCStimC(Node):

    def __init__(self,name):
        super().__init__(name)
        self.outputs = {'spike': 0}
        self.inputs = {'isource': [], 'exc': [], 'inh': []}
        self.vrs = {'Vs':0}
        with open('Stim_C.txt','r') as f:
             spikes=f.readline()
        spikes.strip('/n')
        spikes=spikes.split('\t')
        spikes=[float(x) for x in spikes]
        spikes=[x*1000 for x in spikes]
        self.spiketimes=[int(x) for x in spikes]

    def integrate(self,t,dt):
        if t in self.spiketimes:
           self.outputs['spike'] = 1
        else:
           self.outputs['spike'] = 0

class PCStimE(Node):

    def __init__(self,name):
        super().__init__(name)
        self.outputs = {'spike': 0}
        self.inputs = {'isource': [], 'exc': [], 'inh': []}
        self.vrs = {'Vs':0}
        with open('Stim_E.txt','r') as f:
             spikes=f.readline()
        spikes.strip('/n')
        spikes=spikes.split('\t')
        spikes=[float(x) for x in spikes]
        spikes=[x*1000 for x in spikes]
        self.spiketimes=[int(x) for x in spikes]

    def integrate(self,t,dt):
        if t in self.spiketimes:
           self.outputs['spike'] = 1
        else:
           self.outputs['spike'] = 0

