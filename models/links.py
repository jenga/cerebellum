from rhesus.models.link import Link

class InhSynapse(Link):

    def __init__(self, name, pre, post, g_tau=2.5, a=1.0, f_inc=0.0, d_inc=-0.59,
            f_init=1.00, d_init=1.0, f_tau=1.0, d_tau=77.67):
        super().__init__(name, pre, post, 'spike', 'inh')
        self.g = 0
        self.f_init = f_init
        self.d_init = d_init
        self.f = f_init
        self.d = d_init
        self.g_tau = g_tau
        self.a = a
        self.f_inc = f_inc
        self.d_inc = d_inc
        self.f_tau = f_tau
        self.d_tau = d_tau
    
    def transfer(self, dt):
        s = self.pre.getOutput(self.type)
        self.g = self.g + (-self.g / self.g_tau) * dt + self.a * self.f * self.d * s
        if s > 0.0:
            self.f = self.f + self.f_inc * (1 - self.f)
            self.d = self.d + self.d_inc * self.d
        self.f = self.f + ((self.f_init - self.f) / self.f_tau) * dt
        self.d = self.d + ((self.d_init - self.d) / self.d_tau) * dt
        return self.g*0.5

class ISource(Link):

    def __init__(self, name, pre, post):
        super().__init__(name, pre, post, 'I', 'isource')
    
    def transfer(self, dt):
        return self.pre.getOutput(self.type)

class NoiseLink(Link):

    def __init__(self, name, pre, post):
        super().__init__(name, pre, post, 's', 's')
    
    def transfer(self, dt):
        return self.pre.getOutput(self.type)

class INoise(Link):
    
    def __init__(self, name, pre, post):
        super().__init__(name, pre, post, 'c', 'isource')

    def transfer(self, dt):
        return self.pre.getOutput(self.type)
    
