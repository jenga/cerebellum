from rhesus.models.node import Node

class IClamp(Node):

    def __init__(self, name):
        super().__init__(name)
        self.outputs = {"I": 0}
        self.on = 0.0
        self.dur = 0.0
        self.amp = 0.0
        self.delta = 0.0
        self.freq = 0.0

    def integrate(self, t, dt):
        if (t > self.on):
            self.outputs["I"] = self.amp
        if (t < self.on or t >= self.on + self.dur):
            self.outputs["I"] = 0.0
