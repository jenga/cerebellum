
class Node(object):

    def __init__(self, name):
        self.name = name
        self.inputs = {}
        self.outputs = {}
        self.vrs = {}

    def integrate(self, t, dt):
        '''
        ODEs will be integrated in this function for one timestep of dt. If
        there are inputs, they should be checked and handled within this
        function
        
        '''
        return 0

    def addInput(self, cat, link):
        if cat in self.inputs:
            if self.inputs[cat].count(link) > 0:
                return True
            self.inputs[cat].append(link)
            return True 
        else:
            return False

    def removeInput(self, cat, link):
        if (cat in self.inputs):
            if (self.inputs[cat].count(link) > 0):
                self.inputs[cat].remove(link)
                return True
        return False

    def addOutput(self, cat, link):
        if cat in self.outputs:
            if self.outputs[cat].count(link) > 0:
                return True
            self.outputs[cat].append(link)
            return True
        else:
            return False

    def removeOutput(self, cat, link):
        if cat in self.outputs:
            if self.outputs[cat].count(link) > 0:
                self.outputs[cat].remove(link)
                return True
        else:
            return False

    def getOutput(self, type):
        return self.outputs[type]
