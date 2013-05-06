class Link:

    def __init__(self, name, pre, post, outtype, intype):
        self.name = name
        self.pre = pre;
        self.post = post;
        self.post.addInput(intype, self)
        self.type = outtype

    def delete(self):
        self.pre.removeOutput(self)
        self.post.removeInput(self)

    def transfer(self, dt):
        return self.pre.getOutput(self.type)
