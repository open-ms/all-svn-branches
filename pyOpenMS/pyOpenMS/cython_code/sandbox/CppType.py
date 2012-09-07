class CppType(object):

    def __init__(self, base, targs=None):
        self.base = base
        self.targs = targs if targs is not None else []

    def __str__(self):
        
        if self.targs:
            return self.base + "[%s]" % (", ".join(str(t) for t in self.targs))
        else:
            return self.base

    def identifier(self):
        return str(self).replace("[", "_brl_").replace("]", "_brr_")

