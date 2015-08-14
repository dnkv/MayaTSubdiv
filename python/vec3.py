# from https://sites.google.com/site/williamedwardscoder/vector-class-for-python
class Vec3(tuple):
    def __new__(cls,x,y,z):
        #subclassing a tuple, you have to create it in __new__ since its immutable
        return tuple.__new__(cls,(float(x),float(y),float(z)))
    def __init__(self,x,y,z):
        # nice to tag our own xyz fields onto it; these ought to be immutable,
        # but we don't protect them
        self.x,self.y,self.z = float(x),float(y),float(z)
    def __repr__(self):
        # here's the first example of using ourselves as an iterable
        return "%1.2f,%1.2f,%1.2f"%self
    def dot(self,other):
        return reduce(lambda a,b:a+b,self*other)
    def length(self):
        if not hasattr(self,"_l"): # cache it
            self._l = math.sqrt(self.x**2 + self.y**2 + self.z**2)    
        return self._l
    def angle(self,other):
        return math.acos(self.dot(other)/(self.length()*other.length()))
    def _map(self,other,op):
        try: # missing a nice isiterable() method
            return Vec3(*[op(a,b) for a,b in zip(self,other)])
        except:
            return Vec3(*[op(a,other) for a in self])
    def __add__(self,other):
        return self._map(other,lambda a,b:a+b)
    def __sub__(self,other):              
        return self._map(other,lambda a,b:a-b)
    def __mul__(self,other):              
        return self._map(other,lambda a,b:a*b)
    def __div__(self,other):
        return self._map(other,lambda a,b:a/b)
    def max(self,other):
        return self._map(other,lambda a,b:max(a,b))
    def min(self,other):
        return self._map(other,lambda a,b:min(a,b))
    def abs(self):
        return Vec3(*[abs(x) for x in self])
        
        