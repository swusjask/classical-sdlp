from sage.all import GF, GL, VectorSpace, Parent, UniqueRepresentation, Groups
from sage.structure.element import Element

class SemidirectProductElementEAG(Element):
    def __init__(self, parent, g, x):
        Element.__init__(self, parent)
        self.g = g
        self.x = x
        
    def _repr_(self):
        return f"({self.g}, {self.x})"
        
    def _mul_(self, other):
        g1, x1 = self.g, self.x
        g2, x2 = other.g, other.x
        new_g = g1 + x1 * g2
        new_x = x1 * x2
        
        return self.parent()(new_g, new_x)

    def __hash__(self):
        return hash((self.g, self.x))

    def __eq__(self, other):
        if not isinstance(other, SemidirectProductElementEAG):
            return False
        return self.g == other.g and self.x == other.x
    
    def __invert__(self):
        new_x = self.x^-1
        new_g = -(new_x * self.g)
        
        return self.parent()(new_g, new_x)

class SemidirectProductEAG(Parent, UniqueRepresentation):
    Element = SemidirectProductElementEAG
    
    def __init__(self, p, n):
        """
        Initialize a semidirect product group F_p^n × GL(n, F_p).
        
        Args:
            p: A prime number for the base field
            n: Dimension of the vector space
        """
        self._p = p
        self._n = n
        self._base_ring = GF(p)
        self._V = VectorSpace(self._base_ring, n)
        self._M = GL(n, self._base_ring)
        
        Parent.__init__(self, category=Groups().Finite())
        
    def _element_constructor_(self, g, x):
        return self.element_class(self, g, x)
    
    def base_ring(self):
        return self._base_ring

    def order(self):
        return (self._p^self._n) * self._M.order()
    
    def _repr_(self):
        return f"Semidirect product of Z_{self._p}^{self._n} x Aut(Z_{self._p}^{self._n})"
    
    def random_element(self):
        g = self._V.random_element()
        x = self._M.random_element()
        return self(g, x)
    
    def one(self):
        return self(self._V(0), self._M.one())