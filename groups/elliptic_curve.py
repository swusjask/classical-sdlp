from sage.all import Parent, UniqueRepresentation, Groups
from sage.structure.element import Element

class SemidirectProductElementEC(Element):
    def __init__(self, parent, g, x):
        Element.__init__(self, parent)
        self.g = g
        self.x = x
        
    def _repr_(self):
        return f"({self.g}, {self.x})"
        
    def _mul_(self, other):
        g1, x1 = self.g, self.x
        g2, x2 = other.g, other.x
        new_g = g1 + x1(g2)
        new_x = x1 * x2
        
        return self.parent()(new_g, new_x)

    def __hash__(self):
        return hash((self.g, self.x))

    def __eq__(self, other):
        if not isinstance(other, SemidirectProductElementEC):
            return False
        return self.g == other.g and self.x == other.x
    
    def __invert__(self):
        new_x = self.x^-1
        new_g = -new_x(self.g)
        
        return self.parent()(new_g, new_x)

class SemidirectProductEC(Parent, UniqueRepresentation):
    Element = SemidirectProductElementEC
    
    def __init__(self, E):
        """
        Initialize a semidirect product group E(F_p) × Aut(E(F_p)).
        
        Args:
            E: An elliptic curve over a finite field
        """
        self._E = E
        self._base_ring = E.base_ring()
        self._p = self._base_ring.characteristic()
        
        Parent.__init__(self, category=Groups().Finite())
        
    def _element_constructor_(self, g, x):
        return self.element_class(self, g, x)
    
    def base_ring(self):
        return self._base_ring

    def order(self):
        return (self._E.order()) * len(self._E.automorphisms())
    
    def _repr_(self):
        return f"Semidirect product of E(F_{self._p}) x Aut(E(F_{self._p}))"
    
    def random_element(self):
        from random import choice
        g = self._E.random_element()
        x = choice(self._E.automorphisms())
        return self(g, x)
    
    def one(self):
        return self(self._E(0), self._E.automorphisms()[0])