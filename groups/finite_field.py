from sage.all import IntegerModRing, Zmod, ZZ, gcd, random_prime, euler_phi, Groups, Parent, UniqueRepresentation
from sage.structure.element import Element


class SemidirectProductElementZp(Element):
    def __init__(self, parent, g, x):
        Element.__init__(self, parent)
        self.g = parent._base_ring(g)
        self.x = Zmod(parent._p - 1)(x)
        
    def _repr_(self):
        return f"({self.g}, {self.x})"
        
    def _mul_(self, other):
        g1, x1 = self.g, self.x
        g2, x2 = other.g, other.x
        new_g = ZZ(g1 * g2**x1)
        new_x = (x1 * x2)
        
        return self.parent()(new_g, new_x)
    
    def __hash__(self):
        return hash((self.g, self.x))

    def __eq__(self, other):
        if not isinstance(other, SemidirectProductElementZp):
            return False
        return self.g == other.g and self.x == other.x
    
    def __invert__(self):
        new_x = pow(self.x, -1, self.parent()._p-1)
        new_g = pow(self.g, -new_x, self.parent()._p)
        
        return self.parent()(new_g, new_x)

class SemidirectProductZp(Parent, UniqueRepresentation):
    Element = SemidirectProductElementZp
    
    def __init__(self, p):
        """
        Initialize a semidirect product group F_p* × Aut(F_p*).
        
        Args:
            p: prime modulo p
        """
        self._p = p
        self._base_ring = IntegerModRing(p)
        
        Parent.__init__(self, category=Groups().Finite())
        
    def _element_constructor_(self, g, x):
        return self.element_class(self, g, x)
    
    def base_ring(self):
        return self._base_ring

    def order(self):
        return (self._p - 1) * euler_phi(self._p - 1)
    
    def _repr_(self):
        return f"Semidirect product of Z_{self._p}* x Aut(Z_{self._p}*)"
    
    def random_element(self):
        from random import randrange
        g = randrange(1, self._p)
        while True:
            x = randrange(1, self._p-1)
            if gcd(x, self._p-1) == 1:
                break
        return self(g, x)
    
    def one(self):
        return self(1, 1)

if __name__ == "__main__":
    p = random_prime(2 ** 30)
    G = SemidirectProductZp(p)
    print(f"Created group: {G}\n")
    
    # Generate random element
    a = G.random_element()
    print(f"Random element: {a}\n")
    
    # Test group properties
    print(f"Group order: {G.order()}")
    print(f"Identity element: {G.one()}\n")
    
    # Test invertibility
    a_inv = ~a
    print(f"Inverse of {a} is {a_inv}")
    print(f"a * a^(-1) == 1: {a * a_inv == G.one()}\n")
    
    # Generate another element and test multiplication
    b = G.random_element()
    print(f"Another random element: {b}")
    print(f"a * b = {a * b}\n")
    
    # Verify associativity with a small test
    c = G.random_element()
    print(f"(a * b) * c == a * (b * c): {(a * b) * c == a * (b * c)}\n")

    # Verify lagrange theorem
    print(f"a^{G.order()} == 1: {a**G.order() == G.one()}\n")

    # Verify theorem 3
    rk = a.x.multiplicative_order() * gcd(a.x - 1, G._p - 1)
    print(f"s_g,x({rk}) == 1: {a**(rk) == G.one()}")