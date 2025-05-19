from sage.all import Parent, UniqueRepresentation, Groups
from sage.structure.element import Element

class SemidirectProductElementEC(Element):
    def __init__(self, parent, g, x):
        Element.__init__(self, parent)
        self.g = g
        self.x = x
        
    def _repr_(self):
        return f"({self.g}, {self.x.rational_maps()})"
        
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
        new_x = self.x**-1
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

if __name__ == "__main__":
    from sage.all import EllipticCurve, GF, random_prime

    p = random_prime(2 ** 30)
    E = EllipticCurve(GF(p), [0, 1])  # y^2 = x^3 + 1
    G = SemidirectProductEC(E)
    print(f"Created group: {G}")
    print(f"Elliptic curve: {E}")
    print(f"j-invariant: {E.j_invariant()}")
    print(f"Automorphism group size: {len(E.automorphisms())}")
    print(f"Curve order: {E.order()}\n")
    
    # Generate random element with non identity automorphism
    while True:
        a = G.random_element()
        if a.x != G.one().x:
            break
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
    
    # Verify theorem 9, all possible order of autormophism always divide 24
    rk = 24
    print(f"s_g,x({rk}) == 1: {a**(rk) == G.one()}")

        