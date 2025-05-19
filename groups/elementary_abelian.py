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
        return hash((tuple(self.g), self.x))

    def __eq__(self, other):
        if not isinstance(other, SemidirectProductElementEAG):
            return False
        return self.g == other.g and self.x == other.x
    
    def __invert__(self):
        new_x = self.x**-1
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
        return (self._p**self._n) * self._M.order()
    
    def _repr_(self):
        return f"Semidirect product of F_{self._p}^{self._n} x Aut(F_{self._p}^{self._n})"
    
    def random_element(self):
        g = self._V.random_element()
        x = self._M.random_element()
        return self(g, x)
    
    def one(self):
        return self(self._V(0), self._M.one())

if __name__ == "__main__":
    from sage.all import random_prime, random_matrix, matrix, randint
    
    # Initialize test parameters
    p = random_prime(2 ** 30) 
    n = 7  # Vector space dimension
    G = SemidirectProductEAG(p, n)
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

    # Test case 1: Matrix WITHOUT eigenvalue 1 (general case)
    print("Testing with matrix WITHOUT eigenvalue 1:")
    # Create matrix with eigenvalues > 1
    diag_elements = [randint(2,p-1) for _ in range(n)]
    D = matrix.diagonal(GF(p), diag_elements)
    
    # Apply similarity transformation to make it look random
    while True:
        P = random_matrix(GF(p), n)
        if P.det() != 0:
            break
    A = G._M(P.inverse() * D * P)
    
    # Update our test element with this specific matrix
    a = G(G._V.random_element(), A)
    
    # For matrices without eigenvalue 1, the period divides the matrix order
    rk = a.x.order()
    print(f"Matrix order: {rk}")
    print(f"Matrix eigenvalues: {D.diagonal()}")
    print(f"s_g,x({rk}) == 1: {a**(rk) == G.one()}\n")
    
    # Test case 2: Matrix WITH eigenvalue 1 (worse case bound)
    print("Testing with matrix WITH eigenvalue 1:")
    # Create matrix with first eigenvalue = 1
    diag_elements_with_one = [1] + [randint(2, p-1) for _ in range(n-1)]
    D_with_one = matrix.diagonal(GF(p), diag_elements_with_one)
    
    # Apply similarity transformation
    while True:
        P2 = random_matrix(GF(p), n)
        if P2.det() != 0:
            break
    A_with_one = G._M(P2.inverse() * D_with_one * P2)
    
    # Update our test element with this matrix that has eigenvalue 1
    a = G(G._V.random_element(), A_with_one)
    
    # For matrices with eigenvalue 1, the period can be as large as p * matrix_order
    rk = p * a.x.order()
    print(f"Matrix order: {a.x.order()}")
    print(f"Matrix eigenvalues: {D_with_one.diagonal()}")
    print(f"s_g,x({rk}) == 1: {a**(rk) == G.one()}")
