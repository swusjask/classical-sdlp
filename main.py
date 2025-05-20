#!/usr/bin/env python3
"""
Main script for testing classical algorithms on the Semidirect Discrete Logarithm Problem.
Implements BSGS for DLP and SDLP across different platform groups.
"""

import time
from random import randint
from sage.all import (
    GF,
    EllipticCurve,
    random_prime,
    diagonal_matrix,
    random_matrix,
    gcd,
    matrix,
)

from bsgs import bsgs_dlp, bsgs_sdlp
from period import find_period
from groups.finite_field import SemidirectProductZp
from groups.elementary_abelian import SemidirectProductEA


def test_finite_field(p=None):
    """Test DLP and SDLP on multiplicative group of finite fields."""
    print("=== Finite Field Fp* ===")
    if p is None:
        p = random_prime(2**35, lbound=2**34)
    print(f"Testing with p = {p}")

    # Test standard DLP in Fp*
    print("\n--- Standard DLP ---")
    Fp = GF(p)
    g = Fp.random_element()
    order = g.multiplicative_order()
    print(f"order of g = {order}")
    t = 100 * randint(1, order)
    h = g**t

    start_time = time.time()
    found_t = bsgs_dlp(Fp, h, g, order, operation="*")
    dlp_time = time.time() - start_time
    print(
        f"DLP: Found t = {found_t}, expected = {t % order}"
    )
    assert found_t == t % order
    print(f"DLP time: {dlp_time:.4f}s")

    # Test SDLP in Fp* ⋊ Aut(Fp*)
    print("\n--- SDLP ---")
    G = SemidirectProductZp(p)
    base_elem = G.random_element()

    # Find actual period
    print("Finding actual period...")
    rk = base_elem.x.multiplicative_order() * gcd(base_elem.x - 1, p - 1)
    period = find_period(base_elem, rk)  # Use upper bound from theorem
    print(f"Period r = {period}")
    assert (base_elem**period).g == G.one().g

    # Create SDLP instance
    t = 100 * randint(1, period - 1)
    target = base_elem**t

    # Set up u = (g, σ), v = (1, σ^-1)
    u = base_elem
    v = G(1, pow(base_elem.x, -1, p - 1))
    w = G.one()
    w.g = target.g

    start_time = time.time()
    found_t = bsgs_sdlp(G, w, (u, v), period)
    sdlp_time = time.time() - start_time
    print(
        f"SDLP: Found t = {found_t}, expected = {t % period}"
    )
    assert found_t == t % period
    print(f"SDLP time: {sdlp_time:.4f}s")

def test_elliptic_curve(p=None):
    """Test DLP on elliptic curves (SDLP is O(1) so excluded)."""
    print("\n=== Elliptic Curve E(F_p) ===")
    if p is None:
        p = random_prime(2**30, lbound=2**29)
    E = EllipticCurve(GF(p), [0, 1])  # y^2 = x^3 + 1
    print(f"Testing with p = {p}, curve order = {E.order()}")

    # Test standard DLP
    print("\n--- Standard DLP ---")
    P = E.random_element()
    while P.is_zero():
        P = E.random_element()
    order = P.order()
    print(f"order of P = {order}")
    t = 100 * randint(1, order - 1)
    Q = t * P

    start_time = time.time()
    found_t = bsgs_dlp(E, Q, P, order, operation="+")
    dlp_time = time.time() - start_time
    print(f"DLP: Found t = {found_t}, expected = {t % order}")
    assert found_t == t % order
    print(f"DLP time: {dlp_time:.4f}s")

    print("\n--- SDLP ---")
    print("SDLP: Skipped (O(1) complexity due to small automorphism group)")


def test_elementary_abelian(p=None, n=3):
    """Test SDLP on elementary abelian groups (DLP is O(1) so excluded)."""
    print("\n=== Elementary Abelian Group F_p^n ===")
    if p is None:
        p = random_prime(2**9, lbound=2**8)
    print(f"Testing with p = {p}, n = {n}")

    G = SemidirectProductEA(p, n)
    print("\n--- Standard DLP ---")
    print("DLP: Skipped (O(1) complexity in vector spaces)")

    # Test SDLP
    print("\n--- SDLP ---")

    # Create element with matrix having no eigenvalue 1
    print("\nCase 1: Matrix without eigenvalue 1")
    while True:
        A = G._M.random_element()
        if 1 not in matrix(A).eigenvalues():
            break

    base_elem = G(G._V.random_element(), A)

    # Find period
    matrix_order = A.order()
    period = find_period(base_elem, matrix_order)
    assert (base_elem**period).g == G.one().g
    print(f"Matrix order: {matrix_order}, Period: {period}")

    # Create SDLP instance
    t = 100 * randint(1, period - 1)
    target = base_elem**t

    u = base_elem
    v = G(G._V.zero(), A ** (-1))
    w = G.one()
    w.g = target.g

    start_time = time.time()
    found_t = bsgs_sdlp(G, w, (u, v), period)
    sdlp_time = time.time() - start_time
    print(f"SDLP: Found t = {found_t}, expected = {t % period}")
    assert found_t == t % period
    print(f"SDLP time: {sdlp_time:.4f}s")

    # Test with matrix having eigenvalue 1
    print("\nCase 2: Matrix with eigenvalue 1")
    # Create matrix with eigenvalue 1 in first position
    D = diagonal_matrix(GF(p), [1] + [randint(2, p - 1) for _ in range(n - 1)])
    P_mat = random_matrix(GF(p), n)
    while P_mat.det() == 0:
        P_mat = random_matrix(GF(p), n)
    A_with_one = G._M(P_mat ** (-1) * D * P_mat)

    base_elem = G(G._V.random_element(), A_with_one)

    # Period can be larger due to eigenvalue 1
    matrix_order = A_with_one.order()
    period = find_period(base_elem, p * matrix_order)
    print(f"Matrix order: {matrix_order}, Period: {period}")

    t = 100 * randint(1, period - 1)
    target = base_elem**t

    u = base_elem
    v = G(G._V.zero(), A_with_one ** (-1))
    w = G.one()
    w.g = target.g

    start_time = time.time()
    found_t = bsgs_sdlp(G, w, (u, v), period)
    sdlp_time = time.time() - start_time
    print(f"SDLP: Found t = {found_t}, expected = {t % period}")
    assert found_t == t % period
    print(f"SDLP time: {sdlp_time:.4f}s")


def main():
    """Run all tests comparing DLP and SDLP complexities."""
    print("Classical Algorithms for Semidirect Discrete Logarithm Problem")
    print("=" * 60)

    test_finite_field()
    test_elliptic_curve()
    test_elementary_abelian()

    print("\n" + "=" * 60)
    print("Summary for all bsgs implementation:")
    print("- Finite fields: DLP (O(√p)), SDLP (O(√p))")
    print("- Elliptic curves: DLP (O(√p)), SDLP (O(1))")
    print("- Elementary abelian: DLP (O(1)), SDLP (O(p^(n/2))")


if __name__ == "__main__":
    main()
