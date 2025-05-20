from sage.all import ZZ

def find_period(u, n):
    """
    Finds the period of the cycle X_{g,sigma} for a given element (g, sigma) := u
    in a semidirect product group.
    
    Args:
        g: base element g
        n: An upper bound for the period
        
    Returns:
        The period r of the cycle
    """
    length = ZZ(n)
    for p, e in length.factor():
        while e > 0:
            l = length // p
            if (u**l).g == u.parent().one().g:
                length //= p
                e -= 1
            else:
                break
    
    return length