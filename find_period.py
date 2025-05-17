def find_period(g, n, G):
    """
    Finds the period of the cycle X_{g,sigma} for a given element (g, sigma)
    in a semidirect product group.
    
    Args:
        g: base element g
        n: An upper bound for the period
        
    Returns:
        The period r of the cycle
    """
    u = G.one()
    u.g = g
    length = n
    for p, e in length.factor():
        while e > 0:
            l = length // p
            if u**l == u.parent().one():
                length //= p
                e -= 1
            else:
                break
    
    return length