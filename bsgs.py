def bsgs_dlp(G, h, g, bound, operation="*"):
    """
    Baby-Step Giant-Step algorithm for discrete logarithm problem.
    Given a group G, elements h and g, and a bound, finds t such that:
    - g^t = h (if operation is "*")
    - t*g = h (if operation is "+")
    
    Parameters:
    G - The group
    h - The target element
    g - The base element
    bound - Upper bound for the discrete log
    operation - Either "*" for multiplicative or "+" for additive groups
    """
    from math import ceil, sqrt
    
    m = ceil(sqrt(bound))
    
    T = {}
    
    if operation == "*":
        temp = G.one()
        for i in range(m):
            T[temp] = i
            temp = temp * g
        
        gm_inv = g^(-m)
        temp = h
        for i in range(m):
            if temp in T:
                return T[temp] + i * m
            temp = temp * gm_inv
    
    elif operation == "+":
        temp = G.zero()
        for i in range(m):
            T[temp] = i
            temp = temp + g
        
        gm_neg = -m*g
        temp = h
        for i in range(m):
            if temp in T:
                return T[temp] + i * m
            temp = temp + gm_neg
    
    else:
        raise ValueError("Operation must be either '*' or '+'")
    
    raise ValueError(
        f"Log of {h} to the base {g} does not exist up to the bound {bound}."
    )

def bsgs_sdlp(G, w, base, bound):
    """
    Baby-Step Giant-Step algorithm for semidirect discrete logarithm problem.
    Given a group G, an element w, a base (u,v) and a bound, this function
    returns the discrete logarithm t such that u^t*v^t=w.
    
    Parameters:
    G - The semidirect product group
    w - The target element
    base - A tuple (u, v) where u corresponds to (g, sigma) and v to (1, sigma^-1)
    bound - Upper bound for the discrete log
    """
    from math import ceil, sqrt
    
    u, v = base
    m = ceil(sqrt(bound))
    
    T = {}
    temp = G.one()
    for j in range(m):
        T[temp] = j
        temp = u * temp * v
    
    um_inv = u^(-m)
    vm_inv = v^(-m)
    
    temp = w
    for i in range(m):
        if temp in T:
            # Found a match, return the discrete log
            return i * m + T[temp]
        temp = um_inv * temp * vm_inv
    
    raise ValueError(
        f"Log of {w} to the base {base} does not exist up to bound {bound}."
    )
