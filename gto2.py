import numpy as np
import math

# Evaluation method is disscussed in Taketa, Huzinaga, Ohata, 1966.

#============================================================
# Math Function 
#============================================================
def boys(n, x):
    import scipy.special as sp
    # x -> 0: the denominator goes to 0, zero deviation error will occur.
    # thus, have to return 0-limit value(x->0.)
    if x < 0.: raise
    radius = 0.0001
    if x < radius:  return 1. / (2.0 * n + 1.0)
    else:
        numerator =  sp.gammainc(n+0.5, x)  * sp.gamma(n + 0.5)
        denominator = 2.0 * math.pow(x, n + 0.5)
        return numerator /denominator

def norm2(array):
    return np.linalg.norm(array) ** 2.

def factorial(n):
    if n <= 1:  return 1
    else:       return n * factorial(n - 1)

def factorial2(n):
    if n <= 1:  return 1
    else:       return n * factorial2(n - 2)

def binomial(n, m):
    # return the nCm
    return factorial(n) / factorial(m) / factorial(n-m)

#============================================================
# Gaussian Calculation Function
#============================================================
#XXX These functions don't consider the normalize.
#   Normalize should be done by in the PGTO and CGTO class.

def binomial_prefactor(exponent, pow1, pow2, rx1, rx2):
    #   In 
    #       (x + cons1) ^ pow1 * (x + cons2) ^ pow2 ,
    #   return the coefficient of the [exponent]-degree term
    s = 0.
    for i in xrange(1 + exponent):
        j = exponent - i
        if i <= pow1 and j <= pow2:
            s += binomial(pow1, i) * math.pow(rx1, pow1-i) * binomial(pow2, j) * math.pow(rx2, pow2-j)
    return s

def product_GTO_center(exponent1, center1, exponent2, center2):
    gamma = exponent1 + exponent2
    return (exponent1 * center1 + exponent2 * center2) / gamma

def overlap_1D(exponent1, center1, l1, exponent2, center2, l2):
    s = 0.
    gamma = exponent1 + exponent2
    new_center = (exponent1 * center1 + exponent2 * center2) / gamma
    pa = new_center - center1
    pb = new_center - center2
    for i in xrange(1 + int(math.floor( (l1+l2)/2 )) ):
        s += binomial_prefactor(2*i, l1, l2, pa, pb) * factorial2(2*i-1) / math.pow(2*gamma, i)
    return s

def overlap_3D(exponent1, center1, l1, m1, n1, exponent2, center2, l2, m2, n2):
    # This function does not consider the norm.
    exp1, exp2 = exponent1, exponent2   # just use short name
    distance2 = norm2(center1 - center2)
    gamma = exponent1 + exponent2
    Rp = product_GTO_center(exp1, center1, exp2, center2)
    (Ra, Rb) = (Rp - center1, Rp - center2)
    prefactor = math.pow(math.pi/gamma, 1.5) * math.exp(-exp1*exp2*distance2/gamma)
    sx = overlap_1D(exp1, center1[0], l1, exp2, center2[0], l2)
    sy = overlap_1D(exp1, center1[1], m1, exp2, center2[1], m2)
    sz = overlap_1D(exp1, center1[2], n1, exp2, center2[2], n2)
    return prefactor * sx * sy * sz

#============================================================
# Primitive Gaussian
#============================================================
def product_PGTO_center(lhs, rhs):
    if not isinstance(lhs, primitiveGTO) or not isinstance(rhs, primitiveGTO): raise
    gamma = lhs.exponent + rhs.exponent
    return (lhs.exponent * lhs.center + rhs.exponent * rhs.center) / gamma

def overlap_PGTO(lhs, rhs):
    if not (isinstance(lhs, primitiveGTO) and isinstance(rhs, primitiveGTO)):   raise
    s = overlap_3D(
            lhs.exponent, lhs.center, lhs.l, lhs.m, lhs.n,
            rhs.exponent, rhs.center, rhs.l, rhs.m, rhs.n)
    return s * lhs.norm * rhs.norm

def kinetic_PGTO(pgto1, pgto2):
    if not (isinstance(pgto1, primitiveGTO) and isinstance(pgto2, primitiveGTO)):   raise
    exp1, exp2 = pgto1.exponent, pgto2.exponent
    Ra, Rb = pgto1.center, pgto2.center
    (l1, m1, n1) = pgto1.l, pgto1.m, pgto1.n
    (l2, m2, n2) = pgto2.l, pgto2.m, pgto2.n
    term1 = exp2 * (2*(l2+m2+n2)+3) * overlap_3D(exp1, Ra, l1, m1, n1, exp2, Rb, l2, m2, n2)
    term2 = -2 * math.pow(exp2, 2) * (
            overlap_3D(exp1, Ra, l1, m1, n1, exp2, Rb, l2 + 2, m2, n2) +
            overlap_3D(exp1, Ra, l1, m1, n1, exp2, Rb, l2, m2 + 2, n2) + 
            overlap_3D(exp1, Ra, l1, m1, n1, exp2, Rb, l2, m2, n2 + 2) ) 
    term3 = -0.5 * (
            l2*(l2-1) * overlap_3D(exp1, Ra, l1, m1, n1, exp2, Rb, l2-2, m2, n2) + \
            m2*(m2-1) * overlap_3D(exp1, Ra, l1, m1, n1, exp2, Rb, l2, m2-2, n2) + \
            n2*(n2-1) * overlap_3D(exp1, Ra, l1, m1, n1, exp2, Rb, l2, m2, n2-2) )
    return (term1 + term2 + term3) * pgto1.norm * pgto2.norm

def four_center_PGTO(pgto1, pgto2, pgto3, pgto4):
    #   (12|34) = integrate {1(1)2(1) 1/r 3(2)4(2)} dr1 dr2
    pass

def a_term(i, r, u, l1, l2, PA_x, PB_x, CP_x, gamma):
    # i, r, u: index. see THO paper.
    #  pa_x: product_center -> center1
    #  pb_x: product_center -> center2
    #  pc_x: product_center -> Nuclear
    pre = math.pow(-1, i) * binomial_prefactor(i, l1, l2, PA_x, PB_x)
    numerator = math.pow(-1, u) * factorial(i) * math.pow(CP_x, i - 2*r - 2*u) * math.pow(0.25*gamma, r+u)
    denominator = factorial(r) * factorial(u) * factorial(i - 2*r - 2*u)
    return pre * numerator / denominator

def a_term_reorder(l1, l2, PA_x, PB_x, CP_x, gamma):
    a = [0.]*(1+l1+l2)
    for i in xrange(1 + l1 + l2):
        for r in xrange(1 + int(math.floor(i / 2.0))):
            for u in xrange(1 + int(math.floor((i-2.0*r)/2.0))):
                I = i - 2*r - u # XXX
                a[I] += a_term(i, r, u, l1, l2, PA_x, PB_x, CP_x, gamma)
    return a

def nuclear_attraction_PGTO(pgto1, pgto2, nuclear):
    # < pgto1 | Z/Rzx | pgto2 >
    if not (isinstance(pgto1, primitiveGTO) and (isinstance(pgto2, primitiveGTO))): raise
    (l1, m1, n1) = (pgto1.l, pgto1.m, pgto1.n)
    (l2, m2, n2) = (pgto2.l, pgto2.m, pgto2.n)
    (exp1, exp2) = (pgto1.exponent, pgto2.exponent)
    gamma = exp1 + exp2
    dist2 = norm2(pgto1.center - pgto2.center)
    product_center = product_PGTO_center(pgto1, pgto2)
    pc2 = norm2(nuclear.pos - product_center)
    pre = 2.0 * math.pi / gamma * math.exp(-exp1 * exp2 * dist2 / gamma)

    PA = product_center - pgto1.center
    PB = product_center - pgto2.center
    CP = product_center - nuclear.pos
    a_x = a_term_reorder(l1, l2, PA[0], PB[0], CP[0], gamma)
    a_y = a_term_reorder(m1, m2, PA[1], PB[1], CP[1], gamma)
    a_z = a_term_reorder(n1, n2, PA[2], PB[2], CP[2], gamma)
    s = 0.
    for I in xrange(1+l1+l2):
        for J in xrange(1+m1+m2):
            for K in xrange(1+n1+n2):
                s += a_x[I] * a_y[J] * a_z[K] * boys(I+J+K, gamma*pc2)
    s = s * pre * pgto1.norm * pgto2.norm
    #print "l1: {}, l2: {}, center1: {}, center2: {}, product_center: {}, nuclear: {} => {}".format(l1, l2, pgto1.center, pgto2.center, product_center, nuclear.pos, s)
    return s

def b_term(i1, i2, r1, r2, u, l1, l2, Ax, Bx, Px, gamma1, l3, l4, Cx, Dx, Qx, gamma2):
    PAx = Px - Ax
    PBx = Px - Bx
    QCx = Qx - Cx
    QDx = Qx - Dx
    px  = Qx - Px
    delta = (0.25*gamma1) + (0.25*gamma2)
    term1 = math.pow(-1, i2) * binomial_prefactor(i1, l1, l2, PAx, PBx) * binomial_prefactor(i2, l3, l4, QCx, QDx)
    term2_numerator  = factorial(i1) * factorial(i2) * math.pow(4*gamma1, r1) * math.pow(4*gamma2, r2) * math.pow(delta, r1+r2)
    term2_denominator= math.pow(4*gamma1, i1) * math.pow(4*gamma2, i2) * math.pow(delta, i1+i2) * factorial(r1) * factorial(r2) * factorial(i1-2*r1) * factorial(i2-2*r2)
    term3_numerator = factorial(i1+i2-2*(r1+r2)) * math.pow(-1, u) * math.pow(px, i1+i2-2*(r1+r2)-2*u)*math.pow(delta, u)
    term3_denominator= factorial(u)*factorial(i1+i2-2*(r1+r2)-2*u)
    return term1 * term2_numerator / term2_denominator * term3_numerator / term3_denominator

def b_term_reorder(l1, l2, Ax, Bx, Px, gamma1, l3, l4, Cx, Dx, Qx, gamma2):
    b_term_array = [0.] * (1+l1+l2+l3+l4)
    for i1 in xrange(1+l1+l2):
        for i2 in xrange(1+l3+l4):
            for r1 in xrange(1 + int(math.floor(i1/2.0))):
                for r2 in xrange(1 + int(math.floor(i2/2.0))):
                    for u in xrange(1+int(math.floor((i1+i2)/2.0 - r1-r2))):
                        I = i1+i2-2*(r1+r2)-u
                        b_term_array[I] += b_term(i1, i2, r1, r2, u, l1, l2, Ax, Bx, Px, gamma1, l3, l4, Cx, Dx, Qx, gamma2)
    return b_term_array

def electron_repulsion_PGTO(pgto1, pgto2, pgto3, pgto4):
    P = product_PGTO_center(pgto1, pgto2)
    Q = product_PGTO_center(pgto3, pgto4)
    PQ2 = nom2(Q - P)
    AB2 = norm2(pgto1.center - pgto2.center)
    CD2 = norm2(pgto3.center - pgto4.center)
    gamma1 = pgto1.exponent + pgto2.exponent
    gamma2 = pgto3.exponent + pgto4.exponent
    delta = (0.25*gamma1) + (0.25*gamma2)
    b_array_x = b_term_reorder(pgto1.l, pgto2.l, pgto1.center[0], pgto2.center[0], P[0], gamma1, pgto3.l, pgto4.l, pgto3.center[0], pgto4.center[0], Q[0], gamma2)
    b_array_y = b_term_reorder(pgto1.m, pgto2.m, pgto1.center[1], pgto2.center[1], P[1], gamma1, pgto3.m, pgto4.m, pgto3.center[1], pgto4.center[1], Q[1], gamma2)
    b_array_z = b_term_reorder(pgto1.n, pgto2.n, pgto1.center[2], pgto2.center[2], P[2], gamma1, pgto3.n, pgto4.n, pgto3.center[2], pgto4.center[2], Q[2], gamma2)
    prefactor1= 2 * math.pow(math.pi, 2) / gamma1 / gamma2 * math.sqrt(math.pi / (gamma1 + gamma2))
    prefactor2= math.exp(-(pgto1.exponent*pgto2.exponent*AB2/gamma1) - (pgto3.exponent* pgto4.exponent*CD2/gamma2))
    s = 0.
    for I in (1 + pgto1.l + pgto2.l + pgto3.l + pgto4.l):
        for J in (1 + pgto1.m + pgto2.m + pgto3.m + pgto4.m):
            for K in (1 + pgto1.n + pgto2.n + pgto3.n + pgto4.n):
                s += b_array_x[I] * b_array_y[J] * b_array_z[K] * boys(I+J+K, PQ2/4/delta)
    return prefactor1 * prefactor2 * s * pgto1.norm * pgto2.norm * pgto3.norm * pgto4.norm

class primitiveGTO:
    def __init__(self, exponent, l, m, n, center_x, center_y, center_z):
        self.exponent = exponent
        (self.l, self.m, self.n)  = (l, m, n)
        self.center = np.zeros(3)
        self.center[0] = center_x
        self.center[1] = center_y
        self.center[2] = center_z
        self.norm = 1.
        self.normalize()
    def normalize(self):
        (l,m,n) = (self.l, self.m, self.m)
        alpha = self.exponent
        numerator = math.pow(2, 2*(l+m+n)+1.5) * math.pow(alpha, l+m+n+1.5) 
        denominator = factorial2(2*l-1) * factorial2(2*m-1) * factorial2(2*n-1) * math.pow(math.pi, 1.5)
        self.norm = math.sqrt(numerator / denominator)
    def value(self, x, y, z):
        pos = np.zeros(3)
        pos[0], pos[1], pos[2] = x, y, z
        return math.exp((-1) * self.exponent * norm2(pos - self.center)) * self.norm
    def center(self):
        return self.center

#============================================================
# Contracted Gaussian
#============================================================
def overlap_CGTO(lhs, rhs):
    s = 0.
    for (coeff1, pgto1) in lhs.pgtos:
        for (coeff2, pgto2) in rhs.pgtos:
            s += coeff1 * coeff2 * overlap_PGTO(pgto1, pgto2)
    s *= lhs.norm * rhs.norm
    return s

def kinetic_CGTO(cgto1, cgto2):
    t = 0.
    for (coeff1, pgto1) in cgto1.pgtos:
        for (coeff2, pgto2) in cgto2.pgtos:
            t += coeff1 * coeff2 * kinetic_PGTO(pgto1, pgto2)
    return t * cgto1.norm * cgto2.norm

def nuclear_attraction_CGTO(cgto1, cgto2, nuclear):
    if not (isinstance(cgto1, contractedGTO) and isinstance(cgto2, contractedGTO)): raise
    if not (isinstance(nuclear, Atom)): raise
    t = 0.
    for (coeff1, pgto1) in cgto1.pgtos:
        for (coeff2, pgto2) in cgto2.pgtos:
            t += coeff1 * coeff2 * nuclear_attraction_PGTO(pgto1, pgto2, nuclear)
    return -nuclear.atomic_number * t * cgto1.norm * cgto2.norm 

def electron_repulsion_CGTO(cgto1, cgto2, cgto3, cgto4):
    s = 0.
    for (coeff1, pgto1) in cgto1.pgtos:
        for (coeff2, pgto2) in cgto2.pgtos:
            for(coeff3, pgto3) in cgto3.pgtos:
                for(coeff4, pgto4) in cgto4.pgtos:
                    s += coeff1 * coeff2 * coeff3 * coeff4 * electron_repulsion_PGTO(pgto1, pgto2, pgto3, pgto4)
    return s * cgto1.norm * cgto2.norm * cgto3.norm * cgto4.norm

class contractedGTO:
    def __init__(self, l, m, n, center_x, center_y, center_z):
        self.center = np.zeros(3)
        self.center[0] = center_x
        self.center[1] = center_y
        self.center[2] = center_z
        (self.l, self.m, self.n)  = (l, m, n)
        self.norm = 1.
        self.pgtos = list()
        
    def add_primitiveGTO(self, coeff, exponent):
        pgto = primitiveGTO(exponent, self.l, self.m, self.n, 
                self.center[0], self.center[1], self.center[2] )
        self.pgtos.append( (coeff, pgto) )

    def normalize(self):
        self.norm = 1.
        s = overlap_CGTO(self, self)
        self.norm = 1. / math.sqrt(s)
    
    def value(self, x, y, z):
        s = 0.
        for (coeff, pgto) in self.pgtos:
            s += coeff * pgto.value(x, y, z)
        return s * self.norm

class Atom:
    def __init__(self, pos_x, pos_y, pos_z, atomic_number):
        self.pos = np.zeros(3)
        self.pos[0] = pos_x
        self.pos[1] = pos_y
        self.pos[2] = pos_z
        self.atomic_number = atomic_number

def compute_overlap(bfs):
    dim = len(bfs)
    S = np.zeros( (dim, dim) )
    for i in xrange(dim):
        for j in xrange(dim):
            S[i, j] = overlap_CGTO(bfs[i], bfs[j])
    return S

def compute_T(bfs):
    dim = len(bfs)
    T = np.zeros( (dim, dim) )
    for i in xrange(dim):
        for j in xrange(dim):
            T[i,j] = kinetic_CGTO(bfs[i], bfs[j])
    return T

def compute_K(bfs, atoms):
    dim = len(bfs)
    H = np.zeros( (dim, dim) )
    for atom in atoms:
        v = np.zeros( (dim, dim) )
        #print atom.atomic_number
        for i in xrange(dim):
            for j in xrange(dim):
                v[i,j] = nuclear_attraction_CGTO(bfs[i], bfs[j], atom)
        #print v
        H += v
    return H
    
h1 = contractedGTO(0, 0, 0, 0., 0., 0.)
h1.add_primitiveGTO(0.444635, 0.168856)
h1.add_primitiveGTO(0.535328, 0.623913)
h1.add_primitiveGTO(0.154329, 3.42525)
h1.normalize()

h3 = contractedGTO(0, 0, 0, 1.4, 0., 0.)
h3.add_primitiveGTO(0.444635, 0.168856)
h3.add_primitiveGTO(0.535328, 0.623913)
h3.add_primitiveGTO(0.154329, 3.42525)
h3.normalize()

bfs = [h1, h3]

atoms = list()
atoms.append(Atom(0., 0., 0., 1.0) )
atoms.append(Atom(1.4, 0., 0.0, 1.0) )

S = compute_overlap(bfs)
T = compute_T(bfs)
H =  compute_K(bfs, atoms)

dim = len(bfs)
Cinit = np.zeros( (dim, dim) )
for i in xrange(dim):
    Cinit[i,i] = 1
print Cinit
print S
print T + H
