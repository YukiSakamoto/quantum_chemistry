import numpy as np
import math

# Should move to elsewhere that contain math utils.

def fact(n):
    if n <= 1:
        return 1
    else:
        return n * fact(n - 1)

def fact2(n):
    if n <= 1:
        return 1
    else:
        return n * fact(n - 2)

def binomial(n, m):
    return fact(n) / (fact(m) * fact(n - m))

class PrimitiveGTO:
    def __init__(self, exponent, l, m, n, center):
        self.exponent = exponent
        self.l = l
        self.m = m
        self.n = n
        self.center = center
        self.norm = 1.
        self.normalize()

    def normalize(self):
        (l,m,n) = (self.l, self.m, self.m)
        alpha = self.exponent
        numerator = math.pow(2, 2*(l+m+n)+1.5) * math.pow(alpha, l+m+n+1.5) 
        denominator = fact2(2*l-1) * fact2(2*m-1) * fact2(2*n-1) * math.pow(math.pi, 1.5)
        self.norm = math.sqrt(numerator / denominator)

    def overlap(self, other):
        alpha = self.exponent
        beta = other.exponent
        norm = np.linalg.norm(self.center - other.center)
        v1 = math.pow( math.pi /(alpha + beta), 1.50)
        v2 = math.exp( (-1) * alpha * beta / (alpha + beta) * math.pow(norm, 2.0))
        return v1 * v2 * self.norm * other.norm

    def kinetic(self, other):
        # See Szabo.
        alpha = self.exponent
        beta = other.exponent
        norm = np.linalg.norm(self.center - other.center)
        v1 = alpha * beta / (alpha + beta) * (3 - 2 * alpha * beta / (alpha + beta) * math.pow(norm, 2.0) )
        v2 = math.pow( math.pi / (alpha + beta), 1.5)
        v3 = math.exp( -alpha * beta / (alpha + beta) * pow(norm, 2.0) )
        # Now, only S-S orbitals .
        return v1 * v2 * v3 * self.norm * other.norm
        
    def as_string(self):
        return "P_GTO(exp= -{0:>10}, norm= {1:>10})".format(self.exponent, self.norm)

class ContractedGTO:
    def __init__(self, l, m, n, center_x, center_y, center_z):
        self.center = np.zeros(3)
        self.center[0] = center_x
        self.center[1] = center_y
        self.center[2] = center_z
        self.l = l
        self.m = m
        self.n = n
        self.pgtos = list()
        self.norm = 1.0

    def add_pgto(self, coeff, exponent):
        pgto = PrimitiveGTO(exponent, self.l, self.m, self.n, self.center)
        self.pgtos.append( (coeff, pgto) )

    def overlap(self, other):
        # <self | other>
        S = 0.
        for (self_coeff, self_pgto) in self.pgtos:
            for (other_coeff, other_pgto) in other.pgtos:
                S = S + self_coeff * other_coeff * self_pgto.overlap(other_pgto)
        S = S * self.norm * other.norm
        return S 
    def kinetic(self, other):
        # < self | -1/2 d/dr | other>
        T = 0.
        for(self_coeff, self_pgto) in self.pgtos:
            for (other_coeff, other_pgto) in other.pgtos:
                T = T + self_coeff * other_coeff * self_pgto.kinetic(other_pgto)
        return T * self.norm * other.norm

    def normalize(self):
        S = 0.
        self.norm = 1.
        S = self.overlap(self)
        self.norm = 1. / math.sqrt(S)

    def as_string(self):
        s = ""
        first = True
        s += "Center= {0:>10}, Norm={0:>10}\n".format(self.center, self.norm)
        s += "{0:>15}   {1}\n".format("coefficient", "primitive-GTO")
        for (coeff, pgto) in self.pgtos:
            s += "{0:>15}   {1}\n".format(coeff, pgto.as_string() )
        return s


def compute_overlap(bfs):
    dim = len(bfs)
    S = np.zeros( (dim, dim) )
    for i in xrange(dim):
        for j in xrange(dim):
            S[i, j] = bfs[i].overlap(bfs[j])
    return S

def compute_T(bfs):
    dim = len(bfs)
    T = np.zeros( (dim, dim) )
    for i in xrange(dim):
        for j in xrange(dim):
            T[i,j] = bfs[i].kinetic(bfs[j])
    return T

def compute_coulomb(bfs):
    dim = len(bfs)
    J = np.zeros( (dim, dim, dim, dim) )
    for i in xrange(dim):
        for j in xrange(dim):
            for k in xrange(dim):
                for l in xrange(dim):
                    pass
    return J


h1 = ContractedGTO(0, 0, 0, 0., 0., 0.)
h1.add_pgto(0.444635, 0.168856)
h1.add_pgto(0.535328, 0.623913)
h1.add_pgto(0.154329, 3.42525)
h1.normalize()
print h1.as_string()

h3 = ContractedGTO(0, 0, 0, 1.4, 0., 0.)
h3.add_pgto(0.444635, 0.168856)
h3.add_pgto(0.535328, 0.623913)
h3.add_pgto(0.154329, 3.42525)
h3.normalize()
print h3.as_string()

bfs = [h1, h3]
S = compute_overlap(bfs)
T = compute_T(bfs)
print S
print T


print math.erf(2.56)
