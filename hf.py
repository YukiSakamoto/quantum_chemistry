import numpy as np
import math
# Hartree Fock Method

Debug_enable = True

def debug_output(value):
    if Debug_enable:
        print value
    return value

# SC = FCe


atomic_orbitals = [
    # (Function, center_position)
]

def norm(start, end):
    #dX = end[0] - start[0]
    #dY = end[1] - start[1]
    #dZ = end[2] - start[2]
    #return math.sqrt( dX ** 2 + dY ** 2 + dZ ** 2 )
    displacement = np.array([
        end[0] - start[0],
        end[1] - start[1],
        end[2] - start[2],
    ])
    return np.linalg.norm(displacement)

#============================================================
#   Atomic Orbitals Related
#============================================================
class SlaterTypeFunction:
    def __init__(self, center, exponent, coeff = 1):
        self.center = center
        self.exponent = exponent
        self.coeff = coeff
    def value(self, pos):
        r = norm(pos, self.center)
        return self.coeff * math.exp(self.exponent * r)
class GaussinTypeFunction:
    def __init__(self, center, exponent, coeff = 1):
        pass
    

def gaussian_function_product(orbit_exponent1, centre1, orbit_exponent2, centre2):
    # Just store short named variables.
    Ra = centre1
    Rb = centre2
    a = orbit_exponent1 
    b = orbit_exponent2

    r = a + b
    Rc = (a * Ra + b * Rb) / (a + b)
    k = ((2 / math.pi) ** 2) * ((a * b)**(3/4)) * math.exp( (-a * b/(a + b)) * ((Ra - Rb) ** 2) )
    return (r, Rc, k)

#============================================================
#   Calculate Integrals
#============================================================
def overlap_integral():
    # S_ij = <phi_i | phi_j>
    return None
    
def double_integral():
    # <psi | 1/r | psi>
    return None

def single_integral():
    # <psi | h | psi> = <psi | K + Vne | psi> = <psi | K | psi> + <psi | Vne | psi>
    return None

def initial_guess(dimention):
    return np.zeros((dimention, dimention))

#============================================================
#   SCF Procedure
#============================================================
def rmsdp(D,oldD):
    delta = 0.0
    length = 5
    for i in range (0,length):
        for j in range(0,length):
            delta += (D[i,j] - oldD[i,j]) ** 2
    delta = (delta / (length * length)) ** 0.5
    return delta

def scf():
    threshold = 0.00001
    dimention = 5
    D = initial_guess(dimention)

    S = overlap_integral()
    # F = Hcore + G
    
    while True:
        oldD = D
        update_matrix()
        current_rmsdp = rmsdp(D,oldD)
        if current_rmsdp < threshold:
            break
    return C

a = np.matrix([[1,2,3], [4,5,6], [7,8,9]])
b = np.matrix([[1,2,3], [4,5,6], [7,8,9]])

#print a.dot(b)

#c = np.matrix([[2,1],[1,1]])

#c = np.matrix([[1,0,0], [0,1,0], [0,0,1]])
#print np.linalg.inv(c)

nuclears = [
        ((0.0, 0.0, 0.0), 1.0 ), 
        ((1.0, 0.0, 0.0), 2.0 ),
        ]

f = SlaterTypeFunction( (0.0, 0.0, 0.0), 2.0, 1.0 )
print f.value((1.0, 0.0, 0.0) )

