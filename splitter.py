# Découper par noeud ou par élément un maillage triangulaire en calcul parallèle
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpi4py import MPI

from hilbertcurve import HilbertCurve

OX = 0
OY = 1

def mid_point( p1, p2 ):
    return ( 0.5*(p1[OX]+p2[OX]), 0.5*(p1[OY]+p2[OY]) )

def compute_bounding_box( points : np.array ):
    assert(len(np.shape(points)) == 2)
    
    crd_min = ( np.min(points[:,OX]), np.min(points[:,OY]) )
    crd_max = ( np.max(points[:,OX]), np.max(points[:,OY]) )

    return [ crd_min, crd_max ]

def compute_morton_ordering( vertices : np.array, bbox, N : int ) ->     np.array :
    hilbert_curve = HilbertCurve(N, 2)
    pN : int = (1<<N) - 1
    lgth = [ bbox[1][OX] - bbox[0][OX], bbox[1][OY] - bbox[0][OY]]
    return np.array([ [iVert, 
                       hilbert_curve.distance_from_point([int(pN*(vertices[iVert,OX]-bbox[0][OX])/lgth[OX]),
                                                          int(pN*(vertices[iVert,OY]-bbox[0][OY])/lgth[OY])])]
                    for iVert in range(vertices.shape[0])])

def split_node_mesh():
    return

if __name__ == '__main__':
    import mesh
    import visu_split_mesh as VSM

    