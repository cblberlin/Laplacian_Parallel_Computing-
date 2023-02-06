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

def split_element_mesh( nb_domains : int, elt2vertices : np.array, vertices : np.array ):
    """
    Usage : splitter.split_element_mesh( nbDoms, el2vertices, vertices )
    où nbDoms est le nombre de sous-domaine découpant le domaine initial,
    elt2vertices la connectivité donnant à chaque triangle l'indice des sommets correspondant (dans le sens direct),
    vertices les coordonnées des sommets stockées par point
    """
    nbVerts : int = vertices.shape[0]
    nbElts  : int = elt2vertices.shape[0]
    bbox = compute_bounding_box( vertices )
    # Calcul du barycentre de chaque triangle du maillage :
    bary_coords = np.array( [ (vertices[elt2vertices[iElt,0],:] + vertices[elt2vertices[iElt,1],:] + vertices[elt2vertices[iElt,2],:])/3.
                            for iElt in range(nbElts)])
    morton_numbering = compute_morton_ordering( bary_coords, bbox, 20)
    sort_indices = np.argsort(morton_numbering[:,1])
    morton_numbering = morton_numbering[sort_indices,:]
    nb_loc_elts : int = nbElts//nb_domains
    nb_suppl_elts = nbElts % nb_domains
    splitted_elements = []
    start_indices : int = 0
    for iDom in range(nb_domains):
        nb_elts = nb_loc_elts + (1 if iDom < nb_suppl_elts else 0)
        splitted_elements.append(np.array(morton_numbering[start_indices:start_indices+nb_elts,0]))
        start_indices += nb_elts
    return splitted_elements

if __name__ == '__main__':
    import mesh
    import visu_split_mesh as VSM

    nb_domains = 4
    msh, cl = mesh.read("CarrePetit.msh")
    vert2elts = msh.comp_vertices_to_elements()

    globCom = MPI.COMM_WORLD.Dup()
    nbp     = globCom.size
    rank    = globCom.rank

    print("Dissection du domaine à partir de ses éléments")
    splitted_elt = split_element_mesh(nb_domains, msh.elt2verts, msh.vertices)

    print(f"Domaine {rank} : nombre d'éléments locaux {splitted_elt[rank].shape[0]}")
    print(f"Indices globaux des éléments : {splitted_elt[rank]}")

    


    