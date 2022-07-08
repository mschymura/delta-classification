#########################################################################################################
# the code below is an adaptation of the code for mixed volume classification by Christopher Borger
# https://github.com/christopherborger/mixed_volume_classification/blob/master/volume_classification.sage
#########################################################################################################

load("polytopes.sage")

import logging
import os.path
import sys

# Using the logging package one can conveniently turn off and on the auxiliary messages  

logging.basicConfig(format='%(message)s',stream=sys.stdout,level=logging.INFO)
# After modifying the level from, say, logging.INFO to logging.WARNING , the change will come into force only after _restarting the sage session_ and reloading

# Sandwich is a pair of centrally symmetric lattice polytopes A,B with A being a subset of B.
# For the sake of efficiency, A also comes with its "symmetry-broken" part halfA such that A = halfA \cup -halfA \cup {0}.
# The gap of a sandwich A,B is the difference |B \cap Z^d| - |A \cap Z^d| of the number of integer points in B and A.

# that's the template for names of files, in which we store polytopes
FILE_NAME_DELTA = 'data/dim_%d_delta_%d.txt'
FILE_NAME_DELTA_EXTR = 'data/dim_%d_delta_%d_extremal.txt'


def prepare_sandwiches(m,Delta):
    for basisA in delta_normal_forms(m,Delta):
        # first, we generate A and halfA out of basisA
        mbA = matrix(basisA)
        mA = mbA.augment(-mbA)
        A = Polyhedron(mA.transpose())
        halfA = break_symmetry(A,m)
    
        # second, the outer container B is the centrally symmetric parallelotope spanned by the vectors in basisA
        B = polytopes.parallelotope(mA.transpose())
        
        # B may contain some integral points that are Delta-too-large with respect to A, and so we do:
        B = reduce_sandwich([halfA,A],B,Delta)
        yield [halfA,A],B


def break_symmetry(A,m):
    """
    	takes a centrally symmetric m-dimensional polytope A
    	computes a subset halfA of its vertices I such that I = conv(halfA \cup -halfA)
    """
    vertA = [vector(z) for z in A.vertices_list()]
    halfA = []
    for l in vertA:
    	if (-l in halfA):
    		continue
    	halfA.append(l)
    return halfA


def is_extendable(S,v,Delta):
    """
        Check whether the extension of a set S of vectors by a vector v causes a determinant to exceed Delta.
    """
    m = len(v)
    for C in Combinations(S,m-1):
    	M = matrix(C + [list(v)])
    	if abs(det(M)) > Delta:
    		return false    
    return true


def reduce_sandwich(A,B,Delta):
    """
        For a given sandwich (A,B) and a value of Delta
        the function returns a polytope
        obtained by removing all of the lattice points v of B 
        with the property that if v is added to A, there will be a determinant of absolute value > Delta
    """
    to_be_removed = []
    to_be_kept = []
    for v in B.integral_points():
    	if v in A[1]:
    		continue
    	if (v in to_be_removed or v in to_be_kept):	## this just avoids considering -w in case that w was considered already before
    		continue
    	if is_extendable(A[0],v,Delta):
    		to_be_kept.append(vector(v))
    		to_be_kept.append(-vector(v))
    	else:
    		to_be_removed.append(vector(v))
    		to_be_removed.append(-vector(v))
    Z = [vector(z) for z in B.integral_points()]
    return Polyhedron([z for z in Z if z not in to_be_removed])


def layered_polytope_from_sandwich(A,B):
    """ 3*B is embedded into height 0, two copies of 3*A are embedded into heights 1 and -1.
        Then, one generates a polytope based on these three layers at heights -1,0 and 1
        Note: If A and B are centrally symmetric, then the resulting polytope is centrally symmetric as well.
    """ 
    middleLayer = [tuple(3*vector(v))+(0,) for v in B.vertices()]
    upperLayer = [tuple(3*vector(v))+(1,) for v in A[1].vertices()]
    lowerLayer = [tuple(3*vector(v))+(-1,) for v in A[1].vertices()]
    return Polyhedron(middleLayer+upperLayer+lowerLayer)


def sandwich_normal_form(A,B):
    """
        returns data that allows to distinguish two sandwiches (A,B) 
        (A',B') up to affine unimodular transformations.
    """
#    return affine_normal_form(layered_polytope_from_sandwich(A,B))
#    return layered_polytope_from_sandwich(A,B).lattice_polytope().normal_form(algorithm='palp')
    return layered_polytope_from_sandwich(A,B).lattice_polytope().normal_form(algorithm='palp_native')	# 'palp_native' brings in slower sage implementation
#    return layered_polytope_from_sandwich(A,B).lattice_polytope().normal_form(algorithm='palp_modified')	# 'palp_modified' brings in modified PALP implementation


# Sandwich factory is used to store sandwiches up to affine unimodular transformations.
# A sandwich factory is a dictionary of dictionaries. For each possible gap, a storage
# for sandwiches with this gap is created. The latter storage
# is a dictionary with key,value pairs such that the value is a sandwich and 
# the respective key is the sandwich normal form of this sandwich.


def append_sandwich(sf,A,B):
    """
        If no affine unimodular image of the sandwich (A,B) is in the sandwich factory sf,
        the sandwich (A,B) is appended to sf.
    """
    Gap = B.integral_points_count() - A[1].integral_points_count()
    SNF = sandwich_normal_form(A,B)
    if Gap not in sf.keys():
        sf[Gap] = {}
    if SNF not in sf[Gap].keys():
        sf[Gap][SNF] = [A,B]

            
def new_sandwich_factory(m,Delta):
    sandwich_factory = {}
    for A,B in prepare_sandwiches(m,Delta):
        append_sandwich(sandwich_factory,A,B)
    return sandwich_factory


def sandwich_factory_statistics(sf):
    logging.info("Maximum gap in sandwiches: %d",max(sf.keys()))
    logging.info("Number of sandwiches: %d",sum([len(sf[Gap]) for Gap in sf.keys() if Gap!=0]))
    if 0 in sf.keys():
        logging.info("Number of polytopes found: %d", len(sf[0]))
    logging.info(50*"-")


def delta_classification(m,Delta,extremal):
    """
        runs the sandwich factory algorithm and classifies all centrally symmetric m-dimensional lattice polytopes with largest determinant equal to Delta
        extremal is a Boolean parameter determining whether the whole classification is sought [extremal=false], or only the classification of the extremal examples attaining h(Delta,m) [extremal=true]
    """
    sf = new_sandwich_factory(m,Delta)
    maxGap = max(sf.keys())
    
    # set the known lower bound for h(Delta,m) by Lee et al. 
    if (extremal):
        cmax = m^2 - m + 1 *2*m*Delta
    
    while maxGap > 0:
        
        sandwich_factory_statistics(sf)
        
        for SNF in sf[maxGap].keys():
            A,B = sf[maxGap][SNF]
            
            for v in B.vertices(): # pick a vertex of B which is not in A
                if v not in A[1]:
                    break
            
            blow_up_of_A = Polyhedron(list(A[1].vertices()) + [vector(v)] + [-vector(v)])	## this uses that all points in B are "Delta-ok" for A
            half_of_blow_up_of_A = break_symmetry(blow_up_of_A,m)
            reduction_of_B = Polyhedron([z for z in B.integral_points() if (vector(z) != vector(v) and vector(z) != -vector(v))])
            
            newA = [half_of_blow_up_of_A,blow_up_of_A]
            red_sand = reduce_sandwich(newA,B,Delta)
            if (extremal):
                if (red_sand.integral_points_count() >= cmax):
                    append_sandwich(sf,newA,red_sand)
                    npts_blow_up = blow_up_of_A.integral_points_count()
                    if (npts_blow_up > cmax):
                        cmax = npts_blow_up
                if (reduction_of_B.integral_points_count() >= cmax):
                    append_sandwich(sf,A,reduction_of_B)
            else:
                append_sandwich(sf,newA,red_sand)
                append_sandwich(sf,A,reduction_of_B)
            
        del sf[maxGap]
        maxGap = max(sf.keys())

    sandwich_factory_statistics(sf)
        
    result = []
    for A,B in sf[0].values():
        result.append(A[1])	## only store the polytope in A

    return result


def update_delta_classification_database(m,Delta,extremal):
    # the files storing polytopes are created in the data subfolder
    if not os.path.exists('data'):
        os.mkdir('data')

    # let's see whether the file for the pair (m,Delta) is missing
    if (extremal):
        missingDelta = not os.path.isfile(FILE_NAME_DELTA_EXTR % (m,Delta))
    else:
        missingDelta = not os.path.isfile(FILE_NAME_DELTA % (m,Delta))

    if missingDelta:
        # we should run the delta classification
        
        if (extremal):
            f = open(FILE_NAME_DELTA_EXTR % (m,Delta),'w')
            if (os.path.isfile(FILE_NAME_DELTA % (m,Delta))):
                g = open(FILE_NAME_DELTA % (m,Delta),'r')
                L = eval(g.read().replace('\n',' '))
                g.close()
                hdm = generalized_heller_constant(m,Delta,false)[0]
                result = []
                for P in L:
                    if (Polyhedron(P).integral_points_count() == hdm):
                        result.append(P)
                print([P for P in result],file=f)
                f.close()
            else:
                result = delta_classification(m,Delta,extremal)
                print([[tuple(p) for p in P.vertices()] for P in result],file=f)
                f.close()
        else:
            result = delta_classification(m,Delta,extremal)
            f = open(FILE_NAME_DELTA % (m,Delta),'w')
            print([[tuple(p) for p in P.vertices()] for P in result],file=f)
            f.close()
            
        
def lattice_polytopes_with_given_dimension_and_delta(m,Delta,extremal):
    """
        That's the main function for users of this module. It returns the list of all [extremal=false] or only h(Delta,m)-attaining [extremal=true]
        m-dimensional centrally symmetric lattice polytopes with delta equal to Delta.
    """
    # first, we update the database of lattice polytopes with a given delta
    update_delta_classification_database(m,Delta,extremal)

    # now, we can read the list of polytopes from the corresponding file and return them
    if (extremal):
        f = open(FILE_NAME_DELTA_EXTR % (m,Delta),'r')
    else:
        f = open(FILE_NAME_DELTA % (m,Delta),'r')
    
    L = eval(f.read().replace('\n',' '))
    f.close()
    return [Polyhedron(P) for P in L]


def generalized_heller_constant(m,Delta,extremal):
    """
        Compute the generalized Heller constant h(Delta,m) and a point set attaining it
    """
    
    DeltaPolytopes = lattice_polytopes_with_given_dimension_and_delta(m,Delta,extremal)
    nmax = 0
    for P in DeltaPolytopes:
    	npoints = P.integral_points_count()
    	if npoints > nmax:
    		nmax = npoints
    		Pmax = P
    return nmax , Pmax, len(DeltaPolytopes)
    

