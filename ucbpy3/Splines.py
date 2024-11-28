
#
# Imports 
import numpy as np
import scipy.sparse as sparse

#
# Average inter-knot distances for spherical splines of the form disussed by
# Wang and Dahlen (1995) - here expressed in degrees
AVG_DIST_DEG = np.array([63.4, 31.7, 15.8, 7.9, 3.95, 1.98, 0.99])

#
# ---
#

class CubicBSplines:
    """
    A python class for manipulating and evaluating non-uniformly spaced cubic
    B-splines as described in Megnin and Romanowicz (2000)
    """

    def __init__(self, knots):
        # set the knot locations
        self.knots = knots
        
        # set the number of knots
        self.N = self.knots.size
        
    def evaluate(self, k, r, return_dense = False):
        # initialize the result array b
        bsp = np.zeros(r.shape)

        # establish bounds on summation over basis b-splines
        k_min = max(0,k-2)
        k_max = min(self.N-1,k+2)

        # get spline-knot spacing
        h = np.diff(self.knots[k_min:k_max+1])

        # select which radii in r will evaluate to nonzero b-spline values
        ir_eval = np.nonzero(np.logical_and(r >= self.knots[k_min],
                                            r <= self.knots[k_max]))

        # loop over radii ...
        for ir in ir_eval[0]:
            # find the index of the basis spline for this interval
            below, = np.nonzero(self.knots[k_min:k_max+1] < r[ir])
            if not below.any():
                index_basis = 0
            else:
                index_basis = below[-1]

            # here, we must treat specific splines differently ...
            if k == 0:
                # bottom-most spline
                denom = 3.0 * h[0] * h[0] + 3.0 * h[0] * h[1] + h[1] * h[1]
                if index_basis == 0:
                    a = 4.0 / (h[0] * (h[0] + h[1]) * denom)
                    b = 0.0
                    c = - 12.0 /denom
                    d = 4.0 * (2.0 * h[0] + h[1]) / denom
                elif index_basis == 1:
                    # WARNING a is different than MR00, but same as spl.c
                    a = - 4.0 / ( h[1] * (h[0] + h[1]) * denom)
                    b = 12.0 / ((h[0] + h[1]) * denom)
                    c = - 12.0 * h[1] / ((h[0] + h[1]) * denom)
                    d = 4.0 * h[1] * h[1] / ((h[0] + h[1]) * denom)
            elif k == 1:
                # second to bottom-most spline
                denom = (3.0 * h[0] * h[0] + 4.0 * h[0] * h[1] + h[1] * h[1] + 
                         2.0 * h[0] * h[2] + h[1] * h[2])
                denomsum = sum( h[:3] )
                dd = denomsum * denom
                if index_basis == 0:
                    a = - 4.0 * ((3.0 * h[0] + 2.0 * h[1] + h[2] ) /
                                 (h[0] * (h[0] + h[1]) * dd))
                    b = 0.0
                    c = 12.0 / denom
                    d = 0.0
                elif index_basis == 1:
                    # WARNING a is different than MR00, but same as spl.c
                    a = 4.0 * ((2.0 * h[0] * h[0] + 6.0 * h[0] * h[1] + 
                                3.0 * h[1] * h[1] + 3.0 * h[0] * h[2] +
                                3.0 * h[1] * h[2] + h[2] * h[2]) / 
                               (h[1] * (h[0] + h[1]) * (h[1] + h[2]) * dd))
                    b = - 12.0 * ((3.0 * h[0] + 2.0 * h[1] + h[2]) /
                                  ((h[0] + h[1]) * dd))
                    c = 12.0 * ((- 2.0 * h[0] * h[0] + h[1] * h[1] + h[1] * h[2]) /
                                ((h[0] + h[1]) * dd))
                    d = 4.0 * h[0] * ((4.0 * h[0] * h[1] + 3.0 * h[1] * h[1] +
                                       2.0 * h[0] * h[2] + 3.0 * h[1] * h[2]) /
                                      ((h[0] + h[1]) * dd))
                elif index_basis == 2:
                    dd = dd * (h[1] + h[2])
                    a = - 4.0 * (2.0 * h[0] + h[1]) / (h[2] * dd)
                    b = 12.0 * (2.0 * h[0] + h[1]) / dd
                    c = - 12.0 * (2.0 * h[0] + h[1])* h[2] / dd
                    d = 4.0 * (2.0 * h[0] + h[1]) * h[2] * h[2] / dd
            elif k == self.N - 2:
                # second to top-most spline
                denom = h[0] * h[1] + h[1] * h[1] + 2.0 * h[0] * h[2] + 4.0 * h[1] * h[2] + 3.0 * h[2] * h[2]
                denomsum = sum( h[:3] )
                dd = denomsum * denom
                if index_basis == 0:
                    a = 4.0 * (h[1] + 2.0 * h[2]) / (h[0] * (h[0] + h[1]) * dd)
                    b = 0.0
                    c = 0.0                
                    d = 0.0
                elif index_basis == 1:
                    a = - 4.0 * ((h[0] * h[0] + 3.0 * h[0] * h[1] + 
                                  3.0 * h[1] * h[1] + 3.0 * h[0] * h[2] +
                                  6.0 * h[1] * h[2] + 2.0 * h[2] * h[2]) / 
                                 (h[1] * (h[0] + h[1]) * (h[1] + h[2]) * dd))
                    b = 12.0 * (h[1] + 2.0 * h[2]) / ((h[0] + h[1]) * dd)
                    c = h[0] * b
                    d = h[0] * c / 3.0
                elif index_basis == 2:
                    dd = dd * (h[1] + h[2])
                    a = 4.0 * (h[0] + 2.0 * h[1] + 3.0 * h[2]) / (h[2] * dd)
                    b = - 12.0 * (h[0] + 2.0 * h[1] + 3.0 * h[2]) / dd
                    c = 12.0 * (- h[0] * h[1] - h[1] * h[1] + 2.0 * h[2] * h[2]) / dd
                    d = 4.0 * h[2] * (3.0 * h[0] * h[1] + 3.0 * h[1] * h[1] +
                                      2.0 * h[0] * h[2] + 4.0 * h[1] * h[2]) / dd
            elif k == self.N - 1:
                # top-most spline
                # WARNING denom is different from MR00, but same as spl.c
                denom = (h[0] + h[1]) * (h[0] * h[0] + 3.0 * h[0] * h[1] + 3.0 * h[1] * h[1])
                if index_basis == 0:
                    a = 4.0 / (h[0] * denom)
                    b = 0.0
                    c = 0.0                
                    d = 0.0
                elif index_basis == 1:
                    a = - 4.0 / (h[1] * denom)
                    b = 12.0 / denom
                    c = h[0] * b
                    d = h[0] * c / 3.0
            else:
                # fully internal spline
                denom1 = h[0] + h[1] + h[2] + h[3]
                if index_basis == 0:
                    a = 4.0 / (h[0] * (h[0] + h[1]) * (h[0] + h[1] + h[2]) * denom1)
                    b = 0.0
                    c = 0.0
                    d = 0.0
                elif index_basis == 1:
                    denom2 = (h[0] + h[1]) * (h[0] + h[1] + h[2])
                    denom = denom1 * denom2
                    a = - 4.0 * ((h[0] * h[0] + 3.0 * h[0] * h[1] +
                                  3.0 * h[1] * h[1] + 2.0 * h[0] * h[2] +
                                  4.0 * h[1] * h[2] + h[2] * h[2] + h[0] * h[3] +
                                  2.0 * h[1] * h[3] + h[2] * h[3]) /
                                 (h[1] * (h[1] + h[2]) * (h[1] + h[2] + h[3]) * denom))
                    b = 12.0 / denom
                    c = h[0] * b
                    d = h[0] * c / 3.0
                elif index_basis == 2:
                    denom2 = (h[1] + h[2]) * (h[0] + h[1] + h[2]) * (h[1] + h[2] + h[3])
                    denom = denom1 * denom2
                    a = 4.0 * ((h[0] * h[1] + h[1] * h[1] +
                                2.0 * h[0] * h[2] + 4.0 * h[1] * h[2] + 
                                3.0 * h[2] * h[2] + h[0] * h[3] +
                                2.0 * h[1] * h[3] + 3.0 * h[2] * h[3] +
                                h[3] * h[3]) /
                               (h[2] * (h[2] + h[3]) * denom))
                    b = - 12.0 * (h[0] + 2.0 * h[1] + 2.0 * h[2] + h[3]) / denom
                    c = 12.0 * (- h[0] * h[1] - h[1] * h[1] +
                                h[2] * h[2] + h[2] * h[3]) / denom
                    d = 4.0 * (2.0 * h[0] * h[1] * h[2] + 2.0 * h[1] * h[1] * h[2] +
                               h[0] * h[2] * h[2] + 2.0 * h[1] * h[2] * h[2] +
                               h[0] * h[1] * h[3] + h[1] * h[1] * h[3] +
                               h[0] * h[2] * h[3] + 2.0 * h[1] * h[2] * h[3]) / denom
                elif index_basis == 3:
                    denom2 = (h[2] + h[3]) * (h[1] + h[2] + h[3])
                    denom = denom1 * denom2
                    a = - 4.0 / (h[3] * denom)
                    b = 12.0 / denom
                    c = - h[3] * b
                    d = - h[3] * c / 3.0
            # update the result array
            bsp[ir] = (a * (r[ir] - self.knots[k_min+index_basis]) ** 3 +
                       b * (r[ir] - self.knots[k_min+index_basis]) ** 2 +
                       c * (r[ir] - self.knots[k_min+index_basis]) + d)
        # return the complete result array
        if return_dense:
            return bsp
        else:
            return sparse.csr_matrix(bsp)


class SphericalSplines:
    """
    A python class for manipulating and evaluating spherical splines as
    defined in Wang and Dahlen (1995)
    """
    
    def __init__(self, knots, degrees = True):
        # set spherical spline knot locations, converting to radians
        if degrees:
            self.knots = np.deg2rad(knots[:,0:2])
        else:
            self.knots = knots[:,0:2]

        # set knot levels as integers (for indexing)
        self.level = np.array(knots[:,2], dtype = int)

        # set number of knots
        self.N = self.knots.shape[0]

    def __haversin(self, theta):
        # return the haversine of theta
        return 0.5 * (1.0 - np.cos(theta))
        
    def __dist_haversin(self, points_rad):
        # initialize angular distance delta to zero
        delta = np.zeros((points_rad.shape[0],self.N))

        # loop over spherical spline knots
        for k in range(self.N):
            # compute angular distance using the haversine formula
            delta[:,k] = 2.0 * np.arcsin(
                np.sqrt(self.__haversin(points_rad[:,1] - self.knots[k,1]) +
                        np.cos(points_rad[:,1]) * np.cos(self.knots[k,1]) *
                        self.__haversin(points_rad[:,0] - self.knots[k,0])))
            
        # return the resulting distances
        return delta
        
    def __dist_cos(self, points_rad):
        # initialize angular distance delta to zero
        delta = np.zeros((points_rad.shape[0],self.N))

        # loop over spherical spline knots
        for k in range(self.N):
            # compute angular distance using the spherical law of cosines
            arg = (np.cos(0.5 * np.pi - self.knots[k,1]) * np.cos(0.5 * np.pi - points_rad[:,1]) +
                   np.sin(0.5 * np.pi - self.knots[k,1]) * np.sin(0.5 * np.pi - points_rad[:,1]) *
                   np.cos(points_rad[:,0] - self.knots[k,0]))
            # sanitize
            arg[arg < -1.0] = -1.0
            arg[arg > +1.0] = +1.0
            delta[:,k] = np.arccos(arg)
            
        # return the resulting distances
        return delta

    def evaluate(self, points, degrees = True, return_dense = False, verbose = False,
                 use_fixed_avg_dist = True,
                 use_ucb_expression = True,
                 use_spherical_cosines = True):
        """
        Evaluates the spherical spline knot coefficients at the lon / lat points
        specified in points - can optionally return *sparse* matrix in csr format
        """
        # set sampling locations, converting to radians
        if degrees:
            points_rad = np.deg2rad(points)
        else:
            points_rad = points

        # calculate angular distance from each location to all knots
        if use_spherical_cosines:
            if verbose:
                print(" evaluate: computing distance with the spherical law of cosines")
            delta = self.__dist_cos(points_rad)
        else:
            if verbose:
                print(" evaluate: computing distance with the haversine formula")
            delta = self.__dist_haversin(points_rad)
            
        # output progress
        if verbose:
            print(" evaluate: computed distances")
        
        # extract average inter-neighbor knot distance for each knot, given its level
        if use_fixed_avg_dist:
            avg_dist = np.deg2rad(AVG_DIST_DEG)[self.level - 1]
        else:
            delta_knots = self.__dist(self.knots)
            avg_dist = np.array([ np.mean(np.sort(delta_knots[:,k])[1:6]) for k in range(self.N) ])

        # initialize result array as either zeroed dense or sparse matrix in list-of-lists format
        if return_dense:
            sspl = np.zeros(delta.shape)
        else:
            sspl = sparse.lil_matrix(delta.shape)
        
        # loop over locations 
        for p in range(points.shape[0]):
            # locate and loop over points within one average inter-neighbor distance
            index_delta_1, = np.nonzero(delta[p,:] < avg_dist)
            for d in index_delta_1:
                sspl[p,d] = (  0.75 * ( delta[p,d] / avg_dist[d] ) ** 3
                             - 1.5  * ( delta[p,d] / avg_dist[d] ) ** 2
                             + 1.0)
                
            # locate and loop over points at between one and two averege inter-knot distances
            index_delta_2, = np.nonzero(np.logical_and(delta[p,:] >= avg_dist,
                                                       delta[p,:] <= 2.0 * avg_dist))
            if use_ucb_expression:
                for d in index_delta_2:
                    sspl[p,d] = 0.25 * ( 2.0 - delta[p,d] / avg_dist[d] ) ** 3
            else:
                for d in index_delta_2:
                    sspl[p,d] = (- 0.25 * ( delta[p,d] / avg_dist[d] - 1.0 ) ** 3
                                 + 0.75 * ( delta[p,d] / avg_dist[d] - 1.0 ) ** 2
                                 - 0.75 * ( delta[p,d] / avg_dist[d] - 1.0 )
                                 + 0.25)

            # output progress
            if verbose and p % 10000 == 0:
                print(" evaluate: point %8i" % p)
        
        # return the sampling (possibly as dense matrix)
        if return_dense:
            return sspl
        else:
            return sspl.tocsr()


            
    
