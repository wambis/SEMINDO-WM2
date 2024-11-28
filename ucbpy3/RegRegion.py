from numpy.linalg import inv
import numpy as np

def rot_matrix(Clat=0, Clong=0):
    """
    Build rotation (cartesian) matrix. Input: center latitude and longitude (in deg)
    """
    alpha = 0.0

    Clat = 90.0 - Clat
    if Clong<0.0:
        Clong = 360.0 + Clong

    # Deg2rad
    Clat = np.deg2rad(Clat)
    Clong = np.deg2rad(Clong)
    alpha = np.deg2rad(alpha)

    rot = np.zeros([3,3])
    # Build matrix
    rot[0,0] = np.cos(alpha)*np.cos(Clat)*np.cos(Clong)-np.sin(alpha)*np.sin(Clong)
    rot[0,1] = -np.sin(alpha)*np.cos(Clat)*np.cos(Clong)-np.cos(alpha)*np.sin(Clong)
    rot[0,2] = np.sin(Clat)*np.cos(Clong)
    rot[1,0] = np.cos(alpha)*np.cos(Clat)*np.sin(Clong)+np.sin(alpha)*np.cos(Clong)
    rot[1,1] = -np.sin(alpha)*np.cos(Clat)*np.sin(Clong)+np.cos(alpha)*np.cos(Clong)
    rot[1,2] = np.sin(Clat)*np.sin(Clong)
    rot[2,0] = -np.cos(alpha)*np.sin(Clat)
    rot[2,1] = np.sin(alpha)*np.sin(Clat)
    rot[2,2] = np.cos(Clat)
    return rot

def sph2cart(lat, lon, R=6371):
    colat = 90.0 - lat
    if lon<0.0:
        lon = 360.0 + lon

    colat = np.deg2rad(colat)
    lon = np.deg2rad(lon)

    x = R * np.sin(colat) * np.cos(lon)
    y = R * np.sin(colat) * np.sin(lon)
    z = R * np.cos(colat)
    return np.array([x,y,z])

def cart2sph(x, y, z):
    R = np.sqrt(x**2 + y**2 + z**2)
    lon = np.rad2deg( np.arctan2( y,x ) )
    lat = 90 - np.rad2deg( np.arccos( z / R ) )

    return lat, lon

# Example
"""
# Center coordinates [in deg]
Clat = -20
Clon = -22

# Build rotation matrix
rot = rot_matrix(Clat=Clat, Clong=Clon)

lat = 20
lon = 20
x = np.zeros([3])
x[0] = np.tan( np.deg2rad( lon ) )
x[1] = np.tan( np.deg2rad( lat ) )
D = np.sqrt( 1 + x[0] **2 + x[1] ** 2 )
x[0] = x[0]/D
x[1] = x[1]/D
x[2] = 1./D

x_irot = np.dot(x,rot.T)
lat_c, lon_c = cart2sph(x=x_irot[0], y=x_irot[1], z=x_irot[2])
"""
