#
# Imports
import numpy as np

def geocentric(thetap, degrees=True):
    """
     Calculates geocentric colatitude from geodetic.
     Uses expressions similar to those in:
         /data/25/yuan/src/libsrc/sphere/geocentric.c
     (also similar to that in Dahlen and Tromp (1998) p. 603)

     Arguments:
        thetap : float, possibly array-like
            geodetic colatitude (units depend on `degrees`)
        degrees : bool (optional, default: True)
            assume input is in degrees, return output in degrees
            otherwise, radians
    """
    thetap = np.array(thetap)
    if degrees:
        thetap = np.deg2rad(thetap)
    if np.any(np.sin(thetap) < 0):
        raise ValueError('Supplied colatitude out of acceptable range!')
    FAC = 0.993305621334896
    HPI = 1.570796326794895
    arg = HPI - np.arctan2(FAC * np.cos(thetap), np.sin(thetap))
    if arg.shape == ():
        if arg < 0.0:
            arg += 2.0 * np.pi
        elif arg > np.pi:
            arg -= 2.0 * np.pi
        arg = float(arg)
    else:
        neg = arg < 0
        arg[neg] += 2.0 * np.pi
        arg[np.logical_and(arg > np.pi, np.logical_not(neg))] -= 2.0 * np.pi
    if degrees:
        theta = np.rad2deg(arg)
    else:
        theta = arg
    return theta

def geodetic(theta, degrees=True):
    """
     Calculates geodetic colatitude from geocentric.
     Uses expressions similar to those in:
         /data/25/yuan/src/libsrc/sphere/geocentric.c
     (also similar to that in Dahlen and Tromp (1998) p. 603)

     Arguments:
        theta : float, possibly array-like
            geocentric colatitude (units depend on `degrees`)
        degrees : bool (optional, default: True)
            assume input is in degrees, return output in degrees
    """
    theta = np.array(theta)
    if degrees:
        theta = np.deg2rad(theta)
    if np.any(np.sin(theta) < 0):
        raise ValueError('Supplied colatitude out of acceptable range!')
    FAC = 0.993305621334896
    HPI = 1.570796326794895
    arg = np.arctan2(FAC * np.cos(HPI - theta), np.sin(HPI - theta))
    if arg.shape == ():
        if arg < 0.0:
            arg += 2.0 * np.pi
        elif arg > np.pi:
            arg -= 2.0 * np.pi
        arg = float(arg)
    else:
        neg = arg < 0
        arg[neg] += 2.0 * np.pi
        arg[np.logical_and(arg > np.pi, np.logical_not(neg))] -= 2.0 * np.pi
    if degrees:
        thetap = np.rad2deg(arg)
    else:
        thetap = arg
    return thetap

def delaz(src_loc, sta_loc, degrees=True, delta_only=False):
    """
     Returns source-station epicentral distance and azimuth, the latter 
     measured clockwise from North. Input arguments and return values are 
     assumed in _degrees_ by default, and presumed already to have been
     converted to geocentric coordinates if needed.

     Arguments:
        src_loc : tuple, (float, float)
            location of the source
        sta_loc : tuple, (float, float)
            location of the station
        degrees : bool (optional, default: True)
            assume input is in degrees, return output in degrees
        delta_only : bool (optional, default: False)
            calculate and return only delta (no azimuth)
    """
    # extract source and station locations
    sta_lon, sta_lat = sta_loc
    src_lon, src_lat = src_loc
    # convert to radians
    if degrees:
        sta_lon = np.radians(sta_lon)
        sta_lat = np.radians(sta_lat)
        src_lon = np.radians(src_lon)
        src_lat = np.radians(src_lat)
    # convert to colat
    sta_colat = 0.5 * np.pi - sta_lat
    src_colat = 0.5 * np.pi - src_lat
    # calculate distance
    delta_arg = (np.cos(src_colat) * np.cos(sta_colat) +
                 np.sin(src_colat) * np.sin(sta_colat) * np.cos(sta_lon - src_lon))
    delta_arg = max(-1., min(1., delta_arg))
    delta = np.arccos(delta_arg)
    # calculate azimuth - *not* necessarily CW
    if not delta_only:
        az_arg = ((np.cos(sta_colat) - np.cos(delta) * np.cos(src_colat)) / 
                  (np.sin(delta) * np.sin(src_colat)))
        az_arg = max(-1., min(1., az_arg))
        az = np.arccos(az_arg)
        # correct az if actually CCW
        if np.sin(sta_lon - src_lon) <= 0:
            az = 2 * np.pi - az
    # convert to degrees
    if degrees:
        delta = np.degrees(delta)
        if not delta_only:
            az = np.degrees(az)
    # return
    if delta_only:
        return delta
    else:
        return delta, az

def gc_minor(src_loc, path_spec, dx, path_spec_type='point', degrees=True):
    """
     Returns source-station minor-arc great circle path in increments of dx.
     Input arguments are assumed in _degrees_ by default, and presumed already
     to have been converted to geocentric coordinates if needed.

     Arguments:
        src_loc : tuple, (float, float)
            location of the source
        path_spec : tuple, (float, float)
            specifies path along which to calculate path
            either (sta_lon, sta_lat) or (delta, az)
        dx : float
            path step size
        path_spec_type : str (optional, default: 'point')
            the type of path_spec, either 'point' for a station location or
            'delaz' for a (delta, az) tuple
        degrees : bool (optional, default: True)
            assume input is in degrees, return output in degrees
    """
    # extract source location
    src_lon, src_lat = src_loc
    # convert to radians
    if degrees:
        src_lon = np.radians(src_lon)
        src_lat = np.radians(src_lat)
        dx = np.radians(dx)
    # convert to colat, normalize lon to [0, 360)
    src_colat = 0.5 * np.pi - src_lat
    # parse path specification
    if path_spec_type == 'point':
        sta_lon, sta_lat = path_spec
        # convert to radians
        if degrees:
            sta_lon = np.radians(sta_lon)
            sta_lat = np.radians(sta_lat)
        delta, az = delaz((src_lon, src_lat), (sta_lon, sta_lat), degrees=False)
    elif path_spec_type == 'delaz':
        delta, az = path_spec
        if degrees:
            delta = np.radians(delta)
            az = np.radians(az)
        if az < 0:
            az += 2.0 * np.pi
    else:
        raise ValueError('Unrecognized path_spec type "%s"' % (path_spec_type))
    # flip looking direction to CCW below if necessary
    ccw_fac = 1
    if az > np.pi:
        az = 2 * np.pi - az
        ccw_fac = - 1
    # discretization along the gc path
    deltas = np.arange(0, delta + 1e-5 * dx, dx)
    # recover path colats
    colats_arg = (np.cos(deltas) * np.cos(src_colat) + 
                  np.sin(deltas) * np.sin(src_colat) * np.cos(az))
    colats_arg[colats_arg < -1.] = -1.
    colats_arg[colats_arg >  1.] =  1.
    colats = np.arccos(colats_arg)
    # recover lons
    lons_arg = ((np.cos(deltas) - np.cos(src_colat) * np.cos(colats)) /
                (np.sin(src_colat) * np.sin(colats)))
    lons_arg[lons_arg < -1.] = -1.
    lons_arg[lons_arg >  1.] =  1.
    lons = np.degrees(src_lon + ccw_fac * np.arccos(lons_arg))
    lons = (lons + 360.0) % 360.0
    # convert to desired units
    if degrees:
        lats = 90.0 - np.degrees(colats)
    else:
        lats = 0.5 * np.pi - colats
        lons = np.radians(lons)
    # return
    return lons, lats

def shoot(src_loc, delta, az, degrees=True):
    """
     Returns endpoint of a great-circle path from an anchor location along a
     particular delta / azimuth path. Input arguments are assumed in _degrees_
     by default, and presumed already to have been converted to geocentric
     coordinates if needed.

     Arguments:
        src_loc : tuple, (float, float)
            location of the source
        delta : float
            path distance
        az : float
            path azimuth
        degrees : bool (optional, default: True)
            assume input is in degrees, return output in degrees
    """
    if np.abs(delta) > 0:
        lon, lat = gc_minor(src_loc, (delta, az), delta,
            path_spec_type='delaz', degrees=degrees)
        rval = lon[-1], lat[-1]
    else:
        rval = src_loc
    return rval
