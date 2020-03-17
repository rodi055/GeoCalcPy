import numpy as np
import wgs84
from utils import input_check_Nx3 as _input_check_Nx3
from utils import input_check_Nx3x3 as _input_check_Nx3x3
from utils import input_check_Nx1 as _input_check_Nx1

def earthrad(lat, lat_unit='deg', model='wgs84'):
    """
    Calculate radius of curvature in the prime vertical (East-West) and 
    meridian (North-South) at a given latitude.
    Parameters
    ----------
    lat : {(N,)} array like latitude, unit specified by lat_unit, default in deg
    
    Returns
    -------
    R_N : {(N,)} array like, radius of curvature in the prime vertical (East-West)
    R_M : {(N,)} array like, radius of curvature in the meridian (North-South)
    
    Examples
    --------
    >>> import numpy as np
    >>> from navpy import earthrad
    >>> lat = 0
    >>> Rtransverse, Rmeridian = earthrad(lat)
    >>> Rtransverse
    6378137.0
    >>> Rmeridian
    6335439.3272928288
    >>> lat = [0, np.pi/2]
    >>> Rtransverse, Rmeridian = earthrad(lat,lat_unit='rad')
    >>> Rtransverse
    array([ 6378137.        ,  6399593.62575849])
    >>> Rmeridian
    array([ 6335439.32729283,  6399593.62575849])
    """
    if(lat_unit=='deg'):
        lat = np.deg2rad(lat)
    elif(lat_unit=='rad'):
        pass
    else:
        raise ValueError('Input unit unknown')

    if(model=='wgs84'):
        R_N = wgs84.a/(1-wgs84._ecc_sqrd*np.sin(lat)**2)**0.5
        R_M = wgs84.a*(1-wgs84._ecc_sqrd)/(1-wgs84._ecc_sqrd*np.sin(lat)**2)**1.5
    else:
        raise ValueError('Model unknown')
    
    return R_N, R_M

def lla2ecef(lat, lon, alt, latlon_unit='deg', alt_unit='m', model='wgs84'):
    """
    Convert Latitude, Longitude, Altitude, to ECEF position
    
    Parameters
    ----------
    lat : {(N,)} array like latitude, unit specified by latlon_unit, default in deg
    lon : {(N,)} array like longitude, unit specified by latlon_unit, default in deg
    alt : {(N,)} array like altitude, unit specified by alt_unit, default in m
    
    Returns
    -------
    ecef : {(N,3)} array like ecef position, unit is the same as alt_unit
    """
    lat,N1 = _input_check_Nx1(lat)
    lon,N2 = _input_check_Nx1(lon)
    alt,N3 = _input_check_Nx1(alt)
    
    if( (N1!=N2) or (N2!=N3) or (N1!=N3) ):
        raise ValueError('Inputs are not of the same dimension')
    
    if(model=='wgs84'):
        Rew,Rns = earthrad(lat,lat_unit=latlon_unit)
    else:
        Rew = wgs84.a 
    
    if(latlon_unit=='deg'):
        lat = np.deg2rad(lat)
        lon = np.deg2rad(lon)
    
    x = (Rew + alt)*np.cos(lat)*np.cos(lon)
    y = (Rew + alt)*np.cos(lat)*np.sin(lon)
    z = ( (1-wgs84._ecc_sqrd)*Rew + alt )*np.sin(lat)
    
    ecef = np.vstack((x,y,z)).T

    if(N1==1):
        ecef = ecef.reshape(3)

    return ecef

def ecef2lla(ecef, latlon_unit='deg'):
    """
    Calculate the Latitude, Longitude and Altitude of a point located on earth 
    given the ECEF Coordinates.
    
    References
    ----------
    .. [1] Jekeli, C.,"Inertial Navigation Systems With Geodetic
       Applications", Walter de Gruyter, New York, 2001, pp. 24
    
    Parameters
    ----------
    ecef : {(N,3)} array like input of ECEF coordinate in X, Y, and Z column, unit is meters
    latlon_unit : {('deg','rad')} specifies the output latitude and longitude unit
    
    Returns
    -------
    lat : {(N,)} array like latitude in unit specified by latlon_unit
    lon : {(N,)} array like longitude in unit specified by latlon_unit
    alt : {(N,)} array like altitude in meters
    """
    ecef,N = _input_check_Nx3(ecef)
    if(N>1):
        x = ecef[:,0]; y = ecef[:,1]; z = ecef[:,2]
    else:
        x = ecef[0]; y = ecef[1]; z = ecef[2]

    lon = np.arctan2(y,x)

    # Iteration to get Latitude and Altitude
    p = np.sqrt(x**2 + y**2)
    lat = np.arctan2(z,p*(1-wgs84._ecc_sqrd))

    err = np.ones(N)

    while(np.max(np.abs(err))>1e-10):
        Rew,Rns = earthrad(lat,lat_unit='rad')
        h = p/np.cos(lat) - Rew
        
        err = np.arctan2(z*(1+wgs84._ecc_sqrd*Rew*np.sin(lat)/z),p) - lat
        
        lat = lat + err
    
    if(latlon_unit=='deg'):
        lat = np.rad2deg(lat)
        lon = np.rad2deg(lon)

    return lat, lon, h


def lla2ned(lat, lon, alt, lat_ref, lon_ref, alt_ref, latlon_unit='deg', alt_unit='m', model='wgs84'):
    """
    Convert Latitude, Longitude, Altitude to its resolution in the NED
    coordinate. The center of the NED coordiante is given by lat_ref, lon_ref,
    and alt_ref.
    
    For example, this can be used to convert GPS data to a local NED frame.
    
    Parameters
    ----------
    lat : {(N,)} array like latitude, unit specified by latlon_unit, default in deg
    lon : {(N,)} array like longitude, unit specified by latlon_unit, default in deg
    alt : {(N,)} array like altitude, unit specified by alt_unit, default in m
    
    lat_ref : Reference latitude, unit specified by latlon_unit, default in deg
    lon_ref : Reference longitude, unit specified by latlon_unit, default in deg
    alt : Reference altitude, unit specified by alt_unit, default in m
    
    Returns
    -------
    ned : {(N,3)} array like ecef position, unit is the same as alt_unit        
    """
    ecef  = lla2ecef(lat, lon, alt, latlon_unit=latlon_unit, 
                           alt_unit=alt_unit, model=model)
    ecef0 = lla2ecef(lat_ref, lon_ref, alt_ref,
                           latlon_unit=latlon_unit, 
                           alt_unit=alt_unit, model=model)
    ned  = ecef2ned(ecef-ecef0, lat_ref, lon_ref, alt_ref, 
                          latlon_unit=latlon_unit, alt_unit=alt_unit, model=model)
    return ned

def ned2lla(ned, lat_ref, lon_ref, alt_ref, latlon_unit='deg', alt_unit='m', model='wgs84'):
    """
    Calculate the Latitude, Longitude and Altitude of points given by NED coordinates
    where NED origin given by lat_ref, lon_ref, and alt_ref.
    Parameters
    ----------
    ned : {(N,3)} array like input of NED coordinate in N, E, and D column, unit is meters
    lat_ref : Reference latitude, unit specified by latlon_unit, default in deg
    lon_ref : Reference longitude, unit specified by latlon_unit, default in deg
    alt_ref : Reference altitude, unit specified by alt_unit, default in m
    latlon_unit : {('deg','rad')} specifies the output latitude and longitude unit
    
    Returns
    -------
    lat : {(N,)} array like latitude in unit specified by latlon_unit
    lon : {(N,)} array like longitude in unit specified by latlon_unit
    alt : {(N,)} array like altitude in meters
    Note
    ----
    This method is a wrapper on ned2ecef (add ecef of NED-origin) and ecef2lla.
    """

    ecef = ned2ecef(ned,lat_ref,lon_ref,alt_ref,latlon_unit=latlon_unit,
                                                alt_unit=alt_unit,
                                                model=model)
    # Add vector to ecef representation of NED-origin
    ecef_ref = lla2ecef(lat_ref, lon_ref, alt_ref, latlon_unit=latlon_unit,
                                                   alt_unit=alt_unit,
                                                   model=model)
    ecef += ecef_ref

    lla = ecef2lla(ecef, latlon_unit=latlon_unit)

    return lla


def ned2ecef(ned,lat_ref,lon_ref,alt_ref,latlon_unit='deg',alt_unit='m',model='wgs84'):
    """
    Transform a vector resolved in NED (origin given by lat_ref, lon_ref, and alt_ref)
    coordinates to its ECEF representation. 
    Parameters
    ----------
    ned : {(N,3)} input array, units of meters
    lat_ref : Reference latitude, unit specified by latlon_unit, default in deg
    lon_ref : Reference longitude, unit specified by latlon_unit, default in deg
    alt_ref : Reference altitude, unit specified by alt_unit, default in m
    
    Returns
    -------
    ecef : {(N,3)} array like ned vector, in the ECEF frame, units of meters
    Notes
    -----
    The NED vector is treated as a relative vector, and hence the ECEF representation
    returned is NOT converted into an absolute coordinate.  This means that the 
    magnitude of `ned` and `ecef` will be the same (bar numerical differences).
    
    Examples
    --------
    >>> import navpy
    >>> ned = [0, 0, 1]
    >>> lat_ref, lon_ref, alt_ref = 45.0, -93.0, 250.0 # deg, meters
    >>> ecef = navpy.ned2ecef(ned, lat_ref, lon_ref, alt_ref)
    >>> print("NED:", ned)
    >>> print("ECEF:", ecef)
    >>> print("Notice that 'down' is not same as 'ecef-z' coordinate.")
    """
    lat_ref,N1 = _input_check_Nx1(lat_ref)
    lon_ref,N2 = _input_check_Nx1(lon_ref)
    alt_ref,N3 = _input_check_Nx1(alt_ref)
    
    if( (N1!=1) or (N2!=1) or (N3!=1) ):
        raise ValueError('Reference Location can only be 1')
    
    ned,N = _input_check_Nx3(ned)

    ned = ned.T
    
    C = np.zeros((3,3))

    if(latlon_unit=='deg'):
        lat_ref = np.deg2rad(lat_ref)
        lon_ref = np.deg2rad(lon_ref)
    elif(latlon_unit=='rad'):
        pass
    else:
        raise ValueError('Input unit unknown')

    C[0,0]=-np.sin(lat_ref)*np.cos(lon_ref)
    C[0,1]=-np.sin(lat_ref)*np.sin(lon_ref)
    C[0,2]= np.cos(lat_ref)
    
    C[1,0]=-np.sin(lon_ref)
    C[1,1]= np.cos(lon_ref)
    C[1,2]= 0

    C[2,0]=-np.cos(lat_ref)*np.cos(lon_ref)
    C[2,1]=-np.cos(lat_ref)*np.sin(lon_ref)
    C[2,2]=-np.sin(lat_ref)

    # C defines transoformation: ned = C * ecef.  Hence used transpose.
    ecef = np.dot(C.T,ned)
    ecef = ecef.T

    if(N == 1):
        ecef = ecef.reshape(3)

    return ecef

def ecef2ned(ecef,lat_ref,lon_ref,alt_ref,latlon_unit='deg',alt_unit='m',model='wgs84'):
    """
    Transform a vector resolved in ECEF coordinate to its resolution in the NED
    coordinate. The center of the NED coordiante is given by lat_ref, lon_ref,
    and alt_ref.
    
    Parameters
    ----------
    ecef : {(N,3)} input vector expressed in the ECEF frame
    lat_ref : Reference latitude, unit specified by latlon_unit, default in deg
    lon_ref : Reference longitude, unit specified by latlon_unit, default in deg
    alt : Reference altitude, unit specified by alt_unit, default in m
    
    Returns
    -------
    ned : {(N,3)} array like ecef position, unit is the same as alt_unit
    
    Examples
    --------
    >>> import numpy as np
    >>> from navpy import ecef2ned
    >>> lat 
    """
    lat_ref,N1 = _input_check_Nx1(lat_ref)
    lon_ref,N2 = _input_check_Nx1(lon_ref)
    alt_ref,N3 = _input_check_Nx1(alt_ref)
    
    if( (N1!=1) or (N2!=1) or (N3!=1) ):
        raise ValueError('Reference Location can only be 1')
    
    ecef,N = _input_check_Nx3(ecef)

    ecef = ecef.T
    
    C = np.zeros((3,3))

    if(latlon_unit=='deg'):
        lat_ref = np.deg2rad(lat_ref)
        lon_ref = np.deg2rad(lon_ref)
    elif(latlon_unit=='rad'):
        pass
    else:
        raise ValueError('Input unit unknown')

    C[0,0]=-np.sin(lat_ref)*np.cos(lon_ref)
    C[0,1]=-np.sin(lat_ref)*np.sin(lon_ref)
    C[0,2]= np.cos(lat_ref)
	
    C[1,0]=-np.sin(lon_ref)
    C[1,1]= np.cos(lon_ref)
    C[1,2]= 0

    C[2,0]=-np.cos(lat_ref)*np.cos(lon_ref)
    C[2,1]=-np.cos(lat_ref)*np.sin(lon_ref)
    C[2,2]=-np.sin(lat_ref)

    ned = np.dot(C,ecef)
    ned = ned.T

    if(N == 1):
        ned = ned.reshape(3)

    return ned
