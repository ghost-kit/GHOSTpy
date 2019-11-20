"""
Coordinate System Transformations.  Supported corrdinate systems:

* Geographic
* Solar Ecliptic
* Solar Magnetic (SM)
* Geocentric solar magnetosphere (GSM)

* Transformations provided by cxform:
  * J2000
  * HEEQ
* Transformations provided by geopack:
  * Dipole (MAG)
  * Eqatorial inertial (GEI)
"""
import cxform

try:
    # geopack_2008 has been updated to
    #  * calculate in GSW (rather than GSM) coordinates
    #  * new coefficients to transform 1965 through 2015
    import geopack_08
except:
    try:
        import geopack
    except:
        transformer = 'CXFORM'
    else:
        transformer = 'GEOPACK'        
else:
    transformer = 'GEOPACK_08'




from transforms import *
