"""
This module provides coordinate transformations relevant to geospace modeling.
Many/most of these are wrappers to external C or Fortran libraries, employing
the NumPy vectorize() method to handle ndarray inputs when appropriate.

(x,y,z,
 dx,dy,dz) = SPHtoCAR(phi,theta,rho,dphi,dtheta,drho)
                                  - convert position and directional vectors
                                    from a local spherical to Cartesian grid;
                                    set dphi,dtheta,drho to None, or do not set
                                    to simply obtain position vectors
(phi,theta,rho
 dphi,dtheta,drho) = CARtoSPH(x,y,z,dx,dy,dz)
                                  - convert position and directional vectors 
                                    from a local Cartesian to spherical grid;
                                    set dx,dy,dz to None, or do not set to
                                    simply obtain position vectors
x,y,z = SMtoGSM(x,y,z,dateTime)   - convert from solar magnetic to geocentric
                                    solar magnetospheric coordinates
x,y,z = GSMtoSM(x,y,z,dateTime)   - convert from geocentric solar magnetospheric
                                    to solar magnetic coordinates
x,y,z = GSEtoGSM(x,y,z,dateTime)  - convert from geocentric solar ecliptic to
                                    geocentric magnetospheric coordinates
x,y,z = HEEQtoGSM(x,y,z,dateTime) - convert from heliocentric earth equatorial
                                    to geocentric magnetospheric coordinates
x,y,z = GEOtoMAG(x,y,z,dateTime)  - convert from geographic to dipole (mag) coordinates
x,y,z = MAGtoGEO(x,y,z,dateTime)  - convert from dipole (mag) to geographic coordinates
x,y,z = SMtoMAG(x,y,z,dateTime)   - convert from SM to dipole (mag) coordinates
x,y,z = MAGtoSM(x,y,z,dateTime)   - convert from dipole (mag) to SM coordinates
x,y,z = GEOtoGSM(x,y,z,dateTime)  - convert from geographic to GSM coordinates
x,y,z = GSMtoGEO(x,y,z,dateTime)  - convert from GSM to geographic coordinates
x,y,z = GEOtoSM(x,y,z,dateTime)   - convert from geographic to SM coordinates
x,y,z = SMtoGEO(x,y,z,dateTime)   - convert from SM to geographic coordinates
"""

import ghostpy.transform
import pylab as p
import datetime

def SPHtoCAR(phi,theta,rho,dphi=None,dtheta=None,drho=None):
   """
   convert spherical coordinates (phi,theta,rho) to cartesian (x,y,z)
     - if 3 inputs, convert position vectors
     - if 6 inputs, convert directional vectors
   """
      
   if dphi==None or dtheta==None or drho==None:
      
      if dphi==None and dtheta==None and drho==None:
         
         ## transform position vectors
         
         iSPHtoCAR = 1
         
         ## vectorize geopack_08 call using NumPy's vectorize() method:
         ## Note:  vectorize()'s doc string states that this is not very speed
         ##        optimized...Geopack should be ported to Python, or re-written
         ##        to accept array inputs, if speed is an issue. -EJR 10/2013
         ##
         #output = ghostpy.transform.geopack_08.sphcar_08(rho,theta,phi,0.,0.,0., iSPHtoCAR)
         v_gp08 = p.vectorize(ghostpy.transform.geopack_08.sphcar_08)
         output = v_gp08(rho,theta,phi, 0.,0.,0., iSPHtoCAR)
         
         ## if scalar inputs, generate scalar outputs -EJR 10/2013
         return (tuple(p.asscalar(otmp) for otmp in output[3:6]) if
                 (p.isscalar(phi) and p.isscalar(theta) and p.isscalar(rho)) else
                 (output[3],output[4],output[5]) )
      else:
         raise Exception, 'SPHtoCAR requires exactly 3 or 6 inputs'
         
   else:
      
      iSPHtoCAR = 1
      
      ## transform local directional vectors
         
      ## vectorize geopack_08 call using NumPy's vectorize() method:
      ## Note:  vectorize()'s doc string states that this is not very speed
      ##        optimized...Geopack should be ported to Python, or re-written
      ##        to accept array inputs, if speed is an issue. -EJR 10/2013
      ##
      
      #output = ghostpy.transform.geopack_08.sphcar_08(rho,theta,phi,0.,0.,0., iSPHtoCAR)
      v_gp08a = p.vectorize(ghostpy.transform.geopack_08.sphcar_08)
      _,_,_,x,y,z,_ = v_gp08a(rho,theta,phi, 0.,0.,0., iSPHtoCAR)
      
      #dx,dy,dz = ghostpy.transform.geopack_08.bspcar_08(theta,phi,drho,dtheta,dphi)
      v_gp08b = p.vectorize(ghostpy.transform.geopack_08.bspcar_08)
      dx,dy,dz = v_gp08b(theta,phi, drho,dtheta,dphi)
      
      ## if all scalar inputs, return scalar outputs, otherwise ndarrays -EJR 10/2013
      return (tuple(p.asscalar(otmp) for otmp in (x,y,z,dx,dy,dz)) if 
              (p.isscalar(phi) and p.isscalar(theta) and p.isscalar(rho) and
               p.isscalar(dphi) and p.isscalar(dtheta) and p.isscalar(drho)) else
              (x,y,z,dx,dy,dz) )



def CARtoSPH(x,y,z,dx=None,dy=None,dz=None):
   """
   convert cartesian coordinates (x,y,z) to spherical (phi,theta,rho)
     - if 3 inputs, convert position vectors
     - if 6 inputs, convert directional vectors
   """
   
   
   if dx==None or dy==None or dz==None:
      
      if dx==None and dy==None and dz==None:
         
         ## transform position vectors
         
         iSPHtoCAR = -1
         
         ## vectorize geopack_08 call using NumPy's vectorize() method:
         ## Note:  vectorize()'s doc string states that this is not very speed
         ##        optimized...Geopack should be ported to Python, or re-written
         ##        to accept array inputs, if speed is an issue. -EJR 10/2013
         ##
         #output = ghostpy.transform.geopack_08.sphcar_08(0.,0.,0.,x,y,z, iSPHtoCAR)
         v_gp08 = p.vectorize(ghostpy.transform.geopack_08.sphcar_08)
         output = v_gp08(0.,0.,0.,x,y,z, iSPHtoCAR)
         
         ## if scalar inputs, generate scalar outputs -EJR 10/2013
         return (tuple(p.asscalar(otmp) for otmp in output[2::-1]) if
                 (p.isscalar(x) and p.isscalar(y) and p.isscalar(z)) else
                 (output[2],output[1],output[0]) )
      else:
         raise Exception, 'SPHtoCAR requires exactly 3 or 6 inputs'
   
   else:
         
      iSPHtoCAR = -1
      
      ## transform local directional vectors
         
      ## vectorize geopack_08 call using NumPy's vectorize() method:
      ## Note:  vectorize()'s doc string states that this is not very speed
      ##        optimized...Geopack should be ported to Python, or re-written
      ##        to accept array inputs, if speed is an issue. -EJR 10/2013
      ##
      
      #output = ghostpy.transform.geopack_08.sphcar_08(0.,0.,0.,x,y,z, iSPHtoCAR)
      v_gp08a = p.vectorize(ghostpy.transform.geopack_08.sphcar_08)
      rho,theta,phi,_,_,_,_ = v_gp08a(0.,0.,0.,x,y,z, iSPHtoCAR)
      
      #dx,dy,dz = ghostpy.transform.geopack_08.bspcar_08(theta,phi,drho,dtheta,dphi)
      v_gp08b = p.vectorize(ghostpy.transform.geopack_08.bcarsp_08)
      drho,dtheta,dphi = v_gp08b(x,y,z,dx,dy,dz)
      
      ## if scalar inputs, generate scalar outputs -EJR 10/2013
      return (tuple(p.asscalar(otmp) for otmp in (phi,theta,rho,dphi,dtheta,drho)) if
              (p.isscalar(x) and p.isscalar(y) and p.isscalar(z) and
               p.isscalar(dx) and p.isscalar(dy) and p.isscalar(dz)) else
              (phi,theta,rho,dphi,dtheta,drho) )


def SMtoGSM(x,y,z, dateTime):
    """
    >>> SMtoGSM(1,2,3, datetime.datetime(2009,1,27,0,0,0)) # doctest:+ELLIPSIS
    (-0.126..., 2.0, 3.159...)
    """
    if ghostpy.transform.transformer == 'CXFORM':
        (xGSM,yGSM,zGSM) = ghostpy.transform.cxform.transform('SM','GSM',x,y,z,
                                                            dateTime.year, dateTime.month, dateTime.day,
                                                            dateTime.hour, dateTime.minute, dateTime.second)
        
        return (xGSM,yGSM,zGSM)
    elif ghostpy.transform.transformer == 'GEOPACK':
        ghostpy.transform.geopack.recalc(dateTime.year, dateTime.timetuple().tm_yday,
                                       dateTime.hour, dateTime.minute, dateTime.second)
        iSMtoGSM = 1
        output = ghostpy.transform.geopack.smgsm(x,y,z, 0.,0.,0., iSMtoGSM)
        return (output[3],output[4],output[5])
    elif ghostpy.transform.transformer == 'GEOPACK_08':
        ghostpy.transform.geopack_08.recalc_08(dateTime.year, dateTime.timetuple().tm_yday,
                                             dateTime.hour, dateTime.minute, dateTime.second,
                                             -400.0, 0.0, 0.0)
        iSMtoGSM = 1
        
        ## vectorize geopack_08 call using NumPy's vectorize() method:
        ## Note:  vectorize()'s doc string states that this is not very speed
        ##        optimized...Geopack should be ported to Python, or re-written
        ##        to accept array inputs, if speed is an issue. -EJR 10/2013
        ##
        #output = ghostpy.transform.geopack_08.smgsw_08(x,y,z, 0.,0.,0., iSMtoGSM)
        v_gp08 = p.vectorize(ghostpy.transform.geopack_08.smgsw_08)
        output = v_gp08(x,y,z, 0.,0.,0., iSMtoGSM)
        
        ## if scalar inputs, generate scalar outputs -EJR 10/2013
        return ((output[3],output[4],output[5]) if isinstance(x,p.ndarray) else
                tuple(p.asscalar(otmp) for otmp in output[3:6]))
        
    else:
        raise Exception, 'Coordinate system transformer undefined'

def GSMtoSM(x,y,z, dateTime):
    """
    >>> GSMtoSM(1,2,3, datetime.datetime(2009,1,27,0,0,0)) # doctest:+ELLIPSIS
    (1.997..., 2.0, 2.451...)
    """
    if ghostpy.transform.transformer == 'CXFORM':
        (xSM,ySM,zSM) = ghostpy.transform.cxform.transform('GSM','SM',x,y,z,
                                                         dateTime.year, dateTime.month, dateTime.day,
                                                         dateTime.hour, dateTime.minute, dateTime.second)
        
        return (xSM,ySM,zSM)
    elif ghostpy.transform.transformer == 'GEOPACK':
        ghostpy.transform.geopack.recalc(dateTime.year, dateTime.timetuple().tm_yday,
                                       dateTime.hour, dateTime.minute, dateTime.second)
        iSMtoGSM = -1
        output = ghostpy.transform.geopack.smgsm(0.,0.,0., x,y,z, iSMtoGSM)
        return (output[0],output[1],output[2])
    elif ghostpy.transform.transformer == 'GEOPACK_08':
        ghostpy.transform.geopack_08.recalc_08(dateTime.year, dateTime.timetuple().tm_yday,
                                             dateTime.hour, dateTime.minute, dateTime.second,
                                             -400.0, 0.0, 0.0)
        iSMtoGSM = -1
        
        ## vectorize geopack_08 call using NumPy's vectorize() method:
        ## Note:  vectorize()'s doc string states that this is not very speed
        ##        optimized...Geopack should be ported to Python, or re-written
        ##        to accept array inputs, if speed is an issue. -EJR 10/2013
        ##
        #output = ghostpy.transform.geopack_08.smgsw_08(0.,0.,0., x, y, z, iSMtoGSM)
        v_gp08 = p.vectorize(ghostpy.transform.geopack_08.smgsw_08)
        output = v_gp08(0., 0., 0., x, y, z, iSMtoGSM)
        
        ## if scalar inputs, generate scalar outputs -EJR 10/2013
        return ((output[0],output[1],output[2]) if isinstance(x,p.ndarray) else
                tuple(p.asscalar(otmp) for otmp in output[0:3]))
        
    else:
        raise Exception, 'Coordinate system transformer undefined'

def GSEtoGSM(x,y,z, dateTime):
    """
    >>> GSEtoGSM(1,2,3, datetime.datetime(2009,1,27,0,0,0)) # doctest:+ELLIPSIS
    (0.99..., 0.540..., 3.564...)
    """
    if ghostpy.transform.transformer == 'CXFORM':
        (xGSM,yGSM,zGSM) = ghostpy.transform.cxform.transform('GSE','GSM', x,y,z,
                                                            dateTime.year, dateTime.month, dateTime.day,
                                                            dateTime.hour, dateTime.minute, dateTime.second)
        return (xGSM,yGSM,zGSM)
    elif ghostpy.transform.transformer == 'GEOPACK':
        ghostpy.transform.geopack.recalc(dateTime.year, dateTime.timetuple().tm_yday,
                                       dateTime.hour, dateTime.minute, dateTime.second)
        iGSEtoGSM = -1
        output = ghostpy.transform.geopack.gsmgse(0.,0.,0., x,y,z, iGSEtoGSM)
        return (output[0],output[1],output[2])
    elif ghostpy.transform.transformer == 'GEOPACK_08':
        ghostpy.transform.geopack_08.recalc_08(dateTime.year, dateTime.timetuple().tm_yday,
                                             dateTime.hour, dateTime.minute, dateTime.second,
                                             -400.0, 0.0, 0.0)
        iGSEtoGSM = -1
        
        ## vectorize geopack_08 call using NumPy's vectorize() method:
        ## Note:  vectorize()'s doc string states that this is not very speed
        ##        optimized...Geopack should be ported to Python, or re-written
        ##        to accept array inputs, if speed is an issue. -EJR 10/2013
        ##
        #output = ghostpy.transform.geopack_08.gswgse_08(0.,0.,0., x,y,z, iGSEtoGSM)
        v_gp08 = p.vectorize(ghostpy.transform.geopack_08.gswgse_08)
        output = v_gp08(0., 0., 0., x, y, z, iGSEtoGSM)
        
        ## if scalar inputs, generate scalar outputs -EJR 10/2013
        return ((output[0],output[1],output[2]) if isinstance(x,p.ndarray) else
                tuple(p.asscalar(otmp) for otmp in output[0:3]))
                
    else:
        raise Exception, 'Coordinate system transformer undefined'

def HEEQtoGSM(x,y,z, dateTime):
    """
    >>> HEEQtoGSM(1,2,3, datetime.datetime(2009,1,27,0,0,0)) # doctest:+ELLIPSIS
    (-0.69972..., -2.95..., 2.18...)
    """
    (xGSM,yGSM,zGSM) = ghostpy.transform.cxform.transform('HEEQ', 'GSM', x,y,z,
                                                        dateTime.year, dateTime.month, dateTime.day,
                                                        dateTime.hour, dateTime.minute, dateTime.second)
    return (xGSM,yGSM,zGSM)


def GEOtoMAG(x,y,z, dateTime):
    """
    Converts geographic to dipole (mag) coordinates.
    
    >>> GEOtoMAG(-0.780, -0.613, 0.128, datetime.datetime(2009,1,27,0,0,0)) # doctest:+ELLIPSIS
    ()
    """
    ghostpy.transform.geopack_08.recalc_08(dateTime.year, dateTime.timetuple().tm_yday,
                                             dateTime.hour, dateTime.minute, dateTime.second,
                                             -400.0, 0.0, 0.0)
    
    iGEOtoMAG = 1
    
    ## vectorize geopack_08 call using NumPy's vectorize() method:
    ## Note:  vectorize()'s doc string states that this is not very speed
    ##        optimized...Geopack should be ported to Python, or re-written
    ##        to accept array inputs, if speed is an issue. -EJR 10/2013
    ##
    #output = ghostpy.transform.geopack_08.gswgse_08(x, y, z, 0.,0.,0., iGEOtoMAG)
    v_gp08 = p.vectorize(ghostpy.transform.geopack_08.geomag_08)
    output = v_gp08(x, y, z, 0., 0., 0., iGEOtoMAG)
     
    ## if scalar inputs, generate scalar outputs -EJR 10/2013
    return ((output[3],output[4],output[5]) if isinstance(x,p.ndarray) else
            tuple(p.asscalar(otmp) for otmp in output[3:6]))
    

def MAGtoGEO(x,y,z, dateTime):
    """
    Converts dipole (mag) to geographic coordinates.
    
    >>> MAGtoGEO(0.3133, -0.9313, 0.1857, datetime.datetime(2009,1,27,0,0,0)) # doctest:+ELLIPSIS
    ()
    """
    ghostpy.transform.geopack_08.recalc_08(dateTime.year, dateTime.timetuple().tm_yday,
                                             dateTime.hour, dateTime.minute, dateTime.second,
                                             -400.0, 0.0, 0.0)
    
    iGEOtoMAG = -1
    
    ## vectorize geopack_08 call using NumPy's vectorize() method:
    ## Note:  vectorize()'s doc string states that this is not very speed
    ##        optimized...Geopack should be ported to Python, or re-written
    ##        to accept array inputs, if speed is an issue. -EJR 10/2013
    ##
    #output = ghostpy.transform.geopack_08.gswgse_08(0.,0.,0., x, y, z, iGEOtoMAG)
    v_gp08 = p.vectorize(ghostpy.transform.geopack_08.geomag_08)
    output = v_gp08(0., 0., 0., x, y, z, iGEOtoMAG)
     
    ## if scalar inputs, generate scalar outputs -EJR 10/2013
    return ((output[0],output[1],output[2]) if isinstance(x,p.ndarray) else
            tuple(p.asscalar(otmp) for otmp in output[0:3]))
    


def SMtoMAG(x,y,z, dateTime):
   """
   Converts solar magnetic coordinates to dipole magnetic coordinates
   """
   ghostpy.transform.geopack_08.recalc_08(dateTime.year, dateTime.timetuple().tm_yday,
                                             dateTime.hour, dateTime.minute, dateTime.second,
                                             -400.0, 0.0, 0.0)
   
   iSMtoMAG = -1
   
   ## vectorize geopack_08 call using NumPy's vectorize() method:
   ## Note:  vectorize()'s doc string states that this is not very speed
   ##        optimized...Geopack should be ported to Python, or re-written
   ##        to accept array inputs, if speed is an issue. -EJR 12/2013
   ##
   
   v_gp08 = p.vectorize(ghostpy.transform.geopack_08.magsm_08)
   output =  v_gp08(0., 0., 0., x, y, z, iSMtoMAG)
   
   ## if scalar inputs, generate scalar outputs -EJR 12/2013
   return((output[0],output[1],output[2]) if isinstance(x,p.ndarray) else
            tuple(p.asscalar(otmp) for otmp in output[0:3]))
   

def MAGtoSM(x,y,z, dateTime):
   """
   Converts dipole magnetic to solar magnetic coordinates
   """
   ghostpy.transform.geopack_08.recalc_08(dateTime.year, dateTime.timetuple().tm_yday,
                                             dateTime.hour, dateTime.minute, dateTime.second,
                                             -400.0, 0.0, 0.0)
   
   iMAGtoSM = 1
   
   
   ## vectorize geopack_08 call using NumPy's vectorize() method:
   ## Note:  vectorize()'s doc string states that this is not very speed
   ##        optimized...Geopack should be ported to Python, or re-written
   ##        to accept array inputs, if speed is an issue. -EJR 12/2013
   ##
   v_gp08 = p.vectorize(ghostpy.transform.geopack_08.magsm_08)
   output =  v_gp08(x, y, z, 0., 0., 0., iMAGtoSM)
   
   ## if scalar inputs, generate scalar outputs -EJR 12/2013
   return((output[3],output[4],output[5]) if isinstance(x,p.ndarray) else
            tuple(p.asscalar(otmp) for otmp in output[3:6]))
  


def GEOtoGSM(x,y,z, dateTime):
   """
   Converts geographic to Geocentric Solar Magnetic coordinates
   """
   ghostpy.transform.geopack_08.recalc_08(dateTime.year, dateTime.timetuple().tm_yday,
                                             dateTime.hour, dateTime.minute, dateTime.second,
                                             -400.0, 0.0, 0.0)
   
   iGEOtoGSM = 1
   
   
   ## vectorize geopack_08 call using NumPy's vectorize() method:
   ## Note:  vectorize()'s doc string states that this is not very speed
   ##        optimized...Geopack should be ported to Python, or re-written
   ##        to accept array inputs, if speed is an issue. -EJR 12/2013
   ##
   v_gp08 = p.vectorize(ghostpy.transform.geopack_08.geogsw_08)
   output =  v_gp08(x, y, z, 0., 0., 0., iGEOtoGSM)
   
   ## if scalar inputs, generate scalar outputs -EJR 12/2013
   return((output[3],output[4],output[5]) if isinstance(x,p.ndarray) else
            tuple(p.asscalar(otmp) for otmp in output[3:6]))
  
   
def GSMtoGEO(x,y,z, dateTime):
   """
   Converts Geocentric Solar Magnetic to geographic coordinates
   """
   ghostpy.transform.geopack_08.recalc_08(dateTime.year, dateTime.timetuple().tm_yday,
                                             dateTime.hour, dateTime.minute, dateTime.second,
                                             -400.0, 0.0, 0.0)
   
   iGSMtoGEO = -1
   
   
   ## vectorize geopack_08 call using NumPy's vectorize() method:
   ## Note:  vectorize()'s doc string states that this is not very speed
   ##        optimized...Geopack should be ported to Python, or re-written
   ##        to accept array inputs, if speed is an issue. -EJR 12/2013
   ##
   v_gp08 = p.vectorize(ghostpy.transform.geopack_08.geogsw_08)
   output =  v_gp08(0., 0., 0., x, y, z, iGEOtoGSM)
   
   ## if scalar inputs, generate scalar outputs -EJR 12/2013
   return((output[0],output[1],output[2]) if isinstance(x,p.ndarray) else
            tuple(p.asscalar(otmp) for otmp in output[0:3]))




def GEOtoSM(x,y,z, dateTime):
   """
   Converts geographic to Solar Magnetic coordinates
   """
   ghostpy.transform.geopack_08.recalc_08(dateTime.year, dateTime.timetuple().tm_yday,
                                             dateTime.hour, dateTime.minute, dateTime.second,
                                             -400.0, 0.0, 0.0)
   
   iGEOtoGSM = 1
   iGSMtoSM = -1
   
   ## vectorize geopack_08 call using NumPy's vectorize() method:
   ## Note:  vectorize()'s doc string states that this is not very speed
   ##        optimized...Geopack should be ported to Python, or re-written
   ##        to accept array inputs, if speed is an issue. -EJR 12/2013
   ##
   v_gp08a = p.vectorize(ghostpy.transform.geopack_08.geogsw_08)
   _,_,_,xgsm,ygsm,zgsm =  v_gp08a(x, y, z, 0., 0., 0., iGEOtoGSM)
   
   v_gp08b = p.vectorize(ghostpy.transform.geopack_08.smgsw_08)
   xsm,ysm,zsm,_,_,_ =  v_gp08b(0., 0., 0., xgsm, ygsm, zgsm, iGSMtoSM)
   
   
   ## if scalar inputs, generate scalar outputs -EJR 12/2013
   return((xsm,ysm,zsm) if isinstance(x,p.ndarray) else
            tuple(p.asscalar(otmp) for otmp in (xsm,ysm,zsm)))
  
   
def SMtoGEO(x,y,z, dateTime):
   """
   Converts Solar Magnetic to geographic coordinates
   """
   ghostpy.transform.geopack_08.recalc_08(dateTime.year, dateTime.timetuple().tm_yday,
                                             dateTime.hour, dateTime.minute, dateTime.second,
                                             -400.0, 0.0, 0.0)
   
   iSMtoGSM = 1
   iGSMtoGEO = -1
   
   
   ## vectorize geopack_08 call using NumPy's vectorize() method:
   ## Note:  vectorize()'s doc string states that this is not very speed
   ##        optimized...Geopack should be ported to Python, or re-written
   ##        to accept array inputs, if speed is an issue. -EJR 12/2013
   ##
   v_gp08a = p.vectorize(ghostpy.transform.geopack_08.smgsw_08)
   _,_,_,xgsm,ygsm,zgsm =  v_gp08a(x, y, z, 0., 0., 0., iSMtoGSM)
   
   v_gp08b = p.vectorize(ghostpy.transform.geopack_08.geogsw_08)
   xgeo,ygeo,zgeo,_,_,_ =  v_gp08b(0., 0., 0., xgsm, ygsm, zgsm, iGSMtoGEO)
   
   
   ## if scalar inputs, generate scalar outputs -EJR 12/2013
   return((xgeo,ygeo,zgeo) if isinstance(x,p.ndarray) else
            tuple(p.asscalar(otmp) for otmp in (xgeo,ygeo,zgeo)))
   
