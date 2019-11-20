"""
Exercises the Python bindings to GEOPACK.  You must compile & install
GEOPACK before this test will pass.  See pkgs/geopack-2005 for
details.

Unit test errors suggest that GEOPACK may not have been compiled with
8-byte reals.  See pkgs/geopack-2005/README for details.
"""

import datetime

# Give a more meaningful error to people who don't read instructions:
try:
    import geopack
except ImportError:
    print 'Module "geopack" not found.  This is okay so long as geopack_08 is installed.'
    #raise ImportError, 'No module named "geopack".  Did you remember to compile & install the Python wrapping for GEOPACK?'
else:    
    import unittest
    
    class TestGeopack(unittest.TestCase):
        """
        This class tests the GEOPACK library wrapped via Python & f2py.
        """
        def setUp(self):
            # Inputs inspired by unit tests for my Python bindings to cxform
            self.d = datetime.datetime(year=2005, month=3, day=2,
                                       hour=9, minute=28, second=11)
    
            # Pre-calculated values in a variety of coordinate systems:
            self.geo = ( -664946479.875617,  -640940392.874836,     44869503.669364 )
            self.gse = ( -920803133.4333384,  -70514255.353478357, -46054423.548106655)        
            self.gsm = ( -920803133.4333384,  -58583041.626907118, -60508655.332473055)        
            self.sm =  ( -915530691.724154,   -58583041.626907118, 115513031.39425668 )
    
            # GEOPACK constants.  These should probably be added to the GEOPACK module itself. 
            self.iGSMtoGSE = 1
            self.iGSEtoGSM = -1
    
            self.iSMtoGSM = 1
            self.iGSMtoSM = -1
    
            self.iGEOtoGSM = 1
            self.iGSMtoGEO = -1
    
            self.PERCENT_DIFFERENCE = 0.02
    
            # Must set the date & time prior to carrying out a transformation.
            geopack.recalc(self.d.year,
                           self.d.timetuple().tm_yday,
                           self.d.hour,
                           self.d.minute,
                           self.d.second)

        def __assertPercentDiff(self, x, y):
            p = self.__percentDiff( x, y)
            
            errMesg = "Percent difference between %f and %f = %f exceeds the threshold of %f.  For small differences, this may be attributed to rounding error."
            self.failIf( p > self.PERCENT_DIFFERENCE, msg=( errMesg % (x, y, p, self.PERCENT_DIFFERENCE) ) )
            
    
        def __percentDiff(self, x, y):
            return ( abs((x-y)/x) )
            
    
        def test_GEOGSM(self):
            
            output = geopack.geogsm(self.geo[0], self.geo[1], self.geo[2], 0.0,0.0,0.0,  self.iGEOtoGSM)
    
            self.__assertPercentDiff( self.gsm[0], output[3] )
            self.__assertPercentDiff( self.gsm[1], output[4] )
            self.__assertPercentDiff( self.gsm[2], output[5] )
    
            # Make sure the inverse conversion works        
            output = geopack.geogsm(0.0, 0.0, 0.0, output[3], output[4], output[5], self.iGSMtoGEO)
    
            self.__assertPercentDiff( self.geo[0], output[0] )
            self.__assertPercentDiff( self.geo[1], output[1] )
            self.__assertPercentDiff( self.geo[2], output[2] )
                           
        def test_SMGSM(self):
            
            output = geopack.smgsm(0.0, 0.0, 0.0, self.gsm[0], self.gsm[1], self.gsm[2], self.iGSMtoSM)
    
            self.__assertPercentDiff( self.sm[0], output[0] )
            self.__assertPercentDiff( self.sm[1], output[1] )
            self.__assertPercentDiff( self.sm[2], output[2] )
    
            output = geopack.smgsm(output[0], output[1], output[2], 0.0, 0.0, 0.0, self.iSMtoGSM)
    
            self.__assertPercentDiff( self.gsm[0], output[3] )
            self.__assertPercentDiff( self.gsm[1], output[4] )
            self.__assertPercentDiff( self.gsm[2], output[5] )
    
        def test_GSMGSE(self):
            
            output = geopack.gsmgse(self.gsm[0], self.gsm[1], self.gsm[2], 0.,0.,0., self.iGSMtoGSE)
            self.__assertPercentDiff( self.gse[0], output[3] )
            self.__assertPercentDiff( self.gse[1], output[4] )
            self.__assertPercentDiff( self.gse[2], output[5] )
    
            output = geopack.gsmgse(0.,0.,0., output[3], output[4], output[5], self.iGSEtoGSM)
            self.__assertPercentDiff( self.gsm[0], output[0] )
            self.__assertPercentDiff( self.gsm[1], output[1] )
            self.__assertPercentDiff( self.gsm[2], output[2] )
    
        def test_merkin(self):
            """
            Slava Merkin experienced strange behavior with cxform on
            certain inputs.  Compare those inputs here against GEOPACK.
            """
            d = datetime.datetime(year = 2008, month = 8, day = 1,
                                  hour = 0, minute = 0, second = 0)
    
            geopack.recalc(d.year, d.timetuple().tm_yday,
                           d.hour, d.minute, d.second)
    
            gse = (0.25, -2.63, -0.4)
    
            out = geopack.gsmgse(0.,0.,0., gse[0], gse[1], gse[2], self.iGSEtoGSM)
            self.__assertPercentDiff( 0.25, out[0] )
            self.__assertPercentDiff( -2.6556295971324912, out[1] )
            self.__assertPercentDiff( -0.15662516666844989, out[2] )
    
            out = geopack.gsmgse(out[0],out[1],out[2], 0.,0.,0., self.iGSMtoGSE)
            self.__assertPercentDiff( gse[0], out[3] )
            self.__assertPercentDiff( gse[1], out[4] )
            self.__assertPercentDiff( gse[2], out[5] )
    
        def test_cxform(self):
            d = datetime.datetime(year = 2008, month = 8, day = 1,
                                  hour = 0, minute = 0, second = 0)
    
            geopack.recalc(d.year, d.timetuple().tm_yday,
                           d.hour, d.minute, d.second)
    
            gse = (0.25, -2.63, -0.4)
    
            out = geopack.gsmgse(0.,0.,0., gse[0],gse[1],gse[2], self.iGSEtoGSM)
            out = geopack.smgsm(0.,0.,0., out[0],out[1],out[2], self.iGSMtoSM)
    
            self.__assertPercentDiff( .2830,   out[0] )
            self.__assertPercentDiff( -2.6554, out[1] )
            self.__assertPercentDiff( -0.087301346866137039,  out[2] )
    
            # Check the inverse:
            out = geopack.smgsm(out[0],out[1],out[2], 0.,0.,0., self.iSMtoGSM)
            out = geopack.gsmgse(out[3],out[4],out[5], 0.,0.,0., self.iGSMtoGSE)
    
            self.__assertPercentDiff( gse[0], out[3] )
            self.__assertPercentDiff( gse[1], out[4] )
            self.__assertPercentDiff( gse[2], out[5] )
    
        def test_1990s(self):
            # First day of NHL hockey season in '94-95. Stupid strike!
            d = datetime.datetime(year = 1995, month = 1, day = 20,
                                  hour = 19, minute = 0, second = 0)
    
            geopack.recalc(d.year, d.timetuple().tm_yday,
                           d.hour, d.minute, d.second)
    
            gse = (0.25, -2.63, -0.4)
    
            out = geopack.gsmgse(0.,0.,0., gse[0],gse[1],gse[2], self.iGSEtoGSM)
            out = geopack.smgsm(0.,0.,0., out[0],out[1],out[2], self.iGSMtoSM)
    
            self.__assertPercentDiff( 0.02189803712013394,   out[0] )
            self.__assertPercentDiff( -2.3804487166857329, out[1] )
            self.__assertPercentDiff( -1.2134184699432964,  out[2] )
    
            # Check the inverse:
            out = geopack.smgsm(out[0],out[1],out[2], 0.,0.,0., self.iSMtoGSM)
            out = geopack.gsmgse(out[3],out[4],out[5], 0.,0.,0., self.iGSMtoGSE)
    
            self.__assertPercentDiff( gse[0], out[3] )
            self.__assertPercentDiff( gse[1], out[4] )
            self.__assertPercentDiff( gse[2], out[5] )
    
    if __name__ == "__main__":
        unittest.main()
