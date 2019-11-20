"""
Exercises the Python bindings to CXFORM.  You must compile & install
CXFORM before this test will pass.  See pkgs/cxform-0.71_python for
details.
"""

# Give a more meaningful error to people who don't read instructions:
try:
    import cxform
except ImportError:
    raise ImportError, 'No module named "cxform".  Did you remember to compile & install the Python wrapping for CXFORM?'

import unittest

class TestCXForm(unittest.TestCase):
    """
    This class tests the CX form library wrapped via Python.    
    """
    def setUp(self):
        # Values chosen from test case of cxform distribution
        self.base = 'J2000'
        self.year=2005
        self.month=3
        self.day=2
        self.hour=9
        self.minute=28
        self.second=11
        self.x = -896921337.28302002
        self.y =  220296912.43620300
        self.z =  44419205.01961136

        self.PERCENT_DIFFERENCE = 0.175


    def __assertPercentDiff(self, x, y):
        p = self.__percentDiff( x, y)

        errMesg = "Percent difference between %f and %f = %f exceeds the threshold of %f.  For small differences, this may be attributed to rounding error."
        self.failIf( p > self.PERCENT_DIFFERENCE, msg=( errMesg % (x, y, p, self.PERCENT_DIFFERENCE) ) )
        
    def __percentDiff(self, x, y):
        return ( abs((x-y)/x) )

    def test_2GEO(self):
        out = cxform.transform('J2000', 'GEO',
                               self.x, self.y, self.z,
                               self.year, self.month, self.day, self.hour, self.minute, self.second)

        self.__assertPercentDiff( -664946479.875617, out[0] )
        self.__assertPercentDiff( -640940392.874836, out[1] )        
        self.__assertPercentDiff(   44869503.669364, out[2] )

    def test_2GSE(self):
        out = cxform.transform('J2000', 'GSE',
                               self.x, self.y, self.z,
                               self.year, self.month, self.day, self.hour, self.minute, self.second)
        self.__assertPercentDiff( -920802840.504282, out[0] )
        self.__assertPercentDiff( -64160872.043671027, out[1] )
        self.__assertPercentDiff(  -46046221.090855, out[2] )

    def test_2GSM(self):
        out = cxform.transform('J2000', 'GSM',
                               self.x, self.y, self.z,
                               self.year, self.month, self.day, self.hour, self.minute, self.second)
        self.__assertPercentDiff( -920802840.504282, out[0] )
        self.__assertPercentDiff( -52303890.695484,  out[1] )
        self.__assertPercentDiff( -60503326.877360,  out[2] )

    def test_2SM(self):
        out = cxform.transform('J2000', 'SM',
                               self.x, self.y, self.z,
                               self.year, self.month, self.day, self.hour, self.minute, self.second)
        self.__assertPercentDiff( -915527671.753106, out[0] )
        self.__assertPercentDiff( -52303890.6954846,  out[1] )
        self.__assertPercentDiff( 115531839.327171,  out[2] )

    def test_merkin(self):
        """
        Slava Merkin was experience strange behavior with cxform on
        certain inputs.  Let's try those inputs here & see if we can
        break cxform!
        """
        x=0.25
        y = -2.63
        z = -0.4
        year = 2008
        month = 8
        day = 1
        hour = 0
        minute = 0
        second = 0

        # Transform twice in a row... make sure we get the same thing:
        out = cxform.transform('GSE', 'GSM',
                               x,y,z,
                               year,month,day, hour,minute,second)
        self.__assertPercentDiff( 0.25, out[0] )
        self.__assertPercentDiff( -2.6556295971324912, out[1] )
        self.__assertPercentDiff( -0.15662516666844989, out[2] )

        out = cxform.transform('GSE', 'GSM',
                               x,y,z,
                               year,month,day, hour,minute,second)
        self.__assertPercentDiff( 0.25, out[0] )
        self.__assertPercentDiff( -2.6556295971324912, out[1] )
        self.__assertPercentDiff( -0.15662516666844989, out[2] )

        # Now do the inverse transform & make sure we get the same thing.
        out = cxform.transform('GSM', 'GSE',
                               out[0], out[1], out[2],
                               year,month,day, hour,minute,second)
        self.__assertPercentDiff( 0.25, out[0] )
        self.__assertPercentDiff( -2.6299999999999994, out[1] )
        self.__assertPercentDiff( -0.39999999999999997, out[2] )

    def test_geopack(self):
        """GSE to SM conversion verified by Slava Merkin using GEOPACK."""
        x=0.25
        y = -2.63
        z = -0.4
        year = 2008
        month = 8
        day = 1
        hour = 0
        minute = 0
        second = 0
        out = cxform.transform('GSE', 'SM',
                               x,y,z,
                               year,month,day,hour,minute,second)
        self.__assertPercentDiff( .2830,   out[0] )
        self.__assertPercentDiff( -2.6554, out[1] )
        self.__assertPercentDiff( -.0899,  out[2] )
        
        # Check the inverse:
        out = cxform.transform('SM', 'GSE',
                               out[0], out[1], out[2],
                               year,month,day,hour,minute,second)
        self.__assertPercentDiff( 0.25,  out[0] )
        self.__assertPercentDiff( -2.63, out[1] )
        self.__assertPercentDiff( -0.4,  out[2] )

    def test_geopack(self):
        """GSE to SM conversion verified by Slava Merkin using GEOPACK."""
        x=0.25
        y = -2.63
        z = -0.4
        # First day of NHL hockey season in '94-95. Stupid strike!
        year = 1995
        month = 1
        day = 20
        hour = 19
        minute = 0
        second = 0
        out = cxform.transform('GSE', 'SM',
                               x,y,z,
                               year,month,day,hour,minute,second)
        self.__assertPercentDiff( 0.025938911859381325,   out[0] )
        self.__assertPercentDiff( -2.3804487166857329, out[1] )
        self.__assertPercentDiff( -1.2134184699432964,  out[2] )
        
        # Check the inverse:
        out = cxform.transform('SM', 'GSE',
                               out[0], out[1], out[2],
                               year,month,day,hour,minute,second)
        self.__assertPercentDiff( 0.25,  out[0] )
        self.__assertPercentDiff( -2.63, out[1] )
        self.__assertPercentDiff( -0.4,  out[2] )
 

if __name__ == "__main__":
    unittest.main()
