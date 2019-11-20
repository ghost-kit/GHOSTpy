import cxform

def printCoordinates(system, xyz):
    print 'Input Vector (%s):\t %f\t %f\t %f' % (system, xyz[0], xyz[1], xyz[2])

base = 'J2000'
systems = ['GEO', 'GSE', 'GSM', 'SM']
year=2005
month=3
day=2
hour=9
minute=28
second=11

xyzIn = [-896921337.28302002, 220296912.43620300, 44419205.01961136]
printCoordinates(base, xyzIn)

for trans in systems:
    xyzOut = cxform.transform(base, trans, xyzIn[0], xyzIn[1], xyzIn[2], year,month,day, hour,minute,second)
    printCoordinates(trans, xyzOut)

print "Output exepected:"
print """Input Vector (J2000):\t -896921337.283020\t 220296912.436203\t 44419205.019611
Output Vector (GEO):\t -664946479.875617\t -640940392.874836\t 44869503.669364
Output Vector (GSE):\t -920802840.504282\t -70523436.668648\t -46046221.090855
Output Vector (GSM):\t -920802840.504282\t -58593148.345827\t -60503326.877360
Output Vector (SM):\t -915527671.753106\t -58593148.345827\t 115531839.327171
"""
