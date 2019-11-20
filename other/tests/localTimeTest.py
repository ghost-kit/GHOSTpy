import numpy as np

from prototype import inv_common as ih

# Test localtime conversions.
localtime1 = np.linspace(0, 24, 2000, endpoint=True)
print localtime1
error_count = 0
for lt in localtime1:
    unitVector = ih.get_location_from_localtime(lt)
    calc_localtime = ih.get_local_time_from_location(unitVector)
    if np.isclose(lt, calc_localtime):
        print "Successful Test for localtime: ", lt
    else:
        print "Failure in Test for localtime: ", lt
        print "Given: ", lt
        print "Received: ", calc_localtime
        error_count += 1

print "Number of Errors in reverse localtime lookup: ", error_count