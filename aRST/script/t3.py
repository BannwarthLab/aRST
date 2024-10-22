import math
from collections import Counter


def round_charge(charge,tol):
    # round up to higher int for .tol
    b,a = math.modf(charge)
    if abs(b)>=tol:
        if b<=0:
            out = a-1
        else:
            out = a+1
    else:
        out = a
    return int(out)



print(round_charge(0.3,0.5))