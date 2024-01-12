from numpy import ones,vstack
from numpy.linalg import lstsq

def LineEquation(points = [(1,5),(3,4)]):
    x_coords, y_coords = zip(*points)
    A = vstack([x_coords,ones(len(x_coords))]).T
    m, c = lstsq(A, y_coords,  rcond = None)[0]
    pol = [c, m]
    return pol
    # print("Line Solution is y = {m}x + {c}".format(m=m,c=c))

# mypol = LineEquation()
# print(mypol)