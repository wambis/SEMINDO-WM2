
import numpy as np
from math import acos, pi
"""
This code checks in a point (x, y) in located inside or at the boundary of a
polygone given by [(x1, y1), ..., (xn, yn)].
"""

def check_inside(x, y, poly):
    # check is point is inside
    n = len(poly)
    inside = False

    p1x,p1y = poly[0]
    for i in range(n+1):
       p2x,p2y = poly[i % n]
       if y > min(p1y,p2y):
          if y <= max(p1y,p2y):
             if x <= max(p1x,p2x):
                if p1y != p2y:
                   xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                if p1x == p2x or x <= xints:
                   inside = not inside
       p1x,p1y = p2x,p2y
    if inside:
        return 'IN'

def check_2Dpolygons(x, y, poly):

    n = len(poly)

    # check if point is a vertex
    if (x,y) in poly: return "IN"

    # check if the point is on a boundary
    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x, p2y = poly[i % n]
        if x >= min(p1x, p2x) and x <= max(p1x, p2x):
            if y >= min(p1y, p2y) and y <= max(p1y, p2y):
                if p1x == p2x == x :
                    return "IN"
                else:
                    a = (p2y - p1y)/(p2x - p1x)
                    b = p1y - a*p1x
                    if y == a*x + b: return "IN"
        p1x, p1y = p2x, p2y

    inside = check_inside(x, y, poly)

    return inside

def main():
    import matplotlib.pyplot as plt
    poly2 = [(1,1), (1,5), (5,5), (5,1), (1,1)]
#    poly2 = [(1, 1), (1, 1), (1, 5), (1, 5), (1, 1), (5, 1), (5, 5), (1, 5), (1, 1), (5, 1), (5, 5), (5, 5), (5, 1), (5, 1), (5, 5), (1, 5), (1, 1), (1, 1), (5, 1), (5, 1), (1, 1), (1, 1), (1, 5), (5, 5), (5, 1), (1, 1)]

   # b = []
   # for a in poly2:
   #     if a not in b:
   #         b.append(a)

    x = 3
    y = 3
    if check_2Dpolygons(x, y, poly2) == 'IN':
        print("IN")
    else:
        print("OUT")

    X, Y = zip(*poly2)
    xx ,yy = zip(*poly2)
    plt.scatter(xx, yy, s=300, c='y', marker='*')
    plt.plot(xx, yy, 'k')
    plt.plot(X, Y, 'b')
    plt.plot(x, y, 'ro')
    plt.xlim(0, 10)
    plt.ylim(0, 10)
    plt.show()

if __name__ == '__main__':
    main()
