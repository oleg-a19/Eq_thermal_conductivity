import numpy as np

class Point2D(object):
    def __init__(self,x,y):
        self.x = x
        self.y = y
    def __repr__(self):
        return f'Point ({self.x}, {self.y})'
    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

class Grid:
    def __init__(self, corner_points, K, M):
        self.corner_points = corner_points
        self.M = M
        self.K = K
        self.points = np.empty(K*M, dtype=object)
    
    def __repr__(self):
        return f'Grid({self.points})'
      

    def fill_grid(self):
        hy = (self.corner_points[1].y-self.corner_points[0].y)/(self.M-1)
        hx = (self.corner_points[3].x-self.corner_points[0].x)/(self.K-1)
        # k - индект оси x
        # m - индекс оси y
        for k in range(self.K):
            for m in range(self.M):
                i = m + k*self.M
                self.points[i] = Point2D(self.corner_points[0].x + k*hx, self.corner_points[0].y + m*hy)
