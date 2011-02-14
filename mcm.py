import sys
import math
import random
import numpy
from numpy import linalg

class point:
    def __init__(self, u, v):
        self.u = u
        self.v = v
        self.r = -1
        self.theta = 0
        self.residential_population = 0
        self.commuter_population = 0
        self.max_population = 0
        self.num_towers = 0
        self.points_in_tower_range = 0
        self.towers_required = 0

    def __str__(self):
        string = """(%d,%d) = (%f, %f)
        %f res, %f com, %f max
        %d points in tower range, %d towers in range required
        """ % (self.u, self.v, self.r, self.theta, self.residential_population, self.commuter_population, self.max_population, self.points_in_tower_range, self.towers_required)
        #x = self.r * math.cos(self.theta)
        #y = self.r * math.sin(self.theta)
        #string = "%f %f %f" % (x, y, self.max_population)
        return string

class sqm:
    "square matrix"
    def __init__(self, n):
        "self.n is the square matrix's size, self.a is the matrix"
        self.n = n
        self.a = self.init_zero_array(n,n)
        self.exponent = 1
        self.fill_random()
        self.normalize_rows()

    def assign_matrix(self, a):
        "simply assigns a specifically defined matrix to this instance"
        self.a = a;
        (self.n, b) = self.a.shape

    def init_zero_array(self, m, n):
        "returns an m x n zero array" 
        a = [0] * n
        for i in range(0,n):
            a[i] = [0] * m
        return a

    def fill_random(self):
        "fills the current instance's array with random values in the range [0,1]"
        random.seed()
        for i in range(0,self.n):
            for j in range(0,self.n):
                rand = float(random.randint(0,1000)) / 1000
                self.a[i][j] = rand

    def normalize_rows(self):
        "normalizes the rows of the current matrix so that each row sums to 1"
        for i in range(0,self.n):
            rowsum = 0
            for j in range(0,self.n):
                rowsum += self.a[i][j]
            for j in range(0,self.n):
                self.a[i][j] = self.a[i][j] / rowsum

    def diagonalize(self):
        "finds the eigenvalue and eigenvector matrices for the current matrix to be used to diagonalize self.a"
        (self.eig, self.eigv) = linalg.eig(self.a)
        self.eig = numpy.diag(self.eig)
        self.eigvi = linalg.inv(self.eigv)
        

    def print_diagonalization(self):
        "prints the eigenvalue, eigenvector, and inverse eigenvector inverse matrices, as well as A"
        print "V-1 = ", self.eigvi
        print "L = ", self.eig
        print "V = ", self.eigv
        print "A = ", self.a
        print "V-1AV = ", numpy.dot(numpy.dot(self.eigvi,self.a),self.eigv)

    def iterate(self, k):
        "raises self.a to the kth power"
        self.exponent = k
        #self.ak = numpy.dot(numpy.dot(self.eigv,linalg.matrix_power(self.a,k)),self.eigvi)
        self.ak = linalg.matrix_power(self.a,k)

class LogicError(Exception):
    def __init__(self, value):
        self.parameter = value

    def __str__(self):
        return repr(self.parameter)

class grid:
    "triangular grid on a circular area"
    def __init__(self, gridsize=1, radius=40, population=1000):
        "defines a triangular grid of points with distance between adjacent points equal to `gridsize`, covering a circular area of radiua `radius`"
        self.gridsize = gridsize
        self.radius = radius
        self.points = {}
        self.p = []
        self.population = population
        self.population_vector = []
        self.commuting_vector = []
        self.create_grid()
        self.allocate_random_populations()
        self.transition_matrix = None

    def allocate_random_populations(self):
        "distributes population randomly"
        random.seed()
        total_assigned = 0
        for i in range(0,len(self.p)):
            rand = random.randint(0,1000)
            total_assigned += rand
            self.p[i].residential_population = rand

        for i in range(0,len(self.p)):
            self.p[i].residential_population = float(self.p[i].residential_population) / total_assigned * self.population
            self.population_vector.append(self.p[i].residential_population)


    def assign_population(self, v):
        "sets population vector to v"
        "given a vector whose number of elements corresponds to the number of points on the grid, assigns vector element i as the population of grid point number i"
        if len(v) != len(self.p):
            raise LogicError("Population vector size is not equal to the number of points on the grid.")
        for i in range(0,len(v)):
            self.p[i].residential_population = v[i]
        self.population_vector = v

    def use_transition_matrix(self, tmatrix):
        "calculates commuter and max populations"
        self.transition_matrix = tmatrix
        self.calculate_commuting_population(self.transition_matrix)
        self.calculate_max_populations()

    def calculate_commuting_population(self, tmatrix):
        "multiplies residential population vector by transition matrix"
        self.commuting_vector = numpy.dot(numpy.array(self.population_vector), tmatrix.a)
        for i in range(0,len(self.p)):
            self.p[i].commuter_population = self.commuting_vector[i]

    def calculate_max_populations(self):
        "max of residential and commuter populations"
        for i in range(0,len(self.p)):
            self.p[i].max_population = max(self.p[i].commuter_population,self.p[i].residential_population)

    def create_grid(self):
        "creates the grid for the area using current parameters and saves both a dict and a list of points"
        origin = point(0,0)
        origin.r = 0
        self.points[(0,0)] = origin
        self.p.append(origin)
        
        cur = 1
        while cur * self.gridsize < 2*self.radius/math.sqrt(3):
            num_points = cur * 6
            for i  in range(0,num_points):
                new_point = point(cur,i)
                imod = i % cur
                theta_offset = (i - imod) / cur * 2*math.pi/6
                if imod == 0:
                    new_point.r = cur * self.gridsize
                    new_point.theta = theta_offset
                else:
                    new_point.r = self.law_of_cosines(cur * self.gridsize ,
                                                      imod * self.gridsize ,
                                                      math.pi/3)
                    new_point.theta = theta_offset + math.acos((pow(new_point.r,2) + pow(cur * self.gridsize,2) - pow(imod * self.gridsize,2)) / (2 * cur*self.gridsize * new_point.r))
                if(new_point.r <= self.radius):
                    self.points[(cur,i)] = new_point
                    self.p.append(new_point)
            cur += 1

    def list_points(self):
        "prints all points in the grid"
        #for s in self.points.keys():
        #    print s
        for s in self.p:
            print s

    def list_points_to_file(self, output):
        "prints all points in the grid to file"
        #for s in self.points.keys():
        #    print s
        with open(output,'wt') as fp:
            for s in self.p:
                fp.write("%s\n"%s)

    def distance(self, loc1, loc2):
        "distance between two points"
        x1 = loc1.r * math.cos(loc1.theta)
        y1 = loc1.r * math.sin(loc1.theta)

        x2 = loc2.r * math.cos(loc2.theta)
        y2 = loc2.r * math.sin(loc2.theta)

        return math.sqrt(pow(y2 - y1,2) + pow(x2-x1,2))
        #return self.law_of_cosines_2(loc1.r, loc1.theta, loc2.r, loc2.theta)

    def law_of_cosines(self, r1, r2, gamma):
        "law of cosines given two side lengths and an included angle to find the third side length"
        x = pow(r1,2) + pow(r2,2) - 2 * r1 * r2 * math.cos(gamma)
        return math.sqrt(x)

    def law_of_cosines_2(self, r1, theta1, r2, theta2):
        "law of cosines given two sides separated by the difference theta1 - theta2"
        x = pow(r1,2) + pow(r2,2) - 2 * r1 * r2 * math.cos(delta_theta)
        return math.sqrt(x)

    def calculate_towers_required(self, tower):
        total_tower_points = sum([s.towers_required for s in self.p])
        possible_points_per_tower = tower.possible_point_ranges.keys()
        max_points = max(possible_points_per_tower)
        #print possible_points_per_tower
        items = []
        for s in possible_points_per_tower:
            #item name, value, size
            items.append(("%dpoints"%s,s,s))

        (solution, total) = knapsack01(items,int(total_tower_points))
        points_left = total_tower_points - total
        print solution
        print "%d points left" % points_left
        print "%d towers used" % len(solution)
        
        

class tower:
    def __init__(self, user_capacity=120, radius=20):
        "defines a tower object that can accomodate up to u users within a radius r"
        self.user_capacity = user_capacity
        self.range = radius
        self.counting_radius = self.hex_counting_ceil()
        self.possible_point_ranges = {}

    def calculate_grid_requirements(self, grid):
        """given a grid object where all points have been assigned a maximum population, finds the number of towers
        necessary to provide service for all users at each point, saving the result in the point objects"""
        for s in grid.p:
            s.towers_required = math.ceil(s.max_population / self.user_capacity)

    def hex_counting_ceil(self):
        "uses the tower's range to find the farthest another point can be and still be serviced by this tower"
        return math.ceil(2 * self.range * math.sqrt(3))

    def points_in_range(self, grid, loc):
        "given a grid object and a point on that grid, calculates how many points could be provided service by placing a tower at the given point"
        inrange = [s for s in grid.p if grid.distance(loc, s) <= self.range]

        """print "from (", loc.u, ",", loc.v, ")"
        for s in inrange:
            print "(", s.u, ",", s.v, ")"
        print"""
        
        n = len(inrange)
        self.possible_point_ranges.setdefault(n,0)
        self.possible_point_ranges[n] += 1
        return n

    def assign_ranges(self, grid):
        "runs points_in_range() on every point in a given grid"
        for s in grid.p:
            s.points_in_tower_range = self.points_in_range(grid, s)


def knapsack01(items, capacity):
    """items is a list of format [(item_name, value, size)]

    returns ([knapsack], total_value)"""
    num_items = len(items)
    cost = [0] * (capacity+1)
    best = [0] * (capacity+1)
    knapsack = []
    for j in range(1, num_items+1):
        for i in range(1, capacity+1):
            (name, value, size) = items[j-1]
            if i >= size:
                if cost[i] < cost[i - size] + value:
                    cost[i] = cost[i - size] + value
                    best[i] = j-1

    k = 0
    while capacity - k > 0:
        #print capacity - k, "capacity left: ",knapsack
        (name, value, size) = items[best[capacity - k]]
        knapsack.append((name,value,size))
        k += size
                    

    return (knapsack, cost[capacity])

def test(radius=40, tower_radius=20, pop=10000, output="output.txt"):
    g = grid(radius=40, gridsize=tower_radius/2-1, population=pop)
    a = sqm(len(g.p))
    g.use_transition_matrix(a)

    """
    total_residential = 0
    total_commuter = 0
    total_max = 0
    for i in range(0, len(g.p)):
        print g.p[i]
        print 'Residential Population: ', g.p[i].residential_population
        print 'Commuter Population: ', g.p[i].commuter_population
        print 'Max Population: ', g.p[i].max_population
        print

        total_residential += g.p[i].residential_population
        total_commuter += g.p[i].commuter_population
        total_max += g.p[i].max_population

    print 'Total Residential Population: ', total_residential
    print 'Total Commuter Population: ', total_commuter
    print 'Total Max Population: ', total_max
    """
    t = tower(radius=tower_radius)
    t.calculate_grid_requirements(g)
    t.assign_ranges(g)
    g.list_points()
    g.list_points_to_file(output)
    g.calculate_towers_required(t)
