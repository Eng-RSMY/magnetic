import utils
from utils import R
import numpy as np

class Magnetivity:

    # need to support matrix
    # So the best internal representation is numpy matrix... whereas the external to be r and direction

    def __init__(self, a ,b = None):
        #改了
        if b is not None:
            """
            r = abs;
            theta = theta, in the horizontal plane
            alpha = in vertical
        
            r is an array of abs and so is direction
            """
            self.r = a
            direction = b
            #direction 可能需要改变存储方式为[X,Y,Z]?
            #self.theta = theta
            #self.alpha = alpha
            #方向全部用向量表示
            #print(direction)
            self.direction =  (direction.transpose() / np.linalg.norm(direction, axis = 1))
            
            self.actual = (self.r * self.direction).transpose()
            self.direction = self.direction.transpose()
        else:
            actual = a
            self.actual = actual
            self.r = np.linalg.norm(actual, axis = 1)
            self.direction = (self.actual.transpose() / self.r).transpose()
        
        
    def __add__(self, other):
        
        assert isinstance(other, Magnetivity)
        assert self.actual.shape == other.actual.shape
        
        return Magnetivity(self.actual + other.actual)
    
    def __sub__(self, other):
        """
        self - other
        """
        assert isinstance(other, Magnetivity)
        assert self.actual.shape == other.actual.shape
        
        return Magnetivity(self.actual - other.actual)
        
    def get_angle(self):
        #以北为0度，x正方向,
        #改了
        return utils.vecToAngle(self.direction) #dont need to be transposed
        #horizon_direction = self.get_project(np.array([0,0,1])).actual.transpose()
        #thetas = np.degrees(np.arctan(horizon_direction[1]/horizon_direction[0]))
        #alphas = np.degrees(np.arccos(self.direction.dot(np.array([0,0,1]).transpose())))
        #return thetas,alphas
         
    def get_project(self, plane):
        """
        plane is expressed by a direction vector
        """
        #改了
        assert isinstance(plane, (np.ndarray, np.generic))

        #x_prime = self.direction.dot(plane)
        
        a_primes = self.actual.dot(plane.transpose()).reshape(-1, 1)
        y = np.linalg.norm(plane)
        #x_project_direction = self.direction - plane * x_prime/y
        a_project = self.actual - plane * a_primes /y
        #r = np.linalg.norm(x_project_direction) * self.r
        #x_project_direction = x_project_direction / np.linalg.norm(x_project_direction)
        return Magnetivity(a_project)

    def __setitem__(self, key, value):
        assert key in keys.keys()
        if key == "r":
            self.r = value # value is actually an array
            self.actual = (self.r * self.direction.transpose()).transpose()
        if key == "d":
            self.direction = value
            self.actual = (self.r * self.direction.transpose()).transpose()
        if key == "a":
            self.actual = value
            self.r = np.linalg.norm(value, axis =0)
            self.direction = (self.actual.transpose() / self.r).transpose()
       
    def __mul__(self, scaler):
        assert isinstance(scaler, (float, int))
        return Magnetivity(self.actual * scaler)
        
    def __rmul__(self, scaler):
        return self * scaler 
        
    def __repr__(self):
        return str(self.actual)
        
    def direction_align(self, template_direction):
        """
        主要目的是调转成一个方向，然后好比较大小
        """
        turns = self.direction.dot(template_direction.transpose())
        
        self.r[turns<0] *= -1
        self.direction[turns<0] *= -1
        

class MagneticField:

    #先搞个简单合成再说

    def __init__(self, target_open_loc , target_theta, target_alpha, earth_stable,  earth_wave):
        """
        This class is to return the magnetic intensity of any point and time
        """
        self.target_open_loc = target_open_loc
        self.target_direction = utils.angleToVec(target_theta, target_alpha)
        self.earth_stable = earth_stable
        self.earth_wave = earth_wave

    def earth_magnet(self, t):
        stable = self.earth_stable.reshape(1,3).repeat(len(t), axis=0) # t is a range
        wave = self.earth_wave(t)
        return Magnetivity(stable + wave)
        
    def drill_magnet(self, loc, density = utils.get_target_magnet):

        # should be revised to be used for a sequence
        #This place is not revised yet
        #Loc should already be a 2-d
        relative_lengths = (Magnetivity((self.target_open_loc.transpose() - loc.transpose()).transpose()).get_project(self.target_direction)).r
        #print("r ",relative_length.actual)
        #print("Here: ",self.target_direction)
        return Magnetivity(density(relative_lengths), self.target_direction.copy().reshape(1,-1).repeat(len(relative_lengths),axis=0))
       
    def __getitem__(self, key):
        """
        Calculate the magnetic value and direction at (x,y) and at time range t;
        Thus t could be a period of time
        """
        loc, t = key
        #loc should also be a numpy sequence? 2x X?
       
        assert loc.shape[0] == t.shape[0]

        #主要是因为前面的形状问题没有解决
        return  self.earth_magnet(t) + self.drill_magnet(loc)
    
    
    
if __name__=="__main__":
    earth_stable = np.asarray([29364, 3013, 44520])
    rescue_loc = np.asarray([[100,100,0]])
    target_loc = np.asarray([[0,0,0]])
    rescue_alpha = R(5.5)
    target_alpha = R(5.)
    rescue_theta = R(78)
    target_theta = R(69) #与z轴夹角，地下为正，那么东为y
    
    magneticfield = MagneticField(target_loc, target_theta, target_alpha, earth_stable)
    
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib.pyplot as plt
    
    fig = plt.figure(0)
    ax = fig.gca(projection='3d')
    ax.invert_xaxis() 
    #ax.invert_yaxis() 
    ax.invert_zaxis() 
    xx = np.arange(-100,100,10)
    yy = np.arange(-100,100,10)
    X, Y = np.meshgrid(xx, yy)
    print(X,Y)
    Z = np.zeros(X.shape)
    ax.plot_surface(X,Y,Z,cmap='copper', alpha=0.5)
    ax.xaxis.set_ticks_position('top') 
    ax.yaxis.set_ticks_position('top')  
    points = np.matrix(np.arange(0, 10, 0.1)).transpose()
    print(points, magneticfield.target_direction)
    xyzs = (points * magneticfield.target_direction).transpose()
    figure0 = ax.scatter(rescue_loc[0][0],rescue_loc[0][1],rescue_loc[0][2], c = 'r')
    figure1 = ax.scatter(xyzs[0], xyzs[1], xyzs[2], c = 'b')
    plt.show()
    
    loc = np.asarray([[1,2,3],[2,3,4]])
    t = np.asarray([1 ,2])
    print(magneticfield[loc, t])

