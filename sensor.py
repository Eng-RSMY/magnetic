from magnet import *

class Sensor:
    def __init__(self, kx, ky, kz, magneticfield, kw ,u0 , w ):
        """
        A sensor's properties include u = kB(loc,t) + u0(T) + w
        self.w is the condition of white noise
        """
        self.kx= kx
        self.ky= ky
        self.kz= kz
        self.u0 = u0
        self.w = w
        self.kw = kw
        self.magneticfield = magneticfield
    
    def __getitem__(self, key):
        """
        "t" can be a period of time
        and thus return a numpy array
        the array should be projected to three direction, thus three to be returned
        Questions: X，Y，Z轴分别是朝哪个方向定义的呢？
        温度和白噪声是三个都有还是只有一个呢？
        都不同；
        and returned with a series of theta, and three series of B
        t should be a range of time
        
        XYZ 可以
        
        T也是
        """
        loc, t, T, X, Y, Z = key
        mag = self.magneticfield[loc, t]
        value = mag.r
        direction = mag.direction
        
        x_mag = Magnetivity(self.kx * value + self.u0(T) + self.w(t, self.kw), direction).get_project(X)
        y_mag = Magnetivity(self.ky * value + self.u0(T) + self.w(t, self.kw), direction).get_project(Y)
        z_mag = Magnetivity(self.kz * value + self.u0(T) + self.w(t, self.kw), direction).get_project(Z)
        #magneticfield的输出不对
        # locs to beta is also disturbed, but is in the control module

        return x_mag+y_mag+z_mag
        
    #each has different ?
if __name__=="__main__":
    earth_stable = np.asarray([29364, 3013, 44520])
    rescue_loc = np.asarray([[100,100,0]])
    target_loc = np.asarray([[0,0,0]])
    rescue_alpha = R(5.5)
    target_alpha = R(5.)
    rescue_theta = R(78)
    target_theta = R(69) #与z轴夹角，地下为正，那么东为y
    
    magneticfield = MagneticField(target_loc, target_theta, target_alpha, earth_stable)
    sensor0 = Sensor(0.998, .999, .997, magneticfield)
    
    locs = np.random.randn(100,3)
    t = np.arange(0,100,1)
    print( sensor0[locs, t, 27, utils.x_axis, utils.y_axis, utils.z_axis].__repr__())
    
