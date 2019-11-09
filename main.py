from control import *


if __name__ == "__main__":

    def earth_wave(t):
        #T = several seconds, A = 10nT
        assert isinstance(t , (np.ndarray, np.generic))
        base = np.asarray([[1], [1],[1]])
    
        return (base * np.sin(2*np.pi / 5 * t)).transpose()
        
    def theta_white(length, beta_disturb):
        return np.random.normal(scale = beta_disturb , size = length)
        
    def white(t, k):
        """
        return white noise with amplifier k in time range t
        """
    
        return k * np.random.standard_normal(size=len(t))

    def position_to_temp(T0, locs):
        return T0 + 27/1000 * locs.transpose()[2]

    def temp_volt(T):
        return 0.1*(T-25)
    
    N_th = 5
    earth_stable = np.asarray([29364, 3013, 44520])
    rescue_loc = np.asarray([[-5,5,0]])
    target_loc = np.asarray([[0,0,0]])
    
    rescue_theta, rescue_alpha = utils.vecToAngle(np.array([[-80, 0, 40]]))
    
    #rescue_alpha = R(5.5)
    target_alpha = R(0.)
    #rescue_theta = R(78.)
    target_theta = R(0.)
    #与z轴夹角，地下为正，那么东为y
    #转弯半径，弧长100m，转弯3.5°，此处计算也要写
    #看每次走的距离有多少，再决定能转多少度
    time_step_go = 5.
    velocity = .1
    length = time_step_go * velocity 
    turn_angle = length * R(6.)
    
    magneticfield = MagneticField(target_loc, target_theta, target_alpha, earth_stable, earth_wave)
    sensor0 = Sensor(1.,1.,1., magneticfield,kw =1e-9, u0 = temp_volt, w = white)  #都设成1 都行
    sensor1 = Sensor(1.,1.,1., magneticfield,kw =1e-9, u0 = temp_volt, w = white)
    """
    #系数会非常影响两个传感器求查
    """
    control = Control(rescue_theta, #theta0
        rescue_alpha, #alpha0
        sensor0,
        sensor1,
        0.75, #radius
        velocity, #velocity
        rescue_loc, #current_ctr_loc
        np.pi *2, #  T = 10s #beta0
        0.0, #beta_disturb
        #sampling_time = 30000, #300Hz，测100s
        20, #temperature
        R(0), #sensor0_theta
        R(180), #sensor1_theta
        5., #sensor_distance
        0.1, #ada_p = 
        0.1, #ada_k = 
        R(5), #idea_approaching_angle = 
        utils.enBB*100, #5nT 一下就不要改方向了 #threshold = 
        time_step_go, #time_step_go = 
        10, #time_step_measure = 
        300,
        target_open = target_loc,
        target_direction = utils.angleToVec(target_theta, target_alpha),
        theta_white = theta_white,
        N_th = N_th,
        position_to_temp = position_to_temp,
        turn_angle = turn_angle
        ) #sampling_frequency = 
        
    #control.operate(0)
    #input()
    
    locs = []
    l = 13.7
    r = 14.3
    for i in range(int(l*100),int(r*100)):
        locs.append(np.array([[0, i/100, i/100]]))
    locs = np.asarray(locs) + target_loc
    """
    directions = []
    for i in range(1,20):#theta
        for j in range(1,20):#alpha
            directions.append(utils.angleToVec(R(i/20*360), R(j/20*5)).reshape(1,-1)) 
    """
    #control.experiment_1(locs, l,r)#, #directions)
    control.operate(0)
    #control.experiment_3((4,10), (R(2), R(5)), angles2=(0,0), nums=10)
    input()    

    #传感器的温漂
    
    """
    需要问的问题：
    救援井不知道实际位置
    ①传感器温漂:可线性
    ②转向的函数设计、现有方法可行？Adaptive参数怎么设计？：先设0.5的k
    ③传感器、测角位移的白噪声模拟
    ④滤波方案
    ⑤ 两个传感器的值如果减掉了之后，地磁场可以完全减掉；
    
    地磁场52000nT
    其他：
    ④坐标系检查
    
    控制策略：
    首先磁场是向下的；
    这样就能得到那幅图中的：目标井平行矢量；
    其次beta通过对侧相减，能够得到一个是椭圆形的， 值，那么最大的磁场差值就应该是轴外连线。这样就能够得到那个图中的目标井方向
    这样就能够得到共面的那个面；
    然后理想的靠近方向应该是一个有一定偏转角度的方向，可以设成5°
    然后这个方向在工具面上有个投影；在Zr应该调整到和理想方向一致的方向上。理想方向其实都可以直接得到；
    然后把Zr 即当前的Direction要调整到理想方向上。可以用那个投影来消。也可以直接并过去。
    直接并过去好不好就再说了。
    另外目标井的磁场方向错了。
    
    因此两个传感器有上下差都无所谓。如果是同侧分布的话，应当要把延后180°的拿来减。
    """
    
    
    
    """
    ①二阶系统相频特性，看不同转速的滞后；
    """
    """
    time_step_go也得调
    """
