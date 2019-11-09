import numpy as np
from numpy import deg2rad as R
from numpy import rad2deg as D
import scipy.signal as signal
import matplotlib.pyplot as plt

x_axis = np.asarray([1,0,0])
y_axis = np.asarray([0,1,0])
z_axis = np.asarray([0,0,1])

bw = 20
fH = bw
Kn = 1.22 # 二
BWn = bw*Kn
fL = 1
enBB = 6e-3*np.sqrt(BWn)

power = enBB

def paint(x,i):
    plt.figure(i)
    plt.ion()
    plt.clf()
    plt.plot(x, 'g')
    plt.draw()
    plt.pause(0.01)

def get_target_magnet(xx):
    # xx already an array of number
    #print(xx)
    to_return = 11500*np.exp(1-xx)# - 600
    to_return[to_return<0] = 0
    return to_return

def angleToVec(theta, alpha):
    #north as 0, perpendicular as 0
    # x is north

    return np.asarray([np.sin(alpha)*np.cos(theta), np.sin(alpha)*np.sin(theta), np.cos(alpha)])
    
def vecToAngle(Vec):
    #could be multiple
    vec = Vec.copy()
    assert vec.shape[1]==3
    v_lens = np.linalg.norm(vec, axis=1)
    z_lens = vec.dot(z_axis)
    alphas = np.arccos(z_lens/v_lens)
    
    vec =vec.transpose()    
    thetas = np.arctan(vec[1]/vec[0])
    thetas[vec[0]<0] += np.pi
    return thetas , alphas

    
def get_distance(point, point0, direction0):
    direction = direction0 / np.linalg.norm(direction0)
    return np.linalg.norm((point-point0) - (point-point0).dot(direction)*direction)
    
def get_direction(point, point0, direction0):
    direction = direction0 / np.linalg.norm(direction0)
    d = (point-point0) - (point-point0).dot(direction)*direction
    d /= np.linalg.norm(d)
    return d
    
def lpf2(x, Q, w0, A):
    
    b = np.array([A*w0**2])
    a = np.array([1, w0/Q, w0**2])
    return signal.lfilter(b, a, x)
    
def phase_df(m0, b0s, wn, N_th=3):   #其实就是跟0相位的比较，所以还是只需要一个波形
    #几阶滤波器够用？？
    assert m0.mean() < 1e-7
    bet1 = np.sin(b0s)     #角度不一定从0开始，假设是从0开始的
    bet2 = np.cos(b0s)
    m1 = m0 * bet1
    m2 = m0 * bet2
    #How many times in a round?
    # How frequency is defined in digital?
    """
    Take care here. 归一化截止频率注意算对没有
    """
    bt, at = signal.butter(N_th, wn, 'low')
    m10 = signal.filtfilt(bt, at, m1)
    m20 = signal.filtfilt(bt, at, m2)
    """
    plt.figure(3)
    plt.ion()
    plt.clf()
    plt.plot(m1, 'b')
    plt.plot(m2, 'g')
    plt.plot(m10, 'r')
    plt.plot(m20, 'magenta')
    
    plt.draw()
    plt.pause(0.01)
    """
    m100 = m1.mean()
    m200 = m2.mean()
    
    fai = np.arctan((m20/m10)[len(m10)//2:].mean()) #取后一半，响应时间可能较长
    fai2 = np.arctan(m200/m100)
    #print(D(fai))
    #print("2f: ", D(fai2))
    #fai 是波形相对于sin0偏左的角度
    return fai, fai2
    
    
    """
    假设是二阶低通滤波，然后可以手算
    阻尼比0.707
    得研究转速多快的时候，滞后多少；
    """
