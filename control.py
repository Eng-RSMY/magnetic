from sensor import *
from scipy import signal
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from utils import D

class Control:
    
    def __init__(self,
        theta0,
        alpha0,
        sensor0 ,
        sensor1 ,
        radius ,
        velocity ,
        current_ctr_loc,
        beta0, #  T = 10s
        beta_disturb,
        #sampling_time = 30000, #300Hz，测100s
        temperature,
        sensor0_theta,
        sensor1_theta,
        sensor_distance,
        ada_p,
        ada_k,
        ideal_approaching_angle ,
        threshold, #5nT 一下就不要改方向了
        time_step_go,
        time_step_measure,
        sampling_frequency,
        target_open,
        target_direction,
        theta_white,
        N_th,
        position_to_temp,
        turn_angle
        ):
        
    
        self.theta = theta0
        self.alpha = alpha0
        self.direction = utils.angleToVec(theta0, alpha0).reshape(1,-1)
        self.sensor0 = sensor0
        self.sensor1 = sensor1
        self.radius = radius
        self.velocity = velocity
        self.current_ctr_loc = current_ctr_loc
        self.current_time = 0
        self.beta = beta0 # rad/s
        self.beta_disturb = beta_disturb
        self.sampling_frequency = sampling_frequency
        self.st = sampling_frequency * time_step_measure
        self.temperature = temperature
        self.sensor0_theta = sensor0_theta
        self.sensor1_theta = sensor1_theta
        #self.direction = utils.angleToVec(theta0, alpha0)
        self.sensor_distance = sensor_distance
        self.ada_p = ada_p
        self.ada_k = ada_k
        self.ideal_direction_coefficient = np.sin(ideal_approaching_angle)
        self.threshold = threshold #5nT 一下就不要改方向了
        self.sensor0_loc = self.current_ctr_loc
        self.sensor1_loc = self.current_ctr_loc - self.sensor_distance 
        self.time_step_go = time_step_go
        self.time_step_measure = time_step_measure
        self.target_open = target_open
        self.target_direction =target_direction
        self.theta_white = theta_white
        self.N_th = N_th
        self.position_to_temp = position_to_temp
        self.turn_angle = turn_angle
    
    def __get_locs(self, angles, start_angle, M):
        
        x = np.cos(angles+ start_angle)
        y = np.sin(angles+ start_angle)
        z = np.zeros(len(angles))
        center_to_edge = np.vstack((x,y,z)).transpose()
        real_cte = np.matmul(center_to_edge, M)
        return real_cte

    def __sampling(self, time):
        """
        generate a sequence of location and time.
        get from sensor and merge
        starting from g-high, z=0
        从重力高边开始测；但是传感器的位置可以不在重力高边；对于绝对时间不敏感，只看相对时间
        
        一个采样周期 300-600Hz
        转几圈？0.1-0.2Hz 测量时间一分钟 比如5s一转 转12圈
        截止频率低，要较长时间才能到常数；
        取后面的就行；
        四阶滤波器；
        要把2倍频率滤掉，要转速1/5一下 二倍频的1/10
        """
        #print(np.linalg.norm(self.direction))
        assert np.abs(np.linalg.norm(self.direction, axis=1) - 1) < 1e-7
        angles = self.beta / self.sampling_frequency * (time - time[0]) #real angles #should be rads/(1/300 s)
      
        #计算相对于中心点的向量，因为0和1 号的角位置有差异，所以向量也不同。
        self.g_prll = np.cross(Magnetivity(np.reshape(self.direction, (1,3))).get_project(utils.z_axis).actual, utils.z_axis)
        self.g_high = np.cross(self.direction, self.g_prll)
        M = np.vstack((self.g_high, self.g_prll, self.direction)) * self.radius
      
        rc0 = self.__get_locs(angles, self.sensor0_theta, M)
        rc1 = self.__get_locs(angles, self.sensor1_theta, M)
        
        disturbed_angles =  angles + self.theta_white(len(angles), self.beta_disturb)
        
        #算出来不对
        
        mags0 = self.sensor0[rc0 + self.sensor0_loc, time, self.temperature, utils.x_axis, utils.y_axis, utils.z_axis]
        mags1 = self.sensor1[rc1 + self.sensor1_loc, time, self.temperature, utils.x_axis, utils.y_axis, utils.z_axis]
        """
        plt.figure(3)
        plt.ion()
        plt.clf()
        plt.plot((self.sensor0_loc + rc0).transpose()[0], 'b')
        plt.plot((self.sensor1_loc + rc1).transpose()[0], 'g')
        plt.draw()
        plt.pause(0.01)
        """
        
        """
        fig = plt.figure(3)
        plt.ion()
        fig.clf()
        ax = fig.add_subplot(111,projection='3d')
        ax.scatter((rc0 + self.sensor0_loc).transpose()[0],(rc0 + self.sensor0_loc).transpose()[1],(rc0 + self.sensor0_loc).transpose()[2], c = 'r')
        ax.scatter((rc1 + self.sensor1_loc).transpose()[0],(rc1 + self.sensor1_loc).transpose()[1],(rc1 + self.sensor1_loc).transpose()[2], c = 'b')
        """
        #print(self.sensor0_loc)
        #print(self.sensor1_loc)
        #input()
        #plt.draw()
        #plt.pause(0.01)
        #theta0, theta0_index = self.__merge(da0, mags0)
        #theta1, theta1_index = self.__merge(da1, mags1)
        
        sub_mag = mags0 - mags1 #这个地方只能用r来减 减掉地磁场，但是有一些反向了
        sub_mag.direction_align(self.target_direction) #调转成绝对值
        
        largest_beta = self.__merge(angles, sub_mag.r) # already sutracted phase_difference
        
        rc = self.__get_locs(np.asarray([largest_beta]), 0, M)
        #print(rc)
        rc = Magnetivity(rc + self.direction * 1/2*self.sensor_distance)
        #rc不应该往xy平面上投影！应该往目标井的法平面投影！
        
        #target_direction = mags0.direction.mean(axis = 0)#这个地方疏忽大意了？？？
        target_direction = sub_mag.direction.mean(axis = 0)
        rc = rc.get_project(self.target_direction).direction
        """
        plt.figure(2)
        plt.ion()
        plt.clf()
        plt.plot(mags0.r, 'b')
        plt.plot(mags1.r, 'g')
        plt.draw()
        plt.pause(0.01)
        """
        power = np.sqrt(np.sum(sub_mag.r**2)*1/self.sampling_frequency / self.time_step_measure) #这是均值
        #utils.paint(sub_mag.r,6)
        #print(len(sub_mag.r))
        #input()
        return target_direction, rc, max(sub_mag.r)-sub_mag.r.mean(), power  #实际是以rc0为准， 因为是以rc0为重力高边0度

    def __merge(self, thetas, mag):
        """
        merge to get exact phase difference
        滤波
        get full period
        """
        mags = mag - mag.mean()
        T = 2*np.pi / self.beta
        w0 = self.beta * self.st
        #phase_difference = phase_df(mags, thetas, wn, 2)
        wn = 0.1/T / self.sampling_frequency
        fai,fai2 = utils.phase_df(mags, thetas-thetas[0], wn)
        """
        plt.figure(4)
        plt.title('Phase Difference of Waves')
        plt.ion()
        plt.clf()
        plt.plot((np.sin(thetas)*np.std(mags))[:500], 'b')
        plt.plot(mags[:500], 'g')
        plt.xlabel('Time')
        plt.ylabel('Magnitude')
        plt.draw()
        plt.pause(0.01)
        """
        #严格来讲应该是没有延后的吧，你都不知道本身应该是怎样的波形，又不是跟theta一一对应的
        #to the extrem: this get directly the mean
        #find the largests sin(wn), Fs= 300Hz
        
        #return the largest tool beta
        largest_beta_index = np.where(mag == np.max(mag))
        betas = thetas[largest_beta_index]# - phase_difference
        while betas[betas>2*np.pi].any():
            betas[betas>2*np.pi] -= 2*np.pi
        #得到的应该是一个序列
        largest_beta = betas.mean()
        calculated_beta = np.pi/2 - fai+thetas[0]
        c2beta = np.pi/2 - fai2 + thetas[0]
        #print(largest_beta_index)
        #print("Betas: ",D(calculated_beta), ' ', D(largest_beta),' ', D(c2beta))
        return largest_beta
        
        
    def __turn(self, args):
        """
        adjust theta, alpha, calculate new location, and g-high
        是否使用AdaptiveGradient还需要考虑
        暂时先使用两者的差值向量来做平面修正，用两者的和向量做距离修正
        ？：所以不知道的是井下的具体theta和alpha值？也不知道direction？还是不知道目标井的direction？
        """
        #print(args)
        target_direction, to_target_direction, whether_threshold, power = args
        #print("threshold: ", whether_threshold, ' ', self.threshold ,' ',power, "\n")
        if np.abs(whether_threshold) < self.threshold:
            print("threshold: ", whether_threshold,' st:', self.threshold)
            return
        target_direction /= np.linalg.norm(target_direction)
        #print(to_target_direction)
        #assert 1==0
        to_target_direction /= np.linalg.norm(to_target_direction)
        #print("to: ", to_target_direction)
        #print("target: ", target_direction)
        
        self.direction = to_target_direction * min(self.ideal_direction_coefficient, np.tan(self.turn_angle)) + target_direction 
        self.direction /= np.linalg.norm(self.direction, axis=1)
        
        
    def __refresh_thetaAndAlpha(self):
        #theta and alph
        self.theta, self.alpha = utils.vecToAngle(self.direction)
        #sensor loc
        
    def __refresh_locations(self):
        self.sensor0_loc = self.current_ctr_loc
        self.sensor1_loc = self.current_ctr_loc - self.sensor_distance * self.direction
        self.temperature = self.position_to_temp(self.temperature, self.current_ctr_loc)
        
    def operate(self, start_time):
        """
        automatically operate and print results ( distance to target well)
        only query outside condition
        需要sensor0在上面，否则可能越界
        
        还有个问题，如果测不到，是不是按照原计划走？应该要大于白噪声的nT才能走吧
        """  
        self.current_time=start_time
        self.current_ctr_loc = self.current_ctr_loc + self.direction * self.velocity * self.time_step_go
        self.current_time = start_time + self.time_step_go
        self.__refresh_locations()
        
        fig = plt.figure(0)
        plt.ion()
        fig.clf()
        ax = fig.add_subplot(111,projection='3d')
        ax.invert_xaxis() 
        #ax.invert_yaxis() 
        ax.invert_zaxis() 
        xx = np.arange(-10,10,10)
        yy = np.arange(-10,10,10)
        X, Y = np.meshgrid(xx, yy)
        #print(X,Y)
        Z = np.zeros(X.shape)
        ax.plot_surface(X,Y,Z,cmap='copper', alpha=0.5)
        #平面
        ax.xaxis.set_ticks_position('top') 
        ax.yaxis.set_ticks_position('top')  
        points = np.matrix(np.arange(0, 80, 1)).transpose()
        #print(points, self.target_direction)
        
        f3 = plt.figure(3)
        f3.clf()
        bx = f3.add_subplot(111)
        bx.scatter(self.target_open[0][0], self.target_open[0][1], c='g')
        
        xyzs = (points * self.target_direction).transpose()
        figure0 = ax.scatter(self.current_ctr_loc[0][0],self.current_ctr_loc[0][1],self.current_ctr_loc[0][2], c = 'r')
        
        figure1 = ax.scatter(xyzs[0], xyzs[1], xyzs[2], c = 'b')
        plt.draw()
        #print(self.direction)
        while utils.get_distance(self.current_ctr_loc, self.target_open, self.target_direction) > self.radius:
            #measure
            print("distance: ", utils.get_distance(self.current_ctr_loc, self.target_open, self.target_direction))
            sampling_points = np.arange(0, self.st, 1)
            self.__turn(self.__sampling(sampling_points))
            #go
            self.current_ctr_loc += self.direction * self.velocity * self.time_step_go
            self.current_time = start_time + self.time_step_go + self.time_step_measure
            self.__refresh_locations()
            self.__refresh_thetaAndAlpha()
            figure0 = ax.scatter(self.current_ctr_loc[0][0],self.current_ctr_loc[0][1],self.current_ctr_loc[0][2], c = 'r', alpha=0.7)
            bx.scatter(self.current_ctr_loc[0][0],self.current_ctr_loc[0][1],c='b',alpha=0.7)
            
            plt.draw()
            #print(self.direction)
            plt.pause(0.01)
            
            
            print("CL: ",self.current_ctr_loc , "CT: ",D(self.theta), "CA: ", D(self.alpha))
            
            
    def __sample_exp(self, time):
        assert np.abs(np.linalg.norm(self.direction, axis=1) - 1) < 1e-7
        angles = self.beta / self.sampling_frequency * (time - time[0]) #real angles #should be rads/(1/300 s)
      
        #计算相对于中心点的向量，因为0和1 号的角位置有差异，所以向量也不同。
        self.g_prll = np.cross(Magnetivity(np.reshape(self.direction, (1,3))).get_project(utils.z_axis).actual, utils.z_axis)
        self.g_high = np.cross(self.direction, self.g_prll)
        M = np.vstack((self.g_high, self.g_prll, self.direction)) * self.radius
      
        rc0 = self.__get_locs(angles, self.sensor0_theta, M)
        rc1 = self.__get_locs(angles, self.sensor1_theta, M)
        
        disturbed_angles =  angles + self.theta_white(len(angles), self.beta_disturb)
        
        #算出来不对
        
        mags0 = self.sensor0[rc0 + self.sensor0_loc, time, self.temperature, utils.x_axis, utils.y_axis, utils.z_axis]
        mags1 = self.sensor1[rc1 + self.sensor1_loc, time, self.temperature, utils.x_axis, utils.y_axis, utils.z_axis]
        
        sub_mag = mags0 - mags1 #这个地方只能用r来减 减掉地磁场，但是有一些反向了
        #sub_mag.direction_align(self.target_direction)    
        return disturbed_angles, sub_mag.r        
              
    def experiment_1(self, locs, l, r):
        #10 nT/sqrt(Hz) @1Hz
        noise = utils.power
        # 100nT^2功率的1Hz信号？
        direction = self.target_direction
        sigs = []
        diss = []
        self.__refresh_thetaAndAlpha()
        for loc in locs:
            self.current_ctr_loc = loc
            self.__refresh_locations()
            
            sampling_points = np.arange(0, self.st, 1)
            thetas, datas = self.__sample_exp(sampling_points)
            datass = datas - np.mean(datas)
            power = np.sqrt(np.sum(datass**2)*1/self.sampling_frequency / self.time_step_measure) #这是均值
            SNR = power/noise
            DIS = utils.get_distance(self.current_ctr_loc, self.target_open, self.target_direction)
            #print("Distance: ",utils.get_distance(self.current_ctr_loc, self.target_open, self.target_direction))
            #print("Power_Signal:", power)
            #print("Power_Noise:", noise)
            #print("SNR: ", SNR)
            if np.abs(SNR-1)>0:
                print("Distance: ",DIS)
                print("Power_Signal:", power)
                print("Power_Noise:", noise)
                print("SNR: ", SNR)
            sigs.append(SNR)
            diss.append(np.array(DIS))
            
        plt.figure(5)
        plt.plot(diss, sigs, 'b')
        plt.plot(np.ones(len(diss)), 'r')
        plt.xlim((diss[0], diss[-1]))
        plt.xticks(np.arange(diss[0], diss[-1], 0.1))
        plt.xlabel('Distance')
        plt.ylabel('SNR')
        plt.grid()
        plt.show()
            #self.experiment_2(datas, thetas)
        
    def experiment_2(self, input_wave, thetas, w0 = 20, N_th=2):
        #st = 100
        wn = w0/self.sampling_frequency
        b, a = signal.butter(N_th, wn, 'low')
        filt_wave = signal.filtfilt(b, a, input_wave)
        
        T = 2*np.pi / self.beta
        #phase_difference = phase_df(mags, thetas, wn, 2)
        wn = 0.2/T / self.sampling_frequency
        
        p_d = utils.phase_df(input_wave, thetas, wn, 5)
        print("Phase_Lag: ", utils.D(p_d))
        
    def experiment_3(self, distances, angles, angles2=(0,0),nums=100, start_time=0):
        eligibleD = []
        eligibleA = []
        if angles2 == (0,0):
            #平面
            pass
        diss = (distances[0]+ (distances[1]-distances[0]) * np.random.rand(nums, 3)) / np.sqrt(3) #长度是三倍
        plane_angles = np.random.rand(nums) * (angles[1]-angles[0]) + angles[0]
        base_point = self.target_open + 10. * self.target_direction
        
        k=0
        fig = plt.figure(1)
        plt.ion()
        fig.clf()
        ax = fig.add_subplot(111, projection='3d')
        ax.invert_xaxis() 
        #ax.invert_yaxis() 
        ax.invert_zaxis() 
        xx = np.arange(-10,10,10)
        yy = np.arange(-10,10,10)
        X, Y = np.meshgrid(xx, yy)
        #print(X,Y)
        Z = np.zeros(X.shape)
        ax.plot_surface(X,Y,Z,cmap='copper', alpha=0.5)
        #平面
        ax.xaxis.set_ticks_position('top') 
        ax.yaxis.set_ticks_position('top')  
        points = np.matrix(np.arange(0, 80, 1)).transpose()
        #print(points, self.target_direction)
             
        xyzs = (points * self.target_direction).transpose()
        figure0 = ax.scatter(self.current_ctr_loc[0][0], self.current_ctr_loc[0][1], self.current_ctr_loc[0][2], c = 'b', alpha = 0.6)
                
        figure1 = ax.scatter(xyzs[0], xyzs[1], xyzs[2], c = 'b')
        plt.draw()
        
            
        for dis,angle in zip(diss, plane_angles):
            print(k)
            k+=1
            point = base_point + dis
            self.current_ctr_loc = point
            self.direction = utils.get_direction(point, self.target_open, self.target_direction)*np.tan(angle) + self.target_direction
            self.direction /= np.linalg.norm(self.direction)
            self.__refresh_locations()
            self.__refresh_thetaAndAlpha()           
        
            Dd = utils.get_distance(self.current_ctr_loc, self.target_open, self.target_direction)
            approached = False
            
            
            for i in range(int(np.ceil(Dd)/(self.time_step_go * self.velocity)*10)):
                if Dd< self.radius:
                    approached = True
                    break
                #measure
                #print("distance: ", utils.get_distance(self.current_ctr_loc, self.target_open, self.target_direction))
                sampling_points = np.arange(0, self.st, 1)
                self.__turn(self.__sampling(sampling_points))
                #go
                self.current_ctr_loc += self.direction * self.velocity * self.time_step_go
                self.current_time = start_time + self.time_step_go + self.time_step_measure
                self.__refresh_locations()
                self.__refresh_thetaAndAlpha()
                
                figure0 = ax.scatter(self.current_ctr_loc[0][0],self.current_ctr_loc[0][1],self.current_ctr_loc[0][2], c = 'r', alpha=0.7)
                plt.draw()
                plt.pause(0.01)
                #print(self.direction)
            
                #print("CL: ",self.current_ctr_loc , "CT: ",D(self.theta), "CA: ", D(self.alpha))
            if approached:
                
                eligibleD.append(Dd)
                eligibleA.append(angle)
                
                
        
        print(eligibleD)
        print(eligibleA)
        """
        plt.figure(6)
        plt.scatter(eligibleD, eligibleA, 'b')
#        plt.xlim((8., 10.5))
#        plt.xticks(np.arange(8., 10.5, 0.1))
        plt.grid()
        plt.show()
        """
if __name__=="__main__":
    Control.operate(0)
