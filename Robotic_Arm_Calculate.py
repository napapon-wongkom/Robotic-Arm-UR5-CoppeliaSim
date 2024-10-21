#Code by Napapon Wongkom
import numpy as np
import math as m
import time
import random as rand
from matplotlib import pyplot as plt

d2r = np.pi / 180
r2d = 180 / np.pi

class Robotic_Arm:
    def __init__(self,DH_table):
        self.DH_table = np.array(DH_table)
        pi = np.pi
        self.d2r = pi/180
        self.r2d = 1/self.d2r
        self.Pe = [0, 0, 0, 1]
        self.Ende = [0, 0, 0.19163, 1]
        self.DH_size = len(self.DH_table)
        self.first = 0
        self.last = self.DH_size

    def degarr2radarr(self,deg):
        rad_arr = []
        for i in range(0,len(deg)):
            rad_arr.append(deg[i] * self.d2r)
        rad_arr = np.array(rad_arr)
        return rad_arr

    def DH2T(self,alpha,a,d,theta):
        theta = theta * self.d2r
        alpha = alpha * self.d2r
        Dh = np.array(  [[m.cos(theta), -m.sin(theta), 0, a],
                        [m.sin(theta)*m.cos(alpha), m.cos(theta)*m.cos(alpha), -m.sin(alpha), -m.sin(alpha)*d], 
                        [m.sin(theta)*m.sin(alpha), m.cos(theta)*m.sin(alpha), m.cos(alpha), m.cos(alpha)*d], 
                        [0, 0, 0, 1]])
        return Dh
    
    def Homogeneous(self,first,last):
        Tn = np.identity(4)
        for t in range(first,last):
            alpha = self.DH_table[t,0]
            a = self.DH_table[t,1]
            d = self.DH_table[t,2]
            theta = self.DH_table[t,3]
            T = self.DH2T(alpha,a,d,theta)
            Tn = np.dot(Tn,T)
        return Tn

    def rotation_matrix(self,first,last):
        T = self.Homogeneous(first,last)
        result = T[0:3,0:3]
        return result

    def Euler(self):
        R = self.rotation_matrix(self.first,self.last)
        alpha = np.arctan2(-R[1,2],R[2,2]) * self.r2d
        beta = np.arcsin(R[0,2]) * self.r2d
        gamma = np.arctan2(-R[0,1],R[0,0]) * self.r2d
        euler = np.array([alpha,beta,gamma])
        return euler

    def foward_kinematic(self):
        forward = np.dot(self.Homogeneous(self.first,self.last),self.Pe)
        forward = np.array(forward[:3]).round(6)
        return forward
    
    def Jacobian(self):
        J = []
        for i in range(1,self.DH_size + 1): #1to6
            r = np.dot(self.Homogeneous(i,self.DH_size),self.Ende) #iri = i_6T * 6r6
            R = self.rotation_matrix(0,i) #0iR
            r0 = np.dot(R,r[0:3]) #0ri = 0iR * iri
            k = [0,0,1] #iki
            k0 = np.dot(R,k)  #0ki = 0iR * iki
            J_cross = np.cross(k0,r0) #0ki x 0ri
            Jn = np.concatenate((J_cross,k0))
            J.append(Jn)
        J = np.array(J)
        J = np.transpose(J)
        return J
    
    def X(self):
        Zero = [0,0,0]
        Phi = self.Euler() * self.d2r
        X = np.concatenate((self.foward_kinematic(),Zero))
        return X

    def inverse(self,M):
        M = np.array(M)
        inverse = np.linalg.pinv(M)
        return inverse
    
    def inverse_kinematic(self,desire_pos,alpha,theta):
        desire_pos = np.array(desire_pos)
        Delta_X = np.subtract(desire_pos,self.X())
        invJ = self.inverse(self.Jacobian())
        Delta_theta = np.dot(invJ,Delta_X) #delta theta = inverse Jacobian * delta X
        theta = self.degarr2radarr(theta)
        Ntheta = theta + (alpha * Delta_theta)
        return Ntheta
    
    def MSE(self,desire_pos):
        desire_pos = np.array(desire_pos)
        MSE = np.mean((desire_pos - self.X()) ** 2)
        return MSE
    
    def debug(self):
        print("_____________________________________________")
        print("End point position:    {}".format(self.foward_kinematic().round(4)))
        print("End point orientation: {}".format(self.Euler().round(2)))
        print("_____________________________________________")
        print("jacobian of robot arm : ")
        print(self.Jacobian().round(6))
        print("---------------------------------------------")

class Cubic_polynomial:
    def __init__(self,u0,uf,v0,vf,tf,t):
        self.u0 = u0
        self.uf = uf
        self.v0 = v0
        self.vf = vf
        self.t = t
        self.tf = tf
    def position(self):
        pos = self.u0 + ((3 / self.tf**2) * (self.uf - self.u0) * (self.t ** 2)) - ((2 / self.tf**3) * (self.uf - self.u0) * (self.t ** 3))
        return pos
    
    def velocity(self):
        velocity = ((6 / self.tf ** 2) * (self.uf - self.u0) * self.t) - ((6 / self.tf ** 3) * (self.uf - self.u0) * (self.t ** 2))
        return velocity

    def acceleration(self):
        accelerate = ((6 / self.tf ** 2) * (self.uf - self.u0)) - ((12 / self.tf ** 3) * (self.uf - self.u0) * self.t)
        return accelerate
    
    def debug(self):
         print("time : {} s , Position : {} , Velocity : {} , Accelerate : {} ".format(self.t,self.position(),self.velocity(),self.acceleration()))

class Linear_interpolation:
    def __init__(self,u0,uf,tf,tb,t):
        self.u0 = u0
        self.uf = uf
        self.tf = tf
        self.tb = tb
        self.t = t
        self.vb = (self.u0 - self.uf) / (self.tb - self.tf)
        self.a = self.vb / self.tb

    def position(self):
        while self.t < self.tf:
            if (self.t <= self.tb):
                position = self.u0 + ((1/2) * self.a * (self.t ** 2))
                return position
            elif(self.t <= (self.tf - self.tb)):
                position = self.u0 + ((1/2) * self.a * (self.tb ** 2)) + (self.vb * (self.t - self.tb))
                return position
            elif(self.t <= self.tf):
                position = self.u0 + ((1/2) * self.a * (self.tb ** 2)) + (self.vb * (self.tf - (2 * self.tb))) + (self.vb * (self.t - (self.tf - self.tb))) - ((1/2) * self.a * ((self.t - (self.tf - self.tb)) ** 2))
                return position
    
    def velocity(self):
        while self.t < self.tf:
            if (self.t <= self.tb):
                velocity = self.a * self.t
                return velocity
            elif(self.t <= (self.tf - self.tb)):
                velocity = self.vb
                return velocity
            elif(self.t <= self.tf):
                velocity = self.vb - (self.a * (self.t - (self.tf - self.tb)))
                return velocity
            
    def acceleration(self):
        while self.t < self.tf:
            if (self.t <= self.tb):
                acceleration = self.a
                return acceleration
            elif(self.t <= (self.tf - self.tb)):
                acceleration = 0
                return acceleration
            elif(self.t <= self.tf):
                acceleration = -self.a
                return acceleration
            
    def debug(self):
        print("time : {} s , Position : {} , Velocity : {} , Accelerate : {}".format(self.t,self.position(),self.velocity(),self.acceleration()))

#============================================================================================================================================
def start2stop(th0,thf,deg,joint_n):
    for i in range(joint_n):
        joint_rand0 = rand.randrange(0,deg + 1)
        th0[i] = joint_rand0
        joint_randf = rand.randrange((-deg),0)
        thf[i] = joint_randf

def start2task(th0,d_xyz,deg,joint_n):
    for i in range(joint_n):
        joint_rand0 = rand.randrange(-deg,deg + 1)
        th0[i] = joint_rand0
    for j in range(2):
        xyz_rand = rand.randrange(45,55)
        d_xyz.append(xyz_rand * 0.01)
    d_xyz.append(0.0375)

def set_zero(th0,joint):
    for i in range(joint):
        th0[i] = 0

def normalize_angle(ceil,angle):
    # Normalize multiples of 360 to 0
    if (angle >= 0):
        n_angle = angle % 360
    elif(angle < 0):
        n_angle = angle % -360
    else:
        raise ValueError('Invalid value')
    print(n_angle)
    if n_angle == 0:
        return 0
    # Clamp angles greater than 90 to 90
    elif n_angle > ceil:
        return ceil
    # Clamp angles less than -90 to -90
    elif n_angle < -ceil:
        return -ceil
    # Return the angle if it's within the range (-90, 90)
    else:
        return n_angle

def random_thf(thf,joint):
    for i in range(joint):
        thf[i] = rand.randrange(10,60)

def delay(t):
    init_t = time.time()
    tset = 0
    while tset < t:
      tset = time.time() - init_t

#============================================================================================================================================

if __name__ == "__main__":

    mode = int(input("Please select 0 (Debug) or 1 (Realtime) :"))
    #_____________Set Up Parameter___________________________________________________________________
    t = 0
    tf = 20
    tb = 5
    theta = 0
    th = {}
    th0 = {}
    thf = {}
    d_xyz = [0,0,0]
    set_zero(th0,6)
    random_thf(thf,6)
    #________________________________________________________________________________________________
    if mode == 1:
        print('Find Last Theta')
        while True:
            #-----------------------------------------Change Your DH_Table---------------------------------------------------------
            # Joint Parameter
            UR5_DH_table = [[0,     0,          0.0892,     -90 + th[0]],
                        [90,     0,          0,     90 + th[1]],
                        [0,       0.4251,    0,        th[2]],
                        [0,       0.3922,    0.110,      -90 + th[3]],
                        [-90,       0,    0.0948,        th[4]],
                        [90,       0,    0.07495,        th[5]],
                        [0,       0,    0.19163,        180]
                        ]

            #--------------------------------------------------------initialize robotic arm--------------------------------------------
            UR5 = Robotic_Arm(UR5_DH_table)
            print(UR5.Jacobian())
            d_ori = np.array([0,0,0])
            d_pos = np.concatenate((np.array(d_xyz),d_ori * d2r))
            new_theta = UR5.inverse_kinematic(d_pos,0.0001,thf) * r2d
            thf = dict(enumerate(new_theta))
            error = UR5.MSE(d_pos)
            
            print('Error : {}'.format(error))
            if(error <= 0.001):
                break
        #_________________________________________________________________________________________________________________
        print(thf)
        print("Time : {}".format(tf))
        print("Start")
        t1 = time.time()
        while t < tf:
            for i in range(6):
                Cubic = Cubic_polynomial(th0.get(i),thf.get(i),0,0,tf,t)
                Blend = Linear_interpolation(th0.get(i),thf.get(i),tf,tb,t)
                th[i] = Blend.position()
            #print(th)
            t = time.time() - t1
        print("Stop")   
    #============================================================================================================================================
    elif mode == 0:
        #-----------------------------------------Change Your DH_Table---------------------------------------------------------
        UR5_DH_table = [[0,     0,          0.0892,     -90 +th0.get(0)],
                        [90,     0,          0,     90 + th0.get(1)],
                        [0,       0.4251,    0,        th0.get(2)],
                        [0,       0.3922,    0.110,      -90 + th0.get(3)],
                        [-90,       0,    0.0948,        th0.get(4)],
                        [90,       0,    0.07495+0.19163,        180 + th0.get(5)]
                        ]
        #--------------------------------------------------------initialize robotic arm--------------------------------------------
        UR5 = Robotic_Arm(UR5_DH_table)
        d_ori = UR5.Euler()
        #current_pos = UR5.foward_kinematic()
        #UR5.debug()
        d_pos = np.concatenate((np.array(d_xyz),d_ori * d2r))
        th = UR5.inverse_kinematic(d_pos,0.1,th0) * r2d
        print(circle(0,0,2))
        
        
    

    



    
    


