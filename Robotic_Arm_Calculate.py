import numpy as np
import math as m

class Robotic_Arm:
    def __init__(self,DH_table):
        self.DH_table = np.array(DH_table)
        pi = np.pi
        self.d2r = pi/180
        self.r2d = 1/self.d2r
        self.Pe = [0, 0, 0, 1]
        self.DH_size = len(self.DH_table)
        self.first = 0
        self.last = self.DH_size

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
        euler = np.array([alpha,beta,gamma]).round(4)
        return euler

    def foward_kinematic(self):
        forward = np.dot(self.Homogeneous(self.first,self.last),self.Pe)
        forward = np.array(forward[:3]).round(6)
        return forward
    
    def Jacobian(self):
        J = []
        for i in range(1,self.DH_size + 1):
            r = np.dot(self.Homogeneous(i,self.DH_size),self.Pe)
            R = self.rotation_matrix(0,i)
            r0 = np.dot(R,r[0:3])
            k = [0,0,1]
            k0 = np.dot(R,k)
            J_cross = np.cross(k0,r0)
            Jn = np.concatenate((J_cross,k0))
            J.append(Jn)
        J = np.array(J)
        J = J.transpose()
        return J

#-----------------------------------------Change Your DH_Table---------------------------------------------------------
theta = 0 # Joint Parameter
joint_theta = [theta, theta, theta, theta, theta, theta] 
UR5_DH_table = [[0,     0,          0.0892,     -90 +joint_theta[0]],
            [90,     0,          0,     90 + joint_theta[1]],
            [0,       0.4251,    0,        joint_theta[2]],
            [0,       0.3922,    0.110,      -90 + joint_theta[3]],
            [-90,       0,    0.0948,        joint_theta[4]],
            [90,       0,    0.07495,        0],
            [0,       0,    0.19163,        180]
            ]

#--------------------------------------------------------initialize robotic arm--------------------------------------------
UR5 = Robotic_Arm(UR5_DH_table)
forward = UR5.foward_kinematic()
euler = UR5.Euler()
Jacobian = UR5.Jacobian()
print("_____________________________________________")
print("End point position:    {}".format(forward))
print("End point orientation: {}".format(euler))
print("_____________________________________________")
print("jacobian of robot arm : {}".format(Jacobian))
print("---------------------------------------------")




    
    


