import math
import numpy as np
from neura_dual_quaternions import DualQuaternion, Quaternion

class Line:
        def __init__(self, dq1, dq2, velocity, angular_velocity_max):
                self.p1 = dq1.getPosition().flatten()
                self.p2 = dq2.getPosition().flatten()
                self.q1 = dq1.real
                self.q2 = dq2.real
                
                delta_dq = dq1.inverse()*dq2
                
                self.dist = np.linalg.norm(delta_dq.getPosition())
                self.tangent = np.array([0,0,0])
                if self.dist > 1e-6:
                        self.tangent = (dq2.getPosition() - dq1.getPosition()).flatten()/self.dist

                self.duration = self.dist/velocity
                
                self.theta = abs(delta_dq.real.getAngle())
                e = delta_dq.real.log()
                if e.norm() > 1e-6:
                        e = e.normalize()

                e0 = self.q1*e*self.q1.inverse()
                self.rotation_axis = e0.getVector().flatten()

                if self.duration > 1e-6:
                
                        if self.theta/self.duration > angular_velocity_max:
                                self.duration = self.theta / angular_velocity_max
                                
                else:
                        self.duration = self.theta / angular_velocity_max
        
        def getPosition(self, time):

                translation = self.p1 + self.getVelocity(0)*time

                return translation    
    
        def getAngularVelocity(self):
                return self.rotation_axis*self.theta/self.duration

        def getVelocity(self, time):
                return self.tangent*self.dist/self.duration
    
        def getAcceleration(self, time):
                return np.array([0,0,0])