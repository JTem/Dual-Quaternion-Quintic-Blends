import math
import numpy as np
from neura_dual_quaternions import DualQuaternion, Quaternion

class Line:
        def __init__(self, dq1, dq2, velocity, angular_velocity_max):
                self.p1 = dq1.getPosition().flatten()
                self.p2 = dq2.getPosition().flatten()
                self.q1 = dq1.real
                self.q2 = dq2.real
                
                # compute difference between dual quaternions
                delta_dq = dq1.inverse()*dq2
                
                # compute translational distance between the two key poses
                self.dist = np.linalg.norm(delta_dq.getPosition())
                
                # compute tangent for state functions
                self.tangent = np.array([0,0,0])
                if self.dist > 1e-6:
                        self.tangent = (dq2.getPosition() - dq1.getPosition()).flatten()/self.dist
                
                # compute duration of translational motion
                self.duration = self.dist/velocity
                
                # compute anlge between key orientations
                self.theta = abs(delta_dq.real.getAngle())
                
                # compute and normalize rotation axis
                e = delta_dq.real.log()
                if e.norm() > 1e-6:
                        e = e.normalize()

                e0 = self.q1*e*self.q1.inverse()
                self.rotation_axis = e0.getVector().flatten()
                
                # update duration of angular velocity exceeds angular_velocity_max
                if self.duration > 1e-6:
                
                        if self.theta/self.duration > angular_velocity_max:
                                self.duration = self.theta / angular_velocity_max
                                
                else:
                        self.duration = self.theta / angular_velocity_max
        
        # Calculating position at a given time on the line segement
        def getPosition(self, time):

                translation = self.p1 + self.getVelocity(0)*time

                return translation    
        
        # Calculating angular velocity of the line movement
        def getAngularVelocity(self):
                return self.rotation_axis*self.theta/self.duration
        
        # Calculating velocity at a given time on the line segement
        def getVelocity(self, time):
                return self.tangent*self.dist/self.duration
    
        # Calculating acceleration at a given time on the line segement
        def getAcceleration(self, time):
                return np.array([0,0,0])