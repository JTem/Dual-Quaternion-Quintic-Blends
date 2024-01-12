import math
import numpy as np
from neura_dual_quaternions import Quaternion

class Arc:
        def __init__(self, dq1, dq2, velocity, angular_velocity_max, center, radius, rotation, angle):
                #self.p1 = dq1.getPosition().flatten()
                #self.p2 = dq2.getPosition().flatten()
                self.q1 = dq1.real
                self.q2 = dq2.real
                self.center = center
                self.radius = radius
                self.rotation = rotation
                self.angle = angle
                
                delta_dq = dq1.inverse()*dq2
                
                self.dist = angle * radius
    
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
                point_on_arc = np.array([self.radius*np.cos(time/self.duration*self.angle),
                                         self.radius*np.sin(time/self.duration*self.angle),
                                         0])

                translation = self.center + self.rotation@point_on_arc

                return translation


        def getAngularVelocity(self):
                return self.rotation_axis*self.theta/self.duration


        def getVelocity(self, time):
                velocity = np.array([-self.radius * np.sin(time/self.duration * self.angle) * self.angle/self.duration,
                                     self.radius * np.cos(time/self.duration * self.angle) * self.angle/self.duration,
                                     0])
                
                rotated_velocity = self.rotation @ velocity
                return rotated_velocity

        def getAcceleration(self, time):
                acceleration = np.array([-self.radius * np.cos(time/self.duration * self.angle) * (self.angle/self.duration)**2,
                                         -self.radius * np.sin(time/self.duration * self.angle) * (self.angle/self.duration)**2,
                                         0])
                
                rotated_acceleration = self.rotation @ acceleration
                return rotated_acceleration