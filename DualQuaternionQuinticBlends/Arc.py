import math
import numpy as np
from neura_dual_quaternions import Quaternion

class Arc:
        
        # Initializing the arc with dual quaternions, velocities, center, radius, rotation matrix, and angle
        def __init__(self, dq1, dq2, velocity, angular_velocity_max, center, radius, rotation, angle):
     
                self.q1 = dq1.real
                self.q2 = dq2.real
                self.center = center
                self.radius = radius
                self.rotation = rotation
                self.angle = angle
                
                # Calculating the difference between the initial and final orientations
                delta_q = self.q1.inverse()*self.q2
                
                # calculate the total translational distance
                self.dist = angle * radius
                    
                # calculate the segment duration
                self.duration = self.dist/velocity
                
                # calculate the angle
                self.theta = abs(delta_q.getAngle())
                
                # compute and normalize the rotation axis
                e = delta_q.log()
                if e.norm() > 1e-6:
                        e = e.normalize()

                e0 = self.q1*e*self.q1.inverse()
                self.rotation_axis = e0.getVector().flatten()
                
                # adjust duration of angular velocity exceeds the angular_velocity_max
                if self.duration > 1e-6:
                
                        if self.theta/self.duration > angular_velocity_max:
                                self.duration = self.theta / angular_velocity_max
                                
                else:
                        self.duration = self.theta / angular_velocity_max
            
            
        # Calculating position on the arc at a given time
        def getPosition(self, time):
                point_on_arc = np.array([self.radius*np.cos(time/self.duration*self.angle),
                                         self.radius*np.sin(time/self.duration*self.angle),
                                         0])

                translation = self.center + self.rotation@point_on_arc

                return translation


        # Calculating angular velocity of the arc movement
        def getAngularVelocity(self):
                return self.rotation_axis*self.theta/self.duration


        # Calculating linear velocity at a given time on the arc
        def getVelocity(self, time):
                velocity = np.array([-self.radius * np.sin(time/self.duration * self.angle) * self.angle/self.duration,
                                     self.radius * np.cos(time/self.duration * self.angle) * self.angle/self.duration,
                                     0])
                
                rotated_velocity = self.rotation @ velocity
                return rotated_velocity
        
        # Calculating linear acceleration at a given time on the arc
        def getAcceleration(self, time):
                acceleration = np.array([-self.radius * np.cos(time/self.duration * self.angle) * (self.angle/self.duration)**2,
                                         -self.radius * np.sin(time/self.duration * self.angle) * (self.angle/self.duration)**2,
                                         0])
                
                rotated_acceleration = self.rotation @ acceleration
                return rotated_acceleration