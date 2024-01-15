import math
import numpy as np
from DualQuaternionQuinticBlends.Arc import Arc


class ArcGenerator:
        # Constructor for ArcGenerator class
        def __init__(self):
                pass
        
        # Generates arc segments from a list of Dual Quaternions
        def generateSegments(self, DQ_list, velocity, angular_velocity_max):
                number_of_arcs = len(DQ_list)-2
                self.arc_list = []
                is_last_arc = False
                for i in range(number_of_arcs):

                        if i == number_of_arcs - 1:
                                is_last_arc = True
                                
                        # Generating arcs and appending them to the arc list
                        self.generateArcs(DQ_list[i], DQ_list[i+1], DQ_list[i+2], is_last_arc, velocity, angular_velocity_max)
                
                # Return the list of generated arcs
                return self.arc_list
            
        # Calculates the center and radius of the circle defined by three points
        def calcCenterAndRadius(self, p1, p2, p3):
                
                v = p3-p2
                u = p3-p1
                t = p2-p1

                w = np.cross(t,u)
                
                # Calculate center
                center = p1 + (u*np.dot(t,t) * np.dot(u,v) - t*np.dot(u,u) * np.dot(t,v))*0.5/np.linalg.norm(w)**2
                
                # Calculate radius
                radius = np.linalg.norm((p2 - center))

                return center, radius
        
        # Compute a rotation matrix to align a vector from the center to the start point with another vector from the center to the midpoint
        def computeRotationMatrix(self, v_center_start, v_center_mid):
                z_ = np.cross(v_center_start, v_center_mid)
                z = z_/np.linalg.norm(z_)

                y_ = np.cross(z, v_center_start)

                y = y_/np.linalg.norm(y_)
                
                # Create a rotation matrix
                rotation = np.array([v_center_start, y,  z]).transpose()
                
                # Return the rotation matrix
                return rotation
        
        # Generates two arcs for each triple of Dual Quaternions
        def generateArcs(self, dq1, dq2, dq3, is_last_arc, velocity, angular_velocity_max):
                
                p1 = dq1.getPosition().flatten()
                p2 = dq2.getPosition().flatten()
                p3 = dq3.getPosition().flatten()
                
                center, radius = self.calcCenterAndRadius(p1, p2, p3)
                
                v_center_start_ = p1 - center
                v_center_mid_ = p2 - center
                v_center_end_ = p3 - center
                
                # Generates two arcs for each triple of Dual Quaternions
                v_center_start = v_center_start_/np.linalg.norm(v_center_start_)
                v_center_mid = v_center_mid_/np.linalg.norm(v_center_mid_)
                v_center_end = v_center_end_/np.linalg.norm(v_center_end_)

                rotation = self.computeRotationMatrix(v_center_start, v_center_mid)
  
                angle = self.calculateAngleBetweenVectors(v_center_start, v_center_mid)
                
                # Append the first arc to the list
                self.arc_list.append(Arc(dq1, dq2, velocity, angular_velocity_max, center, radius, rotation, angle))
                
                # Handle the last arc separately to account for all orientations (else last arc would be defined by 3 points instead of two)
                if is_last_arc:
                                
                        rotation = self.computeRotationMatrix(v_center_mid, v_center_end)
                        angle = self.calculateAngleBetweenVectors(v_center_mid, v_center_end)
                        
                        # Append the last arc to the list
                        self.arc_list.append(Arc(dq2, dq3, velocity, angular_velocity_max, center, radius, rotation, angle))
        
        # Calculate the angle between two vectors
        def calculateAngleBetweenVectors(self, v1, v2):
                return math.atan2(np.linalg.norm(np.cross(v1, v2)), np.dot(v1, v2))

        # Getter method to return the list of generated arcs
        def getArcList(self):
                return self.arc_list
    
            
    
