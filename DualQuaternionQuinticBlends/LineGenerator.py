import math
import numpy as np
from DualQuaternionQuinticBlends.Line import Line

class LineGenerator:
        
        # Constructor for LineGenerator class
        def __init__(self):
                pass
        
        # Generates line segments between pairs of Dual Quaternions
        def generateSegments(self, DQ_list, velocity, angular_velocity_max):
                segment_list = []
                for i in range(len(DQ_list)-1):
                        # For each pair of consecutive Dual Quaternions in the list
                        # Create a line segment between them and append it to the segment list
                        segment_list.append(Line(DQ_list[i], DQ_list[i+1], velocity, angular_velocity_max))
                
                # Return the list of generated line segments
                return segment_list

    
