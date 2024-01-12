import math
import numpy as np
from DualQuaternionQuinticBlends.Line import Line

class LineGenerator:
        def __init__(self):
                pass
           
        def generateSegments(self, DQ_list, velocity, angular_velocity_max):
                segment_list = []
                for i in range(len(DQ_list)-1):
                        segment_list.append(Line(DQ_list[i], DQ_list[i+1], velocity, angular_velocity_max))
                        
                return segment_list

    
