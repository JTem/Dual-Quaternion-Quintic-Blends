import numpy as np
import matplotlib.pyplot as plt
from neura_dual_quaternions import DualQuaternion, Quaternion

class DQQBTrajectoryGenerator:
        def __init__(self):
                self.min_blend_duration = 0.002
                
        def generateDynamicTrajectory(self, Segments, a_max, j_max, ang_acc_max, ang_jerk_max):
                self.Segments = Segments

                self.a_max = a_max if a_max >= 0.1 else 0.1
                self.j_max = j_max if j_max >= 1 else 1

                self.ang_acc_max = ang_acc_max if ang_acc_max >= 0.1 else 0.1
                self.ang_jerk_max = ang_jerk_max if ang_jerk_max >= 1 else 1

                for s in self.Segments:
                        if s.duration < self.min_blend_duration:
                                self.Segments.remove(s)

                self.num_segments = len(Segments)
                self.num_blends = len(Segments)+1
                self.duration_blend_list = []

                self.p0_list = []
                self.p1_list = []
                self.v0_list = []
                self.v1_list = []
                self.a0_list = []
                self.a1_list = []
                
                self.q0_list = []
                
                self.ang_v0_list = []
                self.ang_v1_list = []

                self.time_blend_start = []

                self.time_vector = []
                self.generateTrajectoryParameters()
                
                
        def generateTrajectoryParameters(self):
                
                for i in range(0, self.num_blends):
                        self.duration_blend_list.append(self.min_blend_duration)                
                
                
                blend_time = self.min_blend_duration*0.5
                
                #first blend:
                seg = self.Segments[0]
                self.p0_list.append(seg.getPosition(0))
                self.p1_list.append(seg.getPosition(blend_time))

                self.ang_v0_list.append(np.array([0,0,0]))
                self.ang_v1_list.append(seg.getAngularVelocity())

                self.v0_list.append(np.array([0,0,0]))
                self.v1_list.append(seg.getVelocity(blend_time))

                self.a0_list.append(np.array([0,0,0]))
                self.a1_list.append(seg.getAcceleration(blend_time))
                self.q0_list.append(seg.q1)
                
                for i in range(1, self.num_blends-1):
                        first_segment = self.Segments[i-1]
                        second_segment = self.Segments[i]

                        self.p0_list.append(first_segment.getPosition(first_segment.duration - blend_time))
                        self.p1_list.append(second_segment.getPosition(blend_time))

                        self.ang_v0_list.append(first_segment.getAngularVelocity())
                        self.ang_v1_list.append(second_segment.getAngularVelocity())

                        self.v0_list.append(first_segment.getVelocity(first_segment.duration - blend_time))
                        self.v1_list.append(second_segment.getVelocity(blend_time))

                        self.a0_list.append(first_segment.getAcceleration(first_segment.duration - blend_time))
                        self.a1_list.append(second_segment.getAcceleration(blend_time))
                        self.q0_list.append(second_segment.q1)
                        
                #last blend:
                seg = self.Segments[-1]
                self.p0_list.append(seg.getPosition(seg.duration - blend_time))
                self.p1_list.append(seg.getPosition(seg.duration))

                self.ang_v0_list.append(seg.getAngularVelocity())
                self.ang_v1_list.append(np.array([0,0,0]))

                self.v0_list.append(seg.getVelocity(seg.duration - blend_time))
                self.v1_list.append(np.array([0,0,0]))

                self.a0_list.append(seg.getAcceleration(seg.duration - blend_time))
                self.a1_list.append(np.array([0,0,0]))
                self.q0_list.append(seg.q2)
                
                
                cnt = 0
                blend_phase_overlap = True
                max_acc_limit_violated = True
                while (blend_phase_overlap or max_acc_limit_violated) or cnt < 10:
                        
                        cnt += 1
                        print("iteration: ", cnt)
                        
                        if cnt > 150:
                                break
                                
                        max_acc_limit_violated = False
                        for i in range(0, self.num_blends):
                                
                                max_acc, max_jerk = self.maxAccJerkQuintic(self.p0_list[i], self.p1_list[i], self.v0_list[i], self.v1_list[i], self.a0_list[i], self.a1_list[i], self.duration_blend_list[i]) 
                                max_angular_acc, max_angular_jerk = self.maxAccJerkCubic(self.ang_v0_list[i], self.ang_v1_list[i], np.array([0,0,0]), np.array([0,0,0]), self.duration_blend_list[i])
                                
                                a_lim = max([self.a_max, np.linalg.norm(self.a0_list[i]), np.linalg.norm(self.a1_list[i])])
                                
                                lamda_acc = max(max_acc/a_lim, max_angular_acc/self.ang_acc_max)
                                if max_angular_acc > self.ang_acc_max*1.01 or max_acc > a_lim*1.01:
                                        print("acceleration limit violated!")
                                        max_acc_limit_violated = True
                                                
                                self.duration_blend_list[i] *= lamda_acc

                                if self.duration_blend_list[i] < self.min_blend_duration:
                                        self.duration_blend_list[i] = self.min_blend_duration
                                        
                                half_T_blend = self.duration_blend_list[i]*0.5
                                if i == 0:
                                        seg = self.Segments[0]
                                        self.p1_list[i] = (seg.getPosition(half_T_blend))
                                        self.v1_list[i] = (seg.getVelocity(half_T_blend))
                                        self.a1_list[i] = (seg.getAcceleration(half_T_blend))
                                        
                                elif i == self.num_blends-1:
                                        seg = self.Segments[-1]
                                        self.p0_list[i] = (seg.getPosition(seg.duration - half_T_blend))
                                        self.v0_list[i] = (seg.getVelocity(seg.duration - half_T_blend))
                                        self.a0_list[i] = (seg.getAcceleration(seg.duration - half_T_blend))
           
                                else:
                                        seg1 = self.Segments[i-1]
                                        seg2 = self.Segments[i]
                                        self.p0_list[i] = (seg1.getPosition(seg1.duration - half_T_blend))
                                        self.p1_list[i] = (seg2.getPosition(half_T_blend))

                                        self.v0_list[i] = (seg1.getVelocity(seg1.duration - half_T_blend))
                                        self.v1_list[i] = (seg2.getVelocity(half_T_blend))

                                        self.a0_list[i] = (seg1.getAcceleration(seg1.duration - half_T_blend))
                                        self.a1_list[i] = (seg2.getAcceleration(half_T_blend))
                                        
                                        
                                        
                                max_acc, max_jerk = self.maxAccJerkQuintic(self.p0_list[i], self.p1_list[i], self.v0_list[i], self.v1_list[i], self.a0_list[i], self.a1_list[i], self.duration_blend_list[i]) 
                                max_angular_acc, max_angular_jerk = self.maxAccJerkCubic(self.ang_v0_list[i], self.ang_v1_list[i], np.array([0,0,0]), np.array([0,0,0]), self.duration_blend_list[i]) 
                                
                                lamda_jerk = max(np.sqrt(max_angular_jerk/self.ang_jerk_max), np.sqrt(max_jerk/self.j_max))
                                
                                if max_jerk > self.j_max or max_angular_jerk > self.ang_jerk_max:
                                        self.duration_blend_list[i] *= lamda_jerk

                                        if self.duration_blend_list[i] < self.min_blend_duration:
                                                self.duration_blend_list[i] = self.min_blend_duration
                                                
                                        half_T_blend = self.duration_blend_list[i]*0.5
                                        if i == 0:
                                                seg = self.Segments[0]
                                                self.p1_list[i] = (seg.getPosition(half_T_blend))
                                                self.v1_list[i] = (seg.getVelocity(half_T_blend))
                                                self.a1_list[i] = (seg.getAcceleration(half_T_blend))

                                        elif i == self.num_blends-1:
                                                seg = self.Segments[-1]
                                                self.p0_list[i] = (seg.getPosition(seg.duration - half_T_blend))
                                                self.v0_list[i] = (seg.getVelocity(seg.duration - half_T_blend))
                                                self.a0_list[i] = (seg.getAcceleration(seg.duration - half_T_blend))

                                        else:
                                                seg1 = self.Segments[i-1]
                                                seg2 = self.Segments[i]
                                                self.p0_list[i] = (seg1.getPosition(seg1.duration - half_T_blend))
                                                self.p1_list[i] = (seg2.getPosition(half_T_blend))

                                                self.v0_list[i] = (seg1.getVelocity(seg1.duration - half_T_blend))
                                                self.v1_list[i] = (seg2.getVelocity(half_T_blend))

                                                self.a0_list[i] = (seg1.getAcceleration(seg1.duration - half_T_blend))
                                                self.a1_list[i] = (seg2.getAcceleration(half_T_blend))
                                        
                                                
                        blend_phase_overlap = False
                        for i in range(0, self.num_segments):
                                if 0.5*(self.duration_blend_list[i] + self.duration_blend_list[i+1]) > self.Segments[i].duration:

                                        self.Segments[i].duration *= 1.05  
                                        
                                        print("blendphases overlapped")
                                        blend_phase_overlap = True

                                        seg = self.Segments[i]
                                        self.p1_list[i] = (seg.getPosition(self.duration_blend_list[i]*0.5))
                                        self.v1_list[i] = (seg.getVelocity(self.duration_blend_list[i]*0.5))
                                        self.a1_list[i] = (seg.getAcceleration(self.duration_blend_list[i]*0.5))

                                        self.p0_list[i+1] = (seg.getPosition(seg.duration - self.duration_blend_list[i+1]*0.5))
                                        self.v0_list[i+1] = (seg.getVelocity(seg.duration - self.duration_blend_list[i+1]*0.5))
                                        self.a0_list[i+1] = (seg.getAcceleration(seg.duration - self.duration_blend_list[i+1]*0.5))

                                        self.ang_v1_list[i] = (seg.getAngularVelocity())
                                        self.ang_v0_list[i+1] = (seg.getAngularVelocity())
                        

            
                self.time_blend_start.append(0)
                self.time_vector.append(0)
                
                for i in range(0,self.num_segments):
                        duration_sum = self.duration_blend_list[0]*0.5
                        for j in range(i+1):
                                duration_sum += self.Segments[j].duration
                                
                        self.time_vector.append(duration_sum)
                        self.time_blend_start.append(duration_sum - self.duration_blend_list[i+1]*0.5)               
                
                self.time_vector[-1] += self.duration_blend_list[-1]*0.5
                

        def evaluate(self, t):

                
                seg, cnt = self.determineSegmentType(t, self.time_blend_start, self.duration_blend_list)

                if seg == "blend":
                        dt = t - self.time_blend_start[cnt]
                        pos, vel, acc, jerk = self.quinticPolynomial(self.p0_list[cnt], self.p1_list[cnt], self.v0_list[cnt], self.v1_list[cnt], self.a0_list[cnt], self.a1_list[cnt], dt, self.duration_blend_list[cnt])
                        
                        log, angular_velocity, angular_acceleration, angular_jerk = self.cubicPolynomial(self.ang_v0_list[cnt], self.ang_v1_list[cnt], np.array([0,0,0]), np.array([0,0,0]), dt, self.duration_blend_list[cnt])

                        log -= self.ang_v0_list[cnt]*0.5*self.duration_blend_list[cnt]

                        quaternion = Quaternion.exp(Quaternion(0, *log)*0.5)*self.q0_list[cnt]
                
                if seg == "lin":
                        dt = t - self.time_blend_start[cnt] - 0.5*self.duration_blend_list[cnt]
                        jerk = np.array([0,0,0])
                        acc = self.Segments[cnt].getAcceleration(dt)
                        vel = self.Segments[cnt].getVelocity(dt)
                        pos = self.Segments[cnt].getPosition(dt)
                        
                        angular_jerk = np.array([0,0,0])
                        angular_acceleration = np.array([0,0,0])
                        angular_velocity = self.Segments[cnt].getAngularVelocity()

                        log = angular_velocity*dt

                        quaternion = Quaternion.exp(Quaternion(0, *log)*0.5)*self.q0_list[cnt]
                
                
                return pos, vel, acc, jerk, quaternion, angular_velocity, angular_acceleration, angular_jerk
        
        
        def evaluateDQ(self, t):
                seg, cnt = self.determineSegmentType(t, self.time_blend_start, self.duration_blend_list)

                if seg == "blend":
                        dt = t - self.time_blend_start[cnt]
                        pos, vel, acc, jerk = self.quinticPolynomial(self.p0_list[cnt], self.p1_list[cnt], self.v0_list[cnt], self.v1_list[cnt], self.a0_list[cnt], self.a1_list[cnt], dt, self.duration_blend_list[cnt])
                        
                        log, angular_velocity, angular_acceleration, angular_jerk = self.cubicPolynomial(self.ang_v0_list[cnt], self.ang_v1_list[cnt], np.array([0,0,0]), np.array([0,0,0]), dt, self.duration_blend_list[cnt])

                        log -= self.ang_v0_list[cnt]*0.5*self.duration_blend_list[cnt]

                        quaternion = Quaternion.exp(Quaternion(0, *log)*0.5)*self.q0_list[cnt]
                
                if seg == "lin":
                        dt = t - self.time_blend_start[cnt] - 0.5*self.duration_blend_list[cnt]
                        jerk = np.array([0,0,0])
                        acc = self.Segments[cnt].getAcceleration(dt)
                        vel = self.Segments[cnt].getVelocity(dt)
                        pos = self.Segments[cnt].getPosition(dt)
                        
                        angular_jerk = np.array([0,0,0])
                        angular_acceleration = np.array([0,0,0])
                        angular_velocity = self.Segments[cnt].getAngularVelocity()

                        log = angular_velocity*dt

                        quaternion = Quaternion.exp(Quaternion(0, *log)*0.5)*self.q0_list[cnt]
                
                # construct pure quaternions from the computed accelerations and velocities
                w = Quaternion(0, *angular_velocity)
                w_dot = Quaternion(0, *angular_acceleration)
                
                t = Quaternion(0, *pos)
                v = Quaternion(0, *vel)
                a = Quaternion(0, *acc)
                
                # compute and return the respective unit dual quaternions
                dq = DualQuaternion.fromQuatPos(quaternion, pos)
               
                dq_dot = DualQuaternion(0.5*(w*dq.real), 0.5*(v + 0.5*t*w)*dq.real)
                
                dq_ddot = DualQuaternion(0.5*(w_dot + 0.5*w*w)*dq.real, 0.5*(a + 0.5*(v*w + t*w_dot) + 0.5*(v + 0.5*t*w)*w)*dq.real)

                return dq, dq_dot, dq_ddot

        
        def quinticPolynomial(self, pos0, pos1, vel0, vel1, acc0, acc1, t, T):
                h = pos1 - pos0

                p0 = pos0
                p1 = vel0
                p2 = 0.5*acc0
                p3 = 1/(2*T**3)*(20*h - (8*vel1 + 12*vel0)*T - (3*acc0 - acc1)*T**2)
                p4 = 1/(2*T**4)*(-30*h + (14*vel1 + 16*vel0)*T + (3*acc0 - 2*acc1)*T**2)
                p5 = 1/(2*T**5)*(12*h - 6*(vel1 + vel0)*T + (acc1 - acc0)*T**2)


                pos = p0 + p1*t + p2*t**2 + p3*t**3 +p4*t**4 + p5*t**5
                vel = p1 + 2*p2*t + 3*p3*t**2 + 4*p4*t**3 + 5*p5*t**4
                acc = 2*p2 + 6*p3*t + 12*p4*t**2 + 20*p5*t**3
                jerk = 6*p3 + 24*p4*t + 60*p5*t**2

                return pos, vel, acc, jerk
    
        
        def cubicPolynomial(self, vel0, vel1, acc0, acc1, t, T):
                h = vel1 - vel0
                
                p0 = vel0
                p1 = acc0
                p2 = (3.0*h - (2.0*acc0 + acc1)*T)/T**2
                p3 = (-2.0*h + (acc0 + acc1)*T)/T**3
                
                pos = p0*t + 0.5*p1*t**2 + 1.0/3.0*p2*t**3 + 0.25*p3*t**4
                vel = p0 + p1*t + p2*t**2 + p3*t**3
                acc = p1 + 2*p2*t + 3*p3*t**2
                jerk = 2*p2 + 6*p3*t
                
                return pos, vel, acc, jerk
        
        def calcMinJerkTimeCubic(self, vel0, vel1, acc0, acc1, T):
                h = vel1 - vel0
                
                p2 = (3.0*h - (2.0*acc0 + acc1)*T)/T**2
                p3 = (-2.0*h + (acc0 + acc1)*T)/T**3

                time = 0
                if np.linalg.norm(p3) > 1e-6:
                        time = -np.dot(p2, p3)/(3*np.linalg.norm(p3)**2)

                return time
        
        def maxAccJerkQuintic(self, pos0, pos1, vel0, vel1, acc0, acc1, T):
        
                t = self.calcMinJerkTimeCubic(vel0, vel1, acc0, acc1, T)

                for i in range(10):
                        gradient = self.quinticJerkGradient(pos0, pos1, vel0, vel1, acc0, acc1, t, T)
                        t -= 0.02/(i+1)*T*np.sign(gradient)

                _, _, acceleration1, jerk1 = self.quinticPolynomial(pos0, pos1, vel0, vel1, acc0, acc1, 0, T)
                _, _, acceleration2, jerk2 = self.quinticPolynomial(pos0, pos1, vel0, vel1, acc0, acc1, T, T)
                _, _, acceleration3, _ = self.quinticPolynomial(pos0, pos1, vel0, vel1, acc0, acc1, t, T)

                jerk_max = max([np.linalg.norm(jerk1), np.linalg.norm(jerk2)])
                acc_max = max([np.linalg.norm(acceleration1), np.linalg.norm(acceleration2), np.linalg.norm(acceleration3)])

                return acc_max, jerk_max
        
        def quinticJerkGradient(self, pos0, pos1, vel0, vel1, acc0, acc1, t, T):
                eps = 1e-6
                t1 = t + eps
                t2 = t - eps

                _, _, _, jerk1 = self.quinticPolynomial(pos0, pos1, vel0, vel1, acc0, acc1, t1, T)
                _, _, _, jerk2 = self.quinticPolynomial(pos0, pos1, vel0, vel1, acc0, acc1, t2, T)
                gradient = (jerk1 - jerk2)/(2*eps)

                return gradient
        
        def maxAccJerkCubic(self, vel0, vel1, acc0, acc1, T):
                t = self.calcMinJerkTimeCubic(vel0, vel1, acc0, acc1, T)

                _, _, acceleration1, jerk1 = self.cubicPolynomial(vel0, vel1, acc0, acc1, 0, T)
                _, _, acceleration2, jerk2 = self.cubicPolynomial(vel0, vel1, acc0, acc1, T, T)
                _, _, acceleration3, _ = self.cubicPolynomial(vel0, vel1, acc0, acc1, t, T)

                jerk_max = max([np.linalg.norm(jerk1), np.linalg.norm(jerk2)])
                acc_max = max([np.linalg.norm(acceleration1), np.linalg.norm(acceleration2), np.linalg.norm(acceleration3)])

                return acc_max, jerk_max
        
        def determineSegmentType(self, t, time_blend_start, T_blend_list):
                # first find in which area we are
                cnt = 0
                for i in range(0, len(time_blend_start)):
                        if t < time_blend_start[i]:
                                cnt = i-1
                                break;
                
                if t < 0:
                        cnt = 0
                
                if t >= time_blend_start[-1]:
                        return "blend", len(time_blend_start)-1

                if t >= time_blend_start[cnt] and t < time_blend_start[cnt] + T_blend_list[cnt]:
                        return "blend", cnt

                else:
                        return "lin", cnt