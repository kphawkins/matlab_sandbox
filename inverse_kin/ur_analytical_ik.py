import numpy as np
import copy
from numpy import sqrt, sin, cos, arctan2, arccos, arcsin, abs, pi
import roslib
roslib.load_manifest('ur_cart_move')
roslib.load_manifest('hrl_geom')
import rospy
from hrl_geom.pose_converter import PoseConv
from geometry_msgs.msg import PoseStamped
from ur_cart_move.ur_cart_move import RAVEKinematics
from hrl_geom.transformations import rotation_from_matrix as mat_to_ang_axis_point

d1, d2, d3 = 0., 0., 0.
a1, a2, a3 = 1., 0.5, 0.3

h_1 = 0.1273
h_3 = 0.0922
p_1 = 0.0
p_3 = 0.1157

#def inverse_rrr(T04, a=[a1, a2, a3], d=[d1, d2, d3]):
def inverse_rrr(p04x, p04y, x04x, x04y, a=[a1, a2, a3], d=[d1, d2, d3]):
    a1, a2, a3 = a
    d1, d2, d3 = d
    #p04x, p04y = T04[0,3], T04[1,3]
    #x04x = T04[0,0]
    #x04y = T04[1,0]
    qs1 = []
    p13x = p04x - a3*x04x
    p13y = p04y - a3*x04y
    c2 = (p13x**2 + p13y**2 - a1**2 - a2**2) / (2.*a1*a2)
    if abs(c2) > 1. + ZERO_THRESH:
        #print 'low c2'
        #print c2, p13x, p13y
        return []
    for ssign in [1., -1.]:
        q = [0.]*3
        s2 = ssign*np.sqrt(1. - c2**2)
        q[1] = np.arctan2(s2, c2)
        if q[1] < 0.:
            q[1] += 2.*pi
        denom = a1**2 + a2**2 + 2*a1*a2*c2
        c1 = ((a1 + a2*c2) * p13x + a2*s2*p13y) / denom
        s1 = ((a1 + a2*c2) * p13y - a2*s2*p13x) / denom
        q[0] = np.arctan2(s1, c1)
        c12 = np.cos(q[0]+q[1])
        s12 = np.sin(q[0]+q[1])
        c3 = x04x*c12 + x04y*s12
        s3 = c12*x04y - s12*x04x
        q[2] = np.arctan2(s3, c3)
        qs1.append(q)
    #qs2 = []
    #for q in qs1:
    #    qs2.append(q)
    #    if q[0] > 0:
    #        qs2.append([q[0]-2.*pi, q[1], q[2]])
    #    else:
    #        qs2.append([q[0]+2.*pi, q[1], q[2]])
    #qs3 = []
    #for q in qs2:
    #    qs3.append(q)
    #    if q[1] > 0:
    #        qs3.append([q[0], q[1]-2.*pi, q[2]])
    #    else:
    #        qs3.append([q[0], q[1]+2.*pi, q[2]])
    #qs4 = []
    #for q in qs3:
    #    qs4.append(q)
    #    if q[2] > 0:
    #        qs4.append([q[0], q[1], q[2]-2.*pi])
    #    else:
    #        qs4.append([q[0], q[1], q[2]+2.*pi])
    #return qs4
    return qs1

def forward_rrr(q):
    T01 = np.mat([[np.cos(q[0]), -np.sin(q[0]), 0., 0.],
                  [np.sin(q[0]), np.cos(q[0]), 0., 0.],
                  [0., 0., 1., d1],
                  [0., 0., 0., 1.]])
    T12 = np.mat([[np.cos(q[1]), -np.sin(q[1]), 0., a1],
                  [np.sin(q[1]), np.cos(q[1]), 0., 0.],
                  [0., 0., 1., d2],
                  [0., 0., 0., 1.]])
    T23 = np.mat([[np.cos(q[2]), -np.sin(q[2]), 0., a2],
                  [np.sin(q[2]), np.cos(q[2]), 0., 0.],
                  [0., 0., 1., d3],
                  [0., 0., 0., 1.]])
    T34 = np.mat([[1., 0., 0., a3],
                  [0., 1., 0., 0.],
                  [0., 0., 1., 0.],
                  [0., 0., 0., 1.]])
    return T01*T12*T23*T34

def inv_mat(T):
    R = np.mat(np.eye(4))
    R[:3,:3] = T[:3,:3].T
    R[:3,3] = -T[:3,:3].T * T[:3,3]
    return R

Tb0 = np.mat([[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
T6e = np.mat([[0, -1, 0, 0], [0, 0, -1, 0], [1, 0, 0, 0], [0, 0, 0, 1]])

def forward_kin(q, a, d, l):
    Ts = []
    Tlast = Tb0
    for qi, ai, di, li in zip(q, a, d, l):
        cqi, sqi = cos(qi), sin(qi)
        cli, sli = cos(li), sin(li)
        T = np.mat([[cqi, -sqi*cli,  sqi*sli, ai*cqi],
                    [sqi,  cqi*cli, -cqi*sli, ai*sqi],
                    [      0,          sli,          cli,         di],
                    [      0,                0,                0,          1]])
        Tlast = Tlast*T
        Ts.append(T)
    Tlast *= T6e
    return Tlast, Ts

ZERO_THRESH = 0.00001

def inverse_kin(T06, a, d, l, debug=False):
    qs1 = []
    T = inv_mat(Tb0) * T06 * inv_mat(T6e)
    A = d[5]*T[1,2] - T[1,3]
    B = d[5]*T[0,2] - T[0,3]
    R = A*A + B*B
    if abs(A) < ZERO_THRESH:
        print '1: A low'
        return []
    elif abs(B) < ZERO_THRESH:
        print '1: B low'
        return []
    elif d[3]*d[3] > R:
        #print '1: Impossible solution'
        return []
    else:
        for i in range(2):
            qs1.append([0.]*6)
        acos = arccos(d[3] / sqrt(R)) 
        atan = arctan2(B, A)
        pos = acos - atan
        neg = -acos - atan
        if pos >= 0.:
            qs1[0][0] = pos
        else:
            qs1[0][0] = 2.*pi + pos
        if neg >= 0.:
            qs1[1][0] = neg
        else:
            qs1[1][0] = 2.*pi + neg
        #if pos < 0:
        #    qs1[2][0] = pos + 2.*pi
        #else:
        #    qs1[2][0] = pos - 2.*pi
        #if neg < 0:
        #    qs1[3][0] = neg + 2.*pi
        #else:
        #    qs1[3][0] = neg - 2.*pi
    qs2 = []
    for i in range(len(qs1)):
        for j in range(2):
            qs2.append(copy.copy(qs1[i]))
        if debug:
            print 'h', T[0,2]*sin(qs1[i][0]) - T[1,2]*cos(qs1[i][0])
            print 'h2', T
            acos = arccos(T[0,2]*sin(qs1[i][0]) - T[1,2]*cos(qs1[i][0]))
            print 'h3', qs1[i][0]
            print 'h2', (T[0,3]*sin(qs1[i][0]) - T[1,3]*cos(qs1[i][0])-d[3])/d[5]
        acos = arccos((T[0,3]*sin(qs1[i][0]) - T[1,3]*cos(qs1[i][0])-d[3])/d[5])
        if acos >= 0.:
            qs2[i*2+0][4] = acos
            qs2[i*2+1][4] = 2.*pi-acos
        else:
            qs2[i*2+0][4] = -acos
            qs2[i*2+1][4] = 2.*pi+acos

    qs3 = []
    for i in range(len(qs2)):
        for j in range(2):
            qs3.append(copy.copy(qs2[i]))
        s4 = sin(qs2[i][4])
        #print 's4', s4
        #print 'h2', (T[0,0]*sin(qs2[i][0]) - T[1,0]*cos(qs2[i][0]))
        #c1, s1 = cos(qs2[i][0]), sin(qs2[i][0])
        #acos = arctan2(-(T[1,1]*c1-T[0,1]*s1), T[1,0]*c1-T[0,0]*s1)
        #acos = np.arctan(T[2,0]/ T[2,1])
        #print 'k', acos
        #print 'k2', acos
        #acos = ( (-1.)**(i%2+0)* np.sign(T[2,2])**2 *pi/2.
        #        +(-1.)**2* np.sign(T[2,2])**2 *arcsin(T[1,0]) 
        #        +(-1.)**2* np.sign(T[2,2])**2 *qs2[i][0])
        if abs(s4) < ZERO_THRESH:
            #print '6: s4 low'
            qs3[i][5] = 0.
            qs3[i+1][5] = pi
        elif abs(abs(s4) - 1.) < ZERO_THRESH:
            acos = (-1.)**(i%2) * pi/2. + arcsin(T[1,0]) + qs2[i][0]
            if acos >= 0.:
                if T[2,2] >= 0.:
                    qs3[i*2+0][5] = 2.*pi-acos
                    qs3[i*2+1][5] = 2.*pi-acos
                else:
                    qs3[i*2+0][5] = acos
                    qs3[i*2+1][5] = acos
            else:
                if T[2,2] >= 0.:
                    qs3[i*2+0][5] = -acos
                    qs3[i*2+1][5] = -acos
                else:
                    qs3[i*2+0][5] = 2.*pi+acos
                    qs3[i*2+1][5] = 2.*pi+acos
        else:
            acos = arccos((T[0,0]*sin(qs2[i][0]) - T[1,0]*cos(qs2[i][0])) / s4)
        #if abs(cos(acos-qs2[i][0])) - abs(T[0,0]) > ZERO_THRESH:
        #    acos += pi
        #if qs2[0][0] < pi and T[2,2] < 0.:
        #    acos -= pi
            if acos >= 0.:
                #if T[2,2] >= 0.:
                #    qs3[i*1+0][5] = 2.*pi-acos
                #else:
                #    qs3[i*1+0][5] = acos
                qs3[i*2+0][5] = 2.*pi-acos
                qs3[i*2+1][5] = acos
            else:
                #if T[2,2] >= 0.:
                #    qs3[i*1+0][5] = -acos
                #else:
                #    qs3[i*1+0][5] = 2.*pi+acos
                qs3[i*2+0][5] = -acos
                qs3[i*2+1][5] = 2.*pi+acos
        #print 'ssss', s4, qs3[i*1+0][5], qs3[i*1+1][5]
        
        #print '1111111111111', cos(qs3[i][5])*sin(qs3[i][0])*sin(qs3[i][4]) - cos(qs3[i][0])*sin(qs3[i][5]),  cos(qs3[i][5])*sin(qs3[i][0])*sin(qs3[i][4]) + cos(qs3[i][0])*sin(qs3[i][5]), T[0,0], T[1,0]
        for k in [0, 1]:
            if abs(abs(s4) - 1.) < ZERO_THRESH:
                tmp1 = cos(qs3[2*i+k][5])*sin(qs3[2*i+k][0])*sin(qs3[2*i+k][4])
                tmp2 = cos(qs3[2*i+k][0])*sin(qs3[2*i+k][5])
                #print sin(qs3[2*i+k][4])
                if abs(abs(tmp1 - tmp2) - abs(T[0,0])) < ZERO_THRESH:
                    if np.sign(tmp1 - tmp2) != np.sign(T[0,0]):
                        #qs3[2*i+k][5] -= pi
                        #qs3[2*i+k][0] *= -1
                        #qs3[i][5] *= -1
                        if sin(qs3[2*i+k][4]) > 0:
                            qs3[2*i+k][5] = -qs3[2*i+k][5] + 2*qs3[2*i+k][0]
                        else:
                            qs3[2*i+k][5] = -qs3[2*i+k][5] - 2*qs3[2*i+k][0]
                        #print tmp1 - tmp2
                        #print T[0,0]
                        #print 'yo1'
                else:
                    if np.sign(tmp1 + tmp2) != np.sign(T[0,0]):
                        #qs3[i][5] -= pi
                        #qs3[i][0] *= -1
                        #qs3[i][5] *= -1
                        if sin(qs3[2*i+k][4]) < 0:
                            qs3[2*i+k][5] = -qs3[2*i+k][5] + 2*qs3[2*i+k][0]
                        else:
                            qs3[2*i+k][5] = -qs3[2*i+k][5] - 2*qs3[2*i+k][0]
                        #print tmp1 + tmp2
                        #print T[0,0]
                        #print 'yo2'
                while qs3[2*i+k][5] < 0.:
                    qs3[2*i+k][5] += 2.*pi
                while qs3[2*i+k][5] > 2.*pi:
                    qs3[2*i+k][5] -= 2.*pi
        if debug:
            print 'yeh', qs3[i]

        if False:
            print 'wwwwwwwwwwwwwwww', sin(qs3[i][5]+qs3[i][0]), sin(qs3[i][5]-qs3[i][0]), T[0,0]
            print 'qqqqqqqqqqqqqqqq', cos(qs3[i][5]+qs3[i][0]), cos(qs3[i][5]-qs3[i][0]), T[0,1]
            flip_sign_sin, flip_sign_cos, flip_sub_sin, flip_sub_cos = False, False, False, False
            flip_diff = False
            if abs(abs(sin(qs3[i][5]+qs3[i][0])) - abs(T[0,0])) > ZERO_THRESH:
                qs3[i][5] -= 2*qs3[i][0]
                print 'a'
            print 'wwwwwwwwwwwwwwww', sin(qs3[i][5]+qs3[i][0]), sin(qs3[i][5]-qs3[i][0]), T[0,0]

            if abs(sin(qs3[i][5]+qs3[i][0]) - T[0,0]) > ZERO_THRESH:
                flip_sign_sin = True
            if abs(cos(qs3[i][5]+qs3[i][0]) - T[0,1]) > ZERO_THRESH:
                flip_sign_cos = True
            if flip_sign_sin:
                if flip_sign_cos:
                    qs3[i][5] += pi
                    print 'b'
                else:
                    qs3[i][5] = -qs3[i][5] 
                    #qs3[i][5] = -qs3[i][5] - 2*qs3[i][0]
                    qs3[i][0] = -qs3[i][0]
                    print 'c'
            elif flip_sign_cos:
                qs3[i][5] = pi -qs3[i][5]
                #qs3[i][5] = pi -qs3[i][5] - 2*qs3[i][0]
                qs3[i][0] = -qs3[i][0]
                print 'd'
            print 'e'

            print '3333333333333333', sin(qs3[i][5]+qs3[i][0]), sin(qs3[i][5]-qs3[i][0]), T[0,0]
            print '4444444444444444', cos(qs3[i][5]+qs3[i][0]), cos(qs3[i][5]-qs3[i][0]), T[0,1]
            #qs3[i][0] -= pi
            #qs3[i][0] -= 2*acos
        #if T[0,1] >= 0.:
        #    if -T[2,2] >= 0.:
        #        qs3[i][5] -= pi
        #qs3[i*4+2][5] = 2.*pi - acos
        #qs3[i*4+3][5] = -2.*pi + acos
    qs4 = []
    for i in range(len(qs3)):
        c1, s1 = cos(qs3[i][0]), sin(qs3[i][0])
        c5, s5 = cos(qs3[i][4]), sin(qs3[i][4])
        c6, s6 = cos(qs3[i][5]), sin(qs3[i][5])
        x04x = -s5*(T[0,2]*c1 + T[1,2]*s1) - c5*(s6*(T[0,1]*c1 + T[1,1]*s1) - c6*(T[0,0]*c1 + T[1,0]*s1))
        x04y = c5*(T[2,0]*c6 - T[2,1]*s6) - T[2,2]*s5
        p04x = d[4]*(s6*(T[0,0]*c1 + T[1,0]*s1) + c6*(T[0,1]*c1 + T[1,1]*s1)) - d[5]*(T[0,2]*c1 + T[1,2]*s1) + T[0,3]*c1 + T[1,3]*s1
        p04y = T[2,3] - d[0] - d[5]*T[2,2] + d[4]*(T[2,1]*c6 + T[2,0]*s6)
        #_, Ts = forward_kin(qs3[i], a, d, l)
        #T14 = inv_mat(Ts[0]) * T * inv_mat(Ts[5]) * inv_mat(Ts[4])
        #qs_rrr = inverse_rrr(T14, a[1:4], d[1:4])
        if debug:
            print 'lllh', p04x, p04y, x04x, x04y
            print 'kk', c1, s1, c5, s5, c6, s6
        qs_rrr = inverse_rrr(p04x, p04y, x04x, x04y, a[1:4], d[1:4])
        for j in range(len(qs_rrr)):
            qsol = [qs3[i][0], qs_rrr[j][0], qs_rrr[j][1], qs_rrr[j][2], qs3[i][4], qs3[i][5]]
            if abs(-sin(qsol[1] + qsol[2] + qsol[3])*sin(qsol[4]) - T[2,2]) < ZERO_THRESH:
                qs4.append(qsol)
            #Tsol, _ = forward_kin(qsol, a, d, l)
            #print 'yo', qsol
            #print Tsol**-1 * T06
    if False:
        qs4 = np.array(qs4)[np.lexsort(np.mat(qs4).T)[0]]
        unique_sols = []
        qlast = np.array([-999.]*6)
        for i in range(np.size(qs4,0)):
            if np.sum(abs(qlast - qs4[i])) > ZERO_THRESH:
                unique_sols.append(qs4[i])
                qlast = qs4[i]
        return unique_sols
    else:
        return qs4

def main():
    if True:
        #q = [ 4.07545758,  5.71643082, -4.57552159, -2.79061482, -3.17069678, 1.42865389]
        d1, a2, a3, d4, d5, d6 = [0.1273, -0.612, -0.5723, 0.163941, 0.1157, 0.0922]
        a = [0, a2, a3, 0, 0, 0]
        d = [d1, 0, 0, d4, d5, d6]
        l = [pi/2, 0, 0, pi/2, -pi/2, 0]
        kin = RAVEKinematics()
        rospy.init_node("test_ur_ik")
        start_time = rospy.get_time()
        n = 0
        while not rospy.is_shutdown():
            q = (np.random.rand(6)-.5)*4*pi
            x1 = kin.forward(q)
            pos, euler = PoseConv.to_pos_euler(x1)
            m = np.random.randint(-4,5)
            euler = [euler[0], m*np.pi/2 + 0., euler[2]]
            #euler = [euler[0], 0.*np.pi/2 + m*np.pi, euler[2]]
            T = PoseConv.to_homo_mat(pos, euler)
            #q[4] = 0.
            T = kin.forward(q)
            sols = inverse_kin(T,a,d,l)
            print m, len(sols)
            if False and len(sols) == 0:
                print 'wuh', T
                sols = inverse_kin(T,a,d,l,True)
                
            if True:
                for qsol in sols:
                    #Tsol, _ = forward_kin(qsol, a, d, l)
                    Tsol = kin.forward(qsol)
                    #print qsol
                    #print q
                    #print Tsol
                    #print T
                    diff = Tsol**-1 * T
                    ang, _, _ = mat_to_ang_axis_point(diff)
                    dist = np.linalg.norm(diff[:3,3])
                    #print ang, dist
                    if abs(dist) > 1e-5 or abs(ang) > 1e-5:
                        print 'BAD'
                    else:
                        pass
                        #print 'GOOD'
            n += 1
        time_diff = rospy.get_time() - start_time
        print time_diff, n, n/time_diff

    if False:
        #q = [ 4.07545758,  5.71643082, -4.57552159, -2.79061482, -3.17069678, 1.42865389]
        d1, a2, a3, d4, d5, d6 = [0.1273, -0.612, -0.5723, 0.163941, 0.1157, 0.0922]
        a = [0, a2, a3, 0, 0, 0]
        d = [d1, 0, 0, d4, d5, d6]
        l = [pi/2, 0, 0, pi/2, -pi/2, 0]
        kin = RAVEKinematics()
        rospy.init_node("test_ur_ik")
        start_time = rospy.get_time()
        n = 0
        while not rospy.is_shutdown():
            q = (np.random.rand(6)-.5)*4*pi
            T = kin.forward(q)
            sols = inverse_kin(T,a,d,l)
            #print len(sols)
            if False:
                print len(sols)
                for qsol in sols:
                    #Tsol, _ = forward_kin(qsol, a, d, l)
                    Tsol = kin.forward(qsol)
                    diff = Tsol**-1 * T
                    ang, _, _ = mat_to_ang_axis_point(diff)
                    dist = np.linalg.norm(diff[:3,3])
                    #print ang, dist
                    if abs(dist) > 1e-8 or abs(ang) > 1e-8:
                        print 'BAD'
            n += 1
        time_diff = rospy.get_time() - start_time
        print time_diff, n, n/time_diff
        
    if False:
        #q = [ 4.07545758,  5.71643082, -4.57552159, -2.79061482, -3.17069678, 1.42865389]
        q = (np.random.rand(6)-.5)*4*pi
        d1, a2, a3, d4, d5, d6 = [0.1273, -0.612, -0.5723, 0.163941, 0.1157, 0.0922]
        a = [0, a2, a3, 0, 0, 0]
        d = [d1, 0, 0, d4, d5, d6]
        l = [pi/2, 0, 0, pi/2, -pi/2, 0]
        kin = RAVEKinematics()
        T = kin.forward(q)
        print T
        print forward_kin(q,a,d,l)
        print q
        sols = inverse_kin(T,a,d,l)
        for qsol in sols:
            Tsol, _ = forward_kin(qsol, a, d, l)
            diff = Tsol**-1 * T
            ang, _, _ = mat_to_ang_axis_point(diff)
            dist = np.linalg.norm(diff[:3,3])
            if False:
                if abs(dist) > 1e-6:
                    print '-'*80
                else:
                    print '+'*80
                print 'T', T
                print 'qsol', qsol
                print 'q234', np.sum(qsol[1:4])
                print 'q5', qsol[4]
                print '-sin(q2 + q3 + q4)*sin(q5)', -sin(qsol[1] + qsol[2] + qsol[3])*sin(qsol[4])
                print 'z3', T[2,0]
                if abs(dist) > 1e-6:
                    print '-'*80
                else:
                    print '+'*80
            print ang, dist
        print np.sort(sols,0)
        print len(sols)
        #unique_sols = np.array(sols)[np.where(np.hstack(([True], np.sum(np.diff(np.sort(sols,0), 1, 0),1) > 0.000)))[0]]
        #print unique_sols
        #print len(unique_sols)
        #print len(np.hstack(([True], np.sum(np.diff(np.sort(sols,0), 1, 0),1) > 0.000)))
        for qsol in sols:
            #Tsol, _ = forward_kin(qsol, a, d, l)
            Tsol = kin.forward(qsol)
            diff = Tsol**-1 * T
            ang, _, _ = mat_to_ang_axis_point(diff)
            dist = np.linalg.norm(diff[:3,3])
            print ang, dist
            if abs(dist) > 1e-8:
                print 'BAD'
        
        kin.robot.SetDOFValues(np.array([0.]*6))
        rave_sols = kin.manip.FindIKSolutions(T.A, kin.ik_options)
        rave_list = []
        for qsol in rave_sols:
            rave_list.append(np.array(qsol))
            #Tsol, _ = forward_kin(qsol, a, d, l)
            Tsol = kin.forward(qsol)
            diff = Tsol**-1 * T
            ang, _, _ = mat_to_ang_axis_point(diff)
            dist = np.linalg.norm(diff[:3,3])
            print ang, dist
            if abs(dist) > 1e-8:
                print 'BAD'
        print np.sort(rave_list,0)
            #print diff
            #print q
            #print qsol
            #print '-'*80

    if False:
        q = (np.random.rand(3)-0.5)*4.*np.pi
        T04 = forward_rrr(q)
        print T04
        qs = inverse_rrr(T04)
        print qs
        print T04**-1 * forward_rrr(qs[0])
        print T04**-1 * forward_rrr(qs[1])

if __name__ == "__main__":
    import cProfile
    cProfile.run('main()', 'prof')
    #main()
