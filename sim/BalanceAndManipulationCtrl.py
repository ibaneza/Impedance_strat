# author: Aurelien Ibanez

from arboris.homogeneousmatrix import roty
from numpy import array, cos, sin, ones, dot, vstack, zeros, arange, tile, eye, pi, linspace, sqrt, exp, mod, sign
from numpy.linalg import norm, inv
from scipy.interpolate import piecewise_polynomial_interpolate as ppi

from LQPctrl.walk import WalkingCtrl, BalanceCtrl, BalanceCtrlKpAdapt, BalanceTask
from Comp_thread import CompThreadCom

from arboris.homogeneousmatrix import transl, rotx, rotz, roty

from numpy import array, arange, pi, zeros, eye, dot, linspace, ones, sqrt, transpose
from numpy.linalg import norm,inv
import copy
import gc

def _get_matrices(h, dt, hong):
    Px = zeros((h, 3))
    Pu = zeros((h, h))
    for i in range(h):
        Px[i, :] = [1, (i+1)*dt, ((i+1)*dt)**2/2 - hong]
        diag_i = (1 + 3*i + 3*i**2)*dt**3/6 - dt*hong
        Pu[range(i,h), range(h-i)] = diag_i
    return Px, Pu

def _get_matrices_dist(h, dt, M, z, g, F):
    Px = zeros((h, 3))
    Pu = zeros((h, h))
    for i in range(h):
        Px[i, :] = [1, (i+1)*dt, ((i+1)*dt)**2/2 - (M*z)/(M*g-F[2,i])]
        diag_i = (1 + 3*i + 3*i**2)*dt**3/6 - dt*z/g
        Pu[range(i,h), range(h-i)] = diag_i
    return Px, Pu

class BalAndManipCtrl( BalanceCtrlKpAdapt ):
    def __init__(self, comp_thread, com_pos, QonR, horizon, dt, feet, world, hand_box_const, grab_dist_min, task, frame, root_frame, bodies, force_name, box_ctrl, mode, Kp, Kd):
        gc.enable()
        print 'Garbage Collector state = ', gc.isenabled()

    ##A# New stuff
        self.screated = False
        self.world = world
        self.hand_conf = hand_box_const
        self.box_id = 0
        self.good_config = None
        self.grabDistMin = grab_dist_min

        self.task = task
        self.frame = frame
        self.root_frame = root_frame
        self.bodies = bodies

        self.force_name = force_name
        self.box_ctrl= box_ctrl
    ##A# End of New stuff

        self.com_pos = com_pos
        self.QonR = QonR
        self.h = int(horizon/dt)
        box_ctrl.horizon = self.h
        self.dt = dt
        self.feet = feet
        self.root = None

        self.x = zeros((3,))
        self.dx = zeros((3,))
        self.ddx = zeros((3,))
        self.x_des = zeros((3,))
        ##A# Compute total mass
        self.Ma = 0.
        for b in self.bodies:
            self.Ma += b.mass[3,3]

        self.mode = mode
        self._Kp = Kp
        self._Kd = Kd

        self.comp_thread = comp_thread

    def init(self, world, LQP_ctrl, tasks, events, bodies):
        self.comp_thread.init( world, 0.0 )

    ##A# New stuff
        self.count = 0
        ##A# Init factory variables
        self.prevJ = self.frame.jacobian[0:3,:] - self.root_frame.jacobian[0:3,:]
        self.prevdQ = self.world.gvel
        self.iOutput = zeros((3,))
        self.prevdx = zeros((3,))


    ##A# End of New stuff

        from arboris.controllers import WeightController
        self.world = world
        self.LQP_ctrl = LQP_ctrl
        self.tasks = tasks
        self.events = events
        self.bodies = bodies
        self.gravity = 0.
        for c in world._controllers:
            if isinstance(c, WeightController):
                self.gravity = abs(c.gravity)
        self.Hfloor = LQP_ctrl.options['Hfloor']
        self.change_goal(self.com_pos)

    def update(self, pos, vel, counter=0., dt = 0.):
    ##A# New stuff
        print 'Garbage collection'
        gc.collect()
        print 'Garbage collection DONE'
        self.count += 1
        self.CanGrabBox = False
        dist = []
        if self.hand_conf:
            for conf in self.hand_conf:
                dist.extend([self._dist(hc) for hc in conf])
        if min(dist) < self.grabDistMin:
            self.CanGrabBox = True
        # elif min(dist) < 0.1:
            # self.task.ctrl.Kp = 1.


        if not self.CanGrabBox:
            pass
        # print "Dist Hand-Handle: " , min(dist) , self.task.ctrl.Kp

        ##A# Detecting good grasp configuration
        confs = self.hand_conf[self.box_id]
        self.good_config = min([(max([self._dist([c]) for c in conf]), conf) for conf in confs])[1]
        if self.good_config is not None:
             self.has_grabbed = all([c.is_enabled() for c in self.good_config])
             if self.CanGrabBox:
                # print "CAN GRAB BOX!"

                ##A# Kp adaptation
                #self.task.ctrl.Kp = 1.e0

                ##A# \dot{x}^{cmd} computation
                self.iOutput[0:3] += dt*self.task.ctrl.output[0:3]

                # print self.task.ctrl.vel.shape
                self.dx_cmd =  self.iOutput

                ##A# Desired and actual positions
                self.x_des = self.task.ctrl.goal[0:3,3]
                self.x = self.task.ctrl.pos[0:3,3]

                ##A# Actual velocity and acceleration
                self.dx = self.task.ctrl.vel[0:3]
                self.ddx = (self.dx - self.prevdx) / dt
                self.prevdx = self.dx

                ##A# Joint acceleration
                dQ = self.world.gvel
                ddQ = (dQ - self.prevdQ) / dt
                self.prevdQ = dQ

                ##A# Jacobian, Jacobian's derivative, transpose and pseudoinverse
                J = self.frame.jacobian[0:3,:]  - self.root_frame.jacobian[0:3,:]
                self.dJ = (J - self.prevJ)/dt
                self.prevJ = J
                self.Jt = transpose(J)
                self.Ji = dot( self.Jt, inv( dot( J,self.Jt ) ) )

                ##A# Joint Space Mass matrix
                H = self.world.mass

                JHiJt = dot( J, dot(inv(H), self.Jt))
                self.JtiHJi = inv( JHiJt )

                for c in self.good_config:
                    c.enable()
                    #self.force[:,self.count] = c._force[:]
             else:
                #self.force[:,self.count] = zeros(3)[:]
                pass

    ##A# End of New stuff
        #self.change_goal(None)
        com_hat = pos[:, [0, 2]]
        com_h   = pos[0, 1]
        self.com_h = com_h
        self.com_hat = com_hat
        hong = com_h/self.gravity #### WARNING: pb if the ground is not horizontal!!! if n ground != gravity direction

        zmp_ref = self.zmp_ref[counter:counter+self.h]
        zmp_ref = vstack( (zmp_ref, dot(ones((self.h-len(zmp_ref), 1)), self.zmp_ref[-1].reshape((1, 2))) ) )
        #print 'ZMP_REF = ', zmp_ref[0,:]
        self.pref = zmp_ref


        if self.mode == 0:
            Px, Pu = _get_matrices(self.h, self.dt, hong)
            ddV_com = - dot( inv(dot(Pu.T, Pu) + self.QonR*eye(self.h)), dot(Pu.T, dot(Px, com_hat) - zmp_ref) )
        else:
            FDIS = self.box_ctrl.ask_force()
            for i in range(2):
                zmp_ref[:,i] = zmp_ref[:,i] - self.com_h/(self.Ma*self.gravity - FDIS[2,:])*1.*FDIS[i,:] #5.
            #Px, Pu = _get_matrices_dist(self.h, self.dt, self.Ma, self.com_h, self.gravity , self.box_ctrl.ask_force())
            Px, Pu = _get_matrices(self.h, self.dt, hong)
            ddV_com = - dot( inv(dot(Pu.T, Pu) + self.QonR*eye(self.h)), dot(Pu.T, dot(Px, com_hat) - zmp_ref) )



        if self.CanGrabBox:
            if self.mode == 2:
                reduc = 0.1
                self.A = [[1, self.dt, self.dt**2/2],[0,1, self.dt],[0,0,1]]
                self.B = [self.dt**3/6, self.dt**2/2, self.dt]
                if not self.screated:
                    self.values = 10.*ones(2*self.h+1) #(3*self.h)
                    self.values[2*self.h]=100. #:3*self.h] = 100.

                    for i in range(self.h):
                        self.values[i] = ddV_com[i][0]
                        self.values[i+self.h] = ddV_com[i][1]
                    self.screated = True
                boundaries = zeros(2*self.h+1) #3*self.h)
                for (d,f) in [(0,self.h),(self.h,2*self.h)]: #,(2*self.h,3*self.h)]:
                    maxi = 0
                    for i in range(self.h):
                        if abs(self.values[i]) > maxi:
                            maxi = abs(self.values[i])
                    boundaries[d:f] = .5*maxi
                #boundaries[:] = 0.5*abs(self.values[:])
                boundaries[2*self.h] = 20. * 1. / reduc



                self.Kpinit = self.values[2*self.h]
                print 'Optimizing Input with Kp = ', self.Kpinit
                self.comp_thread.set_mode( 3 )
                self.comp_thread.set_kpinit( self.Kpinit )

                FDIS = self.box_ctrl.ask_force()
                self.dJJi = dot( self.dJ,self.Ji )
                self.comp_thread.set_constants( self.Ma, self.com_h, self.gravity, self.dt, self.h )
                self.comp_thread.set_matrices( self.prevJ, self.dJ, self.Ji, self.dJJi, self.JtiHJi )
                self.comp_thread.set_desired_kinematics( self.x_des, transpose(self.pref) )
                self.comp_thread.set_current_kinematics( self.x, self.dx, self.ddx )
                self.comp_thread.set_disturbance( FDIS )
                self.comp_thread.set_simplex_parameters( self.values[0:2*self.h], boundaries[0:2*self.h] )

                self.values[0:2*self.h] = self.comp_thread.update( self.dt )

                self.comp_thread.set_mode( 2 )
                self.comp_thread.set_simplex_parameters( self.values, reduc * boundaries )
                self.values = self.comp_thread.update( self.dt )

                print ddV_com[0]
                ddV_com[0][0] = self.values[0]
                ddV_com[0][1] = self.values[self.h]
                print ddV_com[0]
                self.task.ctrl.Kp = self.values[2*self.h]
                self.task.ctrl.Kd = 2*sqrt(abs(self.task.ctrl.Kp))

                #self.task.ctrl.Kd = values[3*self.h]
            elif self.mode == 1:
                self.task.ctrl.Kp = self._Kp
                self.task.ctrl.Kd = self._Kd
            elif self.mode == 0:
                self.task.ctrl.Kp = self._Kp
                self.task.ctrl.Kd = self._Kd
        self.input = ddV_com[0]
        #print zmp_ref[0,:]
        print self.task.ctrl.Kp, self.task.ctrl.Kd
        return ddV_com[0]