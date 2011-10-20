# author: Aurelien Ibanez

import copy
import socket
import subprocess, shlex
import time
from os.path import abspath
from arboris.observers import SocketCom
from numpy import zeros

class CompThreadCom(SocketCom):
    def __init__(self, Comp_thread_path, host="127.0.0.1", port=5556, options="", verbose=False):
        SocketCom.__init__(self, host, port, timeout = None)
        self.app_call = [abspath(Comp_thread_path), "-socket", self.host, str(self.port)] + shlex.split(options)
        self.verbose = verbose

    def init(self,world,timeline): # timeline and world are ABSOLUTELY useless here, so pass anything you want
        subprocess.Popen(self.app_call)
        SocketCom.init(self, world, timeline)
        self.world = world


    def update(self, dt):
        msg = ""
        msg += self.edit_message(dt)
        try:
            #print "Python>> SENDING"
            self.conn.send(msg)
            #print "Python>> RECEIVING"
            nmsg = self.conn.recv(8192, 0)
            #print "Python>> PARSING"
            result = self.msg_manager(nmsg)
            #print "Python>> RETURNING"
            return result
        except socket.error:
            print "connection lost"
            return "NORES"

    def edit_message(self, dt, prec=4):
        # Send information
        msg = ""
        msg += "Constants {0} {1} {2} {3} {4} {5} {6}\n".format(self.M,self.zc,self.g,self.dt,self.h, self.mode, self.kpinit)
        for (name, matrix) in [("Jacobian",self.J), ("dJacobian",self.dJ),
                                ("JacobianI",self.Ji), ("dJJi", self.dJJi),
                                    ("JtiHJi",self.JtiHJi), ("Pref",self.Pref),
                                        ("FDIS",self.FDIS)]:
            msg += name
            msg += " {0} {1}".format(matrix.shape[0],matrix.shape[1])
            for lin in range(0,matrix.shape[0]):
                for col in range(0,matrix.shape[1]):
                    msg += " {0}".format( matrix[lin,col] )
            msg += "\n"
        for (name, vector) in [("x", self.x), ("dx", self.dx),
                                ("ddx", self.ddx), ("xdes", self.xdes),
                                    ("guess", self.guess), ("increments", self.increments),
                                    ("xc", self.xc), ("dxc", self.dxc),
                                    ("ddxc", self.ddxc)]:
            msg += name
            msg += " {0}".format(vector.shape[0])
            for lin in range(0,vector.shape[0]):
                msg += " {0}".format( vector[lin] )
            msg += "\n"
        return msg

    def msg_manager(self, msg):
        if self.verbose:
            print "read message..."
        if msg == "COMP_THREAD_ERROR":
            print "Computing thread returned with error"
            return
        smsg = msg.split("\n")
        print "Size of message read in ARBORIS/PYTHON = ", int(smsg[0])
        result = zeros(int(smsg[0]))
        for w in range(1,int(smsg[0])+1):
            result[w-1] = (float(smsg[w]))
        return result

    def set_constants(self, M, zc, g, dt, h ):
        self.M = M
        self.zc = zc
        self.g = g
        self.dt = dt
        self.h = h

    def set_matrices(self, J, dJ, Ji, dJJi, JtiHJi ):
        self.J = copy.copy( J )
        self.dJ = copy.copy( dJ )
        self.Ji = copy.copy( Ji )
        self.dJJi = copy.copy( dJJi )
        self.JtiHJi = copy.copy( JtiHJi )

    def set_desired_kinematics(self, xdes, Pref ):
        self.xdes = copy.copy( xdes )
        self.Pref = copy.copy( Pref )

    def set_current_kinematics(self, x, dx, ddx, xc, dxc, ddxc ):
        self.x = copy.copy( x )
        self.dx = copy.copy( dx )
        self.ddx = copy.copy( ddx )
        self.xc = copy.copy( xc )
        self.dxc = copy.copy( dxc )
        self.ddxc = copy.copy( ddxc )

    def set_disturbance(self, FDIS ):
        self.FDIS = copy.copy( FDIS )

    def set_simplex_parameters(self, guess, increments ):
        self.guess = copy.copy( guess )
        self.increments = copy.copy( increments )

    def set_mode( self, mode ):
        self.mode = mode

    def set_kpinit( self, kpinit ):
        self.kpinit = kpinit

    def finish(self):
        SocketCom.finish(self)