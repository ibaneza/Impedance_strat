# author: Aurelien Ibanez

import copy
from arboris.observers import SocketCom

class CompThreadCom(SocketCom):
    def __init__(self, Comp_thread_path, host="127.0.0.1", port=5556, options="", verbose=False):
        SocketCom.__init__(self, host, port)
        self.app_call = [abspath(Comp_thread_path), "-socket", self.host, str(self.port)] + shlex.split(options)
        self.verbose = verbose

    def init(self,world,timeline): # timeline and world are ABSOLUTELY useless here, so pass anything you want
        subprocess.Popen(self.app_call)
        SocketCom.init(self, world, timeline)
        self.world = world


    def update(self, dt):
        rstate = self.LQPctrl._get_robot_state(dt)
        msg = ""
        msg += self.edit_message(rstate, dt)
        try:
            self.conn.send(msg)
            nmsg = self.conn.recv(2048, 0)
            self.msg_manager(nmsg)
        except socket.error:
            print "connection lost"

    def edit_message(self, rstate, dt, prec=4):
        # Send information
        msg = ""
        msg += "Constants {0} {1} {2} {3} {4}\n".format(self.M,self.zc,self.g,self.dt,self.h)
        for (name, matrix) in [("Jacobian",self.J), ("dJacobian",self.dJ),
                                ("JacobianI",self.Ji), ("dJJi", self.dJJi),
                                    ("JtiHJi",self.JtiHJi), ("Pref",self.Pref),
                                        ("FDIS",self.FDIS)]:
            msg += name
            msg += " {0} {1}".format(matrix.shape[0],matrix.shape[1])
            for lin in range(0,matrix.shape[0]):
                for col in range(0,matrix.shape[1]):
                    msg += " {0}".format( matrix(lin,col ) )
            msg += "\n"
        for (name, vector) in [("x", self.x), ("dx", self.dx),
                                ("ddx", self.ddx), ("xdes", self.xdes),
                                    ("guess", self.guess), ("increments", self.increments)]:
            msg += name
            msg += " {0}".format(vector.shape[0])
            for lin in range(0,vector.shape[0]):
                msg += " {0}".format( vector(lin ) )
            msg += "\n"
        return msg

    def msg_manager(self, msg):
        if self.verbose:
            print "read message..."
            print "============================================================"
            print msg
            print "------------------------------------------------------------"
        smsg = msg.split("\n")
        result = []
        for m in smsg:
            if m == "COMP_THREAD_ERROR":
                print "Computing thread returned with error"
                return
            lm = m.split(" ")
            for w in lm:
                results.append(float(w))
        return result

    def set_constants( M, zc, g, dt, h ):
        self.M = M
        self.zc = zc
        self.g = g
        self.dt = dt
        self.host = h

    def set_matrices( J, dJ, Ji, dJJi, JtiHJi ):
        self.J = copy.copy( J )
        self.dJ = copy.copy( dJ )
        self.Ji = copy.copy( Ji )
        self.dJJi = copy.copy( dJJi )
        self.JtiHJi = copy.copy( JtiHJi )

    def set_desired_kinematics( xdes, Pref ):
        self.xdes = copy.copy( xdes )
        self.Pref = copy.copy( Pref )

    def set_current_kinematics( x, dx, ddx ):
        self.x = copy.copy( x )
        self.dx = copy.copy( dx )
        self.ddx = copy.copy( ddx )

    def set_disturbance( FDIS ):
        self.FDIS = copy.copy( FDIS )

    def set_simplex_parameters( guess, increments ):
        self.guess = copy.copy( guess )
        self.increments = copy.copy( increments )

    def finish(self):
        SocketCom.finish(self)