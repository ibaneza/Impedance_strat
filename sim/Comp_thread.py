# author: Aurelien Ibanez

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
        return msg

    def msg_manager(self, msg):
        if self.verbose:
            print "read message..."
            print "============================================================"
            print msg
            print "------------------------------------------------------------"
        smsg = msg.split("\n")
        weight = {}
        for m in smsg:
            lm = m.split(" ")

    def finish(self):
        SocketCom.finish(self)