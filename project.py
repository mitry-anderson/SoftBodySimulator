import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import copy
import math
from PIL import Image

objectFile = "shapes/M.png"
k = 10000. # spring stiffness coefficient
c = 10. # damping coefficient
c_air = 1. # air resistance coefficient
m = 0.2 # mass of a point mass (kg)
epsilon =  1.0 # coefficient of restitution
mu = 0.6 # coefficient of friction

initialAngle = 1. # rad
dropHeight = 0.3# m
objectSpacing = 0.1 # m

h = 0.0 # height of the table (m)

dT = np.sqrt(m/k)/2 # 0.001 # time step (s): bigger k --> this has to be smaller
runTime = 2.0 # total simulation time (s)

f_g = np.array([0.,-9.81])*m # force of gravity (N)


def angle_between_vectors(a,b):
    # 
    if (a[0] == 0 and a[1] == 0) or (b[0] == 0 and b[1] == 0):
        return 0.
    else:
        a_mag = np.linalg.norm(a)
        b_mag = np.linalg.norm(b)
        prod = a@b/(a_mag*b_mag)
    return np.arccos(prod)

def rotation_matrix(angle):
    return np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])

def direction(vector):
    if vector[0] == 0. and vector[1] == 0.:
        return None
    else:
        return vector/np.linalg.norm(vector)




class Point():
    def __init__(self,u,v,x_0,y_0,vx_0,vy_0):
        #self.num = pointNum
        self.u = u
        self.v = v
        self.connectedPoints = []
        
        self.pos_0 = np.array([np.copy(x_0),np.copy(y_0)])
        self.pos = np.array([x_0,y_0])
        self.vel = np.array([vx_0,vy_0])
        self.radius = 0.01
        self.m = m
        self.f = 0
        

    def connect(self,connectedPoints):
        self.connectedPoints = connectedPoints

    def delta_0(self,point):
        return self.pos_0 - point.pos_0

    def delta(self,point):
        return self.pos - point.pos

    def compression(self,point):
        return (np.linalg.norm(self.delta(point)) - np.linalg.norm(self.delta_0(point)) )

    def compression_direction(self,point):
        return -direction(self.pos - point.pos)

    def relative_velocity(self,point):
        return self.vel - point.vel

    def reflect_vel(self,direction):
        vmag = np.linalg.norm(self.vel)
        vel_direction = vmag*direction
        self.vel = self.vel - 2*vel_direction
        return self.vel


    def resolve_collisions(self,points):
        for point in points:
            if point.pos is not self.pos:
                if np.linalg.norm(self.pos - point.pos) <= self.radius + point.radius:
                    point.pos += (self.radius + point.radius)*self.compression_direction(point)
                    point.reflect_vel(self.compression_direction(point))


    def update_force(self):
        f = np.array([0.,0.])
        
        for point in self.connectedPoints:
             
            direction = self.compression_direction(point)

            # spring force  
            f = f + k*self.compression(point)*direction

            # damper force
            f = f - c*(direction@self.vel - direction@point.vel)*direction


        # force of gravity
        f = f + f_g

        # force of air resistance:
        f = f - c_air*self.vel

        self.f = f
        # contact forces yay hybrid system!:
        if self.pos[1] <= h:
            self.handle_contact()
        
        




    def accel(self):
        return self.f/self.m

    def handle_contact(self):
            
        f = self.f
        n_c = np.array([0.,1.])
        n_t = np.array([1.,0.])
        # forces applied when pressing on surface
        # if f[1] < 0:
            
        #     # normal force
        #     f_n = np.array([0.,np.abs(f[1])])
        #     f = f + f_n

        #     # force of friction
        #     sign = -np.sign(self.vel[0])
        #     f = f + sign*mu*np.array([np.linalg.norm(f_n),0.])

        # # impulsive force
        # v_1 = np.array([self.vel[0],-epsilon*self.vel[1]])
        # f_i = m*(v_1 - self.vel)/(dT)
        # f = f + f_i

        if n_c@f < 0:
            # normal force
            f_n = (n_c@(-f))*n_c
            f = f + f_n

            # force of friction
            
            sign = -np.sign((n_t@self.vel)*n_t)
            f = f + sign*mu*(np.linalg.norm(f_n)*n_t)


        # impulsive force
        v_1 = -epsilon*(n_c@self.vel)*n_c + (n_t@self.vel)*n_t
        f_i = m*(v_1 - self.vel)/(dT)
        f = f + f_i

        self.f = f
                


    def update_state(self,points):
        self.vel = self.vel + self.accel()*dT
        self.pos = self.pos + self.vel*dT
        y = 1
        
        


        


def time_step(points):
    for point in points:
        point.update_force()
    for point in points:
        point.update_state(points)
    



def generate_plot_points(points):
    x = []
    y = []
    for point in points:
        x.append(point.pos[0])
        y.append(point.pos[1])
    x.append(points[0].pos[0])
    y.append(points[0].pos[1])
    return (x,y)
        




def load_from_image(filepath,scaleFactor):
    image = Image.open(filepath)
    data = np.asarray(image)
    
    points = []
    
    
    for u in range(data.shape[0]):
        for v in range(data.shape[1]):
            if data[u,v,0] == 0:
                points.append(Point(u,v, u*scaleFactor, v*scaleFactor,0.0,0.0))
            else:
                points.append(None)
                
    
    pointList = []
    for point in points:
        if point is not None:
            u = point.u
            v = point.v
            pointsToConnect = []
            
            if v+1 < data.shape[1] and data[u,v+1,0] == 0:
                pointsToConnect.append(points[v+1+data.shape[0]*u])

            if v-1 > 0 and data[u,v-1,0] == 0:
                pointsToConnect.append(points[v-1+data.shape[0]*u])

            if u+1 < data.shape[0] and data[u+1,v,0] == 0:
                pointsToConnect.append(points[v+data.shape[0]*(u+1)])

            if u-1 > 0 and data[u-1,v,0] == 0:
                pointsToConnect.append(points[v+data.shape[0]*(u-1)])

            if u+1 < data.shape[0] and v+1 < data.shape[1] and data[u+1,v+1,0] == 0:
                pointsToConnect.append(points[v+1+data.shape[0]*(u+1)])

            if u-1 > 0 and v-1 > 0 and data[u-1,v-1,0] == 0:
                pointsToConnect.append(points[v-1+data.shape[0]*(u-1)])

            if u+1 < data.shape[0] and v-1 > 0 and data[u+1,v-1,0] == 0:
                pointsToConnect.append(points[v-1+data.shape[0]*(u+1)])

            if u-1 > 0 and v+1 < data.shape[1] and data[u-1,v+1,0] == 0:
                pointsToConnect.append(points[v+1+data.shape[0]*(u-1)])

            point.connect(pointsToConnect)
            pointList.append(point)
        

    
    print("Created " + str(len(pointList)) + " points.")
    return pointList
    


def reorient_object(points,angle,xOff,yOff):
    translation = np.array([xOff,yOff])
    for point in points:
        point.pos = rotation_matrix(angle)@point.pos
        point.pos_0 = rotation_matrix(angle)@point.pos_0
        point.pos = point.pos + translation
        point.pos_0 = point.pos_0 + translation

    return points






if __name__ == "__main__":

    
    t = 0

    pointRecord = []

    points = load_from_image(objectFile,objectSpacing)
    reorient_object(points, initialAngle,0.0,dropHeight)

    print("(k/m)dt = " + str( ((k/m)*dT) ) )
    print("Beginning Simulation")
    print("dT: " + str(round(dT,6)) + "s. ")
    t_func = 0.0
    while(t<runTime):
        
        printstr = "\rSimulation time: " + str(round(t,3)) + "s | Function time: " + str(round(t_func,3)) + "s "
        print(printstr,end = '')
        t_0 = time.perf_counter()
        time_step(points)
        t_1 = time.perf_counter()
        pointRecord.append(copy.deepcopy(points))
        t += dT
        t_func += t_1 - t_0
    
    timePerCall = t_func/float(len(pointRecord))
    print("\nSimulation Complete with " + str(timePerCall) + "s per iteration.")


    pointRecord = pointRecord[::10]
    def animate(i):
        if i < len(pointRecord):
            pointGen = generate_plot_points(pointRecord[i])
            line.set_xdata(pointGen[0])
            line.set_ydata(pointGen[1])

        return line,
    
    fig, ax = plt.subplots()
    pointGen = generate_plot_points(pointRecord[0])
    line, = ax.plot(pointGen[0],pointGen[1],marker='o',linewidth=0)

    ani = animation.FuncAnimation(fig,animate,interval=0,frames=len(pointRecord),save_count=1000)
        
    plt.xlim(-2,2)
    plt.ylim(-0.3,2)
    plt.gca().set_aspect('equal', adjustable='box')

    # f = "animation.gif" 
    # writergif = animation.PillowWriter(fps=30) 
    # ani.save(f, writer=writergif)
    
    plt.show()
    
   


    
    