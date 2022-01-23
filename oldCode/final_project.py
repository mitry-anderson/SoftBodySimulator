import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import copy
import math

# Globals:
numPoints = 9
k_lin = 1000.
k_max = k_lin*10
c_lin = 1.
m = 0.5
mu = 0.5

dT = 0.01
f_g = np.array([0.,-9.81])*m

def angle_between_vectors(a,b):
    if (a[0] == 0 and a[1] == 0) or (b[0] == 0 and b[1] == 0):
        return 0.
    else:
        a_mag = np.linalg.norm(a)
        b_mag = np.linalg.norm(b)
        prod = a@b/(a_mag*b_mag)
    return np.arccos(prod)

def rotation_matrix(angle):
    return np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])

x_axis = np.array([1,0])
y_axis = np.array([0,1])



class Point():
    def __init__(self,pointNum,x_0,y_0,vx_0,vy_0):
        self.num = pointNum
        self.connectedPoints = []
        
        self.force_applied = False
        self.pos_0 = np.array([np.copy(x_0),np.copy(y_0)])
        self.pos = np.array([x_0,y_0])
        self.vel = np.array([vx_0,vy_0])
        self.pos_prev = np.array([x_0,y_0])
        self.rad = 0.01

    def connect(self,connectedPoints):
        self.connectedPoints = connectedPoints
        

    def mass_vector(self):
        mVector = np.zeros((2,2*numPoints))
        mVector[0,2*self.num] = m
        mVector[1,2*self.num+1] = m
        # for point in self.connectedPoints:
        #     mVector[0,2*point.num] = m/2
        #     mVector[1,2*point.num+1] = m/2
        return mVector
    
    def direction_to_point(self,point):
        delta = point.pos - self.pos
        return delta/np.linalg.norm(delta)

    def force(self):
        f = np.zeros((2,))
        # print("-----------------")
        for point in self.connectedPoints:

            t_0 = self.pos_0 - point.pos_0
            t_now = self.pos - point.pos
        
            delta = np.linalg.norm(t_now) - np.linalg.norm(t_0)
            v_rel_now = (self.vel - point.vel)
            
            if (np.isnan(delta) or np.isnan(t_now[0]) or np.isnan(t_now[1])):
                print("WARNING--undefined behavior")
                d_lin = t_0
            elif delta > 0:
                d_lin = -np.abs(delta)*t_now/np.linalg.norm(t_now)
            elif delta == 0:
                d_lin = 0.
            else:
                d_lin = np.abs(delta)*t_now/np.linalg.norm(t_now)
            # print("t_now")
            # print(t_now/np.linalg.norm(t_now))
            
            # print("d_lin: ")
            # print(d_lin)
            # print("v_rel: ")
            # print(v_rel_now)
            
        
            f += k_lin*d_lin - c_lin*v_rel_now 
        
        if(self.pos[1] <= 0.):
                # print("contact!")
                # impulseTime = 30*dT
                f += np.array([0.,-m*self.vel[1]/dT]) # impulsive force
                
        if(self.pos[1] > 0.):
            f += f_g # gravitational force
        # else:
        #     print(d_lin)


        

        
        # print("force: " + str(f))
        return f

    def collision(self,point):
        if np.linalg.norm(self.pos - point.pos) <= self.rad + point.rad:
            return True
        else:
            return False

    def reflect_vel(self,direction):
        vmag = np.linalg.norm(self.vel)
        vel_direction = vmag*direction
        self.vel = self.vel - 2*vel_direction
        return self.vel
    

    def update(self,a):
        
        self.pos_prev = self.pos
        self.vel = self.vel +  a*dT
        self.pos = self.pos + self.vel*dT
        

        # if(self.pos[1] < 0.):
        #     self.pos[1] = 0.
        #     self.vel[1] = -self.vel[1]
            # self.vel[0] = 0
        # print("pos" + str(self.num))
        # print(self.pos)
        # print("vel" + str(self.num))
        # print(self.vel)

    

def time_step(points):
    mass_matrix = np.zeros((2*numPoints,2*numPoints))
    f = np.zeros((2*numPoints,1))
    a = np.zeros((2*numPoints,1))
    for i,point in enumerate(points):
        # print(mass_matrix[2*i,:])
        massVector = point.mass_vector()
        mass_matrix[2*i,:] = massVector[0,:]
        mass_matrix[2*i+1,:] = massVector[1,:]
        fVec = point.force()
        f[2*i,:] = fVec[0]
        f[2*i+1,:] = fVec[1]
    # print("<><><><><><><><> TIME STEP <><><><><><><><>")
    # print("mass")
    # print(mass_matrix)
    # print(np.linalg.inv(mass_matrix))
    # print("force")
    # print(f)
    a = np.linalg.inv(mass_matrix)@f
    # if np.isnan(a[0]) or np.isnan(a[1]):
    #     print("a isnan")
    # print("accel")
    # print(a)
    

    for i,point in enumerate(points):
        point.update( np.array([a[2*i,0],a[2*i+1,0]]))
        if(point.pos[1] < 0.):
            point.pos[1] = 0.
    for point1 in points:
        for point2 in points:
            if point1 is not point2:
                if point1.collision(point2):
                    direction = point1.direction_to_point(point2)
                    point2.pos += (point1.rad + point2.rad)*direction
                    point2.reflect_vel(direction)
                    



    
def generate_plot_points(points):
    x = []
    y = []
    for point in points:
        x.append(point.pos[0])
        y.append(point.pos[1])
    x.append(points[0].pos[0])
    y.append(points[0].pos[1])
    return (x,y)
        



if __name__ == "__main__":


    p1 = Point(0,0.5,0.5,0,0)
    p1.vel = np.array([0,-1])
    direction = np.array([1,1])
    direction = direction/np.linalg.norm(direction)
    print(p1.reflect_vel(direction))

    point0 = Point(0,0.5,0.5,0,0)
    point1 = Point(1,0.5,0.6,0,0)
    point2 = Point(2,0.6,0.6,0,0)
    point3 = Point(3,0.6,0.5,0,0)

    point4 = Point(4,0.5,0.7,0,0)
    point5 = Point(5,0.6,0.7,0,0)
    point6 = Point(6,0.7,0.7,0,0)
    point7 = Point(7,0.7,0.6,0,0)
    point8 = Point(8,0.7,0.5,0,0)

    
    # superSquare
    point0.connect([point1,point3,point2])
    point1.connect([point0,point4,point2])
    point2.connect([point1,point3,point5,point7,point0,point4,point6])
    point3.connect([point2,point0,point8])

    point4.connect([point1,point5,point2])
    point5.connect([point4,point2,point6])
    point6.connect([point5,point7,point2])
    point7.connect([point2,point6,point8])
    point8.connect([point7,point3,point2])


    points = [point0, point1, point2, point3, point4, point5, point6, point7, point8]

    # angle = 0.0
    # for point in points:
    #     print(point.pos)
    #     point.pos = rotation_matrix(angle)@point.pos
    #     point.pos_0 = rotation_matrix(angle)@point.pos_0
    #     print(point.pos)
    # Triangle:
    # point0.connect([point1,point2])
    # point1.connect([point0,point2])
    # point2.connect([point1,point0])
    # points = [point0, point1, point2]

    
    
    t = 0

    pointRecord = []

    while(t<5):
        #print("time: " + str(t))
        time_step(points)
        pointRecord.append(copy.deepcopy(points))
        t += dT
        # time.sleep(0.01)
       

    def animate(i):
        if i < len(pointRecord):
            pointGen = generate_plot_points(pointRecord[i])
            line.set_xdata(pointGen[0])
            line.set_ydata(pointGen[1])

        return line,
        


    fig, ax = plt.subplots()
    pointGen = generate_plot_points(pointRecord[0])
    line, = ax.plot(pointGen[0],pointGen[1],marker='o',linewidth=0)

    ani = animation.FuncAnimation(fig,animate,interval=dT/10,save_count=50)
        
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    
    
    plt.show()
        



            
            


   

    