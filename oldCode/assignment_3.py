import numpy as np
import matplotlib.pyplot as plt
from final_project_helper import LCPSolve, assignment_3_render

# DEFINE GLOBAL PARAMETERS
L = 0.4
MU = 0.3
EP = 0.5
dt = 0.01
m = 0.3
g = np.array([0., -9.81, 0.])
rg = 1./12. * (2 * L * L)
M = np.array([[m, 0, 0], [0, m, 0], [0, 0, m * rg]])
Mi = np.array([[1./m, 0, 0], [0, 1./m, 0], [0, 0, 1./(m * rg)]])
DELTA = 0.001
T = 150

def get_contacts(q):
    """
        Return jacobian of the lowest corner of the square and distance to contact
        :param q: <np.array> current configuration of the object
        :return: <np.array>, <float> jacobian and distance
    """
    # ------------------------------------------------
    # FILL WITH YOUR CODE
    x = q[0]
    y = q[1]
    theta = q[2]
    nx = 0.0
    ny = 1.0
    l = np.sqrt(2)*L

    points = np.array([[L/2,L/2],[-L/2,L/2],[L/2,-L/2],[-L/2,-L/2]])
    points = np.transpose(points)
    R = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])

    points_rotated = np.matmul(R,points)
    index = np.argmin(points_rotated[1,:])
    ry = points_rotated[1,index]
    rx = points_rotated[0,index]
    
    # nCorner = np.mod(round( ( -theta ) / (np.pi/2) ) , 3) #round(mod(( theta - pi/4 )/(pi/2)), 4)

    # rx = 0.5*l*np.sin(-theta + 0.5*np.pi*nCorner)
    # ry = -0.5*l*np.cos(-theta + 0.5*np.pi*nCorner)

    r = np.array([rx,ry])
    n = np.array([0,1])
    t = np.array([1,0])
    jac = np.array([[ny, nx],[-nx,ny],[np.cross(r,t), np.cross(r,n)]])
    # jac = np.array([ [ny, nx] , [-nx, ny], [-(rx*nx + ry*ny), (rx*ny - ry*nx)]])
    phi = y + ry
    #print("x: " + str(x) + "  y: " + str(y) + "  rx: " + str(rx) + "  ry: " + str(ry) + "  phi: " + str(phi))
    # ------------------------------------------------
    return jac, phi


def form_lcp(jac, v):
    """
        Return LCP matrix and vector for the contact
        :param jac: <np.array> jacobian of the contact point
        :param v: <np.array> velocity of the center of mass
        :return: <np.array>, <np.array> V and p
    """
    # ------------------------------------------------
    # FILL WITH YOUR CODE
    Jt = jac[:,0]
    Jn = jac[:,1]
    Jn_t = np.transpose(Jn)
    Jt_t = np.transpose(Jt)
    #Jhat = np.block([[Jn],[-Jt],[Jt]])

    invM = np.linalg.inv(M)
    fe = m*g
    invMfe = np.matmul(invM,fe)

    #bigBlock = np.matmul(np.matmul(np.transpose(Jhat),invM),Jhat)*dt
    V = np.array([[np.matmul(np.matmul(Jn_t,invM),Jn)*dt , -np.matmul(np.matmul(Jn_t,invM),Jt)*dt, np.matmul(np.matmul(Jn_t,invM),Jt)*dt , 0], \
                  [-np.matmul(np.matmul(Jt_t,invM),Jn)*dt, np.matmul(np.matmul(Jt_t,invM),Jt)*dt , -np.matmul(np.matmul(Jt_t,invM),Jt)*dt, 1], \
                  [np.matmul(np.matmul(Jt_t,invM),Jn)*dt , -np.matmul(np.matmul(Jt_t,invM),Jt)*dt, np.matmul(np.matmul(Jt_t,invM),Jt)*dt , 1], \
                  [MU                                    ,  -1                                   , -1                                    , 0] ])
    
    p = np.array([np.matmul(Jn_t,((1+EP)*v + dt*invMfe) ), \
                  np.matmul(-Jt_t,(v + dt*invMfe))       , \
                  np.matmul(Jt_t,(v + dt*invMfe))        , \
                  0                                         ])
   
    # ------------------------------------------------
    return V, p


def step(q, v):
    """
        predict next config and velocity given the current values
        :param q: <np.array> current configuration of the object
        :param v: <np.array> current velocity of the object
        :return: <np.array>, <np.array> q_next and v_next
    """
    # ------------------------------------------------
    # FILL WITH YOUR CODE
    J,phi = get_contacts(q)
    invM = np.linalg.inv(M)
    fe = m*g
    if(phi >= DELTA):
        
        v_next = v + dt*np.matmul(invM, fe)
        q_next = q + dt*v_next
    else:
        
        V,p = form_lcp(J,v)
        fc = lcp_solve(V,p)
        Jt = J[:,0]
        Jn = J[:,1]
        qp = np.array([0,DELTA,0])
        
        fnet= fe + Jn*fc[0] - Jt*fc[1] + Jt*fc[2]
        v_next = v + dt*np.matmul(invM,fnet)
        q_next = q + dt*v_next + qp
        
        
    
    #q_next = None  # TODO: Replace None with your result
    #v_next = None
    # ------------------------------------------------
    return q_next, v_next


def simulate(q0, v0):
    """
        predict next config and velocity given the current values
        :param q0: <np.array> initial configuration of the object
        :param v0: <np.array> initial velocity of the object
        :return: <np.array>, <np.array> q and v trajectory of the object
    """
    # ------------------------------------------------
    # FILL WITH YOUR CODE

    q = np.zeros((3, T))  
    v = np.zeros((3, T))
    q[:,0] = q0
    v[:,0] = v0
    for i in range(T-1):
        q_next, v_next = step(q[:,i],v[:,i])
        q[:,i+1] = q_next
        v[:,i+1] = v_next
    # ------------------------------------------------
    return q, v


def lcp_solve(V, p):
    """
        DO NOT CHANGE -- solves the LCP
        :param V: <np.array> matrix of the LCP
        :param p: <np.array> vector of the LCP
        :return: renders the trajectory
    """
    sol = LCPSolve(V, p)
    f_r = sol[1][:3]
    return f_r


def render(q):
    """
        DO NOT CHANGE -- renders the trajectory
        :param q: <np.array> configuration trajectory
        :return: renders the trajectory
    """
    assignment_3_render(q)


if __name__ == "__main__":
   
    # to test your final code, use the following initial configs
    q0 = np.array([0.0, 1.5, np.pi / 180. * 30.])
    v0 = np.array([0.,-0.2,0.])#v0 = np.array([0., -0.2, 0.])

    q, v = simulate(q0, v0)
    plt.plot(q[1, :])
    plt.show()

    render(q)
   

    