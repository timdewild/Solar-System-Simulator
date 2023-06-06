import numpy as np
import matplotlib.pyplot as plt
import datetime
import matplotlib.animation as animation
import matplotlib
from astropy.time import Time as time_as
from astroquery.jplhorizons import Horizons as horizons_JPL
matplotlib.rcParams['font.family'] = 'verdana'

class Body:
    """Class that represents bodies (planets, sun) in the solar system."""

    # some data of sun and plantes in solar system
    # grav param mu in AU^3/day^2
    # id is nasa id
    # period T in days
     
    data = {
        'Sun':      {'mu': 0.295912208285591100e-3, 'id': '0', 'T': 0           },
        'Mercury':  {'mu': 0.491248045036476e-10,   'id': '1', 'T': 88.0        },
        'Venus':    {'mu': 0.724345233e-9,          'id': '2', 'T': 224.7       },
        'Earth':    {'mu': 0.888769244512563400e-9, 'id': '3', 'T': 365.2       },
        'Mars':     {'mu': 0.95495486e-10,          'id': '4', 'T': 687.0       },
        'Jupiter':  {'mu': 0.282534584e-06,         'id': '5', 'T': 4331        }, 
        'Saturn':   {'mu': 0.845970607e-07,         'id': '6', 'T': 10755.698   },
        'Uranus':   {'mu': 0.129202482e-07,         'id': '7', 'T': 30685.4     },
        'Neptune':  {'mu': 0.152435734e-07,         'id': '8', 'T': 60189       }
    }

    def __init__(self, name):
        """Initializes body object.

        Args:
            name (str): Name of body: 'Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus' or 'Neptune'
        """

        self.name = name
        self.r = np.zeros(3)
        self.v = np.zeros(3)        
        self.mu = Body.data[self.name]['mu']
        self.id = Body.data[self.name]['id']
        
        #sidereal period body
        self.T = Body.data[self.name]['T']

        # for storing state vectors
        self.orbit = None
        self.vel = None

        # for animation
        self.dot = None
        self.trace = None
        self.traj = None

    def init_containers(self, N_time_steps):
        """Initialize containers for orbit (x,y,z) values at all timesteps and vel (vx, vy, vz) values at all timesteps. 
        Sets orbit and vel to 2d Numpy arrays with shape (N_time_steps, 3). 

        Args: 
            N_time_steps (int): number of timesteps in integration
        """
        self.orbit = np.zeros((N_time_steps, 3))
        self.vel = np.zeros((N_time_steps, 3))

    def get_init_cond(self, start_date):
        """Get initial conditions of body r0 = (x,y,z) and v0 = (vx, vy, vz), units AU and AU/day.
        Accesses NASA's Horizons database using astroquery.jplhorizons 

        Args:
            start_date (str): starting date of integration, format 'YYYY-MM-DD'
        """

        if self.name == 'Sun':
            return

        vec = horizons_JPL(
                    id=self.id, 
                    location='@sun',
                    epochs = time_as(start_date).jd,
        ).vectors()

        r0 = np.array(list(vec['x','y','z'][0]))
        v0 = np.array(list(vec['vx','vy','vz'][0]))

        self.set_r(r0)
        self.set_v(v0)

    def get_r(self):
        """Get Numpy array with current position
        
        Returns:
            self.r (Numpy array): vector containing current position of self

        """
        return self.r

    def set_r(self, r):
        """Set Numpy array with current position"""
        self.r = r

    def get_v(self):
        """Get Numpy array with current velocity
        
        Returns:
            self.v (Numpy array): vector containing current velocity of self

        """
        return self.v

    def set_v(self, v):
        """Set Numpy array with current velocity"""
        self.v = v

    def get_acceleration(self, r, bodies):
        """Get acceleration on self at position r due to bodies in list bodies. 

        Args: 
            r (Numpy array): position to calculate acceleration at
            bodies (list): list with Body objects

        Returns:
            a (Numpy array): return acceleration (ax, ay, az) due to all bodies in bodies except self
        """

        a = np.zeros(3)
        for body in bodies:
            if self.name != body.name: 

                dr = r - body.r
                dr_mag = np.sqrt(dr.dot(dr))

                a += - body.mu * dr / dr_mag ** 3
        
        return a
    
    def get_orbit(self):
        """Get Numpy array with orbit of self
        
        Returns:
            self.orbit (Numpy array): vector containing orbit of self at all timesteps

        """

        return self.orbit
    
    def get_vel(self):
        """Get Numpy array with velocity of self
        
        Returns:
            self.vel (Numpy array): vector containing velocity of self at all timesteps

        """

        return self.vel
    
    
class SolarSystem:
    """Class that represents the solar system."""

    def __init__(self, tmax, dt, start_date):
        """Initializes solar system.

        Args:
            tmax (float): timespan of integration in days
            dt (float): timestep in days
            start_data (str): starting date of integration 'YYYY-MM-DD'
        
        """


        self.tmax = tmax
        self.dt = dt

        # start_date has '%Y-%m-%d' string format
        self.start_date = start_date

        # initilize list with bodies
        self.bodies = []

        # initialize number of timesteps in integrations
        self.N_time_steps = int(self.tmax/self.dt) + 1

        # generate list with dates
        self.date_list = []

        # create datetime.datetime object from the start_date input
        date = datetime.datetime.strptime(self.start_date, '%Y-%m-%d')
        for i in range(self.N_time_steps):
            self.date_list.append(date)
            delta = datetime.timedelta(days = self.dt)
            date += delta 

    def get_date_list(self):
        """Returns date_list"""
        return self.date_list

    def add_body(self, body):
        """Adds body to solar system, gets initial conditions of body from Horizons JPL at start_date. 

        Args: 
            body (Body instance): body to add to solar system
        
        """

        self.bodies.append(body)
        body.get_init_cond(self.start_date)

    def integrate(self, method = 'RK'):
        """
        Runge-Kutta/Forward-Euler integrator for solar system dynamics, calculates orbit and velocity at all supplied times. 

        Args:
            method (str): Keyword to specify integrator, Forward Euler (FE) or Runge-Kutta (RK)
        """

        for body in self.bodies:
            body.init_containers(self.N_time_steps)

            body.orbit[0,:] = body.r
            body.vel[0,:] = body.v

        if method == 'FE':
            for i in range(1,self.N_time_steps):
                for body in self.bodies:

                    # work in rest frame of the sun, do not update sun
                    if body.name == 'Sun':
                        continue
                    
                    a = body.get_acceleration(body.r, self.bodies)

                    v = body.v + a * self.dt
                    r = body.r + body.v * self.dt

                    body.set_v(v)
                    body.set_r(r)

                    body.vel[i,:] = v 
                    body.orbit[i,:] = r

        if method == 'RK':
            for i in range(1,self.N_time_steps):
                for body in self.bodies:

                    # work in rest frame of the sun, do not update sun
                    if body.name == 'Sun':
                        continue
                    
                    k_1r = self.dt * body.v
                    k_1v = self.dt * body.get_acceleration(body.r, self.bodies)

                    k_2r = self.dt * (body.v + 0.5 * k_1v)
                    k_2v = self.dt * body.get_acceleration(body.r + 0.5 * k_1r, self.bodies)

                    k_3r = self.dt * (body.v + 0.5 * k_2v)
                    k_3v = self.dt * body.get_acceleration(body.r + 0.5 * k_2r, self.bodies)

                    k_4r = self.dt * (body.v + k_3v)
                    k_4v = self.dt * body.get_acceleration(body.r + k_3r, self.bodies)

                    v = body.v + 1./6. * (k_1v + 2*k_2v + 2*k_3v + k_4v)
                    r = body.r + 1./6. * (k_1r + 2*k_2r + 2*k_3r + k_4r)

                    body.set_v(v)
                    body.set_r(r)

                    body.vel[i,:] = v 
                    body.orbit[i,:] = r

    def system_energy(self):

        K = np.zeros(self.N_time_steps)
        V = np.zeros(self.N_time_steps)

        # define list with all possible pairs
        pairs = [[body_i, body_j] for body_i in self.bodies for body_j in self.bodies if body_i is not body_j]
        
        for pair in pairs:
            b1, b2 = pair
            V += - b1.mu * b2.mu / np.sqrt( np.sum((b1.get_orbit() - b2.get_orbit())**2, axis=1) )

        for body in self.bodies:
            K += 0.5*body.mu * np.sum(body.get_vel()**2, axis=1)

        E  = K + V

        return K, V, E
    
    def plot_orbits(self):
        """Plot orbits of all bodies in solar system"""

        fig, ax = plt.subplots(figsize = (8,8))
        ax.set_aspect('equal')
        ax.scatter([0],[0], color='k')

        for body in self.bodies:
            orbit = body.get_orbit()
            ax.plot(orbit[:,0], orbit[:,1], lw = 1)
    
    def animate(self, size, fps, file_name):
        """Make Matplotlib animation of integrated orbits. 

        Args:
            size (float): sets bounds of canvas in AU: [xmin, xmax, ymin, ymax] = [-size, +size, -size, +size] 
            fps (int): sets frame rate (per second) of animation
            file_name: sets file_name of animation video

        """

        viz_bodies = {
            'Sun':      {'c': "#feb125", 'ms': 25,  'lms': 1000 },
            'Mercury':  {'c': "#b6a864", 'ms': 9,   'lms': 50   },
            'Venus':    {'c': "#d15127", 'ms': 14,  'lms': 300  },
            'Earth':    {'c': "#b1d6f0", 'ms': 15,  'lms': 300  },
            'Mars':     {'c': "#ce2756", 'ms': 11,  'lms': 100  },
            'Jupiter':  {'c': "#c6b7b3", 'ms': 20,  'lms': 600  }, 
            'Saturn':   {'c': "#954063", 'ms': 20,  'lms': 600  }, 
            'Uranus':   {'c': "#fd6162", 'ms': 20,  'lms': 600  },
            'Neptune':  {'c': "#8ca8cf", 'ms': 20,  'lms': 600  }
        }

        # background color
        c_bkg = "#10112d"

        # color for orbits
        c_orb = "#cbcbcb"

        fig, ax = plt.subplots(figsize = (12,12), constrained_layout = True)
        fig.set_facecolor(c_bkg)
        fig.set_dpi(100)
        ax.set_aspect('equal')
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_facecolor(c_bkg)
        ax.set_frame_on(False)
        ax.set_xlim(-size,size)
        ax.set_ylim(-size,size)

        timestamp = ax.text(.03, .97, 'Date: ', color='w', transform=plt.gca().transAxes, fontsize='x-large')

        for body in self.bodies:
            body.dot, = ax.plot([], [], marker='.', markersize = viz_bodies[body.name]['ms'], color = viz_bodies[body.name]['c'], zorder = 3)
            body.trace, = ax.plot([],[], color = viz_bodies[body.name]['c'], lw = 2)
            body.traj, = ax.plot([],[], color = c_orb, lw = 0.5)

            # for sole purpose of making plotting labels with circles only
            ax.scatter([],[], marker='.', color = viz_bodies[body.name]['c'], label = body.name) 
        
            # if body orbital period is larger than simulated timespan, plot all of (partial) orbit
            if body.T > self.tmax: 
                X, Y = body.get_orbit()[:,0], body.get_orbit()[:,1]
                body.traj.set_data(X,Y)
            
            else:
                n = int(body.T/self.dt) + 1
                X, Y = body.get_orbit()[:n,0], body.get_orbit()[:n,1]
                body.traj.set_data(X,Y)

        # define legend
        lgnd = ax.legend(loc="lower center", ncol = len(self.bodies))

        # remove legend bkg
        lgnd.get_frame().set_alpha(0)

        # legend text color white
        for text in lgnd.get_texts():
            text.set_color("w")

        # set legend handle size 
        for body, handle in zip(self.bodies, lgnd.legendHandles):
            handle._sizes = [viz_bodies[body.name]['lms']]

        def animate_frame(i):

            for body in self.bodies:
                x, y = body.get_orbit()[i,0], body.get_orbit()[i,1] 
                body.dot.set_data(x, y)

            date = self.date_list[i]
            timestamp.set_text('Date: {}'.format(date.strftime('%Y-%m-%d')))


        ani = animation.FuncAnimation(fig, animate_frame, frames=self.N_time_steps, interval=20)
        ani.save(file_name+'.mp4', dpi=150, fps=fps)
