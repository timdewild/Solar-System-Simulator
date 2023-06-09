# Solar-System-Simulator
Python module to simulate and animate the solar system in real time. An example of the animations that can be generated with the code is given below:

https://github.com/timdewild/Solar-System-Simulator/assets/93600756/b6287190-0d8b-46c2-87b6-47cd2a328b51

The notebook `example.ipynb` explains in detail how this is done. You can also download the video [here](inner_solar_system.mp4). 

## Author

Tim de Wild

## Modules required
The modules required are:
```bash
numpy
matplotlib
datetime
astropy
astroquery
```

## How to run
See the example notebook `example.ipynb`. Create a Jupyter notebook in the same directory as where `main.py` is located and import the `Body` and `SolarSystem` classes. Using these classes, the animations can be generated as explained in the example notebook.

## The Solar System as a Gravitational $N$-body System
The solar system is essentially an $N$-body system of point-particles that interact gravitationally according to Newton's law of gravitation. For two bodies $i$ and $j$, separated by $\boldsymbol{r}_{ij}\equiv \boldsymbol{r}_i-\boldsymbol{r}_j$, the force on $i$ due to $j$ is:

```math
\boldsymbol{F}_{ij} = -\frac{Gm_im_j}{r_{ij}^2}\hat{\boldsymbol{r}}_{ij}.
```

The total force on $i$ is obtained by taking the sum over all bodies $j$ except itself:

```math
\boldsymbol{F}_{i} = \sum_{j\neq i}-\frac{Gm_im_j}{r_{ij}^2}\hat{\boldsymbol{r}}_{ij}.
```

The acceleration of body $i$ is obtained by dividing by $m_i$. In celestial mechanics, it is customary to define the gravitational parameter as $\mu_j\equiv Gm_j$, so that the acceleration becomes:

```math
\boldsymbol{a}_i = \ddot{\boldsymbol{r}}_{ij} = \sum_{j\neq i} -\mu_j \frac{\hat{\boldsymbol{r}}_{ij}}{r_{ij}^2}. 
```

To find the orbits of the $N$ bodies, we have to solve this second order differential equation for each of the bodies simultaneously. We will do so using a numerical integration scheme. We extract the initial position and velocity of all the bodies from NASA's JPL Horizons database via the mehod  `get_init_cond()` in the `Body` class. We work in the rest frame of the Sun, so that the sun stays at rest at all times. 

## Integration Schemes
In the code, we have implemented two numerical integration schemes in the `integrate()` method in the `SolarSystem` class: the Forward-Euler scheme and the Runge-Kutta scheme. The former is known for its simplicity, at the expense of accuracy. The latter is known for its accuracy, at the expense of simplicity. Both rely on the descretization of time and relating the state (i.e. velocity and position) of a body at the next timestep to those at previous timesteps. 

### Forward Euler
The Forward Euler integration scheme can be selected via the `method = 'FE'` keyword of `integrate()`. This scheme is given by:
```math
\begin{align}
   \boldsymbol{v}_i^{n+1} & = \boldsymbol{v}_i^{n}+\boldsymbol{a}_i^{n}\;dt, \\
   \boldsymbol{r}_i^{n+1} & = \boldsymbol{r}_i^{n}+\boldsymbol{v}_i^{n}\;dt.
\end{align}
```
The index $i$ labels the bodies, the index $n$ labels the time, which is descritized in steps of $dt$. The acceleration of body $i$ depends solely on its own position and that of the other bodies:

```math
\boldsymbol{a}_i^{n} = \boldsymbol{a}_i(\boldsymbol{r}_i^n, \boldsymbol{r}_j^n). 
```



## References

- [JPL HORIZONS on-line solar system database](https://docs.astropy.org/en/stable/coordinates/solarsystem.html)
