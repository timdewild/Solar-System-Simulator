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
With the required modules installed, simply run the program, e.g. in UNIX command line, with
```bash
python solar_system.py
```
Modify `sim_duration` in the beginning of the code to change the duration of the simulation.

To make a different initial condition, modify and run get_initial_condition.py.

## References

- [JPL HORIZONS on-line solar system database](https://docs.astropy.org/en/stable/coordinates/solarsystem.html)
