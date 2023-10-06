# PySCRDT
A module to calculate the resonance driving terms (RDTs) from the space charge potential, from Twiss tables generated with MAD-X (see [MAD-X documentation](https://mad.web.cern.ch/mad/)) or from X-suite Twiss dictionaries (see [Xsuite documentation](https://xsuite.readthedocs.io/en/latest/). 

Documentation in: https://cds.cern.ch/record/2696190/files/CERN-ACC-NOTE-2019-0046.pdf

For the analytical space charge tune spread estimation: https://github.com/fasvesta/tune-spread

Can be installed with Python environment with `python -m pip install git+https://github.com/fasvesta/PySCRDT.git`.

Contact: `foteini.asvesta@cern.ch` and `elias.walter.waagaard@cern.ch`

### Loading the beam parameters 

PySCRDT either takes MAD-X Twiss tables or X-suite Twiss dictionaries as input. An Xtrack Twiss dictionary can be easily generated from an Xsuite line (from a json file).

```
import xtrack as xt
import json 

with open('some_xsuite_line.json', 'r') as fid:
     input_data = json.load(fid)    
line = xt.Line.from_dict(input_data)
particle_ref = line.particle_ref
line.build_tracker()
twiss_xtrack = line.twiss() 
```

### Initializing the PySCRDT class object 

The PySCRDT object requires beam parameters such as bunch length, bunch intensity, normalized emittances, RMS momentum spread, bunching factor and classical particle radius. 

```
from PySCRDT import PySCRDT

s = PySCRDT()
s.setParameters(
    intensity = bunch_intensity,
    bunchLength = sigma_z,
    emittance_x = nemitt_x,
    emittance_y = nemitt_y, 
    dpp_rms = dpp_rms, 
    bF=None,
    ro = particle_ref.get_classical_particle_radius0() 
)
```
Loading a Twiss dictionary from Xtrack:
```
s.loadTwissFromXsuite(twissTableXsuite=twiss_xtrack)
```
Loading a Twiss table from MAD-X:
```
s.prepareData(twissFile=input_parameters.twiss_file)
```

### Initiate the tune footprint maker

```
from PySCRDT import tune_footprint_maker

Qh = twiss_xtrack['qx']
Qv = twiss_xtrack['qy']
my_tune_footprint = tune_footprint_maker(Qh, Qv)
```

#### Calculating the maximum detuning
```
dQx, dQy = my_tune_footprint.return_max_detuning(s)
```

#### Plot the tune footprint

The plotting parameters such as plotting range for the tune footprint figure - and the plot order and periodicity for the resonance lines - have default values, but can be changed in the constructor of the `tune_footprint_maker` class.

A `matplotlib.pyplot.figure` object is passed to the tune footprint maker. The `axis` object can be returned if `return_fig_axis=True` (default=`False`) to plot additional data in the figure. Two examples for ther PS and SPS ions are given in the folder `examples`.

```
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(8, 8))
my_tune_footprint.generate_tune_footprint(fig, s)
```
![SPS_Pb_tune_footprint](https://github.com/ewaagaard/PySCRDT/assets/68541324/db5e2e45-4d42-4d42-bacb-8045a1c09e72)

