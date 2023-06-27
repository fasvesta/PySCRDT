"""
#  Simple tune footprint maker example with PySCRDT to generate tune footprint
#  with Pb ions in the SPS
#  example by Elias Waagaard 
"""
import matplotlib.pyplot as plt
import xtrack as xt
import json 

from PySCRDT import PySCRDT, tune_footprint_maker

#### Load the sequence ####
fname_line = 'example_data/SPS_2021_Pb_ions_for_tracking.json'
with open(fname_line, 'r') as fid:
     input_data = json.load(fid)    
line = xt.Line.from_dict(input_data)
particle_ref = line.particle_ref
line.build_tracker()
twiss_xtrack = line.twiss() 

#### PS default Pb beam settings #### 
bunch_intensity = 3.5e8 
sigma_z = 0.23
nemitt_x= 1.3e-6
nemitt_y= 0.9e-6 

#### Analytical tune footprint settings ###
Qh = twiss_xtrack['qx']
Qv = twiss_xtrack['qy']
plot_range  =   [[Qh - 0.4, Qh + 0.1], [Qh - 0.4, Qh + 0.1]]   # plot range in Qh & Qv, can be custom-provided to tune_footprint_maker

#### Create instance of PySCRDT, taking beam parameters as input ####
s = PySCRDT()
s.setParameters(
    intensity = bunch_intensity,
    bunchLength = sigma_z,
    emittance_x = nemitt_x,
    emittance_y = nemitt_y, 
    dpp_rms = 1e-3, 
    bF=None,
    ro = particle_ref.get_classical_particle_radius0() 
)
s.loadTwissFromXsuite(twissTableXsuite=twiss_xtrack)

#### Initiate tune footprint maker ####
SPS_tune_footprint = tune_footprint_maker(Qh, Qv, plot_range = plot_range)

#### Print maximum tune shift #### 
dQx, dQy = SPS_tune_footprint.return_max_detuning(s)
print("\ndQx = {:.4f}, dQy = {:.4f}".format(dQx, dQy))

#### Plot the tune footprint ####
fig = plt.figure(figsize=(8, 8))
SPS_tune_footprint.generate_tune_footprint(fig, s)
fig.savefig('output/SPS_Pb_tune_footprint.png', dpi=250)