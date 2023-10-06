"""
#
#   Tune footprint maker based on PySCRDT module
#   - author: Elias Waagaard 
#   - adopted from fasvesta/PySCRDT and fasvesta/tune-spread into one single pip-installable package
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import colors
from .resonance_lines import resonance_lines

SMALL_SIZE = 18
MEDIUM_SIZE = 23
BIGGER_SIZE = 28
plt.rcParams["font.family"] = "serif"
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)   # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)   # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)   # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

class tune_footprint_maker(object):
    """
    Class to generate tune footprints 
    """ 
    def __init__(self,
                 Qh,
                 Qv,
                 plot_range = None,
                 plot_order = 5,
                 periodicity = 16,
                 ):
        self.Qh = Qh
        self.Qv = Qv
        if plot_range is None:
            self.plot_range = [[Qh - 0.3, Qh + 0.1], [Qh - 0.3, Qh + 0.1]]   # range in Qh & Qv for the plot
        else:
            self.plot_range = plot_range
        self.plot_order =  5   # order of resonances to plot
        self.periodicity = periodicity  # periodicity of ring for the colorcode of the plot
        self.detuning_is_calculated = False
        self.p_collection_exists = False 
        
        
    def calculate_detuning_coefficients(self, PySCRDT_object):
        """
        Caclulate detuning coefficients using potentials up to 20th order (needed for up to 3 sigma particles) from PySCRDT class
        """
        detuning=[]
        # order in x
        for i in range(0,21,2):
            # order in y
            for j in range(0,21,2): 
                if (i==0) and (j==0):
                    pass
                elif i+j<21:
                    PySCRDT_object.setOrder([int(i),int(j),'any'])
                    PySCRDT_object.potential()
                    PySCRDT_object.detuning()
                    detuning.append([i,j, PySCRDT_object.getDetuning()])
            
        self.detuning = np.array(detuning)
        self.detuning_is_calculated = True 
    
    
    def return_max_detuning(self, PySCRDT_object):
        """
        Return maximum tune shift 

        Parameters
        ----------
        PySCRDT_object 

        Returns
        -------
        dQx : horizontal maximum detuning 
        dQy : vertical maximum detuning 
        """
        if not self.detuning_is_calculated:
            self.calculate_detuning_coefficients(PySCRDT_object)
    
        detuning_y=[i for i in np.where(self.detuning[:,1]==2)[0] if i in np.where(self.detuning[:,0]==0)[0]][0]
        detuning_x=[i for i in np.where(self.detuning[:,0]==2)[0] if i in np.where(self.detuning[:,1]==0)[0]][0]
        
        dQx = self.detuning[detuning_x, 2]
        dQy = self.detuning[detuning_y, 2]
        
        return dQx, dQy
    

    def initial_xy_polar(self, s_max, s_N, theta_N):
        """
        Generate a grid of polar coordinates 

        Parameters
        ----------
        s_max : maximum value of the radial coordinate
        s_N : number of subdivisions for the radial coordinate
        theta_N : number of subdivisions for the angular coordinate

        Returns
        -------
        numpy array with tuples representing the Cartesian coordinates from the grid of polar coordinates 
        """
        
        return np.array(
            [
                [(s*np.cos(theta), s*np.sin(theta)) for s in np.linspace(0, s_max, s_N+1)]
                for theta in np.linspace(0, np.pi/2., theta_N)
            ])


    def make_polygons_for_footprint(self, PySCRDT_object, s_N=6, s_max=3, theta_N=5):
        """
        Generate polygons for plotting using matplotlib

        Parameters
        ----------
        PySCRDT_object
        s_max : maximum value of the radial coordinate. The default is 6.
        s_N : number of subdivisions for the radial coordinate. The default is 3.
        theta_N : number of subdivisions for the angular coordinate. The default is 5.
        
        Returns
        -------
        None
        """
        if not self.detuning_is_calculated:
            self.calculate_detuning_coefficients(PySCRDT_object)

        #  initialize grid for calculation
        S = self.initial_xy_polar(s_max=s_max, s_N=s_N, theta_N=theta_N)
        
        # Beam parameters and action-angle variables 
        en_x = PySCRDT_object.parameters['emittance_x']
        en_y = PySCRDT_object.parameters['emittance_y']
        beta = PySCRDT_object.parameters['b']
        gamma = PySCRDT_object.parameters['g']
        J_x = S[:,:,0]**2*en_x/2./beta/gamma
        J_y = S[:,:,1]**2*en_y/2./beta/gamma
        
        Qx, Qy = self.Qh, self.Qv 
        
        # Updated tunes from calculated detuning coefficients 
        for x_q, y_q, detuning_coef in self.detuning:
            if x_q:
                Qx+=x_q/2.*detuning_coef*(J_x**(x_q/2.-1))*(J_y**(y_q/2.))
            if y_q:
                Qy+=y_q/2.*detuning_coef*(J_y**(y_q/2.-1))*(J_x**(x_q/2.))
        
        Q = np.dstack(
            (
                [qx.tolist() + [self.Qh] for qx in Qx],
                [qy.tolist() + [self.Qv] for qy in Qy],
            )
        )
        
        Q[:,:,0] += 0.00
        Q[:,:,1] += 0.00
        
        # reshapes the Q array to form four sets of points, p1, p2, p3, and p4, which represent the vertices of the polygons.
        sx = Q.shape[0]-1
        sy = Q.shape[1]-1
        p1 = Q[:-1, :-1, :].reshape(sx*sy, 2)[:, :]
        p2 = Q[1:, :-1, :].reshape(sx*sy, 2)[:]
        p3 = Q[1:, 1:, :].reshape(sx*sy, 2)[:]
        p4 = Q[:-1, 1:, :].reshape(sx*sy, 2)[:]
        
        # colormap for plotting 
        cmap_base = plt.cm.hot
        c_indcs = np.int_(np.linspace(0.1,0.6,s_N+1)*cmap_base.N)
        cmap = colors.ListedColormap([cmap_base(c_indx) for c_indx in c_indcs])
        
        # Stack endpoints to form polygons
        Polygons = np.transpose(np.stack((p1, p2, p3, p4)), (1, 0, 2))
        patches = list(map(matplotlib.patches.Polygon, Polygons))
        self.p_collection = matplotlib.collections.PatchCollection(
            patches, edgecolor='k', linewidth=0.5,
            facecolors=cmap.colors,
            alpha=0.7
        )
        self.p_collection_exists = True 
        
        
    def generate_tune_footprint(self, fig, PySCRDT_object, return_fig_axis=False):
        """
        Generate tune footprint figure from PySCRDT 

        Parameters
        ----------
        fig : matplotlib figure object
        PySCRDT_object 
	Qx: horizontal tune data from particle tracking. Default: None
	Qy: vertical tune data from particle tracking. Default: None
        return_fig_axis : default is False.

        Returns
        -------
        ax : figure object axis for additional plotting, if desirede 
        """
    
        if not self.p_collection_exists:
            self.make_polygons_for_footprint(PySCRDT_object)
       
        # Then also plot this tune footprint with the analytical tune footprint
        ax = fig.add_subplot(1, 1, 1) 
        tune_diagram = resonance_lines(self.plot_range[0],
                    self.plot_range[1], np.arange(1, self.plot_order+1), self.periodicity)
        tune_diagram.plot_resonance(figure_object=fig, interactive=False)
        ax.plot(self.Qh, self.Qv, 'ro', markersize=10.5, alpha=0.8, label='Set tune')
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        ax.add_collection(self.p_collection)
        ax.set_aspect('equal')
        ax.legend(loc = "lower right")
        ax.set_ylabel(r'$Q_{y}$', fontsize=BIGGER_SIZE)
        ax.set_xlabel(r'$Q_{x}$', fontsize=BIGGER_SIZE)
        fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        #plt.show()
        
        if return_fig_axis:
            return ax


    def generate_tune_footprint_with_tracking_data(self, fig, PySCRDT_object, Qx, Qy):
        """
        Generate tune footprint figure from PySCRDT, also including tune data from tracking 

        Parameters
        ----------
        fig : matplotlib figure object
        PySCRDT_object 
	Qx: horizontal tune data from particle tracking. Default: None
	Qy: vertical tune data from particle tracking. Default: None
        return_fig_axis : default is False.

        Returns
        -------
        None
        """
    
        if not self.p_collection_exists:
            self.make_polygons_for_footprint(PySCRDT_object)
       
        # Generate analytical tune footprint
        ax = self.generate_tune_footprint(fig, PySCRDT_object, return_fig_axis=True)

        # Add tune footprint from tracking
        ax.plot(Qx + int(self.Qh), Qy + int(self.Qv), 'go', markersize=4.5, alpha=0.2) 

        plt.show()

