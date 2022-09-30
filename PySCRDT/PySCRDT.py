# - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * 
"""
#
#   PySCRDT
#   A module to calculate the resonance driving terms from the space charge potential
#   
#   Version :   1.1 
#               pre-calculated potentials for faster evaluation
#   Author  : F. Asvesta
#   Contact : fasvesta .at. cern .dot. ch
#
"""
# - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *

from __future__ import absolute_import, division, print_function, unicode_literals
try:
    import numpy as np
except ImportError:
    print("# PySCRDT : numpy module is required. ")
try:
    import sympy as sy
except ImportError:
    print("# PySCRDT : sympy module is required. ")
try:
    import dill
except ImportError:
    print("# PySCRDT : dill module is required. ")


__version   = 1.1
__PyVersion = ["2.7"]
__author    = ["Foteini Asvesta"]
__contact   = ["fasvesta .at. cern .dot. ch"]

class PySCRDT(object):
    """
    class for the calculation of the rdts
    Returns: PySCRDT instance
    """

    def __init__(self, parameters=False, mode=None, twissFile=None, order=None):
        """
        Initialization function
        Input :  parameters : [bool|str]  Parameters needed for the calculations (default=False)
                                          if True the default values of [setParameters] are used 
                                          if str parameters are read in a file
                 mode       : [3|5]       Resonance description mode (default=None)
                 order      : [list]      Resonance order and harmonic (default=None)
                 twissFile  : [str]       MAD-X Twiss file (default=None)
        Returns: void
        """
        self.x, self.y, self.t= sy.symbols('x y t')
        self.a = sy.Symbol('a', positive=True, real=True)
        self.b = sy.Symbol('b', positive=True, real=True)
        self.D = sy.Symbol('D', positive=True, real=True)
        self.fx = sy.Symbol('fx', positive=True, real=True)
        self.fy = sy.Symbol('fy', positive=True, real=True)
        self.V=None
        self.K=None
        self.data=None
        self.factor=None
        self.factor_d=None
        self.rdt=None
        self.rdt_d=None
        self.feed=False
        self.mode=None
        self.order=None
        if type(parameters) is str:
            self.parameters=None
            self.readParameters(parameters)
        else:
            if parameters:
                self.setParameters()
            else:
                self.parameters=None
                print("# PySCRDT : Set parameters with function [setParameters] or read parameters with [readParameters]")
        if twissFile is None:
            print("# PySCRDT : Import Madx twiss file using the function [prepareData]")
        else:
            self.prepareData(twissFile)            
        if (order is None) and (mode is None):
            print("# PySCRDT : Set order in [setOrder]")
        elif order is None:
            print("# PySCRDT : Set order in [setOrder]")
        else:
            if mode is None:
                self.setMode(len(order))
            else:
                self.setMode(mode)
            self.setOrder(order)            
        
    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *
    
    def setMode(self, mode):
        """
        Sets the mode for the characterization of resonances
        Input :  mode  : [3|5] 
        Returns: void
        """
        if mode in [3,5]:
            self.mode=mode
        elif mode is None:
            print("# PySCRDT : Set mode in [setMode]")
        else:
            raise IOError('# PySCRDT::setMode: Mode needs to be 3 or 5 depending on the desired resonance description')
            
    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *
    
    def setOrder(self, args):
        """
        Sets the Resonance Orders
        In "3 mode"
        Input :  m  : [int] order in H
                 n  : [int] order in V
                 l  : [int|'any'] harmonic of the resonance 
        In "5 mode"
        Input :  h  : [int] characteristic of H
                 i  : [int] characteristic of H
                 j  : [int] characteristic of V
                 k  : [int] characteristic of V
                 l  : [int|'any'] harmonic of the resonance 
        Returns: void
        """
        if len(args)==3:
            if self.mode==3 or self.mode==None:
                if type(args[0]) is int:
                    self.m=args[0]
                else:
                    raise IOError('# PySCRDT::setOrder: resonance order needs to be of type int')
                if type(args[1]) is int:
                    self.n=args[1]
                else:
                    raise IOError('# PySCRDT::setOrder: resonance order needs to be of type int')
                if (type(args[2]) is int) or (args[2]=='any'):
                    self.l=args[2]
                else:
                    raise IOError('# PySCRDT::setOrder: harmonic needs to be of type int or set to "any"')
                self.mode=3
            else:
                raise IOError('# PySCRDT::setOrder: You need to define the order using h,i,j,k,l')
        elif len(args)==5:
            if self.mode==5 or self.mode==None:
                if type(args[0]) is int:
                    self.h=args[0]
                else:
                    raise IOError('# PySCRDT::setOrder: resonance order needs to be of type int')
                if type(args[1]) is int:
                    self.i=args[1]
                else:
                    raise IOError('# PySCRDT::setOrder: resonance order needs to be of type int')
                if type(args[2]) is int:
                    self.j=args[2]
                else:
                    raise IOError('# PySCRDT::setOrder: resonance order needs to be of type int')
                if type(args[3]) is int:
                    self.k=args[3]
                else:
                    raise IOError('# PySCRDT::setOrder: resonance order needs to be of type int')
                if (type(args[4]) is int) or (args[4]=='any'):
                    self.l=args[4]
                else:
                    raise IOError('# PySCRDT::setOrder: harmonic needs to be of type int or set to "any"')
                self.mode=5
            else:
                raise IOError('# PySCRDT::setOrder: You need to define the order using m,n,l')
        else:
            if self.mode==3:
                raise IOError('# PySCRDT::setOrder: You need to define the order using m,n,l')
            if self.mode==5:
                raise IOError('# PySCRDT::setOrder: You need to define the order using h,i,j,k,l')
        self.factor=None
        self.factor_d=None
    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *

    def setParameters(self, intensity=41e10, bunchLength=5.96, ro=1.5347e-18, emittance_x=2e-6, emittance_y=1.1e-6, dpp_rms=0.5e-3, dpp=0.0, bF=None, harmonic=1):
        """
        Sets the parameters for the calculation:
        Input :  intensity  : [float] bunch intensity in ppb (Default=41e10)
                 bunchLength: [float] RMS bunch length in m  (Default=5.96)
                 ro         : [float] classical particle radious in m (Default=1.5347e-18 {proton})
                 emittance_x: [float] normalized horizontal emittance in m*rad (Default=2e-6)
                 emittance_y: [float] normalized vertical emittance in m*rad (Default=1.1e-6)
                 dpp_rms    : [float] RMS Dp/p (Default=0.5e-3)
                 dpp        : [float] Single particle Dp/p (Default=0)
                 bF         : [float] Bunching factor (Default=None)
                 harmonic   : [int]   Harmonic number,# of buckets (Default=1)
        Returns: void
        """
        self.parameters={'intensity':intensity, 
                         'bunchLength':bunchLength, 
                         'ro':ro, 
                         'emittance_x':emittance_x, 
                         'emittance_y':emittance_y,
                         'dpp_rms':dpp_rms,
                         'dpp':dpp,
                         'bF':bF,
                         'harmonic':harmonic}
        
    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *
    def readParameters(self, inputFile):
        """
        Reads the parameters from an input file:
        Input Format: intensity   = [float] # bunch intensity in ppb (Default=41e10)
                      bunchLength = [float] # RMS bunch length in m  (Default=5.96)
                      ro          = [float] # classical particle radious in m (Default=1.5347e-18 {proton})
                      emittance_x = [float] # normalized horizontal emittance in m*rad (Default=2e-6)
                      emittance_y = [float] # normalized vertical emittance in m*rad (Default=1.1e-6)
                      dpp_rms     = [float] # RMS Dp/p (Default=0.5e-3)
                      dpp         = [float] # Single particle Dp/p (Default=0.5e-3)
                      bF          = [float] # Bunching factor (Default=None)
                      harmonic    = [int]   # Harmonic number,# of buckets (Default=1)
        Returns: void
        """
        params=np.genfromtxt(inputFile,dtype=str)
        if self.parameters is None:
            self.setParameters()
        if len(np.shape(params))==1:
            if params[0] not in self.parameters.keys():
                raise IOError('# PySCRDT::readParameters: '+ params[0] +' not recognized [checkWriting]')
            else:
                self.parameters[params[0]]=float(params[2])
        else:
            for i in enumerate(params):
                if params[i[0]][0] not in self.parameters.keys():
                    raise IOError('# PySCRDT::readParameters: '+ params[i[0]][0] +' not recognized [checkWriting]') 
                else:
                    self.parameters[params[i[0]][0]]=float(params[i[0]][2])
        if self.data is not None:
            self.beamSize()
            self.ksc()
            
        
    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *
    
    def potential(self, feedDown=False):
        """
        Calculates the space charge potential for the given resonance order
        Inputs : feedDown : [bool] needed when single particle Dp/p non 0 (default=False)
        Returns: Void
        """
        if self.mode==5:
            self.m=self.h+self.i
            self.n=self.j+self.k
        if (self.m is None) or (self.n is None):
            raise IOError('# PySCRDT::potential: You need to define resonance order in [setOrder]')
        if self.m%2!=0 and feedDown==False:
            raise IOError('# PySCRDT::potential: Space charge potential contains only even orders without Dp/p (order given {}), change the order in [setOrder]'.format(str(self.m)))
        if self.n%2!=0:
            raise IOError('# PySCRDT::potential: Space charge potential contains only even orders (order given {}), change the order in [setOrder]'.format(str(self.n)))
        trying=False
        if (self.m+self.n < 21) and feedDown==False:
            trying=True
            try:
                with open(__file__[:__file__.find('PySCRDT.py')]+'potentialsPy3','rb') as f:                                    
                    a=dill.load(f)
                a=np.array(a)
                a=a[np.where(a[:,0]==self.m)[0]]
                self.f=a[np.where(a[:,1]==self.n)][0][2]
            except:
                try:
                    with open(__file__[:__file__.find('PySCRDT.py')]+'potentialsPy2','rb') as f:                                    
                        a=dill.load(f)
                    a=np.array(a)
                    a=a[np.where(a[:,0]==self.m)[0]]
                    self.f=a[np.where(a[:,1]==self.n)][0][2]
                except:
                    trying=False     
                    print('# PySCRDT::potential: Calculating potential')  
        if trying==False:
            V = (-1+sy.exp(-self.x**2/(self.t+2*self.a**2)-self.y**2/(self.t+2*self.b**2)))/sy.sqrt((self.t+2*self.a**2)*(self.t+2*self.b**2))
            if self.m>self.n:
                if feedDown:
                    p1 = sy.series(V, self.x, 0, abs(self.m)+2).removeO()
                else:    
                    p1 = sy.series(V, self.x, 0, abs(self.m)+1).removeO()
                p2 = sy.series(p1, self.y, 0, abs(self.n)+1).removeO()
                termy = sy.collect(p2, self.y, evaluate=False)
                termpowy=termy[self.y**abs(self.n)]
                if feedDown:
                    termpowy=sy.expand(termpowy.subs(self.x,self.x+self.D))
                termx = sy.collect(termpowy, self.x, evaluate=False)
                termpowx=termx[self.x**abs(self.m)]
                sterm=sy.simplify(termpowx)
            else:
                p1 = sy.series(V, self.y, 0, abs(self.n)+1).removeO()
                if feedDown:
                    p2 = sy.series(p1, self.x, 0, abs(self.m)+2).removeO()
                else:    
                    p2 = sy.series(p1, self.x, 0, abs(self.m)+1).removeO()
                termx = sy.collect(p2, self.x, evaluate=False)
                if feedDown:
                    termx=sy.expand(termx.subs(self.x,self.x+self.D))
                termpowx=termx[self.x**abs(self.m)]
                termy = sy.collect(termpowx, self.y, evaluate=False)
                termpowy=termy[self.y**abs(self.n)]
                sterm=sy.simplify(termpowy)
            res = sy.integrate(sterm, (self.t, 0, sy.oo)).doit()
            result=res.doit()
            self.V = sy.simplify(result)
            self.f=sy.lambdify((self.a,self.b,self.D),self.V)
        
    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *

    def ksc(self):
        """
        Calculates the space charge perveance Ksc from the parameters dictionary
        Returns: Void
        """
        if self.parameters is None:
            raise IOError('# PySCRDT::ksc: You need to define parameters in [setParameters]')
        if self.data is None:
            raise IOError('# PySCRDT::ksc: You need to define Madx twiss file in [prepareData]')
        if self.parameters['bF']:
            self.K= 2*self.parameters['intensity']*self.parameters['ro']*(self.parameters['harmonic']/self.parameters['bF'])/(self.parameters['C']*self.parameters['b']**2*self.parameters['g']**3)
        else:
            self.K= 2*self.parameters['intensity']*self.parameters['ro']/(np.sqrt(2*np.pi)*self.parameters['bunchLength']*self.parameters['b']**2*self.parameters['g']**3)
    
    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *

    def beamSize(self):
        """
        Calculates the transverse beam sizes from the parameters dictionary and the twiss file
        Returns: Void
        """
        if self.parameters is None:
            raise IOError('# PySCRDT::beamSize: You need to define parameters in [setParameters]')
        if self.data is None:
            raise IOError('# PySCRDT::ksc: You need to define Madx twiss file in [prepareData]')
        self.sx=np.sqrt(self.parameters['emittance_x']*self.data[:,1]/(self.parameters['b']*self.parameters['g'])+(self.parameters['dpp_rms']*self.data[:,3])**2)
        self.sy=np.sqrt(self.parameters['emittance_y']*self.data[:,2]/(self.parameters['b']*self.parameters['g'])+(self.parameters['dpp_rms']*self.data[:,4])**2)
        
    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *

    def prepareData(self,twissFile): 
        """
        Prepares the data from a MADX Twiss file including at least {s, betx, bety, dx, dy, mux, muy, l}
        Inputs : twissFile : [str] twiss file (default=None)
        Returns: Void
        """
        if twissFile is None:
            raise IOError('# PySCRDT::prepareData: You need to define Madx twiss file in [prepareData]')
        if self.parameters is None:
            raise IOError('# PySCRDT::prepareData: You need to define parameters in [setParameters]')
        params=np.genfromtxt(twissFile,max_rows=40,dtype=str)
        for i in enumerate(params):
            if params[i[0]][1]=='GAMMA':
                self.parameters['g']=float(params[i[0]][3])
                self.parameters['b']=np.sqrt(1-1/self.parameters['g']**2)
            elif params[i[0]][1]=='LENGTH':
                self.parameters['C']=float(params[i[0]][3])
            elif params[i[0]][1]=='Q1':
                self.actualQx=float(params[i[0]][3])
            elif params[i[0]][1]=='Q2':
                self.actualQy=float(params[i[0]][3])
        header=np.genfromtxt(twissFile,skip_header=45,max_rows=1,dtype=str)
        s = np.linspace(0,self.parameters['C'],100000)
        try:
            data=np.loadtxt(twissFile,skiprows=47,usecols=(np.where(header=='S')[0][0]-1,np.where(header=='BETX')[0][0]-1,np.where(header=='BETY')[0][0]-1,np.where(header=='DX')[0][0]-1,np.where(header=='DY')[0][0]-1,np.where(header=='MUX')[0][0]-1,np.where(header=='MUY')[0][0]-1,np.where(header=='L')[0][0]-1,np.where(header=='ALFX')[0][0]-1,np.where(header=='ALFY')[0][0]-1))
            data2=np.zeros((100000,10))
            data2[:,8] = np.interp(s,data[:,0],data[:,8])
            data2[:,9] = np.interp(s,data[:,0],data[:,9])
        except:
            data=np.loadtxt(twissFile,skiprows=47,usecols=(np.where(header=='S')[0][0]-1,np.where(header=='BETX')[0][0]-1,np.where(header=='BETY')[0][0]-1,np.where(header=='DX')[0][0]-1,np.where(header=='DY')[0][0]-1,np.where(header=='MUX')[0][0]-1,np.where(header=='MUY')[0][0]-1,np.where(header=='L')[0][0]-1))
            data2=np.zeros((100000,8))
        data2[:,1] = np.square(np.interp(s,data[:,0],np.sqrt(data[:,1])))
        data2[:,2] = np.square(np.interp(s,data[:,0],np.sqrt(data[:,2])))
        data2[:,3] = np.interp(s,data[:,0],self.parameters['b']*data[:,3])
        data2[:,4] = np.interp(s,data[:,0],self.parameters['b']*data[:,4])
        data2[:,5] = np.interp(s,data[:,0],data[:,5])
        data2[:,6] = np.interp(s,data[:,0],data[:,6])
        data2[:,7] += self.parameters['C']/len(s)
        data2[:,0] = s
        self.data=data2
        self.beamSize()
        self.ksc()
        
    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *

    def reIndexing(self,factor):
        """
        Auxiliary Function
        Returns: [dict]
        """
        self.dictionary={}
        for i in factor.keys():
            if len(i.args)==0:
                self.dictionary[i]=factor[i]
            else:
                self.dictionary[sy.exp(i.args[0]/1.)]=factor[i]
        return self.dictionary
    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *

    def calculateFactor(self,Detuning=False):
        """
        Auxiliary Function
        Returns: [float]
        """
        if Detuning==True:
            if self.m==0:
                det1=sy.cos(self.fy)**abs(self.n)
                det3=det1.rewrite(sy.exp)
                det2=sy.expand(det3)
                self.factor_d=float(sy.collect(det2,sy.exp(1j*self.fy),evaluate=False)[1])
            elif self.n==0:
                det1=sy.cos(self.fx)**abs(self.m)
                det3=det1.rewrite(sy.exp)
                det2=sy.expand(det3)
                self.factor_d=float(sy.collect(det2,sy.exp(1j*self.fx),evaluate=False)[1])
            else:
                det1=(sy.cos(self.fx)**abs(self.m)*sy.cos(self.fy)**abs(self.n))
                det3=det1.rewrite(sy.exp)
                det2=sy.expand(det3)
                factor1=sy.collect(det2,sy.exp(1j*self.fx),evaluate=False)[1]
                self.factor_d=float(sy.collect(factor1,sy.exp(1j*self.fy),evaluate=False)[1])
            return self.factor_d
        else:
            if self.mode==5:
                self.m=self.h-self.i
                self.n=self.j-self.k
            if self.m==0:
                det1=sy.cos(self.fy)**abs(self.n)
                det3=det1.rewrite(sy.exp)
                det2=sy.expand(det3)
                factor=sy.collect(det2,sy.exp(1j*self.fy),evaluate=False)
                dictionary=self.reIndexing(factor)
                self.factor=float(2.*dictionary[sy.exp(abs(self.n)*1j*self.fy)])
            elif self.n==0:
                det1=sy.cos(self.fx)**abs(self.m)
                det3=det1.rewrite(sy.exp)
                det2=sy.expand(det3)
                factor=sy.collect(det2,sy.exp(1j*self.fx),evaluate=False)
                dictionary=self.reIndexing(factor)
                self.factor=float(2.*dictionary[sy.exp(abs(self.m)*1j*self.fx)])
            else:
                det1=(sy.cos(self.fx)**abs(self.m)*sy.cos(self.fy)**abs(self.n))
                det3=det1.rewrite(sy.exp)
                det2=sy.expand(det3)
                factor1=sy.collect(det2,sy.exp(1j*self.fx),evaluate=False)
                dictionary=self.reIndexing(factor1)
                factor1=dictionary[sy.exp(abs(self.m)*1j*self.fx)]
                factor2=sy.collect(factor1,sy.exp(1j*self.fy),evaluate=False)
                dictionary=self.reIndexing(factor2)
                self.factor=float(2.*dictionary[sy.exp(abs(self.n)*1j*self.fy)])
            return self.factor
    
    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *

    def resonanceDrivingTerms(self,feedDown=False):
        """
        Calculates the resonance driving terms for the resonance requested
        Returns: Void
        """
        self.feed=feedDown
        if self.V is None:
            self.potential(feedDown=self.feed)
        if self.data is None:
            raise IOError('# PySCRDT::resonanceDrivingTerms: You need to run [prepareData] first')
        if self.K is None:
            self.ksc()
        if self.factor is None:
            self.calculateFactor()
        if self.mode==3:
            if self.l=='any':
                self.rdt_s=self.factor*self.data[:,7]*self.K/2./(2*np.pi)*(np.sqrt(2*self.data[:,1])**abs(self.m))*(np.sqrt(2*self.data[:,2])**abs(self.n))*self.f(self.sx,self.sy,self.parameters['dpp']*self.data[:,3])*np.exp(1j*2*np.pi*(self.m*self.data[:,5]+self.n*self.data[:,6]))
            else:
                self.rdt_s=self.data[:,7]*self.factor*self.K/2./(2*np.pi)*(np.sqrt(2*self.data[:,1])**abs(self.m))*(np.sqrt(2*self.data[:,2])**abs(self.n))*self.f(self.sx,self.sy,self.parameters['dpp']*self.data[:,3])*np.exp(1j*(self.m*2*np.pi*self.data[:,5]+self.n*2*np.pi*self.data[:,6]+(self.l-self.m*self.actualQx-self.n*self.actualQy)*2*np.pi*self.data[:,0]/self.parameters['C']))
        else:
            if self.l=='any':
                self.rdt_s=self.factor*self.data[:,7]*self.K/2./(2*np.pi)*(np.sqrt(2*self.data[:,1])**abs(self.h+self.i))*(np.sqrt(2*self.data[:,2])**abs(self.j+self.k))*self.f(self.sx,self.sy,self.parameters['dpp']*self.data[:,3])*np.exp(1j*2*np.pi*((self.h-self.i)*self.data[:,5]+(self.j-self.k)*self.data[:,6]))
            else:
                self.rdt_s=self.data[:,7]*self.factor*self.K/2./(2*np.pi)*(np.sqrt(2*self.data[:,1])**abs(self.h+self.i))*(np.sqrt(2*self.data[:,2])**abs(self.j+self.k))*self.f(self.sx,self.sy,self.parameters['dpp']*self.data[:,3])*np.exp(1j*((self.h-self.i)*2*np.pi*self.data[:,5]+(self.j-self.k)*2*np.pi*self.data[:,6]+(self.l-(self.h-self.i)*self.actualQx-(self.j-self.k)*self.actualQy)*2*np.pi*self.data[:,0]/self.parameters['C']))
        self.rdt=sum(self.rdt_s)

    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *
    
    def detuning(self):
        """
        Calculates the non linear detuning terms
        Returns: Void
        """
        if self.V is None:
            self.potential()
        if self.data is None:
            raise IOError('# PySCRDT::resonanceDrivingTerms: You need to run [prepareData] first')
        if self.K is None:
            self.ksc()
        if self.factor_d is None:
            self.calculateFactor(Detuning=True)
        if self.mode==3:
            self.rdt_s_d=self.factor_d*self.data[:,7]*self.K/2./(2*np.pi)*(np.sqrt(2*self.data[:,1])**abs(self.m))*(np.sqrt(2*self.data[:,2])**abs(self.n))*self.f(self.sx,self.sy,self.parameters['dpp']*self.data[:,3])
        else:
            self.rdt_s_d=self.factor_d*self.data[:,7]*self.K/2./(2*np.pi)*(np.sqrt(2*self.data[:,1])**abs(self.h+self.i))*(np.sqrt(2*self.data[:,2])**abs(self.j+self.k))*self.f(self.sx,self.sy,self.parameters['dpp']*self.data[:,3])
        self.rdt_d=sum(self.rdt_s_d)

    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *

    def updateParameters(self,**kwargs):
        """
        Updates the parameter dictionary
        Input :  any of  'intensity' 
                         'bunchLength' 
                         'ro'
                         'emittance_x'
                         'emittance_y'
                         'dpp_rms'
                         'dpp'
                         'b'
                         'g'
                         'bF'
                         'harmonic'
        Returns: void
        """
        if kwargs is not None:
            for key, value in kwargs.iteritems():
                if key not in self.parameters.keys():
                    raise IOError('# PySCRDT::updateParameters: '+key+' not recognized [checkWriting]')
                else:
                    self.parameters[key]=value
            if self.data is not None:
                self.beamSize()
                self.ksc()
                
    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *

    def getParameters(self):
        """
        Returns the parameters dictionary
        Returns: [dict] the parameters dictionary
        """
        if self.parameters is None:
            raise IOError('# PySCRDT::getParameters: You need to define parameters in [setParameters]|[readParameters]')
        else:
            return self.parameters 

    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *
    
    def getWorkingPoint(self):
        """
        Returns the tunes (Qx, Qy)
        Returns: [tuple] 
        """
        if self.data is None:
            raise IOError('# PySCRDT::getWorkingPoint: You need to define Madx twiss file in [prepareData]')
        else:
            return self.actualQx, self.actualQy

    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *

    def getOrder(self):
        """
        Returns the orders of the resonances
        Returns: [tuple] 
        """
        if self.mode is None:
            raise IOError('# PySCRDT::getOrder: You need to define resonance mode description in [setMode]')
        elif self.mode==3:
            if self.m is None:
                raise IOError('# PySCRDT::getOrder: You need to define the order in [setOrder]')
            else:
                return self.m, self.n, self.l
        elif self.mode==5:
            if self.h is None:
                raise IOError('# PySCRDT::getOrder: You need to define the order in [setOrder]')
            else:
                return self.h, self.i, self.j, self.k, self.l

    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *

    def getMode(self):
        """
        Returns the resonance mode description
        Returns: [int] 
        """
        if self.mode is None:
            raise IOError('# PySCRDT::getMode: You need to define resonance mode description in [setMode]')
        else:
            return self.mode
    
    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *
    
    def getKsc(self):
        """
        Returns the space charge perveance Ksc
        Returns: [float] 
        """
        if self.K is None:
            self.ksc()
        return self.K

    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *
   
    def getPotential(self):
        """
        Returns the potential V
        Returns: [sympy expression] 
        """
        if self.V is None:
            self.potential()
        return self.V

    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *

    def getResonanceDrivingTerms(self,feedDown=False):
        """
        Returns the RDTs
        Inputs : feedDown : [bool] needed only if [resonanceDrivingTerms] has not been already used (default=False)
        Returns: [dict] the resonance driving terms 
        """
        self.feed=feedDown
        if self.rdt is None:
            self.resonanceDrivingTerms(feedDown=self.feed)
        self.RDT={'RDT':self.rdt, 'Amplitude': abs(self.rdt), 'Phase': np.angle(self.rdt)}
        return self.RDT

    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *

    def getDetuning(self):
        """
        Returns the RDTs
        Returns: [dict] the resonance driving terms
        """
        if self.rdt_d is None:
            self.detuning()
        return self.rdt_d

    # - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *
    
    def checkWriting(self):
        """
        Returns the correct writing format for setting or updating parameters
        Returns: [dict]
        """
        return {'Set & Update':['intensity', 'bunchLength', 'emittance_x', 'emittance_y', 'dpp_rms', 'dpp', 'ro'], 'Update only': ['b', 'g']}

