# Potencial.py
# generates a potential ramp #
# *** warning supresion
import warnings
warnings.filterwarnings("ignore")
import time

# *** numeric libraries
try:
	import numpy as np #pip3 install numpy
	from functools import reduce
	import scipy.io
	from scipy.interpolate import interp1d
except: 
	print('WARNNING :: Potencial.py :: can NOT correctly import numerical libraries')
	print('Install by: ( pip3 install scipy )')
	
# *** graph libraries
try:
	import matplotlib.pyplot as plt #pip3 install matplotlib
	from matplotlib.animation import FuncAnimation
	import ipywidgets as widgets
except: 
	print('WARNNING :: Potencial.py :: can NOT correctly import graph libraries')
	print('Install by: ( pip3 install matplotlib )')
	
class Capacitance(object):
	def __init__(self, 	Er=None, E0=None, Xh=None, 
						z=None,  e=None,  Kb=None, 
						T=None,  ni=None, 
						C=None,  Xdl=None,
						K1=None, K2=None, K3=None, ):
		# == Physical parameters == # 
		self.Er  	= Er # Solution permittivity	| 	Float 	|
		self.E0  	= E0 # vacuum permittivity		|	Float 	|
		self.Xh 	= Xh #
		self.z 		= z  # 1Q
		self.e 		= e  #
		self.Kb		= Kb #
		self.T 		= T  # 298K | temperatura 	| 	Float 	|
		self.ni		= ni #
		self.initialize(Er, E0, Xh, z, e, Kb, T, ni)
		
		self.K1		= K1 
		self.K2		= K2
		self.K3		= K3
		self.evaluate_K()

		self._C = None
		self.Ch = None
		self.Cgc= None

		self._U = None
		self.Xdl 	= Xdl # 0.3nm  in water + 1M

	@property
	def C(self):
		return self._C
		
	@C.setter
	def C(self, U=None, K1=None, K2=None, K3=None):
		K1 = K1 if type(K1) != type(None) else self.K1
		K2 = K2 if type(K2) != type(None) else self.K2
		K3 = K3 if type(K3) != type(None) else self.K3
		U  = U  if type(U)  != type(None) else self.U

		Ch = np.repeat(K1, U.shape)
		K1, K2, K3 = 1,1,1 
		Cgc = K2*np.cosh(K3*U) 
		self._C = Ch*Cgc/(Ch+Cgc) 
		self.Ch  = Ch
		self.Cgc = Cgc

	@property
	def U(self):
		return self._U
		
	@U.setter
	def U(self, u):
		self.C = u
		self._U = u

	def initialize(self, Er=None, E0=None, Xh=None, 
						 z =None, e =None, Kb=None, 
						 T =None, ni=None, save=True):
		# == Physical parameters == # 
		E0 = E0 if type(E0) != type(None) else 500 # vacuum permittivity		|	Float 	|
		Er = Er if type(Er) != type(None) else E0*80.2 # Solution permittivity	| 	Float 	|
		Xh = Xh if type(Xh) != type(None) else 1
		z  = z  if type(z)  != type(None) else 1 # 1Q
		e  = e  if type(e)  != type(None) else 1
		Kb = Kb if type(Kb) != type(None) else 1
		T  = T  if type(T)  != type(None) else 298 # 298K | temperatura 	| 	Float 	|
		ni = ni if type(ni) != type(None) else 1
		
		if save:
			self.Er, self.E0, self.Xh = Er, E0, Xh
			self.Kb, self.T,  self.ni = Kb, T , ni
			self.z , self.e = z, e

		return Er, E0, Xh, z, e, Kb, T, ni

	def print(self, ):
		if type(self.U) == np.ndarray:
			print(f'')

	def timer(func):
		def wrapper(*args, **kwargs):
			before = time.time()
			func(*args, **kwargs)
			print(f'{func} {time.time()-before}s ')
		
		return wrapper

	def evaluate_K(self, 	Er=None, E0=None, Xh=None, 
							z=None,  e=None,  Kb=None, 
							T=None, ni=None,  save=True):

		Er = Er if type(Er) != type(None) else self.Er
		E0 = E0 if type(E0) != type(None) else self.E0
		Xh = Xh if type(Xh) != type(None) else self.Xh
		z  = z  if type(z)  != type(None) else self.z
		e  = e  if type(e)  != type(None) else self.e
		Kb = Kb if type(Kb) != type(None) else self.Kb
		T  = T  if type(T)  != type(None) else self.T
		ni = ni if type(ni) != type(None) else self.ni

		K1 = Er*E0/Xh 
		K2 = (2*Er*E0*z**2*e**2*ni/(Kb*T))**0.5 
		K3 = z*e/(2*Kb*T)

		if save:
			self.K1, self.K2, self.K3 = K1, K2, K3
			self.Er, self.E0, self.Xh = Er, E0, Xh
			self.Kb, self.T,  self.ni = Kb, T , ni
			self.z ,self.e = z, e

		return K1, K2, K3

	def GC_thickness(self, 	Er=None, E0=None, Xh=None, 
							z=None,  e=None,  Kb=None, 
							T=None, ni=None,  save=True):
		Er = Er if type(Er) != None else self.Er
		E0 = E0 if type(E0) != None else self.E0
		Xh = Xh if type(Xh) != None else self.Xh
		z  = z  if type(z)  != None else self.z
		e  = e  if type(e)  != None else self.e
		Kb = Kb if type(Kb) != None else self.Kb
		T  = T  if type(T)  != None else self.T
		ni = ni if type(ni) != None else self.ni

		Xdl = (E0*Er*Kb*T/ (2*ni*z**2*e**2) )**0.5

		if save:
			self.Xdl = Xdl
			self.Er, self.Er, self.Xh = Er, E0, Xh
			self.Kb, self.T,  self.ni = Kb, T , ni
			self.z ,self.e = z, e
		return self.Xdl

	@timer
	def plot(self, ax=None):
		# === Initialize variables === #
		if type(ax) == type(None): 	fig, ax = plt.subplots(1,1)

		# === PLOT === #
		ax.plot(self.U, self.C, 					lw=2, ls='-', color=(0.7,0.4,0.4), alpha=0.9, label='Gouy–Chapman–Stern (GCS) capacitance model' )
		ax.plot(self.U, self.Ch, 					lw=2, ls='-', color=(0.7,0.7,0.4), alpha=0.9, label='Helmholtz Capacitance' )
		ax.plot(self.U, self.Cgc , 	lw=2, ls='-', color=(0.7,0.4,0.7), alpha=0.9, label='Gouy–Chapman Capacitance' )

		ax.set_ylim(np.min(self.C)*0.9, np.min(self.Cgc)*1.7)
		ax.set_xlabel('Potential ')
		ax.set_ylabel('Captacitance')
		ax.set_title('Double layer - Gouy–Chapman–Stern (GCS) capacitance model')

		# === LABEL hansdler === #
		handles, labels = ax.get_legend_handles_labels()

		# reverse the order
		ax.legend(handles[::-1], labels[::-1])

		# or sort them by labels
		import operator
		hl = sorted(zip(handles, labels),
		            key=operator.itemgetter(1))
		handles2, labels2 = zip(*hl)

		ax.legend(handles2, labels2)

	@timer
	def cookbook(self, ):
		'''
		# @property
		# @prop.setter
		# @deleter
		U = Potential()
		U.generate(duration=5, dt=0.1, ciclos=6, slope=2)
		U.plot()
		plt.show()
		'''
		self.initialize()
		self.U = np.arange(-20,20,0.1)
		self.plot()
		plt.show()

'''
C = Capacitance()
C.cookbook()
'''