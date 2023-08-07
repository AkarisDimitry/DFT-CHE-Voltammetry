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
	
class Potential(object):
	def __init__(self, 	U=None, time=None, dt=None, duration=None, ciclos=None, 
						smooth=None, atenuation=None, Umax=None):
		# == Fisical parameters == # 
		self.U 		= U  				# Potentail 				|  ARRAY  |  shape = duration/dt
		self.time 	= time 				# Time variables  			|  ARRAY  |  shape = duration/dt 
		self.dt 	= dt 				# Time partition  			|  FLOAT  |  
		self.duration = duration  		# one cicle duration 		|  FLOAT  |
		self.ciclos 	= ciclos		# Total number of cicles 	|  INT    |   

		# == Curve shape control == #
		self.smooth 	= smooth  		# potential ramp smoothness |  FLOAT  | 
		self.atenuation = atenuation  	# atenuation of the ramp    |  FLOAT  | 
		self.Umax       = Umax   		# Maximun value of the ramp |  FLOAT  | 

	def print(self, ):
		if type(self.U) == np.ndarray:
			print(f'')

	def timer(func):
		def wrapper(*args, **kwargs):
			before = time.time()
			func(*args, **kwargs)
			print(f'{func} {time.time()-before}s ')
		
		return wrapper

	@timer
	def generate(self, 	duration=None 		, dt=None 		, ciclos=None 	  , 
						shape='triangular'	, smooth=None 	, atenuation=None , Umax=None, slope=None,
						save=True, ):
		
		def segment(	duration=None , dt=None    , t0=None, value=None, slope=None,
						time=None 	  , signal=None, ):
			return {'time'	:	np.concatenate((time  , np.linspace(t0 	, t0+duration 	, num=int(duration/dt)) 			)), 
					'signal':	np.concatenate((signal, np.linspace(0 	, duration  	, num=int(duration/dt))*slope+value )) }

		def softened():	pass

		time 	= np.array([])
		signal 	= np.array([])
		if shape == 'triangular':
			(_, time), (_, signal) = segment(	duration=duration, dt=dt, t0=0, value=0, slope=0,
												time=time, signal=signal).items()
			for c in range(0, ciclos):
				(_, time), (_, signal) = segment(	duration=duration, dt=dt, t0=duration*(c+1), value=duration*slope*(c%2), slope=-slope*(c%2*2-1),
													time=time, signal=signal).items()

			(_, time), (_, signal) = segment(	duration=duration, dt=dt, t0=duration*(c+2), value=-duration*slope*(c%2-1), slope=0,
												time=time, signal=signal).items()

		if save:
			self.time 		= time
			self.U 			= signal 

			self.duration 	= duration
			self.dt 		= dt
			self.ciclos 	= ciclos
			self.slope 		= slope

	@timer
	def plot(self, ax=None):
		# === Initialize variables === #
		if type(ax) == type(None): 	fig, ax = plt.subplots(1,1)
		ciclos   = self.ciclos   if type(self.ciclos)   == int 			else 0
		duration = self.duration if type(self.duration) == int 			else 0
		U 		 = self.U 		 if type(self.U) 		== np.ndarray 	else 0

		# === PLOT === #
		for c in range(1, 2+ciclos):
			ax.plot( [c*duration, c*duration], [np.min(U), np.max(U)], lw=1, ls='--', color=(0.4,0.4,0.4), alpha=0.4 )

		ax.plot(self.time, self.U, lw=2, ls='-', color=(0.7,0.4,0.4), alpha=0.9, label='Applied Potential' )
		ax.set_xlabel('Time (s)')
		ax.set_ylabel('Potential (V)')
		ax.set_title('Potential Ramp')

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
		self.generate(duration=5, dt=0.1, ciclos=6, slope=2)
		self.plot()
