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
						smooth=None, atenuation=None, Umin=None, Umax=None):
		# == Fisical parameters == # 
		self._U		= U  				# Potentail 				|  ARRAY  |  shape = duration/dt
		self._dU	= None  			# Potentail derivate		|  ARRAY  |  shape = duration/dt - 1
		
		self.time_dU = None 				# dU time
		self.time 	 = time 				# Time variables  			|  ARRAY  |  shape = duration/dt 
		self.dtime   = None 				# Time derivate
		
		self.dt 	 = dt 				# Time partition  			|  FLOAT  |  
		self.duration = duration  		# one cicle duration 		|  FLOAT  |
		self.ciclos 	= ciclos		# Total number of cicles 	|  INT    |   

		# == Curve shape control == #
		self.smooth 	= smooth  		# potential ramp smoothness |  FLOAT  | 
		self.atenuation = atenuation  	# atenuation of the ramp    |  FLOAT  | 
		self.Umin       = Umin   		# Maximun value of the ramp |  FLOAT  | 
		self.Umax       = Umax   		# Maximun value of the ramp |  FLOAT  | 

	@property
	def dU(self):
		return self._dU
	
	@dU.setter
	def dU(self, U):

		if type(U) != type(None):
			self.dtime = self.time[1:]-self.time[:-1]
			self._dU = (U[1:]-U[:-1]) / self.dtime
			self.time_dU = (self.time[1:]+self.time[:-1])/2
		else: self._dU = U

	@property
	def U(self):
		return self._U
	
	@U.setter
	def U(self, U):
		self._U = U
		self.dU = U

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
	def generate(self, 	duration=None 		,dt=None 		,ciclos=None 	  , 
						shape='triangular'	,smooth=None 	,atenuation=None, Umin=None, Umax=None, slope=None,
						save=True, ):
		
		def segment(	duration=None , dt=None    , t0=None, value=None, slope=None,
						time=None 	  , signal=None, ):
			return {'time'	:	np.concatenate((time  , np.linspace(t0 	, t0+duration 	, num=int(duration/dt), endpoint=False) 			)), 
					'signal':	np.concatenate((signal, np.linspace(0 	, duration  	, num=int(duration/dt), endpoint=False)*slope+value )) }

		def atenuation_peaks():
			signal_peak = np.zeros(ciclos-1)
			time_peak 	= np.zeros(ciclos-1)
			for c in range(ciclos-1):
				signal_peak[c] = np.max(signal)*((c-1)%2)*(1-atenuation)
				time_peak[c] = int(duration/dt)*(c+2)*dt

			return time_peak, signal_peak

		def softened():
			cicle_lenght = int(duration/dt)
			cicle_sparce_lenght = int(int(int(duration/dt)*(1-smooth))/2)*2
			cicle_sparce_loss = cicle_lenght - cicle_sparce_lenght

			filter_smooth = np.ones( cicle_lenght*(ciclos+2) ).astype(int)
			for c in range(1, ciclos+2):
				filter_smooth[int(c*cicle_lenght-cicle_sparce_loss/2):int(c*cicle_lenght+cicle_sparce_loss/2)] = 0

			time_sparce_smooth, signal_sparce_smooth = np.extract(filter_smooth, time), np.extract(filter_smooth, signal)
			if type(atenuation) != type(None):	time_sparce_peak, 	signal_sparce_peak   = atenuation_peaks()
			else: time_sparce_peak, signal_sparce_peak = np.array([]), np.array([])

			time_sparce 	= np.concatenate((time_sparce_smooth,   time_sparce_peak  ))
			signal_sparce  	= np.concatenate((signal_sparce_smooth, signal_sparce_peak))

			fsp = interp1d(time_sparce, signal_sparce, kind='cubic')
			signal_sp = fsp(time)

			time_sp = time

			return time_sp, signal_sp

		smooth 		= smooth 	 if type(smooth)   	 != type(None) else self.smooth 	if type(self.smooth) 	 != type(None) else 0
		atenuation 	= atenuation if type(atenuation) != type(None) else self.atenuation 
		Umax  		= Umax 		 if type(Umax) 	 	 != type(None) else self.Umax   	if type(self.Umax)   	 != type(None) else np.inf
		Umin 		= Umin  	 if type(Umin)  	 != type(None) else self.Umin   	if type(self.Umin)   	 != type(None) else -np.inf

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

		time, signal = softened()

		# == Set max signal value
		signal[ signal>Umax ] = Umax
		signal[ signal<Umin ] = Umin

		if save:
			self.time 		= time
			self.U 			= signal 

			self.duration 	= duration
			self.dt 		= dt
			self.ciclos 	= ciclos
			self.slope 		= slope

	@timer
	def plot(self, ax=None, step=-1):
		# === Initialize variables === #
		if type(ax) == type(None): 	fig, ax = plt.subplots(1,1)
		ciclos   = self.ciclos   if type(self.ciclos)   == int 			else 0
		duration = self.duration if type(self.duration) == int 			else 0
		U 		 = self.U 		 if type(self.U) 		== np.ndarray 	else 0

		# === PLOT === #
		for c in range(1, 2+ciclos):
			ax.plot( [c*duration, c*duration], [np.min(U), np.max(U)], lw=1, ls='--', color=(0.4,0.4,0.4), alpha=0.4 )

		ax.plot(self.time_dU[:step],  self.dU[:step], 'o', ms=2, lw=2, ls='-', color=(0.4,0.7,0.4), alpha=0.6, label='Applied Potential change (dU/dt)'  )
		ax.plot(self.time[:step],     self.U[:step],  'o', ms=2, lw=2, ls='-', color=(0.7,0.4,0.4), alpha=0.6, label='Applied Potential (U)' )
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
		self.generate(duration=5, dt=0.1, ciclos=6, slope=2, smooth=0.5, atenuation=0.1, Umin=0)
		self.plot()
		plt.show()

'''
U = Potential()
U.cookbook()
'''