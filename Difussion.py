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
	
class Difussion(object):
	def __init__(self, 	C=None, D=None, F=None, reactans=None,
						duration=None, dt=None,
						x=None, dx=None, ):
		# == Fisical parameters == # 
		self._C 		= None
		self.D 			= D
		self.F  		= F 
		self.reactans 	= reactans 

		self.time 		= None
		self.duration 	= duration
		self.dt 		= dt

		self.dx 		= dx
		self.x 			= x

		self.steps_t 	= None
		self.steps_x 	= None
		self.C_x_t   	= C

		self._dCdx 		= None
		self._dCdt 		= None
		self._dCdt_t 	= None

		self.CCbulk 	= False
		self.CCsurface 	= False

		self.color 	= [ 	
			'#DC143C', # 	crimson 			#DC143C 	(220,20,60)
			'#ADFF2F', #	green yellow 		#ADFF2F 	(173,255,47)
			'#40E0D0', #	turquoise 			#40E0D0 	(64,224,208)
			'#FF8C00', #  	dark orange 		#FF8C00 	(255,140,0)
			'#BA55D3', #	medium orchid 		#BA55D3 	(186,85,211)
			'#1E90FF', #	dodger blue 		#1E90FF 	(30,144,255)
			'#FF1493', #	deep pink 			#FF1493 	(255,20,147)
			'#8B4513', #	saddle brown 		#8B4513 	(139,69,19)
			'#FFD700', #	gold 				#FFD700 	(255,215,0)
			'#808000', #	Olive 				#808000 	(128,128,0)
			'#808080', #	Gray 				#808080 	(128,128,128)
			'#FF00FF', #	Magenta / Fuchsia 	#FF00FF 	(255,0,255)
			'#00FFFF', #	Cyan / Aqua 		#00FFFF 	(0,255,255)
			'#000000', #	Black 				#000000 	(0,0,0)
							# ------- REPEAT -------- #
			'#DC143C', # 	crimson 			#DC143C 	(220,20,60)
			'#ADFF2F', #	green yellow 		#ADFF2F 	(173,255,47)
			'#40E0D0', #	turquoise 			#40E0D0 	(64,224,208)
			'#FF8C00', #  	dark orange 		#FF8C00 	(255,140,0)
			'#BA55D3', #	medium orchid 		#BA55D3 	(186,85,211)
			'#1E90FF', #	dodger blue 		#1E90FF 	(30,144,255)
			'#FF1493', #	deep pink 			#FF1493 	(255,20,147)
			'#8B4513', #	saddle brown 		#8B4513 	(139,69,19)
			'#FFD700', #	gold 				#FFD700 	(255,215,0)
			'#808000', #	Olive 				#808000 	(128,128,0)
			'#808080', #	Gray 				#808080 	(128,128,128)
			'#FF00FF', #	Magenta / Fuchsia 	#FF00FF 	(255,0,255)
			'#00FFFF', #	Cyan / Aqua 		#00FFFF 	(0,255,255)
			'#000000', #	Black 				#000000 	(0,0,0)
			]

	@property
	def C(self):
		return self.C_x_t

	@C.setter
	def C(self, C):
		self.C_x_t = C


	@property
	def dCdx(self):
		if type(self._dCdx) == type(None):
			self._dCdx = (self.C[:,1:,:] - self.C[:,:-1,:]) / self.dx 
		return self._dCdx

	@dCdx.setter
	def dCdx(self, dCdx):
		self._dCdx = dCdx

	@property
	def dCdt(self):
		if type(self._dCdt) == type(None):
			self._dCdt = (self.C[1:,:,:] - self.C[:-1,:,:]) / self.dt 
		return self._dCdt

	@dCdt.setter
	def dCdt(self, dCdx):
		self._dCdt = dCdt

	def print(self, ):
		if type(self.U) == np.ndarray:
			print(f'')

	def timer(func):
		def wrapper(*args, **kwargs):
			before = time.time()
			r = func(*args, **kwargs)
			print(f'{func} {time.time()-before}s ')
			return r
		return wrapper

	def difussion_step(self, 	C_x_t=None,   D=None, 
								step=None,    dt=None, 
								CCsurface=None, CCbulk=None, 
								save=None):
		'''
		References
			1 Quarteroni A., Sacco R., Saleri F. (2007) Numerical Mathematics (Texts in Applied Mathematics). New York: Springer.
			2 Durran D. R. (1999) Numerical Methods for Wave Equations in Geophysical Fluid Dynamics. New York: Springer.
    		3 Fornberg B. (1988) Generation of Finite Difference Formulas on Arbitrarily Spaced Grids, Mathematics of Computation 51, no. 184 : 699-706. PDF.
    	'''
		C_x_t  		= C_x_t 	if type(C_x_t) 		!= type(None) else self.C_x_t
		D      		= D     	if type(D)     		!= type(None) else self.D
		dt 	   		= dt    	if type(dt)    		!= type(None) else self.dt
		CCsurface 	= CCsurface if type(CCsurface)  != type(None) else self.CCsurface
		CCbulk	  	= CCbulk 	if type(CCbulk)  	!= type(None) else self.CCbulk

		#C_x_t1 = C_x_t[step, 1:, :] -  (C_x_t[step, 1:, :]-C_x_t[step, :-1, :]) * D * dt
		#C_x_t2 = C_x_t[step, :-1, :] - (C_x_t[step, 1:, :]-C_x_t[step, :-1, :]) * D * dt
		#C_x_t1 = C_x_t1[:-1]
		gradient_rigth 			=  C_x_t[step, 1:, :] - C_x_t[step, :-1, :]
		gradient_left  			=  C_x_t[step, :-1, :]- C_x_t[step, 1:, :]
		C_x_t[step+1, 1:-1, :]  = (C_x_t[step, :-2, :]+ C_x_t[step, 2:, :])/2 + (gradient_rigth[:-1]+gradient_left[1:])*D*dt

		#grad = np.gradient(C_x_t[step, :, :], axis=0, edge_order=2)
		#C_x_t3 = C_x_t[step, :, :] - (grad) * D * dt 
		#C_x_t[step+1, 1:-1, :] = C_x_t3[1:-1, :]

		if not CCsurface:		C_x_t[step+1, 0, :]  = C_x_t1[0,:]
		if not CCbulk:			C_x_t[step+1, -1, :] = C_x_t1[-1,:]
		else:					
			C_x_t[step+1, -1, :] = C_x_t[step, -1, :]
			C_x_t[step+1,  0, :] = C_x_t[step,  0, :] + gradient_rigth[0]*D*dt
			
		C_x_t[C_x_t<0] = 0

		return	C_x_t[step+1, :, :]

	@timer
	def difussion_system(self, steps_t=None, save=None):
		steps_t = steps_t if type(steps_t) != type(None) else self.steps_t
		
		for step in range(steps_t):
			self.C_x_t[step+1, :, :] = self.difussion_step(step=step)


	def set_Cxi(self, step_x, C=None, reactans=None):
		'''
		=== Edit Difusssion conditions ===
		step_x 		| 	INT 	|   edit a particular distance from electrode
		C 	 		| 	ARRAY 	|   concentration array [time, reactant] or [reactant] 
		reactans	| 	list 	|   list of reactants used to sync data 

		'''
		if type(reactans) == type(None): 
			# a-sync array as input 
			self.C_x_t[:, step_x, :] = C
		else: 
			# sincronize by reactans name
			if len(C.shape) == 2: 
				# matrix [time, reactant] as input 
				# CC time dependant 
				for i, r in enumerate(self.reactans):
					if r in reactans:
						self.C_x_t[:, step_x, i] = C[:, reactans.index(r)]
			
			elif len(C.shape) == 1: 
				# matrix [reactant] as input 
				# CC time INdependant 
				for i, r in enumerate(self.reactans):
					if r in reactans:
						self.C_x_t[:, step_x, i] = C[reactans.index(r)]

		return self.C_x_t

	def set_Cx0(self, C=None, reactans=None):
		self.CCsurface = True
		C_x_t = self.set_Cxi(step_x=0, C=C, reactans=reactans)
		self.C_x_t = C_x_t
		return self.C_x_t

	def set_Cxf(self, C=None, reactans=None):
		self.CCbulk = True
		C_x_t = self.set_Cxi(step_x=-1, C=C, reactans=reactans)
		self.C_x_t = C_x_t
		return self.C_x_t

	def set_Cx0t0(self, step_x=0, step_t=0, C=None):
		self.C_x_t[step_t, step_x, :] = C0
		return self.C_x_t

	def set_Ct0(self, step_x=0, C=None, reactans=None, CCbulk=None, CCsurface=None):
		if type(CCsurface)  != type(None): self.CCsurface 	= CCsurface
		if type(CCbulk) 	!= type(None): self.CCbulk 		= CCbulk 

		if type(reactans) == type(None): 
			# a-sync array as input 
			self.C_x_t[:, 0, :] = C
		else: 
			# sincronize by reactans name
			if len(C.shape) == 2: 
				# matrix [time, reactant] as input 
				# CC time dependant 
				for i, r in enumerate(self.reactans):
					if r in reactans:
						self.C_x_t[0, :, i] = C[:, reactans.index(r)]
			
			elif len(C.shape) == 1: 
				# matrix [reactant] as input 
				# CC time INdependant 
				for i, r in enumerate(self.reactans):
					if r in reactans:
						self.C_x_t[0, :, i] = C[reactans.index(r)]

		return self.C_x_t

	def set_Cx0ti(self, step_t=0, C=None, reactans=None, CCbulk=None, CCsurface=None):
		if type(CCsurface)  != type(None): self.CCsurface 	= CCsurface

		if type(reactans) == type(None): 
			# a-sync array as input 
			self.C_x_t[step_t, 0, :] = C
		else: 
			# sincronize by reactans name
			if len(C.shape) == 2: 
				# matrix [time, reactant] as input 
				# CC time dependant 
				for i, r in enumerate(self.reactans):
					if r in reactans:
						self.C_x_t[step_t, 0, i] = C[:, reactans.index(r)]
			
			elif len(C.shape) == 1: 
				# matrix [reactant] as input 
				# CC time INdependant 
				for i, r in enumerate(self.reactans):
					if r in reactans:
						self.C_x_t[step_t, 0, i] = C[reactans.index(r)]

		return self.C_x_t

	def allocate_mem(self, 	duration=None, dt=None, 
							x=None, dx=None,
							reactans=None, D=None,
							save=True):
		duration= duration  if type(duration)  	!= type(None) else self.duration
		dt 		= dt 	 	if type(dt) 		!= type(None) else self.dt
		x  		= x  	 	if type(x)  		!= type(None) else self.x
		dx 		= dx 	 	if type(dx) 		!= type(None) else self.dx
		reactans= reactans 	if type(reactans)	!= type(None) else self.reactans

		F = len(reactans)
		steps_t = int(duration/dt)
		steps_x = int(x/dx)
		time = np.arange(0, duration, dt)

		C_x_t = np.zeros((steps_t+1, steps_x, F))

		if save:
			self.D 			= np.array(D)
			self.duration 	= duration
			self.dt 		= dt
			self.time 		= time
			
			self.x  		= x
			self.dx 		= dx
			self.steps_t 	= steps_t
			self.steps_x 	= steps_x
			self.C_x_t   	= C_x_t

			self.F    	 	= len(reactans)
			self.reactans 	= reactans

		return True

	@timer
	def plot(self,  ax=None, steps_t=None, steps_x=None,
					x=None, C_x_t=None, reactans=None):
		# === Initialize variables === #
		if type(ax) == type(None): 	fig, ax = plt.subplots(1,1)
		x  		= x  	 	if type(x)  	 != type(None) else self.x
		C_x_t 	= C_x_t 	if type(C_x_t) 	 != type(None) else self.C_x_t
		steps_t = steps_t 	if type(steps_t) != type(None) else self.steps_t
		steps_t = steps_t   if not type(steps_t) == None   else 0
		steps_x = steps_x 	if type(steps_x) != type(None) else self.steps_x


		# === PLOT === #
		for f, r in enumerate(self.reactans):
			if type(reactans) == type(None) or r in reactans:
					ax.plot(	np.linspace(0, x, steps_x), C_x_t[steps_t, :, f], lw=2, ls='-', 
								color=self.color[f], alpha=0.7, label=f'{r}' )


		ax.set_xlabel('X (mm)')
		ax.set_ylabel('CC (mM)')
		ax.set_title('Difussion (fick law)')
		ax.set_ylim(np.min(self.C_x_t), 	np.max(self.C_x_t))

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
		self.allocate_mem(t=10, dt=0.01,x=5, dx=0.1, names=['H2O', 'O2'], D=[1,1] )
		self.set_Ct0(C0=[0,0.5])
		self.set_Cx0t0(step_x=0, step_t=0, C0=[0,1])
		self.difussion_system()
		self.plot(steps_t=-1)
		plt.show()


'''
difussion = Difussion()
difussion.cookbook()
'''
