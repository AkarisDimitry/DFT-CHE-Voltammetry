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
	
class Reaction(object):
	def __init__(self, 	time=None, dt=None, step=None, steps=None,
						k=None, t=None, C=None, C0=None,
						react=None, react_user=None, reactans=None):
		# == Fisical parameters == # 
		self.time 	= time 				# Time variables  			|  ARRAY  |  shape = duration/dt 
		self.dt  	= dt
		self.step 	= step
		self.steps  = steps

		self.t  	= t
		self.C 		= C  				# Potentail 				|  ARRAY  |  shape = duration/dt
		self.C0 	= C0
		self._dCdt_c  = None
		self._dCdt_t  = None

		self.k 		= k

		self.react   	= react 		# [[reactivos], [productos], k, [stekiometria]]
		self.react_user = react_user

		self.reactans   = reactans

	@property
	def dCdt_t(self):
		if type(self._dCdt_t) == type(None):
			self._dCdt_t = (self.time[1:] + self.time[:-1])/2
		return self._dCdt_t

	@dCdt_t.setter
	def dCdt_t(self, dCdt_t):
		self.dCdt_t = dCdt_t

	@property
	def dCdt(self):
		if type(self._dCdt) == type(None):
			self._dCdt = (self.C[1:] - self.C[:-1]) / self.dt
		return self._dCdt

	@dCdt.setter
	def dCdt(self, dCdt):
		self.dCdt = dCdt
		
	def print(self, ):
		pass

	def timer(func):
		def wrapper(*args, **kwargs):
			before = time.time()
			r = func(*args, **kwargs)
			print(f'{func} {time.time()-before}s ')
			return r

		return wrapper

	def set_Cti(self, step_t=0, C=None, reactans=None):
		if type(reactans) == type(None): 
			# a-sync array as input 
			self.C_x_t[step_t, :] = C
		else: 
			# sincronize by reactans name
			if len(C.shape) == 2: 
				# matrix [time, reactant] as input 
				# CC time dependant 
				for i, r in enumerate(self.reactans):
					if r in reactans:
						self.C[step_t, i] = C[:, reactans.index(r)]
			
			elif len(C.shape) == 1: 
				# matrix [reactant] as input 
				# CC time INdependant 
				for i, r in enumerate(self.reactans):
					if r in reactans:
						self.C[step_t, i] = C[reactans.index(r)]

		return self.C

	def reaction_step(self, react=None, C=None, step=None, dt=None):
		react 	= react if type(react) 	!= type(None) else self.react
		C 		= C 	if type(C) 		!= type(None) else self.C
		step 	= step 	if type(step) 	!= type(None) else self.step
		dt 		= dt 	if type(dt) 	!= type(None) else self.dt
		
		for i, c in enumerate(C[step,:]):
			C[step+1,i] = c

		for (reactans, products, rate, coef_react, coef_prod) in react:
			v 					 = np.prod(C[step,reactans]**coef_react) * rate[step]*dt
			v = 0.1 if v>0.1 else -0.1 if v<-0.1 else v

			C[step+1,reactans] 	-= v*coef_react
			C[step+1,products] 	+= v*coef_prod

		C[step+1, C[step+1,:]<0] = 0

		#print(self.reactans) 
		#C[step+1, self.reactans.index('O2')] = 0.25 # !!!!!
		return C[step+1,:]

	@timer
	def reaction_evaluate(self, react=None, C=None, dt=None, time=None,
								save=True, v=True ):
		react 	= react if type(react) 	!= type(None) else self.react
		C 		= C 	if type(C) 		!= type(None) else self.C
		dt 		= dt 	if type(dt) 	!= type(None) else self.dt
		time	= time 	if type(time) 	!= type(None) else self.time

		for step, t in enumerate(time):
			C[step+1,:] = self.reaction_step(react=react, C=C, step=step, dt=dt)

		if save:
			self.C 		= C
			self.react 	= react
			self.step 	= step
			self.dt 	= dt

		if v: self.summary()

		return C

	def allocate_mem(self, 	react_user=None, dt=None, time=None, IC=None,
							reactans=None, save=True	):
		reactans 	= reactans 	if not type(reactans) == type(None) else self.reactans if not type(self.reactans) == type(None) else []
		time_N = int(time/dt)
		react = []
		time = np.array([ n*dt for n in range(time_N) ])

		for (key,r) in react_user.items():
			for (key,value) in r['reactans'].items():
				if not key in reactans: 	reactans.append(key)

			for (key,value) in r['products'].items():
				if not key in reactans: 	reactans.append(key)

		C = np.zeros( (time_N+1, len(reactans)) )

		for (key,r) in react_user.items():
			reactans_bin= np.array([ True if r_o in r['reactans'] else False for r_o in reactans ] )
			products_bin= np.array([ True if r_o in r['products'] else False for r_o in reactans ] )
			rate   		= np.repeat(r['rate'] ,time_N) if type(r['rate']) == float else r['rate']
			coef_react	= np.array([ value for (key,value) in r['reactans'].items() ] )
			coef_prod 	= np.array([ value for (key,value) in r['products'].items() ] )

			react.append([reactans_bin, products_bin, rate, coef_react, coef_prod])


		for i, r in enumerate(reactans):
			if r in IC:
				C[0,i] = IC[r]

		if save:
			self.time_N 	= time_N
			self.time  		= time
			self.react 		= react
			self.C 	 		= C
			self.dt 		= dt
			self.reactans 	= reactans
			self.react_user = react_user

		return react
		# react_user={'Reaction':{'reactans':['H2O':1],'products':['O2':1],'rate':1 } }, dt=0.1, time=10,
		# [[reactivos], [productos], k, [stekiometria]]

	def summary(self, ):
		print('='*5+' Reaction '+'='*5)
		print(f'\t Steps             : {self.step}')
		print(f'\t dt                : {self.dt}s')
		print(f'\t Reaction duration : {self.time_N*self.dt}s')
		print(f'\t Reactions         : {len(self.react)}')

		for i, (key,r) in enumerate(self.react_user.items()):
			reac = ' + '.join([f'{int(r1)} {key1}' if r1!=1 else f'{key1}' for (key1,r1) in r['reactans'].items() ])
			prod = ' + '.join([f'{int(r1)} {key1}' if r1!=1 else f'{key1}' for (key1,r1) in r['products'].items() ])
			print( f'     {reac} --(k{i+1})--> {prod}' )


	@timer
	def plot(self, ax=None, reactans=None):
		# === Initialize variables === #
		if type(ax) == type(None): 	fig, ax = plt.subplots(1,1)

		# === PLOT === #
		for c in range(self.C.shape[1]):
			if type(reactans) == type(None) or self.reactans[c] in reactans:
				ax.plot( self.time, self.C[1:,c], lw=1, ls='--', alpha=0.8, label=self.reactans[c])

		ax.set_xlabel('Time (s)')
		ax.set_ylabel('CC')
		ax.set_title('Reaction plot')

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
		self.allocate_mem( {'Reaction':{'reactans':{'H2O':1},'products':{'O2':1},'rate':0.1 } }, 
						IC={'H2O':0.1, 'O2':0.1}, dt=0.1, time=100 )
		self.reaction_evaluate()
		self.plot()
		plt.show()


'''
# == Eg. == #
R = Reaction()
R.cookbook()
'''