# Potencial.py
# generates a potential ramp #

# *** warning supresion
import warnings
warnings.filterwarnings("ignore")
import time, operator

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

# *** Voltammetry libraries *** #
try: 	import Difussion
except: print('WARNNING :: Voltagram.py :: can NOT correctly import voltametry libraries :: Difussion.py pack')

try: 	import Reaction
except: print('WARNNING :: Voltagram.py :: can NOT correctly import voltametry libraries :: Reaction.py pack')

try: 	import Potential
except: print('WARNNING :: Voltagram.py :: can NOT correctly import voltametry libraries :: Potential.py pack')

try: 	import Capacitance
except: print('WARNNING :: Voltagram.py :: can NOT correctly import voltametry libraries :: Capacitance.py pack')

try: 	import ORR
except: print('WARNNING :: Voltagram.py :: can NOT correctly import voltametry libraries :: ORR.py pack')

class Voltagram(object):
	def __init__(self, 	difussion=None, reaction=None, potential=None, capacitance=None,
						time=None):
		# == Fisical parameters == # 
		self.difussion	= difussion
		self.D = None

		self.reaction	= reaction
		self.R = None
		
		self.potential	= potential
		self.P = None
		
		self.capacitance= capacitance
		self.C = None
		
		self.ORR = None
		
		self.time 	= time 

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

	@timer
	def initialize(self,):

 		# multiple peaks : 
 		# DOI: 10.1039/C9SC01545K
 		# DOI: 10.1021/acs.jchemed.7b00361.
		# === Make/Setting potential model === # 
		duration = 10 	# seg
		dt = 0.1		# seg
		ciclos =  3
		slope = 0.05 	# V/seg

		self.potential   = Potential.Potential()
		self.potential.generate(duration=duration, dt=dt, ciclos=ciclos, slope=slope, smooth=0.5, atenuation=0.1, Umin=0)
		self.potential.U *= -1 
		self.potential.U += 0.5

		# === Make/Setting ORR model === # 
		self.ORR = ORR.OxigenReaction() 
		self.ORR.calculate(sys={'E':-805.522,'ZPE':0.0,'S':0.0}, sys_O={'E':-811.362,'ZPE':0.07,'S':0.0}, sys_OH={'E':-815.785,'ZPE':0.35,'S':0.0}, sys_OOH={'E':-820.329,'ZPE':0.43,'S':0.0}, sys_O2=None, 
								H2O={'E':-14.213,'ZPE':0.56,'S':0.67}, H2={'E':-6.76,'ZPE':0.27,'S':0.41}, Gb=-0.0*4)
		self.ORR.summary()
		# estimate k #

		self.ORR.G_U(U = self.potential.U )
		k1a, k2a, k3a, k4a = self.ORR.G2K(kmax=1)
		k1b, k2b, k3b, k4b = 1-self.ORR.ORR['ki_U_ORR']
		# *** SIMPLIFICATION one limiting step #
		print( np.max((k1a, k2a, k3a, k4a), axis=0) )
		k1a, k2a, k3a, k4a = [np.min((k1a, k2a, k3a, k4a), axis=0)]*4
		k1b, k2b, k3b, k4b = [np.max((k1b, k2b, k3b, k4b), axis=0)]*4
		print(k1a, k2a, k3a, k4a)
	
		# === === ===Make/Setting reaction model === === === # 
		self.reation     = Reaction.Reaction()
		self.reation.allocate_mem( 
						{
						 'R1a':{'reactans':{'A':1, 'O2':1},'products':{'OH-':1, 'B':1},'rate':k1a }, 
						 'R1b':{'reactans':{'OH-':1, 'B':1},'products':{'A':1, 'O2':1},'rate':k1a }, 

						 'R2a':{'reactans':{'B':1},'products':{'OH-':1, 'C':1},'rate':k2a }, 
						 'R2b':{'reactans':{'OH-':1, 'C':1},'products':{'B':1},'rate':k2a }, 

						 'R3a':{'reactans':{'C':1},'products':{'OH-':1, 'D':1},'rate':k3a }, 
						 'R3b':{'reactans':{'OH-':1, 'D':1},'products':{'C':1},'rate':k3b }, 

						 'R4a':{'reactans':{'D':1},'products':{'OH-':1, 'A':1},'rate':k4a }, 
						 'R4b':{'reactans':{'OH-':1, 'A':1},'products':{'D':1},'rate':k4a }, 
						 }, 
						IC={'A':0.1, 'B':0.0, 'C':0.0, 'D':0.0, 'OH-':0.00005, 'O2':0.25}, dt=dt, time=(ciclos+2)*duration )
		self.reation.summary()		
		#self.reation.reaction_evaluate()

		# === === === Make/Setting difussion model === === === # 
		self.difussion   = Difussion.Difussion()
		self.difussion.allocate_mem(duration=(ciclos+2)*duration, dt=dt,x=3, dx=0.01, 
									reactans=['O2', 'OH-'], D=[0.1,0.3] )
		# SET electrode conditions #
		self.difussion.set_Cx0(C=self.reation.C, reactans=self.reation.reactans)
		# SET bulk conditions #
		self.difussion.set_Cxf(C=self.reation.C[0,:], reactans=self.reation.reactans)
		self.difussion.set_Ct0(C=self.reation.C[0,:], reactans=self.reation.reactans)

		#self.difussion.difussion_system()
		self.coupled_system(time=self.potential.time)

		# === === === Make/Setting capacitance model === === === # 
		self.capacitance = Capacitance.Capacitance()
		self.capacitance.U = self.potential.U

		self.current_evaluate()

		#self.ORR.plot_k()
		#self.potential.plot()
		#self.capacitance.plot()
		#self.difussion.plot(steps_t=1200, reactans=['O2', 'OH-'])
		#self.reation.plot( )#reactans=['O2', 'OH-'] )

		self.plot()
		plt.show()

		return True

	@timer
	def coupled_system(self, dt=None, time=None):

		def coupled_step(step):
			# === reacction step === #
			self.reation.C[step+1,:] = self.reation.reaction_step(step=step)

			# === Difussion step === #
			self.difussion.set_Cx0ti(	step_t=step, 
										C=self.reation.C[step+1,:], 
										reactans=self.reation.reactans)

			C = self.difussion.difussion_step(step=step)
			self.reation.set_Cti(step_t=step+1, C=C[0,:], reactans=self.difussion.reactans)

			return C

		for step, t in enumerate(time):
			coupled_step(step)


	@timer
	def current_faradic_evaluate(self, ):
		# -z*F*D * dCC[0,0] * dt *0.1)
		# Effects of Voltage-Dependence of the Constant Phase Element and Ohmic Parameters in the Modeling and Simulation of Cyclic Voltammograms
		# DOI:10.1149/1945-7111/abbe5d
		
		# Determination of Constant Phase Element Parameters under Cyclic Voltammetry Conditions Using a Semi-theoretical Equation
		# https://doi.org/10.5796/electrochemistry.18-00082
		f_dcdx = interp1d(self.difussion.time , self.difussion.dCdx[1:,0, 1], kind='cubic')
		dCdx = f_dcdx( self.potential.time_dU )
		
		#Ic = 1 * np.e**( (-1*self.potential.U)/(1) )
		#plt.figure(10), plt.plot(self.potential.U, Ic)
		#plt.show()
		print(self.reation.reactans)
		#self.If = -dCdx + (self.reation.C[1:-1,2]-self.reation.C[:-2,2])/self.potential.dt
		self.If =  (self.reation.C[1:-1,2]-self.reation.C[:-2,2])/self.potential.dt


	@timer
	def current_nonfaradic_evaluate(self, ):
		#Inf.append(-z*F*D * dCC[0,0] * dt *0.1)
		f_ct = interp1d(self.potential.time , self.capacitance.C, kind='cubic')
		Ct = f_ct( self.potential.time_dU )

		#f_difdt = interp1d( If_t , If, kind='cubic') 
		#difdt = f_difdt( self.potential.time_dU )
		
		self.Inf = Ct * ( (self.potential.dU) + 0.01) 

		return self.Inf

	@timer
	def current_evaluate(self, ):
		self.current_faradic_evaluate()
		self.current_nonfaradic_evaluate()

		self.I = -self.If + self.Inf

		return self.I

	@timer
	def time_step(self,):
		pass
			
	def allocate_mem(self, 	t=None, dt=None):
		self.If  = np.zeros( int(self.t/self.dt) )
		self.Inf = np.zeros( int(self.t/self.dt) )
		self.I   = np.zeros( int(self.t/self.dt) )

	@timer
	def plot(self, ax=None, time_step=10):
		def animation_frame(i):
			step = i*time_step
			if i%10==0: print(f'step {i}')
			# === PLOT potential curve === # 
			#ax[0,0].clear()
			#self.potential.plot(ax[0,0], step=step)
			plot_POTENTIAl0.set_xdata( self.potential.time[:step] )
			plot_POTENTIAl0.set_ydata( self.potential.U[:step] )
			plot_POTENTIAl1.set_xdata( self.potential.time_dU[:step],  )
			plot_POTENTIAl1.set_ydata( self.potential.dU[:step] )
			
			# === PLOT reaction rate === # 
			#ax[0,1].clear()
			#self.ORR.plot_k(ax=ax[0,1], step=step)
			plot_RATE0.set_xdata( self.potential.time_dU[:step],  )
			plot_RATE0.set_ydata( self.ORR.ORR['ki_U_ORR'][0][:step] )
			plot_RATE1.set_xdata( self.potential.time_dU[:step],  )
			plot_RATE1.set_ydata( self.ORR.ORR['ki_U_ORR'][1][:step] )
			plot_RATE2.set_xdata( self.potential.time_dU[:step],  )
			plot_RATE2.set_ydata( self.ORR.ORR['ki_U_ORR'][2][:step] )
			plot_RATE3.set_xdata( self.potential.time_dU[:step],  )
			plot_RATE3.set_ydata( self.ORR.ORR['ki_U_ORR'][3][:step] )
			
			# === PLOT reaction profile === # 
			#ax[0,2].clear()						
			#self.ORR.plot2(ax=ax[0,2], 
			#	G1=self.ORR.ORR['G1_U_ORR'][step],
			#	G2=self.ORR.ORR['G2_U_ORR'][step],
			#	G3=self.ORR.ORR['G3_U_ORR'][step],
			#	G4=self.ORR.ORR['G4_U_ORR'][step])
			G1 = self.ORR.ORR['Gi_U_ORR'][0][step]
			G2 = self.ORR.ORR['Gi_U_ORR'][1][step]+G1
			G3 = self.ORR.ORR['Gi_U_ORR'][2][step]+G1+G2
			G4 = self.ORR.ORR['Gi_U_ORR'][3][step]+G1+G2+G3
			G1f = self.ORR.ORR['Gi_U_ORR'][3][step]+G1+G2+G3+G4
			plot_G1.set_xdata( [0.0, 1.0],  )
			plot_G1.set_ydata( [G1,G1]  )
			plot_G1f.set_xdata( [6.0, 7.0],  )
			plot_G1f.set_ydata( [G1f,G1f]  )
			plot_G2.set_xdata( [1.5, 2.5],  )
			plot_G2.set_ydata( [G2,G2]  )
			plot_G3.set_xdata( [3.0, 4.0],  )
			plot_G3.set_ydata( [G3,G3]  )
			plot_G4.set_xdata( [4.5, 5.5],  )
			plot_G4.set_ydata( [G4,G4]  )
			
			# === PLOT Difussion profile === # 
			#ax[1,0].clear()						
			#self.difussion.plot(ax=ax[1,0], steps_t=step)
			plot_D1.set_xdata( np.linspace( 0, self.difussion.x, int(self.difussion.steps_x) )  )
			plot_D2.set_xdata( np.linspace( 0, self.difussion.x, int(self.difussion.steps_x) ) )
			for i, r in enumerate(self.difussion.reactans):
				if r == 'O2': 			plot_D1.set_ydata( self.difussion.C_x_t[step,:,i] )
				if r == 'OH-': 			plot_D2.set_ydata( self.difussion.C_x_t[step,:,i] )
			
			# === PLOT CC evolution === # 
			#ax[1,1].clear()						
			#self.reation.plot(ax=ax[1,1])
			plot_R1.set_xdata( self.potential.time[:step] )
			plot_R2.set_xdata( self.potential.time[:step] )

			for i, r in enumerate(self.reation.reactans):
				if r == 'OH-': 			plot_R2.set_ydata( self.reation.C[:step,i] )
				if r == 'O2': 			plot_R1.set_ydata( self.reation.C[:step,i] )

			# === PLOT CC evolution === # 
			#ax[1,1].clear()						
			#self.reation.plot(ax=ax[1,1])
			plot_V.set_xdata( self.potential.U[:step] )
			plot_V.set_ydata( self.I[:step] )


		fig, ax = plt.subplots(2,3)

		# === parameters === #
		Umin, Umax = np.min(self.potential.U), np.max(self.potential.U)
		Usig = (Umax-Umin)
		Tmin, Tmax = np.min(self.potential.time), np.max(self.potential.time)
		Tsig = (Tmax-Tmin)

		Kmin, Kmax = np.min(self.ORR.ORR['ki_U_ORR']), np.max(self.ORR.ORR['ki_U_ORR'])
		Ksig = (Kmax-Kmin)
		Gmin, Gmax = np.min(self.ORR.ORR['Gi_U_ORR']), np.max(self.ORR.ORR['Gi_U_ORR'])
		Gsig = (Gmax-Gmin)

		Cmin, Cmax = np.min(self.difussion.C_x_t), np.max(self.difussion.C_x_t)
		Csig = (Cmax-Cmin)
		Imin, Imax = np.min(self.I), np.max(self.I)
		Isig = (Imax-Imin)

		# === PLOT potential curve === #  # === PLOT potential curve === # 		
		plot_POTENTIAl0, = ax[0,0].plot(	0,0, '-o', ms=2, lw=2, ls='-', 
										color=(0.7,0.4,0.4), alpha=0.6, 
										label='Applied Potential (U)')
		plot_POTENTIAl1, = ax[0,0].plot(	0,0, '-o', ms=2, lw=2, ls='-', 
										color=(0.4,0.7,0.4), alpha=0.6, 
										label='Applied Potential change (dU/dt)')

		ax[0,0].set_xlabel('Time (s)')
		ax[0,0].set_ylabel('Potential (V)')
		ax[0,0].set_title('Potential Ramp')
		ax[0,0].set_xlim( (Tmin-0.1*Tsig, Tmax+0.1*Tsig) )
		ax[0,0].set_ylim( (Umin-0.1*Usig, Umax+0.1*Usig) )

		# === PLOT Reaction rate === # 	# === PLOT Reaction rate === # 	
		plot_RATE0, = ax[0,1].plot(	0,0, '-o', ms=2, lw=2, ls='-', 
									color=self.color[0], alpha=0.6, 
									label='Step 1 (k1)')
		plot_RATE1, = ax[0,1].plot(	0,0, '-o', ms=2, lw=2, ls='-', 
									color=self.color[1], alpha=0.6, 
									label='Step 2 (k2)')
		plot_RATE2, = ax[0,1].plot(	0,0, '-o', ms=2, lw=2, ls='-', 
									color=self.color[2], alpha=0.6, 
									label='Step 3 (k3)')
		plot_RATE3, = ax[0,1].plot(	0,0, '-o', ms=2, lw=2, ls='-', 
									color=self.color[3], alpha=0.6, 
									label='Step 4 (k4)')
		ax[0,1].set_xlabel('Time (s)')
		ax[0,1].set_ylabel('Kinetic constant')
		ax[0,1].set_title('Kinetic constant Norkov model')
		ax[0,1].set_xlim( (Tmin-0.1*Tsig, Tmax+0.1*Tsig) )
		ax[0,1].set_ylim( (Kmin-0.1*Ksig, Kmax+0.1*Ksig) )

		# === PLOT Reaction rate === # 	# === PLOT Reaction rate === # 	
		plot_G1, = ax[0,2].plot(	0,0, '-o', ms=2, lw=2, ls='-', 
									color=self.color[0], alpha=0.6, 
									label='G1 (eV)')
		plot_G1f, = ax[0,2].plot(	0,0, '-o', ms=2, lw=2, ls='-', 
									color=self.color[0], alpha=0.6, 
									label='G1 (eV)')
		plot_G2, = ax[0,2].plot(	0,0, '-o', ms=2, lw=2, ls='-', 
									color=self.color[1], alpha=0.6, 
									label='G2 (eV)')
		plot_G3, = ax[0,2].plot(	0,0, '-o', ms=2, lw=2, ls='-', 
									color=self.color[2], alpha=0.6, 
									label='G3 (eV)')
		plot_G4, = ax[0,2].plot(	0,0, '-o', ms=2, lw=2, ls='-', 
									color=self.color[3], alpha=0.6, 
									label='G4 (eV)')
		ax[0,2].set_xlabel('Step')
		ax[0,2].set_ylabel('Free energy')
		ax[0,2].set_title('Reaction profile | Potential dependant')
		ax[0,2].set_xlim( (-0.5, 8.0) )
		ax[0,2].set_ylim( (Gmin-4.1*Gsig, Gmax+2.1*Gsig) )


		# === PLOT Difussion === # 	# === PLOT Difussion === # 	
		plot_D1, = ax[1,0].plot(	0,0, '-o', ms=2, lw=2, ls='-', 
									color=self.color[0], alpha=0.6, 
									label='O2 CC.')
		plot_D2, = ax[1,0].plot(	0,0, '-o', ms=2, lw=2, ls='-', 
									color=self.color[1], alpha=0.6, 
									label='OH- CC.')

		ax[1,0].set_xlabel('Distance to electrode')
		ax[1,0].set_ylabel('CC')
		ax[1,0].set_title('Difussion profile | Fick law')
		ax[1,0].set_xlim( (-0.1, self.difussion.x +0.1) )
		ax[1,0].set_ylim( (Cmin-Csig*0.1, Cmax+.1*Csig) )

		# === PLOT reaction evolution === # 	# === PLOT reaction evolution === # 	
		plot_R1, = ax[1,1].plot(	0,0, '-o', ms=2, lw=2, ls='-', 
									color=self.color[0], alpha=0.6, 
									label='O2 CC.')
		plot_R2, = ax[1,1].plot(	0,0, '-o', ms=2, lw=2, ls='-', 
									color=self.color[1], alpha=0.6, 
									label='OH- CC.')

		ax[1,1].set_xlabel('Time (s)')
		ax[1,1].set_ylabel('CC')
		ax[1,1].set_title('CC on the electrode surface')
		ax[1,1].set_xlim( (Tmin-0.1*Tsig, Tmax+0.1*Tsig) )
		ax[1,1].set_ylim( (Cmin-Csig*0.1, Cmax+.1*Csig) )

		# === PLOT reaction evolution === # 	# === PLOT reaction evolution === # 	
		plot_V , = ax[1,2].plot(	0,0, '-o', ms=2, lw=2, ls='-', 
									color=(0.7,0.4,0.4), alpha=0.6, 
									label='Voltagram')
		ax[1,2].set_xlabel('Potential(U)')
		ax[1,2].set_ylabel('I')
		ax[1,2].set_title('Voltagram')
		ax[1,2].set_xlim( (Umin-0.1*Usig, Umax+0.1*Usig) )
		ax[1,2].set_ylim( (Imin-Isig*0.1, Imax+.1*Isig) )

		# === LABEL hansdler === #
		for n1 in range(2):
			for n2 in range(3):
				handles, labels = ax[n1,n2].get_legend_handles_labels()
				ax[n1,n2].legend(handles[::-1], labels[::-1])
				hl = sorted(zip(handles, labels),
				            key=operator.itemgetter(1))
				handles2, labels2 = zip(*hl)
				ax[n1,n2].legend(handles2, labels2)

		animation = FuncAnimation(fig, func=animation_frame, 
						frames=int(self.potential.time.shape[0]/time_step), interval=0)
		plt.show()


	def plot(self, ):
		import numpy as np
		import matplotlib.pyplot as plt
		from IPython.display import Image

		# Asumiendo que tienes tus 6 matrices de datos como m1, m2, m3, m4, m5, m6
		m1 = np.random.rand(10, 10)
		m2 = np.random.rand(10, 10)
		m3 = np.random.rand(10, 10)
		m4 = np.random.rand(10, 10)
		m5 = np.random.rand(10, 10)
		m6 = np.random.rand(10, 10)

		data = [m1, m2, m3, m4, m5, m6]  # Reemplaza esto con tus propios datos

		fig, ax = plt.subplots(2,3)  # Ajusta esto según el número de paneles que necesitas

		def animate(i):
		    axs[0, 0].imshow(data[i % 6])  # Actualiza el panel 1
		    axs[0, 1].imshow(data[(i+1) % 6])  # Actualiza el panel 2
		    axs[0, 2].imshow(data[(i+2) % 6])  # Actualiza el panel 3
		    axs[1, 0].imshow(data[(i+3) % 6])  # Actualiza el panel 4
		    axs[1, 1].imshow(data[(i+4) % 6])  # Actualiza el panel 5
		    axs[1, 2].imshow(data[(i+5) % 6])  # Actualiza el panel 6

		ani = FuncAnimation(fig, animate, frames=20, interval=200)

		ani.save('animation.gif', writer='pillow', fps=60)

		plt.close()


	@timer
	def cookbook(self, ):
		pass

V = Voltagram()
V.initialize()


