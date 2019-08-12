# -*- coding: iso-8859-1 -*-
def twostr_func(wavelength, F_s, solarzenithangle,albedo_dif, 
			albedo_dir, temp_ground, w_0, g, tau_n, temp_c):
	"""
	This function is an implementation of the two-stream formalism. 
	It is based on Toon et al (1989). 
	This version updates deviates from the Toon et al formalism in that 
	it allows for separate direct and diffuse albedos 
	(i.e. allowing for dependence of the albedo on zenith angle).

	#Inputs:
	wavelength: wavelength center of bin, in nm.
	F_s: solar input at TOA, in erg/cm2/s/nm; np.pi*F_s=the incident solar flux at the TOA. This is equivalent to F_s from Toon et al, discussed just after Equation 3. 
	solarzenithangle: Angle of direct beam delative to the vertical, in radians. Rugheimer et al set the solar zenith angle to 60 degrees (np.pi/3).
	albedo_dif: diffuse radiation albedo.
	albedo_dir: direct radiation albedo.
	temp_ground: temperature of surface in K. Should be continuous with temperature of bottom layer of atmosphere.
	
	w_0[n]: single-scattering albedo, for a layer n. If 0, all extinction due to absorption. If 1, all extinction due to scattering.
	g[n]: scattering asymmetry parameter, for a given layer. If 0, particle scatters symmetrically (e.g. Rayleigh). 
	tau_n[n]: optical depth of each layer
	temp_c[n]: temperature at the upper edge of layer n. So, temp_c[0] is the temperature at the upper edge of the atmosphere, while temp_c[N] ends up corresponding to the temperature of the surface of the planet (or, assuming continuity, the bottom edge of the atmosphere).


	#Outputs:
	F_plus_tau0[n]=F_plus evaluated at tau=0 for every layer n
	F_plus_taumax[n]=F_plus evaluated at tau=tau_n[n] for every layer n
	F_minus_tau0[n]=F_minus evaluated at tau=0 for every layer n
	F_minus_taumax[n]=F_minus evaluated at tau=tau_n[n] for every layer n
	
	F_net[n]: net flux at the bottom of layer n. 
	AMEAN[n]: 4*np.pi*J_n, where J_n is the mean intensity at the bottom of layer n. 
	surface_intensity: an  estimate of the total (hemispherically-integrated) intensity received by a point at the surface of the planet. It is equal to the direct intensity plus half the diffuse intensity.

	#Calling example:
	F_plus_tau0, F_plus_taumax, F_minus_tau0, F_minus_taumax, F_net, AMEAN, surface_intensity=twostr_func.twostream_toon_fun(wavelength, F_s, solarzenithangle,albedo_diffuse, albedo_direct,  w_0, g, tau_n, temp_c)

	Notes:
	1. This model includes both solar (direct beam) and blackbody emission, for completeness. In practice, BB flux is not relevant at non-solar wavelengths, and typical codes (twostr.f included) have separate codes for the solar and terrestrial emission wavelength regimes.
	2. This model does not track energy deposition and hence cannot be used for climate calculations. This will come later.
	3. This function uses a subfunction (defined at the bottom of the file) to calculate the Planck flux. 
	"""
	
	########################
	###Import useful libraries
	########################
	import numpy as np
	import pdb
	import scipy.linalg




	########################
	###Define model parameters
	########################
	#Properties of the ground
	emissivity_ground=1.-albedo_dif #emissivity of ground. 1=perfect BB emitter.

	#Optical depth structure
	Nlayer=len(tau_n) #number of layers in the atmospheric model.
	
	tau_c=np.zeros(Nlayer+1)# tau_c[n] is the cumulative optical depth at the upper edge of layer n. So tau_c[0]=0, and tau_c[N] is the maximum possible.
	for n in range(0, Nlayer):
		tau_c[n+1]=tau_c[n]+tau_n[n] 

	#In the Toon formalism, j=0 corresponds to space, and j=N+1 corresponds to the planet surface.
	#These points in wavelength space define the edges of the bins in tau space. 
	#Other terminology:
	#	tau_c=cumulative optical depth of layers *above* layer n. 
	#	tau_n=total optical depth of the layer n
	#	tau=total optical depth at any point within a layer n, hence satisfying 0<tau<tau_n

	mu_0=np.cos(solarzenithangle) #"incident direction of solar beam"


	########################
	###Determine the two-stream approximation coefficients.
	########################
	#Eddington and quadrature are good at solar wavelengths (i.e., not thermal blackbody dominated). delta scalings of Joseph et al (1976) recommended to replace w_0, g, tau in this case. However, when dominated by internal isotropic sources like the Planck function, hemispheric mean approximation is preferable. When w_0=0, quadrature case has problems. This happens esp at thermal wavelengths. Again this favors using hemispheric mean at these wavelengths
	
	#We use quadrature because 1) we are at solar wavelengths for this UV work and 2) that's what twostr.f does (which is our comparison case)
	gamma_1= np.sqrt(3.)*(2.-w_0*(1.+g))/2. #consistent with Toon et al; consistent with Pierrehumbert gamma_1
	gamma_2=np.sqrt(3.)*w_0*(1.-g)/2. #consistent with Toon et al; consistent with Pierrehumbert gamma_2
	gamma_3=(1.-np.sqrt(3.)*g*mu_0)/2. #consistent with Toon et al; equal to the Pierrehumbert gamma_plus/w_0
	gamma_4=1.-gamma_3 #consistent with Toon et al; equal to the Pierrehumbert gamma_minus/w_0
	mu_1=1./np.sqrt(3.)+np.zeros(np.shape(gamma_1))#In Toon paper (eqn 18), this is given by: (1.-w_0)/(gamma_1-gamma_2). For the quadrature approximation, it is 1./np.sqrt(3.). Given its use, it seems to relate most closely to gamma_B from Pierrehumbert (see eqs 5.27, 5.30)

	##Eddington
	#gamma_1= (7.-w_0*(4.+3.*g))/4.
	#gamma_2=-1.*(1.-w_0*(4.-3.*g))/4.
	#gamma_3=(2.-3.*g*mu_0)/4.
	#gamma_4=1.-gamma_3 #consistent with Toon et al; equal to the Pierrehumbert gamma_minus/w_0
	#mu_1=1./2.+np.zeros(np.shape(gamma_1))#In Toon paper (eqn 18), this is given by: (1.-w_0)/(gamma_1-gamma_2). For the quadrature approximation, it is 1./np.sqrt(3.). Given its use, it seems to relate most closely to gamma_B from Pierrehumbert (see eqs 5.27, 5.30)

	alambda=np.sqrt(np.abs(gamma_1*gamma_1-gamma_2*gamma_2)) #this is the lower-case lambda, from eqn 21 of Toon et al
								 #The absolute value was added based on the code Toon just sent us. This corresponds to his AK(L,J) parameter. But it should not matter since gamma_1>gamma_2 for w_0<1.
	clambda=(gamma_1-alambda)/(gamma_2) #this is the upper-case lambda, from eqn 22 of Toon et al

	EMLT=np.exp(-alambda*tau_n) #this appears to be a prefactor used to facilitate computation of eqn 44 of Toon et al
	e1=1.+clambda*EMLT
	e2=1.-clambda*EMLT
	e3=clambda+EMLT
	e4=clambda-EMLT

	########################
	###Set up calculation
	########################
	"""
	The fundamental equation we are solving is of form:
	A_{l}*Y_{l-1}+B_{l}*Y_{l}+D{l+1}=E_{l}  (equation 39 of Toon et al)
	Here, A_l, B_l, D_l, E_l are quantities we determine, and the Y_l is what we solve for.
	Hence, we can summarize that we are solving a matrix equation that takes form:
	PY=E
	where Y[l]=Y_l
	      E[l]=E_l
	      P[l, l-1]=A_l [row, column]
	      P[l, l]=B_l
	      P[l, l+1]=D_l
	      P[i,j]=0 else
	Toon et al use 1-indexing. Hence n runs from 1 to N, l runs from 1 to 2N, where N is the number of layers, and they have:
	Y_l=Y_{1n} for l=1,3,5,...2n-1...2N-1
	Y_l=Y_{2n} for l=2,4,6,...2n...2N

	However, we use Python, which has 0-indexing. Hence *we* choose that n runs from 0 to N-1, l runs from 0 to 2N-1, and:
	Y_l=Y_{1n} for l=0,2,4...2n...2N-2
	Y_l=Y_{2n} for l=1,3,5...2n+1...2N-1

	The Y_{1n} and Y_{2n} are related to F^+_n and F^-_n via equations 31 and 32 of Toon et al.
	This parametrization has been done to remove exponentials with positive operands (ie ones that could grow large and lead to numerical instabilities) from the matrix.

	Note: The mapping of this PQ=R to the F+ and F- is unclear because of 1) this parametrization in terms of Y_l (done to eliminate numerical instabilities) and 2)further linear combinations done to convert a pentagiagonal matrix to an even simpler tridiagonal matrix. Hence intuitive checks are hard.
	"""

	########################
	###Set up surface flux
	########################
	S_sfc=albedo_dir*mu_0*np.exp(-tau_c[-1]/mu_0)*np.pi*F_s+emissivity_ground*np.pi*Planck(temp_ground, wavelength)
	#Surface emission. Formed by adding blackbody emission from the ground to the reflected energy from the direct beam. The direct beam's reflected energy is assumed to be purely diffuse. This corresponds to equations 37 and 38 of Toon et al. Note that this does NOT match equation 5.31 of Pierrehumbert because it does not include the reflected diffuse radiation. So, this implicitly assumes the diffuse albedo to be 0. 

	########################
	###Set up C-values
	########################
	#In the reshuffled set of parameters used in this formalism, these seem analagous to the forcing term in Pierrehumbert. All the added radiation is contained in here.

	def C_plus(n, tau): #implementation of superposition of eqns 23 and 27 from Toon et al
		solarrad_denominator=alambda[n]**2.-1./mu_0**2.
		solarrad_prefactor=w_0[n]*F_s*np.pi
		solarrad_exponential=np.exp(-1.*(tau_c[n]+tau)/mu_0)
		solarrad_factor=((gamma_1[n]-1./mu_0)*gamma_3[n]+gamma_4[n]*gamma_2[n])
		solarrad=solarrad_prefactor*solarrad_factor*solarrad_exponential/solarrad_denominator #units of flux: erg/s/cm2/nm
		
		blackbody_prefactor=2*np.pi*mu_1[n]
		B0n=Planck(temp_c[n], wavelength)
		B1n=(Planck(temp_c[n+1], wavelength)-B0n)/tau_n[n] #this is effectively a slope
		blackbody_factor=B0n+B1n*(tau+1./(gamma_1[n]+gamma_2[n]))
		blackbody=blackbody_prefactor*blackbody_factor #start with units of the Planck function, which are: erg/s/cm2/nm/sr. But multiplying by 2pi sr restores the units of flux. So can safely add them. 
		
		result=solarrad+blackbody
		return result

	def C_minus(n, tau): #implementation of superposition of eqns 24 and 27 from Toon et al
		solarrad_denominator=alambda[n]**2.-1./mu_0**2.
		solarrad_prefactor=w_0[n]*F_s*np.pi
		solarrad_exponential=np.exp(-1.*(tau_c[n]+tau)/mu_0)
		solarrad_factor=((gamma_1[n]+1./mu_0)*gamma_4[n]+gamma_3[n]*gamma_2[n])
		solarrad=solarrad_prefactor*solarrad_factor*solarrad_exponential/solarrad_denominator #units of flux: erg/s/cm2/nm
		
		blackbody_prefactor=2*np.pi*mu_1[n]
		B0n=Planck(temp_c[n], wavelength)
		B1n=(Planck(temp_c[n+1], wavelength)-B0n)/tau_n[n] #this is effectively a slope
		blackbody_factor=B0n+B1n*(tau-1./(gamma_1[n]+gamma_2[n]))
		blackbody=blackbody_prefactor*blackbody_factor #start with units of the Planck function, which are: erg/s/cm2/nm/sr. But multiplying by 2pi sr restores the units of flux. So can safely add them. 
		
		result=solarrad+blackbody
		return result

	########################
	###Calculate matrix coefficients
	#########################
	#initialize the A, B, D, and E.
	A=np.zeros(Nlayer*2)
	B=np.zeros(np.shape(A))
	D=np.zeros(np.shape(A))
	E=np.zeros(np.shape(A))


	#For l=0 (n=0) we have the boundary condition that the downward diffuse flux at the top of the first layer is equal to any incident diffuse downward flux. We set this to be zero.
	A[0]=0.
	B[0]=e1[0]
	D[0]=-1.*e2[0]
	E[0]=0.-1*C_minus(0,0) #This is really F_minus[0,0], i.e. we are assuming there is no downward diffuse flux from the top of the atmosphere.

	#for l=2N-1 (n=N-1), we have the boundary condition that the upward flux at the surface is the sume of the reflected downward diffuse flux and energy from any other sources (e.g. reflected direct beam, BB emission of the ground)/np.sqrt(3.)
	A[2*Nlayer-1]=e1[Nlayer-1]-albedo_dif*e3[Nlayer-1]
	B[2*Nlayer-1]=e2[Nlayer-1]-albedo_dif*e4[Nlayer-1]
	D[2*Nlayer-1]=0.
	E[2*Nlayer-1]=S_sfc-C_plus(Nlayer-1, tau_n[Nlayer-1])+albedo_dif*C_minus(Nlayer-1, tau_n[Nlayer-1])

	#There is a problem in the Toon paper. As written, the l=2n depends on e_n+1, running over the array edge. twostr.f resolves this by adopting a different mapping: their definition reduces to defining l=2(n+1) and running n from 0 to N-1. In this case, l=2 (The third value in the list of ls) depends on n=0 and n=1. This eliminates the overflow problem. We have implemented this below.
	
	##For n=1,2,3...N-1, l=2,4,6,...2N-2:
	for n in range(0, Nlayer-1):
		l=2*(n+1)
		A[l]=e2[n]*e3[n]-e4[n]*e1[n]
		B[l]=e1[n]*e1[n+1]-e3[n]*e3[n+1]
		D[l]=e3[n]*e4[n+1]-e1[n]*e2[n+1]
		
		E[l]=e3[n]*(C_plus(n+1, 0.)-C_plus(n, tau_n[n]))+e1[n]*(C_minus(n,tau_n[n])-C_minus(n+1,0.))


	#For n=0...N-2, l=1,3...2N-3:
	for n in range(0, Nlayer-1):
		l=2*n+1
		A[l]=e2[n+1]*e1[n]-e3[n]*e4[n+1]
		B[l]=e2[n]*e2[n+1]-e4[n]*e4[n+1]
		D[l]=e1[n+1]*e4[n+1]-e2[n+1]*e3[n+1]
		
		E[l]=e2[n+1]*(C_plus(n+1, 0.)-C_plus(n, tau_n[n]))-e4[n+1]*(C_minus(n+1, 0)-C_minus(n, tau_n[n])) #twostr.f has a -1*e_{4,n+1}. We have applied the same even though this is NOT what is written in the Toon et al paper. We have done this because Toon told us (6/26/2015) that there are some sign errors in the coefficients, and we currently trust the validated CLIMA code over the paper we know has errors in it. EDIT: Looking at the code Toon shared with us, he does the same. 


	########################
	###Assemble matrix equation components
	#########################
	P=np.zeros([Nlayer*2,Nlayer*2])

	#l=0: no "A" coefficient b/c l-1 has no meaning
	P[0,0]=B[0]
	P[0,1]=D[0]

	#l=2N-1: no "D" coefficient b/c l+1 has no meaning
	P[2*Nlayer-1,2*Nlayer-1-1]=A[2*Nlayer-1]
	P[2*Nlayer-1,2*Nlayer-1]=B[2*Nlayer-1]

	for l in range(1, Nlayer*2-1): #This populates the matrix P in PY=E. 
		P[l, l-1]=A[l]
		P[l,l]=B[l]
		P[l,l+1]=D[l]

	########################
	###Invert matrix
	#########################
	#Y=np.linalg.solve(P, E) #this is the Y_l
	
	#try using a specialized solver
	ab=np.zeros([3,2*Nlayer])
	ab[0,:]=np.append(0.0, np.diag(P, k=1))
	ab[1,:]=np.diag(P, k=0)
	ab[2,:]=np.append(np.diag(P, k=-1),0.0)
	#pdb.set_trace()
	Y=scipy.linalg.solve_banded((1,1), ab, E) #this is the Y_l


	########################
	###Convert from Y_l to Y_1n, Y_2n
	#########################
	#The Y_1n as defined in Toon et al correspond to l=1,3, 5...2N-1. Adjusting for the zero-indexing of Python as we have done, they instead correspond to l=0,2,...2N-2
	#The Y_2n as defined in Toon et al correspond to l=2,4,6...2N. Adjusting for Python zero-indexing as we have done, they instead correspond to l=1,3,5...2N-1.
	#For detail, see eq. 40.
	Y_1=np.zeros(Nlayer)
	Y_2=np.zeros(Nlayer)
	for n in range(0, Nlayer):
		Y_1[n]=Y[2*n]
		Y_2[n]=Y[2*n+1] 
		#last number called is Nlayer-1=N-1, so is consistent.
	
	########################
	###Convert from Y_1n, Y_2n to F_plus, F_minus
	#########################
	def F_plus(n,tau): #defined from Eqn 31 of Toon et al.
		term1=Y_1[n]*(np.exp(-alambda[n]*(tau_n[n]-tau))+clambda[n]*np.exp(-alambda[n]*tau))
		term2=Y_2[n]*(np.exp(-alambda[n]*(tau_n[n]-tau))-clambda[n]*np.exp(-alambda[n]*tau))
		term3=C_plus(n,tau)
		
		result=term1+term2+term3
		return result

	def F_minus(n, tau): #defined from Eqn 32 of Toon et al.
		term1=Y_1[n]*(clambda[n]*np.exp(-alambda[n]*(tau_n[n]-tau))+np.exp(-alambda[n]*tau))
		term2=Y_2[n]*(clambda[n]*np.exp(-alambda[n]*(tau_n[n]-tau))-np.exp(-alambda[n]*tau))
		term3=C_minus(n,tau)
		
		result=term1+term2+term3
		return result
	
	########################
	###Evaluate F_plus, F_minus at boundary edges
	#########################
	F_plus_tau0=np.zeros(np.shape(tau_n))
	F_plus_taumax=np.zeros(np.shape(tau_n))
	F_minus_tau0=np.zeros(np.shape(tau_n))
	F_minus_taumax=np.zeros(np.shape(tau_n))

	for n in range(0, Nlayer):
		F_plus_tau0[n]=F_plus(n, 0.)
		F_plus_taumax[n]=F_plus(n, tau_n[n])
		F_minus_tau0[n]=F_minus(n, 0.)
		F_minus_taumax[n]=F_minus(n, tau_n[n])


	########################
	###Convert from Y_1n, Y_2n to F_net, mean intensity.
	#########################
	#test if diffuse flux dominates over direct flux. If direct flux dominant, instead set mu_1=mu_0
	
	#if F_minus_taumax[-1]<mu_0*np.pi*F_s*np.exp(-tau_c[-1]/mu_0):
		#mu_1=np.zeros(np.shape(mu_1))+mu_0
	#mu_1=np.zeros(np.shape(mu_1))+mu_0
	
	F_net=np.zeros(np.shape(tau_n)) #defined from Eqn 48 of Toon et al. This quantity is the net flux at the BASE of layer n.
	for n in range(0, Nlayer):
		direct=mu_0*np.pi*F_s*np.exp(-(tau_c[n]+tau_n[n])/mu_0) #eqn 50 of Toon et al

		term1=Y_1[n]*(e1[n]-e3[n])
		term2=Y_2[n]*(e2[n]-e4[n])
		term3=C_plus(n, tau_n[n])-C_minus(n, tau_n[n])
		
		F_net[n]=term1+term2+term3 -direct

	AMEAN=np.zeros(np.shape(tau_n)) #defined from Eqn 49 of Toon et al. This is the equivalent of the quantity AMEAN in the twostr.f code. It is equal to 4*np.pi*J_n, where J_n is the mean intensity at the base of layer n. Hence this quantity AMEAN should be equal to the total intensity received by a point at the base of layer n. 
	for n in range(0, Nlayer):
		direct=mu_0*np.pi*F_s*np.exp(-(tau_c[n]+tau_n[n])/mu_0) #eqn 50 of Toon et al
	
		term1=Y_1[n]*(e1[n]+e3[n])
		term2=Y_2[n]*(e2[n]+e4[n])
		term3=C_plus(n, tau_n[n])+C_minus(n, tau_n[n])
		
		#AMEAN[n]=(1./mu_1[n])*(term1+term2+term3)+direct/mu_0	
		AMEAN[n]=(1./mu_1[n])*(F_plus_taumax[n]+F_minus_taumax[n])+direct/mu_0	
	
	########################
	###Compute "surface intensity"
	#########################	
	#"Surface intensity" refers to the total intensity that would be intercepted by a particle at the surface of the planet. Whereas the total intensity is equal to (F_plus[-1]+F_minus[-1])/mu_1+direct[-1]/mu_0, the surface intensity is instead equal to (F_minus[-1])/mu_1+direct[-1]/mu_0, i.e. the downwelling diffuse intensity (since the bottom intensity is cut out due to there being a planet there) plus the direct intensity
	
	surface_intensity=(F_minus_taumax[-1]/mu_1[-1])+(np.pi*F_s)*np.exp(-(tau_c[-1])/mu_0)
	
	########################
	###Return Result
	#########################
	#F_minus_tau0
	#np.max(np.abs((F_minus_taumax[:-1]-F_minus_tau0[1:]))/F_minus_tau0[1:])
	#np.max(np.abs((F_plus_taumax[:-1]-F_plus_tau0[1:]))/F_plus_tau0[1:])
	
	return (F_plus_tau0, F_plus_taumax, F_minus_tau0, F_minus_taumax, F_net, AMEAN, surface_intensity)





########################
###Define blackbody function
#########################
def Planck(T, wav):
	"""
	Implementation of the Planck function for brightness temperature B
	wav=wavelength in nm
	T=blackbody temp in K
	"""
	wav_cm=wav*1.e-7 #convert wavelengths from nm to cm.
	c=2.99792e10 #speed of light, in cm/s
	h=6.62607e-27#Planck constant, in erg*s
	kb=1.38065e-16#Boltzmann constant, in erg/K
	
	import numpy as np
	result_cm=(2.*h*c**2./wav_cm**5.)*1./(np.exp(h*c/(wav_cm*kb*T))-1) #ergs/cm^3/s/steradian 
	#Will return RunTime warnings for extremal values, which occur at these wavelengths. 
	result=result_cm*1.e-7 #convert to units of ergs/cm^2/nm/s/steradian 
	return result #result is in units of 


if __name__ == "__main__":
	
	wavelength 		 = 
	F_s 			 = 
	solarzenithangle = 
	albedo_dif 		 = 
	albedo_dir 		 = 
	temp_ground 	 = 
	w_0 			 = 
	g 				 = 
	tau_n 			 = 
	temp_c 			 = 
	
	output = twostr_func(wavelength, F_s, solarzenithangle,albedo_dif, 
						 albedo_dir, temp_ground, w_0, g, tau_n, temp_c)





