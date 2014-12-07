
'''A class which wraps data from a molecular dynamics simulation.'''

try:
	import MDAnalysis
	from amx.tools import *
	from numpy import *
except: print 'warning: import failure in simset'

class SimSet:

	def __init__(self):

		#---MDAnalysis universe
		self.universe = None
		#---time
		self.nframes = 0
		self.vecs,self.vecs_index = [],[]
		self.time_list = []
		#---files
		self.universe_structfile = ''
		self.universe_trajfile = ''
		
	def load_trajectory(self,files,resolution=None,lenscale=None):

		'''
		Load a molecular dynamics trajectory into the MembraneSet instance.

		This function takes a structure and trajectory file and creates a Universe instance from the
		MDAnalysis library which serves as a member of MembraneSet. The function also populates meta-data 
		inside the MembraneSet instance which keep track of the frames and timining of the trajectory.

		Args:
			files (tuple): A tuple containing the path to the structure (usually a GRO file) and trajectory 
			for the target simulation.

		Kwargs:
			resolution (str): Specify either "CGMD" or "AAMD" to denote the "graining" of simulation. This is
			only used in a few functions, for bookkeeping purposes.
		   
			lenscale (float): Set the natural lengthscale of the incoming data. Since this function uses the
			MDAnalysis package, the default value is 10, corresponding to 10 Angstroms, so that all 
			downstream units are in nanometers.
		   
		'''
		
		#---default lengthscale is 10 because many MD programs use Angstroms and we compute in nm
		self.lenscale = lenscale if lenscale != None else 10.
		if resolution != None: self.resolution = resolution
		if self.nframes != 0:
			status('status: clearing trajectory')
			self.vecs,self.vecs_index = [],[]
			self.nframes = 0
		self.universe = MDAnalysis.Universe(files[0],files[1])
		self.universe_structfile = files[0]
		self.universe_trajfile = files[1]
		self.nframes = len(self.universe.universe.trajectory)
		if hasattr(self.universe.trajectory[0],'time'): self.time_start = self.universe.trajectory[0].time
		if hasattr(self.universe.trajectory,'totaltime'): self.time_total = self.universe.trajectory.totaltime
		if hasattr(self.universe.trajectory,'dt'): self.time_dt = self.universe.trajectory.dt
		#---if trajectory is available check the timestamps
		if hasattr(self.universe,'trajectory'):
			self.time_list = [self.universe.trajectory[i].time 
				for i in range(len(self.universe.trajectory)) 
				if hasattr(self.universe.trajectory[i],'time')]
		status('status: the trajectory file has '+str(self.nframes)+' frames')

	def gotoframe(self,frameno):
		
		'''Iterate to another frame without delay.'''
		
		if len(self.universe.trajectory) == 1:
			print 'warning: only one frame is available'
		elif frameno == 0:
			frame = self.universe.trajectory[frameno]
		elif self.universe.trajectory.frame-1 < frameno:
			[self.universe.trajectory.next() for i in range(frameno-(self.universe.trajectory.frame-1))]
		elif self.universe.trajectory.frame-1 > frameno:
			frame = self.universe.trajectory[frameno]

	def vec(self,frameno):

		'''Returns box vectors and stores them for quicker retrieval.'''

		if frameno not in self.vecs_index:
			if self.universe.trajectory.frame-1 != frameno: self.gotoframe(frameno)
			vec = self.universe.dimensions[0:3]
			self.vecs.append(vec)
			self.vecs_index.append(frameno)
			return vec
		else: return self.vecs[self.vecs_index.index(frameno)]

	def identify_monolayers(self,atomdirectors,startframeno=0):
		'''General monolayer identifier function. Needs: names of outer, inner atoms on lipids.'''
		status('status: identifying monolayers')
		if 0: status('status: moving to frame '+str(startframeno))
		self.gotoframe(startframeno)
		pointouts = self.universe.selectAtoms(atomdirectors[0])
		pointins = [self.universe.selectAtoms(atomdirectors[j]).coordinates() 
			for j in range(1,len(atomdirectors))]
		whichlayer = [0 if i > 0.0 else 1 for i in pointouts.coordinates()[:,2] - mean(pointins,axis=0)[:,2]]
		monos = []
		#---monos separates the lipids by absolute index into monolayers from index zero
		monos.append([pointouts.resids()[i]-1 for i in range(len(whichlayer)) if whichlayer[i] == 0])
		monos.append([pointouts.resids()[i]-1 for i in range(len(whichlayer)) if whichlayer[i] == 1])
		#---monolayer rerack hack if some residue IDs are missing
		#---Nb this may affect the tilter, mesher, identify_residues, and batch_gr functions so beware
		if (max(monos[0]+monos[1])-min(monos[0]+monos[1])) != len(monos[0]+monos[1])-1:
			print 'warning: resorting the monolayer indices because there is a mismatch'
			#---reracker is a sorted list of all of the absolute indices
			reracker = list(sort(monos[0]+monos[1]))
			#---monos_rerack is a copy of reracker in relative indices which is separated into monolayers 
			monos_rerack = [[reracker.index(i) for i in monos[m]] for m in range(2)]
			self.monolayer_residues_abs = monos
			#---resids_reracker is in absolute units
			self.resids_reracker = reracker
			#---monolayer_residues is in relative units
			self.monolayer_residues = monos_rerack
		else: self.monolayer_residues = monos
		if len(monos[0]) != len(monos[1]):
			print 'warning: there is a difference in the number of lipids per monolayer'

	def identify_residues(self,selector):
		'''General monolayer identifier function. Needs: names of outer, inner atoms on lipids.'''
		self.monolayer_by_resid = []
		self.resnames = []
		self.resids = []
		#---selector holds the proposed resnames
		for sel in selector:
			selection = self.universe.selectAtoms('resname '+sel)
			if len([i-1 for i in selection.resids()]):
				#---resids holds the absolute residue numbers
				self.resids.append([i-1 for i in selection.resids()])
				self.resnames.append(sel)
		#---when re-racker has been defined by identify_monolayers then provide distinct residue numberings
		if not hasattr(self,'resids_reracker'): self.resids_reracker = []
		if self.resids_reracker != []: 
			self.resids_abs = self.resids
			newresids = [[self.resids_reracker.index(i) for i in j] for j in self.resids]
			self.resids = newresids
		else:
			self.resids_abs = self.resids
		#---monolayer_residues is in relative indices
		for monolayer in self.monolayer_residues:
			monolayer_resids = []
			#for resids in residue_ids:
			#---self.resids is in relative indices
			for resids in self.resids:
				monolayer_resids.append(list(set.intersection(set(monolayer),set(resids))))
			#---monolayer_by_resid is in relative indices
			self.monolayer_by_resid.append(monolayer_resids)
		if self.resids_reracker == []:
			self.monolayer_by_resid_abs = self.monolayer_by_resid
		else:
			#---pack absolute residue indices into monolayer_by_resid_abs
			self.monolayer_by_resid_abs = \
				[[[self.resids_reracker[r] for r in restype] 
				for restype in mono] for mono in self.monolayer_by_resid]

