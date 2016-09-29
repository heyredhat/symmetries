# SYMMETRIES
# usage: python symmetries.py GROUP_NAME POLYPHONY
# right now we have: dihedral1, dihedral2, diherdral3, dihedral4, dihedral5, and qubit
# POLYPHONY is an integer. beware of setting it too high!
# also, stupidly, it will crash if POLYPHONY is > than # of transformations in ensemble
############################################################################################
from pyo import *
import time
import math
import cmath
import random
import itertools
import numpy as np
import numpy.linalg as la
import sys
np.set_printoptions(precision=2)

############################################################################################
# Gives the parity of a permutation
def perm_parity(lst):
    parity = 1
    for i in range(0,len(lst)-1):
        if lst[i] != i:
            parity *= -1
            mn = min(range(i,len(lst)), key=lst.__getitem__)
            lst[i],lst[mn] = lst[mn],lst[i]
    return parity    

def powerset(iterable):
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))

def sigmoid(x):
  return 1 / (1 + math.exp(-x))

class Ensemble:
	# A state is some random n vector
	# Transformations are (names, n x n matrices, midi_codes, and images)
	# And you can specify what order of commutator
	# to calculate up to (aka how much polyphony is allowed)
	def __init__(self, state, transformations, commutator=3):
		self.n = len(transformations)
		self.state = state
		basic_names = []
		basic_matrices = []
		basic_midi_codes = []
		basic_images = []
		for t in transformations:
			basic_names.append(t[1])
			basic_matrices.append(t[0])
			basic_midi_codes.append(t[2])
			basic_images.append(t[3])
		self.transformation_names = []
		self.transformation_matrices = []
		self.transformation_midi_codes = []
		self.transformation_eigs = []
		self.transformation_images = []
		# So we unwrap the transformations
		# but now we need to find the commutators
		# between the transformations
		# and add those as transformations
		# because that's how we handle chords aka
		# playing more than one note at a time
		ops = range(len(basic_matrices))
		cchords = list(powerset(ops))
		m = 0
		for i in range(len(cchords)):
			if len(cchords[i]) > commutator:
				m = i
				break
		for chord in cchords[:m]:
			if chord != ():
				chord = list(chord)
				C = np.matrix(np.zeros((len(state), len(state))))
				ch = range(0, len(chord))
				for permutation in itertools.permutations(ch):
					if len(permutation) > 0:
						M = np.identity(len(state))
						for i in permutation:
							M = np.dot(M, basic_matrices[chord[i]])
						p = list(permutation)
						M = M*perm_parity(p[:])
						C = C + M
				self.transformation_names.append("".join([basic_names[i] for i in chord]))
				self.transformation_matrices.append(np.matrix(C))
				self.transformation_midi_codes.append([basic_midi_codes[i] for i in chord])
				self.transformation_eigs.append(la.eig(C))
				self.transformation_images.append([basic_images[i] for i in chord])

	def apply(self, midi_keys, verbose=True):
		current_state = "THE STATE: %s " % str(self.state)
		i = self.transformation_midi_codes.index(sorted(midi_keys))
		name = "TRANSFORMATION: %s" % self.transformation_names[i]
		current_name = name
		matrix = self.transformation_matrices[i]
		current_matrix = str(matrix)
		eigenvalues, eigenvectors = self.transformation_eigs[i]
		if verbose:
			print "APPLYING %s" % (name)
		if verbose:
			print "STATE: %s" % (self.state)
		if verbose:
			print "MATRIX:\n%s" % (matrix)
		notes = []
		current_projections = ""
		# we use the provided midi keys pressed
		# to get the associated transformation matrices
		# then we go through all the eigenvalues and eigenvectors of the transformation matrix
		# and we find the projection of the current state onto each eigenvector
		# finally, we actually update the state by applying the matrix
		# and then package it up and send it off
		for i, l in enumerate(eigenvalues):
			u = eigenvectors[:,i]
			p = np.inner(u.T, np.conjugate(self.state))
			p =  np.inner(p,np.conjugate(p))[0,0]
			notes.append((complex(l), complex(p)))
			if verbose:
				l_r, l_th = cmath.polar(l)
				p_r, p_th = cmath.polar(p)
				print "\t--> L: (%.2f, %.2f) | P: (%.2f, %.2f) | U: %s" % (l_r, l_th, p_r, p_th, u.T)
				current_projections += "(%.2f, %.2f) %s\n->(%.2f, %.2f)\n" % (l_r, l_th, u.T, p_r, p_th)
		self.state = self.state*matrix
		if verbose:
			print "TRANSFORMED STATE: %s" % (self.state)
			print
		current_transformed_state = "TRANSFORMED STATE:\n%s" % (self.state)
		return notes, [current_state, current_name, current_matrix, current_projections, current_transformed_state]

	def __str__(self):
		s = "{STATE: %s\nTRANSFORMATIONS:\n\n" % (self.state.T)
		for i in range(len(self.transformation_matrices)):
			s += "NAME: %s\nMATRIX:\n%s\nMIDI CODES: %s\n" % (self.transformation_names[i], self.transformation_matrices[i], self.transformation_midi_codes[i])
			eigenvalues, eigenvectors = self.transformation_eigs[i]
			for j, l in enumerate(eigenvalues):
				s += "\teigenvalue %s | eigenvector %s\n" % ("({0.real:.2f} + {0.imag:.2f}i)".format(l), eigenvectors[:,j].T)
			s += "\n"
		s += "}\n"
		return s

# Generates an ensemble whose transformation matrices
# are given by the dihedral group of order n
# dihedral 3, for example, is the symmetries of an equilaterial triangle
def dihedral_ensemble(n, commutator=3):
	transformations = []
	rotations = []
	reflections = []
	k = 60
	for i in range(n):
		c = 2.*math.pi*float(i)/float(n)
		r = np.matrix([[math.cos(c), -1*math.sin(c)],[math.sin(c), math.cos(c)]])
		s = np.matrix([[math.cos(c), math.sin(c)],[math.sin(c), -1*math.cos(c)]])
		transformations.append((r, "R%d" % i, k, "R%d.jpg" % i))
		k += 1
		transformations.append((s, "S%d" % i, k, "S%d.jpg" % i))
		k += 1
	return Ensemble(np.array([complex(1./math.sqrt(2)), complex(1./math.sqrt(2))]), transformations, commutator=commutator)

# Generates an ensemble whose transformation matrices
# are given by certain symmetries of a qubit:
# the pauli spin matrices (aka X, Y, Z rotations)
# the hadamard transform
# and some number r of phase rotation transforms
def qubit_ensemble(r=12):
	SI = np.matrix([[complex(1),complex(0)],[complex(0),complex(1)]])
	SX = np.matrix([[complex(0),complex(1)],[complex(1),complex(0)]])
	SY = np.matrix([[complex(0),complex(0,-1)],[complex(0,1), complex(0)]])
	SZ = np.matrix([[complex(1),complex(0)],[complex(0),complex(-1)]])
	H = (1./math.sqrt(2))*np.matrix([[complex(1), complex(1)],[complex(1), complex(-1)]])
	transformations = [(SI, "SI", 60, "I.jpg"), (SX, "SX", 61, "SX.jpg"), (SY, "SY", 62, "SY.jpg"), (SZ, "SZ", 63, "SZ.jpg"), (H, "H", 64, "H.jpg")]
	for i in range(1,r):
		R = np.matrix([[complex(1), complex(0)],[complex(0), cmath.exp(complex(0,2*math.pi*float(i)/float(r)))]])
		transformations.append((R, "R%d" % i, 64+i, "R.jpg"))
	return Ensemble(np.array([complex(1./math.sqrt(2)), complex(1./math.sqrt(2))]), transformations)

############################################################################################
# Starts up the pyo midi server
server = Server()
server.setMidiInputDevice(2) # ?
server.boot()
server.start()

# Constructs ensemble of your choice
print "CONSTRUCTING ENSEMBLE..."
ensemble, GROUP_NAME = None, None
if __name__ == "__main__":
	if len(sys.argv) != 3:
		ensemble, GROUP_NAME  = dihedral_ensemble(3, commutator=3), "dihedral3"
		print "Defaulting to dihedral3"
	else:
   		if sys.argv[1] == "dihedral1":
   			ensemble, GROUP_NAME = dihedral_ensemble(1, commutator=int(sys.argv[2])), "dihedral1"
   		elif sys.argv[1] == "dihedral2":
   			ensemble, GROUP_NAME  = dihedral_ensemble(2, commutator=int(sys.argv[2])), "dihedral2"
   		elif sys.argv[1] == "dihedral3":
   			ensemble, GROUP_NAME  = dihedral_ensemble(3, commutator=int(sys.argv[2])), "dihedral3"
   		elif sys.argv[1] == "dihedral4":
   			ensemble, GROUP_NAME  = dihedral_ensemble(4, commutator=int(sys.argv[2])), "dihedral4"
   		elif sys.argv[1] == "dihedral5":
   			ensemble, GROUP_NAME  = dihedral_ensemble(5, commutator=int(sys.argv[2])), "dihedral5"
   		elif sys.argv[1] == "qubit":
   			ensemble, GROUP_NAME  = qubit_ensemble(r=12), "qubit"

# Audio setup
midi = Notein(ensemble.n)
FREQ_BASE = 220
sines = []
active = []

############################################################################################
# Here lies the gui
import wx, os

class TestFrame(wx.Frame):
    def __init__(self, *args, **kwargs):
    	# Here we build up the interface
    	# loading up the images associated with
    	# the different transformations
    	# arranging the textboxes that display the data
    	# beware: images must be in ./images/group_name!
    	global ensemble
        wx.Frame.__init__(self, *args, **kwargs)
        self.SetBackgroundColour('white')
        self.Maximize(True)
        self.jpgs = GetJpgList("./images/"+GROUP_NAME)
        self.MaxImageSize = 400
        self.images = {}
        for i in range(len(self.jpgs)):
        	self.images[self.jpgs[i]] = self.load_image(self.jpgs[i])        
        bitmap = wx.EmptyBitmap(self.MaxImageSize, self.MaxImageSize)
        dc = wx.MemoryDC( bitmap )
        myColor = 'white'
        dc.SetBackground( wx.Brush( myColor ) )
        dc.Clear()
        dc.SelectObject( wx.NullBitmap )
        self.Image = wx.StaticBitmap(self, bitmap=bitmap)
        font = wx.Font(16, wx.MODERN, wx.NORMAL, wx.BOLD, underline=False, faceName="", encoding=wx.FONTENCODING_DEFAULT)
        self.CurrentState = wx.StaticText(self, wx.ID_ANY, label="", style=wx.ALIGN_LEFT)
        self.CurrentState.SetFont(font)
        self.CurrentName = wx.StaticText(self, wx.ID_ANY, label="", style=wx.ALIGN_LEFT)
        self.CurrentName.SetFont(font)
        self.CurrentMatrix = wx.StaticText(self, wx.ID_ANY, label="", style=wx.ALIGN_RIGHT)
        self.CurrentMatrix.SetFont(font)
        self.CurrentProjections = wx.StaticText(self, wx.ID_ANY, label="", style=wx.ALIGN_RIGHT)
        self.CurrentProjections.SetFont(font)
        self.CurrentTransformedState = wx.StaticText(self, wx.ID_ANY, label="", style=wx.ALIGN_RIGHT)
        self.CurrentTransformedState.SetFont(font)
        self.RandomState = wx.StaticText(self, wx.ID_ANY, label="", style=wx.ALIGN_RIGHT)
        self.RandomState.SetFont(font)
        box = wx.BoxSizer(wx.HORIZONTAL)
        lbox = wx.BoxSizer(wx.VERTICAL)
        lbox.Add(self.CurrentState, 1, wx.ALIGN_LEFT, 10)
        lbox.Add(self.CurrentName, 1, wx.ALIGN_LEFT, 10)
        lbox.Add(self.CurrentMatrix, 1, wx.ALIGN_LEFT, 10)
        box.Add(lbox, 1, wx.ALIGN_LEFT, 10)
        box.Add(self.Image, 1, wx.ALIGN_CENTER, 10)
        rbox = wx.BoxSizer(wx.VERTICAL)
        rbox.Add(self.CurrentProjections , 1, wx.ALIGN_RIGHT, 10)
        rbox.Add(self.CurrentTransformedState, 1, wx.ALIGN_RIGHT, 10)
        rbox.Add(self.RandomState, 1, wx.ALIGN_RIGHT, 10)
        box.Add(rbox, 1, wx.ALIGN_RIGHT, 10)
        self.SetSizerAndFit(box)
        wx.EVT_CLOSE(self, self.OnCloseWindow)

    def load_image(self, file):
    	Img = wx.Image(file, wx.BITMAP_TYPE_JPEG)
    	W = Img.GetWidth()
        H = Img.GetHeight()
        if W > H:
            NewW = self.MaxImageSize
            NewH = self.MaxImageSize * H / W
        else:
            NewH = self.MaxImageSize
            NewW = self.MaxImageSize * W / H
        Img = Img.Scale(NewW,NewH)
        return wx.BitmapFromImage(Img)

    def DisplayNext(self, event=None, chord=None, text=[]):
    	global ensemble
    	# updates the display given the current state of the ensemble
    	# does some bitmap magic
    	# to create images for commutator transformations
    	# by overlaying the images of transformations in the commutator
    	if chord != None:
    		filenames = ensemble.transformation_images[ensemble.transformation_midi_codes.index(sorted(chord))]
        	combo = wx.EmptyBitmap(self.MaxImageSize, self.MaxImageSize)
	        dc = wx.MemoryDC(combo)
	        myColor = 'white'
	        dc.SetBackground(wx.Brush(myColor))
	        dc.Clear()
	        parity = 1
	        for nex in filenames:
	        	nex_i = self.images["./images/"+GROUP_NAME+"/"+nex]
	        	nex_dc = wx.MemoryDC()
	        	nex_dc.SelectObject(nex_i)
	        	dc.Blit(0,0, self.MaxImageSize, self.MaxImageSize, nex_dc, 0, 0, rop=wx.XOR)
	        	nex_dc.SelectObject(wx.NullBitmap)
	        	parity = parity*-1
	        dc.SelectObject(wx.NullBitmap)
	        white = wx.EmptyBitmap(self.MaxImageSize, self.MaxImageSize)
	        dc = wx.MemoryDC(white)
	        myColor = 'white'
	        dc.SetBackground(wx.Brush(myColor))
	        dc.Clear()
	        nex_dc = wx.MemoryDC()
	        nex_dc.SelectObject(combo)
	        imm = wx.ImageFromBitmap(combo)
	        if parity == -1:
	        	nex_dc.Blit(0,0, self.MaxImageSize, self.MaxImageSize, dc, 0, 0, rop=wx.XOR)
	        	nex_dc.SelectObject(wx.NullBitmap)
	        dc.SelectObject(wx.NullBitmap)
	        self.Image.SetBitmap(combo)
	        if text != []:
	        	self.CurrentState.SetLabel(text[0])
	        	self.CurrentName.SetLabel(text[1])
	        	self.CurrentMatrix.SetLabel(text[2])
	        	self.CurrentProjections.SetLabel("EIGVALS, EIGVECS, PROJECTS:\n" + text[3])
	        	self.CurrentTransformedState.SetLabel("\n"+text[4])
	        	self.RandomState.SetLabel(text[5])
	        self.Layout()
	        self.Refresh()

    def OnCloseWindow(self, event):
        self.Destroy()

def GetJpgList(dir):
    jpgs = [f for f in os.listdir(dir) if f[-4:] == ".jpg"]
    return [os.path.join(dir, f) for f in jpgs]

class App(wx.App):
    def OnInit(self):
        frame = TestFrame(None, id=-1, title="SYMMETRIES", pos=(0,25),size=(600,600))
        self.SetTopWindow(frame)
        frame.Show(True)
        self.frame = frame
        return True

app = App(0)

############################################################################################
# Given a "note" to play (aka a pair of complex numbers)
# create the sine waves, play them, and add them to the bin of current sinewaves
# complex number A in polar form is a magnitude and an angle
# you hear two notes: the first is the magnitude*a base frequency
# the second is the first times the angle, squished between 0 and 1
# so the angle of the complex number gives the musical interval
# as for complex number B, its magnitude gives the amplitude of the two sine waves
# times a base amplitude, and its angle gives the length of time the note spends
# fading in and out
def note_handler(notes):
	global sines
	global active
	for l, p in notes:
		freq = float(abs(l))*FREQ_BASE
		phase = float(cmath.phase(l))/(2*math.pi)
		if phase < 0:
			phase = 1 - abs(phase)
		if phase == 0:
			phase = 1
		phase +=1
		amp =  float(abs(p))
		fade = (float(abs(cmath.phase(p)))/(2*math.pi))
		if fade < 0:
			fade = 1 - abs(fade)
		envelope_a = Fader(fadein=fade+0.1, fadeout=fade+0.1, dur=0.5, mul=(1.-1./(amp+1))*0.85, add=0)
		a = Sine(freq=freq, phase=0, mul=envelope_a, add=0).out(chnl=0)
		envelope_b = Fader(fadein=fade+0.1, fadeout=fade+0.1, dur=0.5, mul=(1.-1./(amp+1))*0.85, add=0)
		b = Sine(freq=freq*phase, phase=0, mul=envelope_b, add=0).out(chnl=1)
		envelope_a.play()
		envelope_b.play()
		sines.append((None, envelope_a, a))
		sines.append((None, envelope_b, b))

# Removes a sine wave from the bin of sines
def stop(pitch):
	global sines
	to_remove = []
	for i in range(len(sines)):
		if sines[i][0] == pitch:
			sines[i][1].stop()
			to_remove.append(i)
	for i in to_remove[::-1]:
		del sines[i]

# Call back for midi note press
# takes currently active keys
# and applies the corresponding transformation in the ensemble
# if it exists
# then updates the display
# NOTE: If you press key #59, aka the one just below the one we start from
# when assigning midi notes
# it generates a random state and replaces the current one with that
def on_note_on(voice):
	global app
	global sines
	global active
	active = list(set(active))
	pitch = int(midi.get("pitch", True)[voice])
	current_random = ""
	if pitch == 59:
		print "*RANDOM STATE*"
		ensemble.state = np.array([[complex(random.random(), random.random()) for i in range(len(ensemble.transformation_matrices[0]))]])
		print ensemble.state
		print
		current_random = "*RANDOM STATE*\n%s" % (ensemble.state)
		app.frame.RandomState.SetLabel(current_random)
		app.frame.Layout()
		app.frame.Refresh()
	if [pitch] in ensemble.transformation_midi_codes+[[59]]:
		if [pitch] != [59]:
			active.append(pitch)
		active = list(set(active))
		if len(active) > 0 and sorted(active) in ensemble.transformation_midi_codes:
			notes, text = ensemble.apply(active)
			if notes != None:
				note_handler(notes)
			app.frame.DisplayNext(chord=active, text=text+[current_random])

# Call back for midi note unpress
# basically the same as above, but removes a note
# from the currently active list
def on_note_off(voice):
	global sines
	global active
	pitch = int(midi.get("pitch", True)[voice])
	if [pitch] in ensemble.transformation_midi_codes:
		if pitch in active:
			active.remove(pitch)
		active = list(set(active))
		if len(active) > 0 and sorted(active) in ensemble.transformation_midi_codes:
			notes, text = ensemble.apply(active)
			if notes != None:
				note_handler(notes)
			app.frame.DisplayNext(chord=active, text=text+[""])

# Set up the call backs
note_on = TrigFunc(midi["trigon"], on_note_on, arg=range(ensemble.n))
note_off = TrigFunc(midi["trigoff"], on_note_off, arg=range(ensemble.n))

# And we're off!
print "GO!"
app.MainLoop()
