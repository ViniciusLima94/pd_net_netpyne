'''
	Potjans & Diesmann cortical network.
	The cell-type specific cortical microcircuit: relating structure and activity in a full-scale spiking network model.
	https://www.ncbi.nlm.nih.gov/pubmed/23203991
'''
from netpyne import specs, sim
import numpy as np 

netParams = specs.NetParams()


###############################################################################
# Network parameters
###############################################################################
'''
	Population size per layer
                    2/3e   2/3i   4e    4i    5e    5i    6e     6i 
''' 

table = np.array([[0.101,  0.169, 0.044, 0.082, 0.032, 0.,     0.008, 0. ],
                  [0.135,  0.137, 0.032, 0.052, 0.075, 0.,     0.004, 0.],
                  [0.008,  0.006, 0.050, 0.135, 0.007, 0.0003, 0.045, 0.],
                  [0.069,  0.003, 0.079, 0.160, 0.003, 0.,     0.106, 0.],
                  [0.100,  0.062, 0.051, 0.006, 0.083, 0.373,  0.020, 0.],
                  [0.055,  0.027, 0.026, 0.002, 0.060, 0.316,  0.009, 0.],
                  [0.016,  0.007, 0.021, 0.017, 0.057, 0.020,  0.040, 0.225],
                  [0.036,  0.001, 0.003, 0.001, 0.028, 0.008,  0.066, 0.144]])

n_layer = (np.array([20683, 5834, 21915, 5479, 4850, 1065, 14395, 2948])/10.).astype(int)
layer_labels = np.array(['L23e', 'L23i', 'L4e', 'L4i', 'L5e', 'L5i', 'L6e', 'L6i'])
bg_layer_labels = np.array(['L23e_', 'L23i_', 'L4e_', 'L4i_', 'L5e_', 'L5i_', 'L6e_', 'L6i_'])
bg_layer = (np.array([1600, 1500 ,2100, 1900, 2000, 1900, 2900, 2100])/10.).astype(int)
tau_syn = 0.5
'''
	Creanting cortical layers
'''
for i in range(len(n_layer)):
	netParams.popParams[layer_labels[i]] = {'cellType': 'RS', 'numCells': n_layer[i], 'cellModel': 'Izhi'}
	netParams.stimSourceParams[bg_layer_labels[i]] = {'type': 'NetStim', 'rate': 8*bg_layer[i], 'start': 0, 'noise': 1}

'''
	Cell mechanisms
'''
cellRule = {'conds': {'cellType': 'RS', 'cellModel': 'Izhi'},  'secs': {}} 						# cell rule dict
cellRule['secs']['soma'] = {'geom': {}, 'pointps': {}}  												# soma params dict
cellRule['secs']['soma']['geom'] = {'diam': 10.0, 'L': 10.0, 'cm': 31.831}  							# soma geometry
cellRule['secs']['soma']['pointps']['Izhi'] = {'mod':'Izhi2007b', 'C':1, 'k':0.7, 
	'vr':-60, 'vt':-40, 'vpeak':35, 'a':0.03, 'b':-2, 'c':-50, 'd':100, 'celltype':1}  					# soma hh mechanisms
netParams.cellParams['RSrule'] = cellRule  		
'''
	External Poisson input
'''
netParams.synMechParams['exc'] = {'mod': 'ExpSyn', 'tau': tau_syn, 'e': 0, }
netParams.synMechParams['inh'] = {'mod': 'ExpSyn', 'tau': tau_syn, 'e': -80} 

for i in range(len(n_layer)):
	conn_dir = bg_layer_labels[i] + '->' + layer_labels[i]
	netParams.stimTargetParams[conn_dir] = {
        	'source': bg_layer_labels[i],
        	'sec':'soma',
        	'loc': 0.5,
    		'weight': 1.5e-2,
    		'delay': .1,
        	'conds': {'pop':layer_labels[i] }
        }

'''
	Connect neurons
'''

for c in range(len(layer_labels)):
	for r in range(len(layer_labels)):
		nsyn = int(np.log(1.0-table[r][c])/np.log(1.0 - (1.0/float(n_layer[c]*n_layer[r]))))
		conn_dir = layer_labels[c] + '->' + layer_labels[r]
		if layer_labels[c][-1] == 'e':
			if layer_labels[c] == 'L4e' and layer_labels[r] == 'L23e':
				netParams.connParams[conn_dir] = { 
				        'preConds': {'pop': layer_labels[c]}, 
				        'postConds': {'pop': layer_labels[r]}, 
				        #'convergence': nsyn,
					'probability': nsyn/(n_layer[c]*n_layer[r]),
				        'weight': '2*max( 0, normal(1.5e-2,1.5e-3) )',                        
				        'delay':  'max( 0.1, normal(1.5, 0.75) )',   
				        #'threshold': -10,                                
				        'synMech': 'exc'}
			else:
				netParams.connParams[conn_dir] = { 
			        'preConds': {'pop': layer_labels[c]}, 
			        'postConds': {'pop': layer_labels[r]}, 
			        #'convergence': nsyn,
				'probability': nsyn/(n_layer[c]*n_layer[r]),    
			        'weight': 'max( 0, normal(1.5e-2, 1.5e-3) )',                        
			        'delay':  'max( 0.1, normal(1.5, 0.75) )',   
			        #'threshold':-10,                                
			        'synMech': 'exc'}
		if layer_labels[c][-1] == 'i':
			netParams.connParams[conn_dir] = { 
			        'preConds': {'pop': layer_labels[c]}, 
			        'postConds': {'pop': layer_labels[r]}, 
			        #'convergence': nsyn,
				'probability': nsyn/(n_layer[c]*n_layer[r]),            
			        'weight': '4*max( 0, normal(1.5e-2, 1.5e-3) )',                        
			        'delay':  'max( 0.1, normal(0.8, 0.40) )',      
			        #'threshold': -10,                             
			        'synMech': 'inh'}       

simConfig = specs.SimConfig()           # object of class SimConfig to store simulation configuration
simConfig.saveCellSecs=0 
simConfig.saveCellConns=0 
simConfig.gatherOnlySimData=0
#simConfig.printPopAvgRates = True
simConfig.duration = 1000                     # Duration of the simulation, in ms
simConfig.dt = 0.025                            # Internal integration timestep to use
simConfig.verbose = False                       # Show detailed messages
simConfig.recordCells = []
simConfig.recordTraces = {}#{'V_soma':{'sec':'soma','loc':0.5,'var':'v'}}  # Dict with traces to record
simConfig.recordStep = 0.1                      # Step size in ms to save data (e.g. V traces, LFP, etc)
simConfig.filename = 'model_output'  # Set file output name
simConfig.savePickle = False            # Save params, network and sim output to pickle file
simConfig.hParams = {'celsius': 36, 'v_init': -65.0, 'clamp_resist': 0.001}
simConfig.analysis['plotRaster'] = {'include': layer_labels, 'timeRange': [0,1000],'popRates': True, 'saveFig': 'PYR_raster.png', 'showFig': False}  
simConfig.analysis['plotTraces'] = False#{'include': [0, 100, 200]}                     # Plot recorded traces for this list of cells
sim.createSimulateAnalyze(netParams, simConfig)
