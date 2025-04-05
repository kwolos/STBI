# %% import essential libraries
import numpy as np
import pandas as pd
from itertools import product
import sys, copy, os
from scipy.optimize import least_squares, differential_evolution
sys.path.append('./additionalFunctions/')
import PythonAdditionalFunctions as PAF

# %% turn off warnings occuring during loading .xlsx files
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='openpyxl')

# %% follow logs
import logging
logging.basicConfig(filename="./logsANG.log", filemode="a", level=logging.DEBUG,
                    format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')
logging.getLogger('matplotlib.font_manager').disabled = True
# %% reading and building tree
tree = {}
tree['treeStructure'], tree['branchingPoints'] = PAF.buildTreeCircle('./PythonPWAExtension_library/Data/biffurcationTree.xlsx')

# scaling arterial tree to patient height
# calculating approximate distance from the heart to the foot 
pathLeg = np.array([14,                         # common carotid  
                    17, 26, 27, 34, 36, 38, 40, # from aortic arch to illiac 
                    41, 42, 44, 47])            # from iliac to tibial

# gathering data
L = 0 
for i in range(len(pathLeg)):
    L = L + tree['treeStructure'][1, pathLeg[i]]

# scaling to 175 cm
L = 175

Patients = np.load('./data/dataANGE.npy', allow_pickle=True).item()
PatientsID = Patients.keys()

method = "PSO"
add_method = "LSQ"
n = 6

IfVasoData = pd.read_csv('./data/DataPatients.csv', delimiter=',', on_bad_lines='skip',)

# %% loop with computing 
for PatientID in PatientsID: 
    if not os.path.exists('./results'):
        os.makedirs('results')

    folderToSave = f"results/fitResult_{PatientID.replace('CH_', '')}.npy"
    mss0 = f"considered patient: {PatientID}"
   
    ID = int(PatientID.split('_')[1])
    print(f"Patient ID: {ID}, is empty: {IfVasoData[IfVasoData['ID'] == ID]['Vdosing']}")
    if (IfVasoData[IfVasoData['ID'] == ID]['Vdosing'].empty) | (IfVasoData[IfVasoData['ID'] == ID]['Vdosing'].isna().all()):
        mss0a = f"considered patient {ID} is not within the patietns with followed vasopressor dosage"
        logging.info(mss0a)
        print(mss0a)

    logging.info(mss0)
    print(mss0)
    GivenPatient = Patients[PatientID]

    # load measurements from arms and ankles
    CH = GivenPatient['PATIENT_MEASUREMENTS']
    ch_ = {i+1: 0 if isinstance(CH['CUFF_WAVEFORMS'][i+1], int) else 1 for i in range(4)}
    
    # I: check if we have at least one arm and one leg
    # II: only one arm and one leg
    if (isinstance(CH['CUFF_WAVEFORMS'][1], int) and isinstance(CH['CUFF_WAVEFORMS'][2], int))\
        or (isinstance(CH['CUFF_WAVEFORMS'][3], int) and isinstance(CH['CUFF_WAVEFORMS'][4], int)): 
        mss = f"Patient {PatientID} has no measurements from the two arms or legs\n=====================================\n\n"
        logging.info(mss)
        continue

    # if we don't have two arms - stop
    if ch_[1]==0 & ch_[2]==0:
        mss1 = f"Patient {PatientID} has no measurements in two arms"
        logging.info(mss1)
        continue

    # # II: if we have to arms, take only left
    # if ch_[1]==1 & ch_[2]==1:
    #     ch_[2] = 0 
    #     mss1a = f"Patient {PatientID} has measurements in two arms; we take measurement from left"
    #     logging.info(mss1a)
    
    # if we have two legs, take only left
    # if ch_[3] & ch_[4]:
    #     ch_[4] = 0 
    mss2 = f"Considered cuffs: {ch_}"
    logging.info(mss2)
    print(mss2)

    PatientCharacteristics = GivenPatient['PATIENT_DATA']
    data = PAF.loadDataANGE(PatientCharacteristics, CH, L)

    data['model']   = 2
    data['message'] = False 
    
    SPa = GivenPatient['PATIENT_DATA']['SP']
    DPa = GivenPatient['PATIENT_DATA']['DP']

    badSP = [k for k, v in SPa.items() if (v == -2147483648) | (v == 0)]
    badDP = [k for k, v in DPa.items() if (v == -2147483648) | (v == 0)]

    for k in badSP:
        SPa[k] = np.mean([v for k,v in SPa.items() if k not in badSP]) 
    for k in badDP: 
        DPa[k] = np.mean([v for k,v in DPa.items() if k not in badDP]) 

    data['SP'] = SPa
    data['DP'] = DPa

    # Scaling tree according to the scale factor
    treeP = copy.deepcopy(tree)
    treeP['treeStructure'][1:4,:] = treeP['treeStructure'][1:4,:]*data['ScaleFactor']
    
    # we need to round L to given precision
    aux = np.round(treeP['treeStructure'][1,:],2)
    treeP['treeStructure'][1,:] = aux

    # making sure that where to take measurement variable is ok
    indx = treeP['treeStructure'][8,:]>treeP['treeStructure'][1,:]
    treeP['treeStructure'][8,indx] = np.floor(treeP['treeStructure'][1,indx]/2)

    params = PAF.initializeParamsVaso(data, treeP)
    paramsToFit = ['tm', 'scaleRes', 'Emax', 'scaleComp', 'k3']

    simSettings = {}
    simSettings['secondsToSimulate'] = 8.0  # number of seconds to simulate 
    simSettings['pps'] = 1000               # Hz, number of points per second in the output 
    simSettings['dt']  = 2e-04

    x0, lb, ub = PAF.startingPointVaso(params, paramsToFit)
    lb, ub = np.array(lb), np.array(ub)
    method = "PSO"
    add_method = "LSQ"

    if method == "PSO":
        print("========== PSO fitting ============\n")
        print(f"== Fitting: {paramsToFit} ===")
        print("====================================\n")
        PSOoptions = {'n_particles': 50,
                      'iters': 1}
        joint_vars, res = PAF.Optimize_PSO(paramsToFit, PSOoptions, lb, ub,
                                           treeP, simSettings, data, nrm=True,
                                           n=n, ange=ch_, init=None)
    
    if method == "DIFF_EV":
        print("===== differential evolution fitting ======")
        print(f"== Fitting: {paramsToFit} ===")
        print("===========================================")
        
        bounds_real = list(zip(lb, ub))
        bounds_norm = (np.zeros_like(lb), np.ones_like(ub))
        res = differential_evolution(
                                    PAF.opt_funcFFT_scipy,
                                    bounds=bounds_norm, 
                                    strategy='best1bin',
                                    disp=True, 
                                    x0=x0,
                                    args=(
                                        treeP,
                                        simSettings,
                                        data,
                                        paramsToFit,
                                        bounds_real, 
                                        n,
                                        None, 
                                        ch_
                                        ),
                                    )
    
    if add_method == "LSQ":
        print("========== Additional method: LSQ ==========\n")
        bounds_real = list(zip(lb, ub))
        bounds_norm = (np.zeros_like(lb), np.ones_like(ub))

        try: 
            x0 = res.x
        except AttributeError:
            x0 = res

        res_lsq = least_squares(PAF.opt_funcFFT_scipy, 
                                x0 = x0,
                                diff_step=0.05,
                                verbose=2,
                                ftol=1e-6,
                                xtol=1e-6,
                                gtol=1e-6,
                                method='lm',
                                args=(
                                    treeP,
                                    simSettings,
                                    data,
                                    paramsToFit,
                                    bounds_real, 
                                    n,
                                    None, 
                                    ch_, 
                                    None
                                ),
                                )
        joint_vars = [res_lsq.x[i] * (bounds_real[i][1] - bounds_real[i][0]) + bounds_real[i][0] for i in range(len(res_lsq.x))]    

    err, pulse, P, Q = PAF.F_ANGE_FFT(x=joint_vars, tree=treeP, simSettings=simSettings, data=data, 
                                      paramsToFit=paramsToFit, CH=ch_, n=n)
    mss3 = f"Computed error: {np.sum([np.sum([e**2 for e in v]) for _, v in err.items()])}\n====================================================\n\n"
    logging.info(mss3)
    print(f"Computed error: {np.sum([np.sum([e**2 for e in v]) for _, v in err.items()])}")

    toSave = {}
    toSave['P']           = P 
    toSave['treeP']       = treeP
    toSave['dataP_1']     = data
    toSave['simSettings'] = simSettings
    toSave['paramsToFit'] = paramsToFit
    toSave['pulse']       = pulse
    toSave['params']      = joint_vars
    toSave['Q_inlet']     = Q[0, :]
    toSave['error']       = err

    np.save(folderToSave, toSave)