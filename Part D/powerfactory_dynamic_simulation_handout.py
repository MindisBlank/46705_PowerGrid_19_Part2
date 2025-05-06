# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 13:16:32 2025

@author: Benjamin Vilmann, PhD student at Technical University of Denmark 
"""
# ======================== CONFIGURATION ========================
pf_python_path = 'C:\\Program Files\\DIgSILENT\\PowerFactory 2023 SP4A\\Python'
project_name = 'SM_fault_analysis'
export_path = r'[YOUR]\\[PATH]\\[HERE]'

# ======================== SETUP ========================
# --- --- --- YOU SHOULD NOT TOUCH THE CODE HERE SETUP --- --- ---
# Import the powerfactory module
import sys
# The presumed python folder location in PowerFactory
# Read the current used python version
python_version = sys.version.split('.')[1] 
# Version requirements for the PowerFactory connection
assert int(python_version) in [7,8,9,10,11], "You must run python "
# Add the powerfactory python module path to your temporary system path 
sys.path.insert(0,f'{pf_python_path}\\3.{python_version}')
# Now you can import the powerfactory module
import powerfactory
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# ============= PowerFactory Automation Wrapper =============
class PF_Wrapper:
    def __init__(self,project_name, export_path):
        try:
            self.app = app = powerfactory.GetApplicationExt()
        except powerfactory.ExitError as error:
            print(error)
            print('error.code = %d' % error.code)
        
        # Show application window
        app.Show() 
        
        # activate project
        app.ActivateProject(project_name)
        
        # Get folders and windows
        self.study = app.GetProjectFolder('study')
        self.netdat = app.GetProjectFolder('netdat')
        self.export_path = export_path

        self.variations = app.GetProjectFolder('scheme')
        self.operations = app.GetProjectFolder('scen')

        # 
        self.study_case_name = None
        self.events = {}

        return

    def create_study_case(self,name:str,folder:str=None):
        if folder is not None:
            # Check if folder exists already
            if folder not in [obj.loc_name for obj in self.study.GetContents()]:
                folder_obj = self.study.CreateObject("IntFolder", folder)
            else:
                folder_obj = [obj for obj in self.study.GetContents() if obj.loc_name == folder][0]
            if name not in [obj.loc_name for obj in folder_obj.GetContents()]:
                self.study_case = study_case = folder_obj.CreateObject("IntCase", name)
            else:
                self.study_case = study_case = [obj for obj in folder_obj.GetContents() if obj.loc_name == name][0]
                
            self.study_case_name = study_case_name = f'{folder}_{name}' 
        else:
            self.study_case = study_case = self.study.CreateObject("IntCase", name)
            self.study_case_name = name 

        study_case.Activate()

        # Assign folders to the study case
        self.ldf = self.app.GetFromStudyCase('ComLdf') # Load flow
        self.inc = self.app.GetFromStudyCase('ComInc') # Initial conditions
        self.sim = self.app.GetFromStudyCase('ComSim') # Simulation
        self.evt = self.app.GetFromStudyCase('IntEvt') # Event handler
        self.res = self.app.GetFromStudyCase('ComRes') # Results
        
        # Get grid NB: Assuming only ONE grid in the the powerfactory model!
        grid = [obj for obj in pf.app.GetProjectFolder('netdat').GetContents() if obj.GetClassName() == 'ElmNet'][0]
        grid.Activate()

        # Fetch result object depending on already existence
        if len(self.app.GetCalcRelevantObjects(f'{self.study_case_name}.ElmRes', 0)) >= 1:
            self.results = self.app.GetCalcRelevantObjects(f'{self.study_case_name}.ElmRes', 0)[0]
        else:
            self.results = study_case.CreateObject("ElmRes", self.study_case_name)

        return

    def setup_test(self,name):
        self.setup_result()
        
        return

    def define_initial_conditions(self,**kwargs):        
        for key, val in kwargs.items():
            try:
                setattr(self.inc, key, val)
            except TypeError as e:
                setattr(self.inc, key, float(val))
        
        return
    
    def define_result_settings(self,**kwargs):

        for key, val in kwargs.items():
            try:
                setattr(self.res, key, val)
            except TypeError as e:
                setattr(self.res, key, float(val))
        
        self.res.f_name = f'{self.export_path}\\{self.study_case_name}.csv'
        self.res.pResult = self.results

        return

    def set_params(self,params:dict):
        for elm, vals in params.items():
            PF_object = self.app.GetCalcRelevantObjects(elm, 0)[0]         
            for key,val in vals.items():
                setattr(PF_object,key,eval(f'{type(val).__name__}({val})'))
        return

    def define_result_variables(self,result_variables:dict):
        for key, vals in result_variables.items():
            PF_object = self.app.GetCalcRelevantObjects(key, 0)[0]         
            for val in vals:
                self.results.AddVariable(PF_object, val)
        return

    def add_event(self,evt_type:str,name:str,time:float,shc_type=None,Rf:float=0,Xf:float=0,var=None,val=None):
        if evt_type in ['shc']:
            assert shc_type is not None, "shc_type must be in ['1phg', '2ph', '2phg', '3phg', 'clear']"
            self.events[name] = self.evt.CreateObject('EvtShc', f'{name}_{shc_type}')
            self.events[name].p_target = self.app.GetCalcRelevantObjects(name, 0)[0]  
            self.events[name].i_shc = ['3phg', '2ph', '1phg', '2phg', 'clear'].index(shc_type)
            self.events[name].time = time
            self.events[name].ZfaultInp = 0 # fault representation 0: Z(Rf, Xf)
            self.events[name].R_f = Rf
            self.events[name].X_f = Xf
        elif evt_type in ['param']:
            assert var is not None and val is not None, "Define 'var' and 'val' for parameter event"
            self.events[name] = self.evt.CreateObject('EvtParam', f'{name}_{var}={val}')
            self.events[name].p_target = self.app.GetCalcRelevantObjects(name, 0)[0]  
            self.events[name].variable = var
            self.events[name].value = str(val)
            self.events[name].time = time        

        return

    def simulate(self,sim_time:float,verbose:bool=True):
        if verbose: print(f'Simulating: {self.study_case_name}')
        # Get initial conditions (run load flow for dynamic simulations)
        self.inc.Execute()

        # Run test
        self.sim.tstop = sim_time
        self.sim.Execute()

        # Export results
        self.res.Execute()
        return

# ======================== LOAD PROJECT ========================
# --- --- --- YOU ARE WELCOME TO EDIT THE CODE BELOW --- --- ---
# Fetch the powerfactory wrapper
pf = PF_Wrapper(project_name,export_path)

# Define result variables
result_variables = {'B1.ElmTerm':['m:u0','m:u1','m:u2','m:phiu0','m:phiu1','m:phiu2'],
                    'B3.ElmTerm':['m:u0','m:u1','m:u2','m:phiu0','m:phiu1','m:phiu2'],
                    'Generator.ElmSym':['m:I:bus1:A','m:phii:bus1:A','m:I:bus1:B','m:phii:bus1:B','m:I:bus1:C','m:phii:bus1:C'],
                    }
params = {'Z.ElmSind':{'xrea':0.0,'rrea':0.0}}

# Create a study case
for study_case_name in ['3phg', '2ph', '1phg', '2phg']:
    pf.create_study_case(study_case_name,folder='faults')
    
    # Define initial conditions
    pf.define_initial_conditions(iopt_sim='ins',iopt_net = 'rst',dtemt=0.0001,iopt_fastchk=0,iopt_fastout=0,iopt_fastsol=0)
    pf.define_result_settings(iopt_exp=6,iopt_csel = 0, iopt_vars = 0,ciopt_locn=1,ciopt_head=1)
    pf.define_result_variables(result_variables)
    pf.set_params(params)
    
    # Add events
    pf.add_event('shc','B3.ElmTerm',.1,shc_type=study_case_name)
    
    # Simulate
    pf.simulate(10) # hereby calculate initial conditions, run simulation, and export results

#%%
# ======================== INVESTIGATE RESULTS ========================
files = [f for f in os.listdir(export_path) if '.csv' in f]
for file in files:

    print(f'Reading file: {file}')
    df = pd.read_csv(os.path.join(export_path,file),sep=';',header=[0,1]).set_index((file.split('.')[0],'b:tnow in s'))
    df.index.name = 't'

    # Change multiindex pandas dataframe into a dictionary of element related values
    dfs = {}
    for elm in ['B1','B3','Generator']:
        dfs[elm] = df.xs(elm,axis=1,level=0)
    
    # Plot generator values
    fig, ax = plt.subplots(1,1,dpi=200,sharex=True)
    fig.suptitle(file)
    for col in dfs['Generator'].columns:       
        ax.plot(dfs['Generator'][col],alpha=0.5,label='$I_A$ (PowerFactory)')

