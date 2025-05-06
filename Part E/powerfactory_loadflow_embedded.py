import powerfactory as pf
import pandas as pd

# Get the PowerFactory application object
app = pf.GetApplication()

# Get the active project
prj = app.GetActiveProject()

# Get the active project
script = app.GetCurrentScript()

# ========== LOAD FLOW ==========
# Get the load flow (check the documentation for other studies)
ldf = app.GetFromStudyCase("ComLdf")

# LOAD FLOW CONFIGURATION (see load flow config window for details)
ldf.iopt_net = script.GetAttribute('load_flow_type')

# ... Configure the load flow for your 

# Execute the load flow
ldf.Execute()

# ========== GET DATA ==========
app.PrintInfo("Acquiring data from relevant objects")

#######################################################
# ASSESS NODE DATA
#######################################################

# Get node info
app.PrintInfo("Bus results")

# Get all terminal objects in PowerFactory
nodes = app.GetCalcRelevantObjects("*.ElmTerm")

# Initialzing object to obtain data
data_node = {'bus': [], 'v': [], 'Pgen': [], 'Qgen': [], 'Pload': [], 'Qload': []}

# Iterate through each terminal to store desired values
for node in nodes:
    data_node['bus'].append(node.loc_name)
    data_node['v'].append(node.GetAttribute("m:u"))
    data_node['Pgen'].append(node.GetAttribute("m:Pgen"))
    data_node['Qgen'].append(node.GetAttribute("m:Qgen"))
    data_node['Pload'].append(node.GetAttribute("m:Pload"))
    data_node['Qload'].append(node.GetAttribute("m:Qload"))

# Convert data object into a pandas DataFrame
data_node = pd.DataFrame(data_node).set_index('bus')

# Print table to console
app.PrintInfo(data_node)

#######################################################
# ASSESS LINE DATA
#######################################################
# Similar as before 
app.PrintInfo("Line data")

# Get all the line elements in PowerFactory
lines = app.GetCalcRelevantObjects(...)

# Initializing data object
data_line = {'line': [], ...}

# Iterate through lines
for line in lines:
    data_line['line'].append(line.loc_name)
    ...

# Convert data object into a pandas DataFrame
data_line = pd.DataFrame(data_line).set_index('line')

# Print table to console
app.PrintInfo(data_line)

#######################################################
# ASSESS TRANSFORMER DATA
#######################################################
# Implement same data acquisition for transformers...

# ========== SAVE DATA ==========
# save the pandas dataframes to excel
# ....