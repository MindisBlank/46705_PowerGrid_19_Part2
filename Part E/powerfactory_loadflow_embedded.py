import powerfactory as pf
import pandas as pd
import traceback

app = pf.GetApplication()
prj = app.GetActiveProject()
script = app.GetCurrentScript()

try:
    # ================= 1) INITIALIZE =================
    app.PrintInfo("Initialization: fetching PF objects")
    # (app, prj, script already fetched above)
except Exception as e:
    app.PrintError("Initialization error: " + str(e))
    app.PrintError(traceback.format_exc())
    raise

try:
    # ================= 2) BASE-CASE LOAD FLOW =================
    app.PrintInfo("Running base-case load flow")
    ldf = app.GetFromStudyCase("ComLdf")
    ldf.iopt_net = script.GetAttribute('load_flow_type')
    # … add any other ldf.* config calls here …
    if ldf.Execute() != 0:
        raise RuntimeError("Base-case load flow did not converge")
    app.PrintInfo("Base-case load flow succeeded")
except Exception as e:
    app.PrintError("Base-case load flow error: " + str(e))
    app.PrintError(traceback.format_exc())
    raise

# ================= 3) EXTRACT BASE-CASE RESULTS =================
try:
    app.PrintInfo("=== Extracting base-case node results ===")
    nodes = app.GetCalcRelevantObjects("*.ElmTerm")
    data_node = {
        'bus': [], 'v_pu': [], 'Pgen_MW': [], 'Qgen_MVar': [],
        'Pload_MW': [], 'Qload_MVar': []
    }
    for nd in nodes:
        try:
            data_node['bus'].append(nd.loc_name)
            data_node['v_pu'].append(nd.GetAttribute("m:u"))
            data_node['Pgen_MW'].append(nd.GetAttribute("m:Pgen"))
            data_node['Qgen_MVar'].append(nd.GetAttribute("m:Qgen"))
            data_node['Pload_MW'].append(nd.GetAttribute("m:Pload"))
            data_node['Qload_MVar'].append(nd.GetAttribute("m:Qload"))
        except Exception as e:
            app.PrintError(f"  Node data error ({nd.loc_name}): {e}")
            app.PrintError(traceback.format_exc())
    df_node = pd.DataFrame(data_node).set_index('bus')
except Exception as e:
    app.PrintError("Error extracting base-case node results: " + str(e))
    app.PrintError(traceback.format_exc())
    raise

try:
    app.PrintInfo("=== Extracting base-case line results ===")
    lines = app.GetCalcRelevantObjects("*.ElmLne")
    data_line = {'line': [], 'from_bus': [], 'to_bus': [], 'loading_%': []}
    for ln in lines:
        try:
            data_line['line'].append(ln.loc_name)
            data_line['from_bus'].append(ln.bus1.loc_name)
            data_line['to_bus'].append(  ln.bus2.loc_name)
            data_line['loading_%'].append( ln.GetAttribute("m:Loading") * 100 )
        except Exception as e:
            app.PrintError(f"  Line data error ({ln.loc_name}): {e}")
            app.PrintError(traceback.format_exc())
    df_line = pd.DataFrame(data_line).set_index('line')
except Exception as e:
    app.PrintError("Error extracting base-case line results: " + str(e))
    app.PrintError(traceback.format_exc())
    raise

try:
    app.PrintInfo("=== Extracting base-case transformer results ===")
    xfmrs = app.GetCalcRelevantObjects("*.ElmTr2")
    data_xf = {'xfmr': [], 'from_bus': [], 'to_bus': [], 'loading_%': []}
    for xf in xfmrs:
        try:
            data_xf['xfmr'].append(xf.loc_name)
            data_xf['from_bus'].append(xf.bus1.loc_name)
            data_xf['to_bus'].append(  xf.bus2.loc_name)
            data_xf['loading_%'].append(xf.GetAttribute("m:Loading") * 100)
        except Exception as e:
            app.PrintError(f"  Transformer data error ({xf.loc_name}): {e}")
            app.PrintError(traceback.format_exc())
    df_xf = pd.DataFrame(data_xf).set_index('xfmr')
except Exception as e:
    app.PrintError("Error extracting base-case transformer results: " + str(e))
    app.PrintError(traceback.format_exc())
    raise

try:
    # ================= 4) SAVE BASE-CASE TO EXCEL =================
    app.PrintInfo("Saving base-case results to Excel")
    with pd.ExcelWriter("nordic32_basecase_results.xlsx") as writer:
        df_node.to_excel(writer, sheet_name="Nodes")
        df_line.to_excel(writer, sheet_name="Lines")
        df_xf.to_excel(writer, sheet_name="Transformers")
    app.PrintInfo("Base-case results saved: nordic32_basecase_results.xlsx")
except Exception as e:
    app.PrintError("Error saving base-case results: " + str(e))
    app.PrintError(traceback.format_exc())
    raise

# ================= 5) CONTINGENCY (N-1) ANALYSIS =================
app.PrintInfo("=== Starting contingency (N-1) analysis ===")
cont_elems = app.GetCalcRelevantObjects("*.ElmLne")  # change/add classes as needed
cont_results = []

for elem in cont_elems:
    name = elem.loc_name
    orig_status = elem.outserv
    app.PrintInfo(f"Processing contingency: {name}")
    try:
        # take it out of service
        elem.outserv = 1
        # run LF
        if ldf.Execute() != 0:
            raise RuntimeError("Load-flow failed to converge")
        # collect metrics
        lines_act = app.GetCalcRelevantObjects("*.ElmLne")
        terms_act = app.GetCalcRelevantObjects("*.ElmTerm")
        max_loading = max(l.GetAttribute("m:Loading")*100 for l in lines_act)
        min_voltage = min(t.GetAttribute("m:u") for t in terms_act)
        cont_results.append({
            'contingency':    name,
            'status':         'OK',
            'max_loading_%':  max_loading,
            'min_voltage_pu': min_voltage
        })
    except Exception as e:
        app.PrintError(f"  Error on {name}: {e}")
        app.PrintError(traceback.format_exc())
        cont_results.append({
            'contingency': name,
            'status':      'ERROR',
            'max_loading_%': None,
            'min_voltage_pu': None
        })
    finally:
        # restore service
        elem.outserv = orig_status

try:
    # build DataFrame and save
    df_cont = pd.DataFrame(cont_results).set_index('contingency')
    df_cont.to_excel("nordic32_contingency_results.xlsx")
    app.PrintInfo("Contingency results saved: nordic32_contingency_results.xlsx")
except Exception as e:
    app.PrintError("Error saving contingency results: " + str(e))
    app.PrintError(traceback.format_exc())
    raise

app.PrintInfo("Script completed successfully")
