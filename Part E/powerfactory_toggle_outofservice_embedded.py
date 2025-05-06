import powerfactory 

# Get application
app = powerfactory.GetApplication()

# Get relevant PowerFactory element
trafo = app.GetCalcRelevantObjects("T2031-4031.ElmTr2")[0]

# Get "out of service"
x_service = trafo.GetAttribute('outserv')

# What is going on here?
x_service_inv = not bool(x_service)
trafo.SetAttribute('outserv',x_service_inv)
