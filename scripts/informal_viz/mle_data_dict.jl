using Turing, DynamicHMC, PythonPlot, JLD2, KernelDensity, DataFrames, PythonCall, CellFitElectrolyte, CellFitElectrolyte.OCV
np = pyimport("numpy")
#pygui(true)

FOLDERNAME = "results/outputs0117_elec/"



cycle_vec = []
mle_vec = []
cellgeometry = CellFitElectrolyte.cell_geometry()
cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
sym = :εₑ⁻
CELL = "VAH12"
for cycle in keys(data_dict[CELL]["distributions"])

    U = data_dict[CELL]["distributions"][cycle][sym].kde

    mle = U.x[argmax(U.density)]
    push!(cycle_vec, cycle)
    push!(mle_vec, mle)
end

figure(1)
clf()
scatter(cycle_vec, mle_vec)

