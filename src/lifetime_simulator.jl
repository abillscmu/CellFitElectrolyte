function simulate_standard_loadstep!(integrator, end_time, cache, cellgeometry, cathodeocv, anodeocv, input_type, input_value)
    while integrator.t < end_time
        p = integrator.p
        #if CCCV, have to do some weird control stuff.
        if input_type == 5
            Voltage,_,_ = CellFitElectrolyte.calc_voltage_new(integrator.u,@view(p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
            if (Voltage < p.p_phys.vfull * 0.8)
                integrator.opts.dtmax = 10.0
            else
                integrator.opts.dtmax = 5.0
            end
            step!(integrator)
            if ((Voltage >= p.p_phys.vfull) & ((p.p_phys.cccv_switch != true)))
                integrator.p.p_phys.cccv_switch = true
            end
            if ((((integrator.u[8] >= p.p_phys.ifull) & ((p.p_phys.cccv_switch == true)) & (input_type ==5))) & (p.p_phys. cccv_switch_2 != true))
                u = copy(integrator.u)
                #going to rest, have to reset integrator val :(
                u[8] = 0
                OrdinaryDiffEq.set_u!(integrator, u)
                integrator.p.p_phys.cccv_switch_2 = true
            end
        else
            integrator.opts.dtmax = end_time - integrator.t
            step!(integrator)
        end
        if integrator.sol.retcode != :Default
            return false
        end
    end
    return true
end

function simulate_normal_cycle!(integrator, cycle_array, cache, cellgeometry, cathodeocv, anodeocv, num_steps)
    #Normal array
    times = cycle_array[2:num_steps+1]
    times .+= integrator.t
    types = cycle_array[num_steps+2:num_steps+1+num_steps]
    values = cycle_array[num_steps+2+num_steps:num_steps+1+2*num_steps]
    #by definition, times[1] == integrator.t
    for t in (@view times[2:end])
        add_tstop!(integrator,t)
    end
    #Step through the cycle
    for step::Int in 1:num_steps-1
        p = integrator.p
        Temp = 320
        @pack! p.p_phys = Temp
        input_type = types[step]
        input_value = values[step]
        end_time = times[step+1]
        @pack! p.p_phys = input_type,input_value
        #make sure we can deal with the initial condition
        if (input_type == 0.0) | (input_type == 5.0)
            u = copy(integrator.u)
            u[8] = input_value
            OrdinaryDiffEq.set_u!(integrator, u)
        elseif input_type == 1.0
            u = copy(integrator.u)
            Voltage, _, _, _, _ = CellFitElectrolyte.calc_voltage_new(integrator.u,@view(integrator.p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
            u[8] = input_value/Voltage
            OrdinaryDiffEq.set_u!(integrator, u)
        end
        success = CellFitElectrolyte.simulate_standard_loadstep!(integrator, end_time, cache, cellgeometry, cathodeocv, anodeocv, input_type, input_value)
        if !(success)
            return success
        end
    end
    return true
end


function simulate_rpt!(integrator, cycle_array, cache, cellgeometry, cathodeocv, anodeocv)
    integrator.opts.dtmax = 10.0
    #RPT
    input_type = 3.0
    input_value = cycle_array[2]
    @pack! integrator.p.p_phys = input_type, input_value
    Voltage, _, _, _, _ = CellFitElectrolyte.calc_voltage_new(integrator.u,@view(integrator.p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
    u = copy(integrator.u)
    u[8] = input_value
    OrdinaryDiffEq.set_u!(integrator, u)
    while Voltage >= 2.5
        step!(integrator)
        if integrator.sol.retcode != :Default
            error("integrator failed")
        end
        Voltage, _, _, _, _ = CellFitElectrolyte.calc_voltage_new(integrator.u,@view(integrator.p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
    end
    input_type = 0.0
    input_value = 0.0
    @pack! integrator.p.p_phys = input_type, input_value
    end_time = integrator.t + cycle_array[4]
    add_tstop!(integrator, end_time)
    while integrator.t < end_time
        step!(integrator)
        if integrator.sol.retcode != :Default
            return error("integrator failed")
        end
        Voltage, _, _, _, _ = CellFitElectrolyte.calc_voltage_new(integrator.u,@view(integrator.p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
    end
    input_type = 5.0
    input_value = cycle_array[3]
    @pack! integrator.p.p_phys = input_type, input_value
    u = copy(integrator.u)
    u[8] = input_value
    Voltage, _, _, _, _ = CellFitElectrolyte.calc_voltage_new(integrator.u,integrator.p.p_phys,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
    OrdinaryDiffEq.set_u!(integrator, u)
    while integrator.u[8] < -0.01
        step!(integrator)
        Voltage, _, _, _, _ = CellFitElectrolyte.calc_voltage_new(integrator.u,@view(integrator.p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
        if integrator.sol.retcode != :Default
            error("integrator failed")
        end
        if ((Voltage >= integrator.p.p_phys.vfull) & ((integrator.p.p_phys.cccv_switch != true)))
            integrator.p.p_phys.cccv_switch = true
        end
    end
    end_time = integrator.t + cycle_array[5]
    add_tstop!(integrator, end_time)
    input_type = 0.0
    input_value = 0.0
    u = deepcopy(integrator.u)
    u[8] = input_value
    OrdinaryDiffEq.set_u!(integrator, u)
    @pack! integrator.p.p_phys = input_type, input_value
    while integrator.t < end_time
        step!(integrator)
        if integrator.sol.retcode != :Default
            false
        end
    end
    return true
end


function construct_odeproblem(f, u, p)
    vars = ones(length(u))
    vars[8] = 0
    mm = diagm(vars)
    func = ODEFunction(f, mass_matrix=mm)
    prob = ODEProblem(func,u,(0.0,Inf),p)
    integrator = init(prob,QNDF(autodiff=false),save_everystep=false, verbose=false)
    integrator.opts.maxiters = 1e7
    return integrator
end


