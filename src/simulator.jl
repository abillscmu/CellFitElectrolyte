function CC_discharge_untildead!(integrator, p, cathode_ocv, anode_ocv, iapp_discharge, v_bn, v_sn, T)
    ## Discharge
    u_curr::Array{T,1} = integrator.u
    input_type = 3
    input_value = iapp_discharge
    @pack! integrator.p = input_type, input_value

    #Trying to take big steps during the early part of discharge
    voltage = voltage_from_solution_scalar(p, integrator.u, cathode_ocv, anode_ocv)
    while voltage > 3.0#slightly hacky but works and is close enough
        set_proposed_dt!(integrator, 1000.0)
        step!(integrator)
        t2 = integrator.t
        voltage = voltage_from_solution_scalar(p, integrator.u, cathode_ocv, anode_ocv)
        u_curr = integrator.u
        q_max = u_curr[10]
        qbn = u_curr[4]
        q_max_bn = q_max * v_bn / (v_bn + v_sn)
        xbn = qbn / q_max_bn
    end

    #Now slow down as we get into the "curvy" part of the ocv
    while (xbn > 0.05) && (voltage > 2.5)
        set_proposed_dt!(integrator, 1.0)
        step!(integrator)
        voltage = voltage_from_solution_scalar(p, integrator.u, cathode_ocv, anode_ocv)
        u_curr = integrator.u
        q_max = u_curr[10]
        qbn = u_curr[4]
        q_max_bn = q_max * v_bn / (v_bn + v_sn)
        xbn = qbn / q_max_bn
    end
end

function rest_fortime!(integrator, rest_time)
    ## Rest 1
    input_type = 0
    input_value = 0.0
    @pack! integrator.p = input_type, input_value
    endtime = integrator.t + rest_time

    #Should be pretty quick -- rest time is easy
    while integrator.t <= endtime
        set_proposed_dt!(integrator, endtime - integrator.t)
        step!(integrator)
    end
end

function CCCV!(integrator,p,cathode_ocv,anode_ocv,iapp_charge,V_full)
    ## CCCV
    voltage = voltage_from_solution_scalar(p, integrator.u, cathode_ocv, anode_ocv)
    input_value = iapp_charge
    input_type = 3
    @pack! integrator.p = input_type, input_value
    #Fast as we can through CC Charge     
    while voltage < (2.5+0.9*(V_full - 2.5))
        step!(integrator)
        voltage = voltage_from_solution_scalar(p, integrator.u, cathode_ocv, anode_ocv)
    end
    while voltage < (V_full)
        set_proposed_dt!(integrator, 1.0)
        step!(integrator)
        voltage = voltage_from_solution_scalar(p, integrator.u, cathode_ocv, anode_ocv)
    end
    input_type = 2
    input_value = V_full
    iapp = integrator.u[end]
    @pack! integrator.p = input_type, input_value

    #Fast through the begining of CV
    while iapp < -0.2
        #set_proposed_dt!(integrator)
        step!(integrator)
        iapp = integrator.u[end]
    end
    #Slow Down at the end
    while iapp < -0.01
        set_proposed_dt!(integrator, 1)
        step!(integrator)
        iapp = integrator.u[end]
    end
end

function normal_step!(integrator,end_time)
    ## Normal cycle discharge and rest steps
    DT = end_time - integrator.t
    dt_first = sqrt(DT) 

    while integrator.t < end_time - dt_first
        set_proposed_dt!(integrator, dt_first)
        step!(integrator)
    end
    while integrator.t < end_time
        set_proposed_dt!(integrator, 1.0)

        step!(integrator)
    end
end