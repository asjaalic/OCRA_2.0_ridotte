# SOLVE OPTIMIZATION PROBLEM

function solveRollingPlanning(InputParameters::InputParam, Battery::BatteryParam, ResultsOpt::Results, vec_prices)

    @unpack (NYears, NMonths, NStages, Big, NHoursStep, disc, Hours_rolling, Hours_saved) = InputParameters;                #NSteps, NHoursStage
    @unpack (min_SOC, max_SOC, Eff_charge, Eff_discharge, min_P, max_P, max_SOH, min_SOH, Nfull,fix ) = Battery;
    @unpack (soc, cap, soc_quad) = ResultsOpt;

    println("Solving Optimization Problem Case 2")
    k_end = min_SOH/(2*Nfull)

    charge_end = []
    discharge_end = []
    state_end = []
    degradation_end = []
    capacity_end = []

    stateQuad_end = []
    x_end = []
    y_end = []
    z_end = []
    u_end = []
    w_xx_end = []
    w_yy_end = []
    w_zz_end = []
    w_xy_end = []
    w_xz_end = []
    w_zy_end = []
    w_uu_end = []
    w_xu_end = []
    w_yu_end = []
    w_zu_end = []

    push!(capacity_end, cap[end])
    push!(state_end, soc[end])
    push!(stateQuad_end, soc_quad[end])

    subset_price = zeros(Hours_rolling)
    plan=1824
    new_vec = []

    problem_end = BuildRolling(InputParameters, SolverParameters, Battery)

    for stage=1:plan        #plan

        subset_price[:] = vec_prices[Hours_saved*(stage-1)+1:Hours_saved*(stage+1)]

        if stage==1
            for i=1:Hours_rolling
            push!(new_vec,subset_price[i])
            end
        else
            for i=(Hours_saved+1):Hours_rolling
            push!(new_vec,subset_price[i])
            end
        end


            for j=1:Hours_rolling               # aggiorno i 48 valori di prezzo
                
                set_objective_coefficient(
                    problem_end.M_f,
                    problem_end.discharge_f[j],
                    NHoursStep*subset_price[j]
                )

                set_objective_coefficient(
                    problem_end.M_f,
                    problem_end.charge_f[j],
                    -NHoursStep*subset_price[j]
                )  

            end

            #AGGIORNARE LO STATE OF CHARGE E LA CAPACITY
            if stage==1
                JuMP.set_normalized_rhs(
                    problem_end.State_initial[1],
                    soc[end],        #input da problema ottimizzazione precedente
                )

                JuMP.set_normalized_rhs(
                    problem_end.initial_capacity_f[1],
                    cap[end],        #input da problema ottimizzazione precedente
                )
            else
                JuMP.set_normalized_rhs(
                    problem_end.State_initial[1],
                    state_end[end],
                )

                JuMP.set_normalized_rhs(
                    problem_end.initial_capacity_f[1],
                    capacity_end[end],
                )

            end


        @timeit to "Solve optimization" optimize!(problem_end.M_f)
        
        if termination_status(problem_end.M_f) != MOI.OPTIMAL
            println("NOT OPTIMAL: ", termination_status(problem.M_f))
        else
            println("Optimization finished")
        end
        

        @timeit to "Collecting results" begin

            if stage==plan
                for iStep=1:Hours_rolling
                    push!(state_end, JuMP.value(problem_end.soc_f[iStep+1]))
                    push!(charge_end, JuMP.value(problem_end.charge_f[iStep]))
                    push!(discharge_end, JuMP.value(problem_end.discharge_f[iStep]))
                    push!(degradation_end,JuMP.value(problem_end.deg_f[iStep]))     #*k_end
                    push!(capacity_end, JuMP.value(problem_end.capacity_f[iStep+1]))
                    push!(stateQuad_end, JuMP.value(problem_end.soc_quad_f[iStep+1]))
                end
            else
                for iStep=1:Hours_saved
                    push!(state_end, JuMP.value(problem_end.soc_f[iStep+1]))
                    push!(charge_end, JuMP.value(problem_end.charge_f[iStep]))
                    push!(discharge_end, JuMP.value(problem_end.discharge_f[iStep]))
                    push!(degradation_end,JuMP.value(problem_end.deg_f[iStep]))     #*k_end
                    push!(capacity_end, JuMP.value(problem_end.capacity_f[iStep+1]))
                    push!(stateQuad_end, JuMP.value(problem_end.soc_quad_f[iStep+1]))
                end
            end

        if capacity_end[end] <= 0.21*min_SOH
            break
        end

    end
    
end

    println("Collected results Rolling Case 2")

    return ResultsEndLife(
        problem_end,
        charge_end,
        discharge_end,
        state_end,
        stateQuad_end,
        degradation_end,
        capacity_end,
        x_end,
        y_end,
        z_end,
        u_end,
        w_xx_end,
        w_yy_end,
        w_zz_end,
        w_uu_end,
        w_xy_end,
        w_xz_end,
        w_zy_end,
        w_xu_end,
        w_yu_end,
        w_zu_end,
        new_vec,
    )

end





function solveCase1(InputParameters::InputParam, Battery::BatteryParam, ResultsOpt::Results, vec_p, Tot_steps)

    @unpack (NYears, NMonths, NStages, Big, NHoursStep, disc,) = InputParameters;                #NSteps, NHoursStage
    @unpack (min_SOC, max_SOC, Eff_charge, Eff_discharge, min_P, max_P, max_SOH, min_SOH, Nfull,fix ) = Battery;

    println("Solving Optimization Problem Case 1")

    K_e = min_SOH/(2*Nfull)

    charge_e1 = zeros(Tot_steps)
    discharge_e1 = zeros(Tot_steps)
    soc_e1 = zeros(Tot_steps+1)
    deg_e1 = zeros(Tot_steps)
    cap_e1 = zeros(Tot_steps+1)

    soc_quad_e1 = zeros(Tot_steps+1)
    x_e1 = zeros(Tot_steps+1)
    y_e1 = zeros(Tot_steps+1)
    z_e1 = zeros(Tot_steps+1)
    u_e1 = zeros(Tot_steps+1)
    w_xx_e1 = zeros(Tot_steps+1)
    w_yy_e1 = zeros(Tot_steps+1)
    w_zz_e1 = zeros(Tot_steps+1)
    w_xy_e1 = zeros(Tot_steps+1)
    w_xz_e1 = zeros(Tot_steps+1)
    w_zy_e1 = zeros(Tot_steps+1)
    w_uu_e1 = zeros(Tot_steps+1)
    w_xu_e1 = zeros(Tot_steps+1)
    w_yu_e1 = zeros(Tot_steps+1)
    w_zu_e1 = zeros(Tot_steps+1)

    problem_e1 = endlife(InputParameters, Battery, ResultsOpt, vec_p, Tot_steps)

    @timeit to "Solve optimization" optimize!(problem_e1.M_e)

    if termination_status(problem_e1.M_e) != MOI.OPTIMAL
        println("NOT OPTIMAL: ", termination_status(problem_e1.M_e))
    else
        println("Optimization finished")
    end

    @timeit to "Collecting results" begin

        
        for iStep=1:Tot_steps
            soc_e1[iStep] = JuMP.value(problem_e1.soc_e[iStep])
            charge_e1[iStep] = JuMP.value(problem_e1.charge_e[iStep])
            discharge_e1[iStep] = JuMP.value(problem_e1.discharge_e[iStep])
            deg_e1[iStep] = JuMP.value(problem_e1.deg_e[iStep])*K_e
            cap_e1[iStep] = JuMP.value(problem_e1.capacity_e[iStep])

            soc_quad_e1[iStep] = JuMP.value(problem_e1.soc_quad_e[iStep])
            x_e1[iStep] = JuMP.value(problem_e1.x_e[iStep])
            y_e1[iStep] = JuMP.value(problem_e1.y_e[iStep])
            z_e1[iStep] = JuMP.value(problem_e1.z_e[iStep])
            u_e1[iStep] = JuMP.value(problem_e1.u_e[iStep])
            w_xx_e1[iStep] = JuMP.value(problem_e1.w_xx_e[iStep])
            w_yy_e1[iStep] = JuMP.value(problem_e1.w_yy_e[iStep])
            w_zz_e1[iStep] = JuMP.value(problem_e1.w_zz_e[iStep])
            w_xy_e1[iStep] = JuMP.value(problem_e1.w_xy_e[iStep])
            w_xz_e1[iStep] = JuMP.value(problem_e1.w_xz_e[iStep])
            w_zy_e1[iStep] = JuMP.value(problem_e1.w_zy_e[iStep])

            w_uu_e1[iStep] = JuMP.value(problem_e1.w_uu_e[iStep])
            w_xu_e1[iStep] = JuMP.value(problem_e1.w_xu_e[iStep])
            w_yu_e1[iStep] = JuMP.value(problem_e1.w_yu_e[iStep])
            w_zu_e1[iStep] = JuMP.value(problem_e1.w_zu_e[iStep])

        end

        soc_e1[end] = JuMP.value(problem_e1.soc_e[end])
        soc_quad_e1[end] = JuMP.value(problem_e1.soc_quad_e[end])
        x_e1[end] = JuMP.value(problem_e1.x_e[end])
        y_e1[end] = JuMP.value(problem_e1.y_e[end])
        z_e1[end] = JuMP.value(problem_e1.z_e[end])
        u_e1[end] = JuMP.value(problem_e1.u_e[end])
        w_xx_e1[end] = JuMP.value(problem_e1.w_xx_e[end])
        w_yy_e1[end] = JuMP.value(problem_e1.w_yy_e[end])
        w_zz_e1[end] = JuMP.value(problem_e1.w_zz_e[end])
        w_xy_e1[end] = JuMP.value(problem_e1.w_xy_e[end])
        w_xz_e1[end] = JuMP.value(problem_e1.w_xz_e[end])
        w_zy_e1[end] = JuMP.value(problem_e1.w_zy_e[end])

        w_uu_e1[end] = JuMP.value(problem_e1.w_uu_e[end])
        w_xu_e1[end] = JuMP.value(problem_e1.w_xu_e[end])
        w_yu_e1[end] = JuMP.value(problem_e1.w_yu_e[end])
        w_zu_e1[end] = JuMP.value(problem_e1.w_zu_e[end])

        cap_e1[end] = JuMP.value(problem_e1.capacity_e[end])
    
        

    end
    
    println("Collected results Case 1")

    return ResultsCase1(
        soc_e1,
        charge_e1,
        discharge_e1,
        deg_e1,
        soc_quad_e1,
        x_e1,
        y_e1,
        z_e1,
        u_e1,
        w_xx_e1,
        w_yy_e1,
        w_zz_e1,
        w_uu_e1,
        w_xy_e1,
        w_xz_e1,
        w_zy_e1,
        w_xu_e1,
        w_yu_e1,
        w_zu_e1,
        cap_e1,  
    )

end

function solveNewDegradation(InputParameters::InputParam, Battery::BatteryParam, ResultsOpt::Results, prices)

    @unpack (NYears, NMonths, NStages, Big, NHoursStep, disc, Hours_rolling, Hours_saved) = InputParameters;                #NSteps, NHoursStage
    @unpack (min_SOC, max_SOC, Eff_charge, Eff_discharge, min_P, max_P, max_SOH, min_SOH, Nfull,fix ) = Battery;
    @unpack (soc, cap, soc_quad) = ResultsOpt;

    println("Solving Optimization Problem Case 3")

    charge_new_deg = []
    discharge_new_deg = []
    state_new_deg = []
    degradation_new_deg = []
    capacity_new_deg = []

    stateQuad_new_deg = []
    x_new_deg = []
    y_new_deg = []
    z_new_deg = []
    u_new_deg = []
    w_xx_new_deg = []
    w_yy_new_deg = []
    w_zz_new_deg = []
    w_xy_new_deg = []
    w_xz_new_deg = []
    w_zy_new_deg = []
    w_uu_new_deg = []
    w_xu_new_deg = []
    w_yu_new_deg = []
    w_zu_new_deg = []

    push!(capacity_new_deg, cap[end])
    push!(state_new_deg, soc[end])
    push!(stateQuad_new_deg, soc_quad[end])

    subset_price_new = zeros(Hours_rolling)
    plan=1824
    new_vec_deg = []

    problem_new = New_deg_formulation(InputParameters, SolverParameters, Battery)

    for stage=1:plan        #plan

        subset_price_new[:] = prices[Hours_saved*(stage-1)+1:Hours_saved*(stage+1)]

        if stage==1
            for i=1:Hours_rolling
            push!(new_vec_deg,subset_price_new[i])
            end
        else
            for i=(Hours_saved+1):Hours_rolling
            push!(new_vec_deg,subset_price_new[i])
            end
        end


            for j=1:Hours_rolling               # aggiorno i 48 valori di prezzo
                
                set_objective_coefficient(
                    problem_new.M_new,
                    problem_new.discharge_new[j],
                    NHoursStep*subset_price_new[j]
                )

                set_objective_coefficient(
                    problem_new.M_new,
                    problem_new.charge_new[j],
                    -NHoursStep*subset_price[j]
                )  

            end

            #AGGIORNARE LO STATE OF CHARGE E LA CAPACITY
            if stage==1
                JuMP.set_normalized_rhs(
                    problem_new.State_initial_new[1],
                    soc[end],        #input da problema ottimizzazione precedente
                )

                JuMP.set_normalized_rhs(
                    problem_new.initial_capacity_new[1],
                    cap[end],        #input da problema ottimizzazione precedente
                )
            else
                JuMP.set_normalized_rhs(
                    problem_new.State_initial_new[1],
                    state_new_deg[end],
                )

                JuMP.set_normalized_rhs(
                    problem_new.initial_capacity_new[1],
                    capacity_new_deg[end],
                )

                for a=1:Hours_rolling

                    JuMP.set_normalized_coefficient(
                        problem_new.deg_new_1[a],
                        soc_quad_new[a],
                        -capacity_new_deg[end]/(2*Nfull*state_new_deg[end]^2)
                    )

                    JuMP.set_normalized_coefficient(
                        problem_new.deg_new_1[a],
                        soc_quad_new[a+1],
                        +capacity_new_deg[end]/(2*Nfull*state_new_deg[end]^2)
                    )

                    JuMP.set_normalized_coefficient(
                        problem_new.deg_new_1[a],
                        soc_new[a+1],
                        -capacity_new_deg[end]/(Nfull*state_new_deg[end])
                    )

                    JuMP.set_normalized_coefficient(
                        problem_new.deg_new_1[a],
                        soc_new[a],
                        +capacity_new_deg[end]/(Nfull*state_new_deg[end])
                    )

                    JuMP.set_normalized_coefficient(
                        problem_new.deg_new_2[a],
                        soc_quad_new[a+1],
                        -capacity_new_deg[end]/(2*Nfull*state_new_deg[end]^2)
                    )

                    JuMP.set_normalized_coefficient(
                        problem_new.deg_new_2[a],
                        soc_quad_new[a],
                        +capacity_new_deg[end]/(2*Nfull*state_new_deg[end]^2)
                    )

                    JuMP.set_normalized_coefficient(
                        problem_new.deg_new_2[a],
                        soc_new[a],
                        -capacity_new_deg[end]/(Nfull*state_new_deg[end])
                    )

                    JuMP.set_normalized_coefficient(
                        problem_new.deg_new_2[a],
                        soc_new[a+1],
                        +capacity_new_deg[end]/(Nfull*state_new_deg[end])
                    )

                end

            end


        @timeit to "Solve optimization" optimize!(problem_new.M_new)
        
        if termination_status(problem_new.M_new) != MOI.OPTIMAL
            println("NOT OPTIMAL: ", termination_status(problem_new.M_new))
        else
            println("Optimization finished")
        end
        

        @timeit to "Collecting results" begin

            if stage==plan
                for iStep=1:Hours_rolling
                    push!(state_new_deg, JuMP.value(problem_new.soc_new[iStep+1]))
                    push!(charge_new_deg, JuMP.value(problem_new.charge_new[iStep]))
                    push!(discharge_new_deg, JuMP.value(problem_new.discharge_new[iStep]))
                    push!(degradation_new_deg,JuMP.value(problem_new.deg_new[iStep]))     #*k_end
                    push!(capacity_new_deg, JuMP.value(problem_new.capacity_new[iStep+1]))
                    push!(stateQuad_new_deg, JuMP.value(problem_new.soc_quad_new[iStep+1]))
                end
            else
                for iStep=1:Hours_saved
                    push!(state_new_deg, JuMP.value(problem_new.soc_new[iStep+1]))
                    push!(charge_new_deg, JuMP.value(problem_new.charge_new[iStep]))
                    push!(discharge_new_deg, JuMP.value(problem_new.discharge_new[iStep]))
                    push!(degradation_new_deg,JuMP.value(problem_new.deg_new[iStep]))     #*k_end
                    push!(capacity_new_deg, JuMP.value(problem_new.capacity_new[iStep+1]))
                    push!(stateQuad_new_deg, JuMP.value(problem_new.soc_quad_new[iStep+1]))
                end
            end

        if capacity_new_deg[end] <= 0.21*min_SOH
            break
        end

    end
    
end

    println("Collected results Case 3")

    return ResultsNewDeg(
        problem_new_deg,
        charge_new_deg,
        discharge_new_deg,
        state_new_deg,
        stateQuad_new_deg,
        degradation_new_deg,
        capacity_new_deg,
        x_new_deg,
        y_new_deg,
        z_new_deg,
        u_new_deg,
        w_xx_new_deg,
        w_yy_new_deg,
        w_zz_new_deg,
        w_uu_new_deg,
        w_xy_new_deg,
        w_xz_new_deg,
        w_zy_new_deg,
        w_xu_new_deg,
        w_yu_new_deg,
        w_zu_new_deg,
        new_new_deg,
    )

end
