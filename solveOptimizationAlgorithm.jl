# SOLVE OPTIMIZATION PROBLEM

function solveOptimizationProblem_3(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam)

    @unpack (NYears, NMonths, NStages, Big, NHoursStep, bin) = InputParameters;                #NSteps, NHoursStage
    @unpack (min_SOC, max_SOC, Eff_charge, Eff_discharge, min_P, max_P, max_SOH, min_SOH, Nfull,fix ) = Battery;

    println("Solving Optimization Problem")

    k = min_SOH/(2*Nfull)

    objective = 0                   
    #revenues_per_stage = zeros(NStages)
    gain_stage = zeros(NStages)
    cost_rev = zeros(NStages)
    deg_stage = zeros(NStages)

    charge = zeros(NSteps)
    discharge = zeros(NSteps)
    soc = zeros(NSteps+1)
    deg = zeros(NSteps)
    rev= zeros(NStages)
    cap = zeros(NSteps+1)
    e = zeros(NStages)
    #bin = zeros(NSteps)

    #rev_vendita = zeros(NStages-1)
    rev_vendita = zeros(NStages)
    rev_acquisto = zeros(NStages)

    soc_quad = zeros(NSteps+1)
    x = zeros(NSteps+1)
    y = zeros(NSteps+1)
    z = zeros(NSteps+1)
    xy = zeros(NSteps+1)
    xz = zeros(NSteps+1)
    yz = zeros(NSteps+1)
    xyz = zeros(NSteps+1)

    problem = BuildStageProblem_3(InputParameters, SolverParameters, Battery)

    @timeit to "Solve optimization" optimize!(problem.M)

    if termination_status(problem.M) != MOI.OPTIMAL
        println("NOT OPTIMAL: ", termination_status(problem.M))
    else
        println("Optimization finished")
    end

    @timeit to "Collecting results" begin
        objective = JuMP.objective_value(problem.M)
        
        for iStep=1:NSteps
            soc[iStep] = JuMP.value(problem.soc[iStep])
            charge[iStep] = JuMP.value(problem.charge[iStep])
            discharge[iStep] = JuMP.value(problem.discharge[iStep])
            deg[iStep] = JuMP.value(problem.deg[iStep])*k
            cap[iStep] = JuMP.value(problem.capacity[iStep])
            #bin[iStep] = JuMP.value(problem.binary[iStep])

            soc_quad[iStep] = JuMP.value(problem.soc_quad[iStep])
            x[iStep] = JuMP.value(problem.x[iStep])
            y[iStep] = JuMP.value(problem.y[iStep])
            z[iStep] = JuMP.value(problem.z[iStep])
            xy[iStep] = JuMP.value(problem.xy[iStep])
            xz[iStep] = JuMP.value(problem.xz[iStep])
            yz[iStep] = JuMP.value(problem.yz[iStep])
            xyz[iStep] = JuMP.value(problem.xyz[iStep])

        end

        soc[end] = JuMP.value(problem.soc[end])
        soc_quad[end] = JuMP.value(problem.soc_quad[end])
        x[end] = JuMP.value(problem.x[end])
        y[end] = JuMP.value(problem.y[end])
        z[end] = JuMP.value(problem.z[end])
        xy[end] = JuMP.value(problem.xy[end])
        xz[end] = JuMP.value(problem.xz[end])
        yz[end] = JuMP.value(problem.yz[end])
        xyz[end] = JuMP.value(problem.xyz[end])

        cap[end] = JuMP.value(problem.capacity[end])
     
        for iStage=1:NStages
            rev[iStage] = JuMP.value(problem.revamping[iStage])
            e[iStage] = JuMP.value(problem.e[iStage])
            deg_stage[iStage] = sum(deg[iStep] for iStep=(Steps_stages[iStage]+1):(Steps_stages[iStage+1]))
            rev_acquisto[iStage] = JuMP.value(problem.rev_acquisto[iStage])
            rev_vendita[iStage] = JuMP.value(problem.rev_vendita[iStage])
        end
        
        for iStage=2:(NStages-1)
            #revenues_per_stage[iStage] = sum(Power_prices[iStep]*NHoursStep*(discharge[iStep]-charge[iStep]) for iStep=(Steps_stages[iStage]+1):(Steps_stages[iStage+1])) + Battery_price_sale[iStage]*(cap[Steps_stages[iStage-1]+1]+rev[iStage-1]-rev_vendita[iStage-1]) - Battery_price_purchase[iStage]*(cap[Steps_stages[iStage]+1]+rev[iStage]+rev_acquisto[iStage]) - e[iStage]*fix
            #revenues_per_stage[iStage] = sum(Power_prices[iStep]*NHoursStep*(discharge[iStep]-charge[iStep]) for iStep=(Steps_stages[iStage]+1):(Steps_stages[iStage+1])) - Battery_price_purchase[iStage]*rev[iStage] - e[iStage]*fix
            gain_stage[iStage] = sum(Power_prices[iStep]*NHoursStep*(discharge[iStep]-charge[iStep]) for iStep=(Steps_stages[iStage]+1):(Steps_stages[iStage+1]))
            #cost_rev[iStage] = Battery_price_purchase[iStage]*(cap[Steps_stages[iStage]+2]) - Battery_price_sale[iStage]*(cap[Steps_stages[iStage]+1]-rev_vendita[iStage]) + e[iStage]*fix
            cost_rev[iStage] = Battery_price_purchase[iStage]*(cap[Steps_stages[iStage]+2]+rev_acquisto[iStage]) - Battery_price_sale[iStage]*(cap[Steps_stages[iStage]+1]-rev_vendita[iStage]) + e[iStage]*fix
            
           
        end

        #revenues_per_stage[1] = sum(Power_prices[iStep]*NHoursStep*(discharge[iStep]-charge[iStep]) for iStep=(Steps_stages[1]+1):(Steps_stages[2])) - Battery_price_purchase[1]*rev[1] - fix*e[1]
        gain_stage[1] = sum(Power_prices[iStep]*NHoursStep*(discharge[iStep]-charge[iStep]) for iStep=(Steps_stages[1]+1):(Steps_stages[2]))
        cost_rev[1] = Battery_price_purchase[1]*rev[1] + fix*e[1]

        #revenues_per_stage[NStages] = sum(Power_prices[iStep]*NHoursStep*(discharge[iStep]-charge[iStep]) for iStep=(Steps_stages[NStages]+1):(Steps_stages[NStages+1])) + Battery_price_sale[NStages+1]*(cap[end]-min_SOH)- fix*e[NStages] - Battery_price_purchase[NStages]*rev[NStages]
        gain_stage[NStages]= sum(Power_prices[iStep]*NHoursStep*(discharge[iStep]-charge[iStep]) for iStep=(Steps_stages[NStages]+1):(Steps_stages[NStages+1]))
        cost_rev[NStages] = e[NStages]*fix + Battery_price_purchase[NStages]*(cap[Steps_stages[NStages]+2]+rev_acquisto[NStages]) - Battery_price_sale[NStages+1]*(cap[end]-min_SOH) -Battery_price_sale[NStages]*(cap[Steps_stages[NStages]+1]-rev_vendita[NStages])

        #= revenues_per_stage[NStages] = sum(Power_prices[iStep]*NHoursStep*(discharge[iStep]-charge[iStep]) for iStep=(Steps_stages[NStages]+1):(Steps_stages[NStages+1])) - Battery_price_purchase[NStages]*capacity[Steps_stages[iStage]+2] + Battery_price_sale[NStages]*capacity[Steps_stages[NStages]+1] - e[NStages]*fix
        cost_rev[iStage] = Battery_price_purchase[iStage]*capacity[Steps_stages[iStage]+2] - Battery_price_sale[iStage]*capacity[Steps_stages[iStage]+1] + e[iStage]*fix
        =#
        

    end
    
    println("Collected results")

    return Results_3(
        objective,
        #revenues_per_stage,
        gain_stage,
        cost_rev,
        deg_stage,
        soc,
        charge,
        discharge,
        #bin,
        deg,
        soc_quad,
        x,
        y,
        z,
        xy,
        xz,
        yz,
        xyz,
        rev,
        cap,  
        e,
        rev_vendita,
        rev_acquisto,
    )

end


function solveOptimizationProblem_4(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam)

    @unpack (NYears, NMonths, NStages, Big, NHoursStep, bin) = InputParameters;                #NSteps, NHoursStage
    @unpack (min_SOC, max_SOC, Eff_charge, Eff_discharge, min_P, max_P, max_SOH, min_SOH, Nfull,fix ) = Battery;

    println("Solving Optimization Problem")

    k = min_SOH/(2*Nfull)

    objective = 0                   
    #revenues_per_stage = zeros(NStages)
    gain_stage = zeros(NStages)
    cost_rev = zeros(NStages)
    deg_stage = zeros(NStages)

    charge = zeros(NSteps)
    discharge = zeros(NSteps)
    soc = zeros(NSteps+1)
    deg = zeros(NSteps)
    rev= zeros(NStages)
    cap = zeros(NSteps+1)
    e = zeros(NStages)
    #bin = zeros(NSteps)

    #rev_vendita = zeros(NStages-1)
    rev_vendita = zeros(NStages)
    rev_acquisto = zeros(NStages)

    soc_quad = zeros(NSteps+1)
    x = zeros(NSteps+1)
    y = zeros(NSteps+1)
    z = zeros(NSteps+1)
    u = zeros(NSteps+1)
    xy = zeros(NSteps+1)
    xz = zeros(NSteps+1)
    xu = zeros(NSteps+1)
    yz = zeros(NSteps+1)
    yu = zeros(NSteps+1)
    zu = zeros(NSteps+1)
    xyzu = zeros(NSteps+1)

    problem = BuildStageProblem_4(InputParameters, SolverParameters, Battery)

    @timeit to "Solve optimization" optimize!(problem.M)

    if termination_status(problem.M) != MOI.OPTIMAL
        println("NOT OPTIMAL: ", termination_status(problem.M))
    else
        println("Optimization finished")
    end

    @timeit to "Collecting results" begin
        objective = JuMP.objective_value(problem.M)
        
        for iStep=1:NSteps
            soc[iStep] = JuMP.value(problem.soc[iStep])
            charge[iStep] = JuMP.value(problem.charge[iStep])
            discharge[iStep] = JuMP.value(problem.discharge[iStep])
            deg[iStep] = JuMP.value(problem.deg[iStep])*k
            cap[iStep] = JuMP.value(problem.capacity[iStep])
            #bin[iStep] = JuMP.value(problem.binary[iStep])

            soc_quad[iStep] = JuMP.value(problem.soc_quad[iStep])
            x[iStep] = JuMP.value(problem.x[iStep])
            y[iStep] = JuMP.value(problem.y[iStep])
            z[iStep] = JuMP.value(problem.z[iStep])
            u[iStep] = JuMP.value(problem.u[iStep])
            xy[iStep] = JuMP.value(problem.xy[iStep])
            xz[iStep] = JuMP.value(problem.xz[iStep])
            xu[iStep] = JuMP.value(problem.xu[iStep])
            yz[iStep] = JuMP.value(problem.yz[iStep])
            yu[iStep] = JuMP.value(problem.yu[iStep])
            zu[iStep] = JuMP.value(problem.zu[iStep])
            xyzu[iStep] = JuMP.value(problem.xyzu[iStep])

        end

        soc[end] = JuMP.value(problem.soc[end])
        soc_quad[end] = JuMP.value(problem.soc_quad[end])
        x[end] = JuMP.value(problem.x[end])
        y[end] = JuMP.value(problem.y[end])
        z[end] = JuMP.value(problem.z[end])
        u[end] = JuMP.value(problem.u[end])
        xy[end] = JuMP.value(problem.xy[end])
        xz[end] = JuMP.value(problem.xz[end])
        xu[end] = JuMP.value(problem.xu[end])
        yz[end] = JuMP.value(problem.yz[end])
        yu[end] = JuMP.value(problem.yu[end])
        zu[end] = JuMP.value(problem.zu[end])
        xyzu[end] = JuMP.value(problem.xyzu[end])

        cap[end] = JuMP.value(problem.capacity[end])
     
        for iStage=1:NStages
            rev[iStage] = JuMP.value(problem.revamping[iStage])
            e[iStage] = JuMP.value(problem.e[iStage])
            deg_stage[iStage] = sum(deg[iStep] for iStep=(Steps_stages[iStage]+1):(Steps_stages[iStage+1]))
            rev_acquisto[iStage] = JuMP.value(problem.rev_acquisto[iStage])
            rev_vendita[iStage] = JuMP.value(problem.rev_vendita[iStage])
        end
        
        for iStage=2:(NStages-1)
            #revenues_per_stage[iStage] = sum(Power_prices[iStep]*NHoursStep*(discharge[iStep]-charge[iStep]) for iStep=(Steps_stages[iStage]+1):(Steps_stages[iStage+1])) + Battery_price_sale[iStage]*(cap[Steps_stages[iStage-1]+1]+rev[iStage-1]-rev_vendita[iStage-1]) - Battery_price_purchase[iStage]*(cap[Steps_stages[iStage]+1]+rev[iStage]+rev_acquisto[iStage]) - e[iStage]*fix
            #revenues_per_stage[iStage] = sum(Power_prices[iStep]*NHoursStep*(discharge[iStep]-charge[iStep]) for iStep=(Steps_stages[iStage]+1):(Steps_stages[iStage+1])) - Battery_price_purchase[iStage]*rev[iStage] - e[iStage]*fix
            gain_stage[iStage] = sum(Power_prices[iStep]*NHoursStep*(discharge[iStep]-charge[iStep]) for iStep=(Steps_stages[iStage]+1):(Steps_stages[iStage+1]))
            #cost_rev[iStage] = Battery_price_purchase[iStage]*(cap[Steps_stages[iStage]+2]) - Battery_price_sale[iStage]*(cap[Steps_stages[iStage]+1]-rev_vendita[iStage]) + e[iStage]*fix
            cost_rev[iStage] = Battery_price_purchase[iStage]*(cap[Steps_stages[iStage]+2]+rev_acquisto[iStage]) - Battery_price_sale[iStage]*(cap[Steps_stages[iStage]+1]-rev_vendita[iStage]) + e[iStage]*fix
        end

        #revenues_per_stage[1] = sum(Power_prices[iStep]*NHoursStep*(discharge[iStep]-charge[iStep]) for iStep=(Steps_stages[1]+1):(Steps_stages[2])) - Battery_price_purchase[1]*rev[1] - fix*e[1]
        gain_stage[1] = sum(Power_prices[iStep]*NHoursStep*(discharge[iStep]-charge[iStep]) for iStep=(Steps_stages[1]+1):(Steps_stages[2]))
        cost_rev[1] = Battery_price_purchase[1]*rev[1] + fix*e[1]

        #revenues_per_stage[NStages] = sum(Power_prices[iStep]*NHoursStep*(discharge[iStep]-charge[iStep]) for iStep=(Steps_stages[NStages]+1):(Steps_stages[NStages+1])) + Battery_price_sale[NStages+1]*(cap[end]-min_SOH)- fix*e[NStages] - Battery_price_purchase[NStages]*rev[NStages]
        gain_stage[NStages]= sum(Power_prices[iStep]*NHoursStep*(discharge[iStep]-charge[iStep]) for iStep=(Steps_stages[NStages]+1):(Steps_stages[NStages+1]))
        cost_rev[NStages] = e[NStages]*fix + Battery_price_purchase[NStages]*(cap[Steps_stages[NStages]+2]+rev_acquisto[NStages]) - Battery_price_sale[NStages+1]*(cap[end]-min_SOH) -Battery_price_sale[NStages]*(cap[Steps_stages[NStages]+1]-rev_vendita[NStages])

        #= revenues_per_stage[NStages] = sum(Power_prices[iStep]*NHoursStep*(discharge[iStep]-charge[iStep]) for iStep=(Steps_stages[NStages]+1):(Steps_stages[NStages+1])) - Battery_price_purchase[NStages]*capacity[Steps_stages[iStage]+2] + Battery_price_sale[NStages]*capacity[Steps_stages[NStages]+1] - e[NStages]*fix
        cost_rev[iStage] = Battery_price_purchase[iStage]*capacity[Steps_stages[iStage]+2] - Battery_price_sale[iStage]*capacity[Steps_stages[iStage]+1] + e[iStage]*fix
        =#
        

    end
    
    println("Collected results")

    return Results_4(
        objective,
        #revenues_per_stage,
        gain_stage,
        cost_rev,
        deg_stage,
        soc,
        charge,
        discharge,
        #bin,
        deg,
        soc_quad,
        x,
        y,
        z,
        u,
        xy,
        xz,
        xu,
        yz,
        yu,
        zu,
        xyzu,
        rev,
        cap,  
        e,
        rev_vendita,
        rev_acquisto,
    )

end