# STAGE MAXIMIZATION PROBLEM FORMULATION

function BuildStageProblem(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam)       

    @unpack (MIPGap, MIPFocus, Method, Cuts, Heuristics) = SolverParameters;
  
    @unpack (NYears, NMonths, NStages, Big, NHoursStep, conv, bin) = InputParameters;     #NSteps,NHoursStage
    @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull, fix) = Battery ;         
    disc = 7

    k = min_SOH/(2*Nfull)
    Small = 0.64

    M = Model(Gurobi.Optimizer)
    set_optimizer_attribute(M, "MIPGap", 0.05)

    # DEFINE VARIABLES

    @variable(M, min_SOC <= soc[iStep=1:NSteps+1] <= max_SOC, base_name = "Energy")                # MWh   energy_Capacity NSteps
    @variable(M, min_SOC^2 <= soc_quad[iStep=1:NSteps+1] <= max_SOC^2, base_name = "Square energy")

    @variable(M, min_P <= charge[iStep=1:NSteps] <= max_P, base_name= "Charge")      #max_disc   0<=discharge<=1
    @variable(M, min_P <= discharge[iStep=1:NSteps] <= max_P, base_name= "Discharge")
    #@variable(M, binary[iStep=1:NSteps], Bin, base_name = "Binary_operation")
    
    @variable(M, 0 <= deg[iStep=1:NSteps] <= Small, base_name = "Degradation")

    @variable(M, 0 <= revamping[iStage=1:NStages] <= (max_SOH-min_SOH), base_name = "Revamping")
    @variable(M, min_SOH <= capacity[iStep=1:NSteps+1] <= max_SOH, base_name = "Energy_Capacity")        #energy_Capacity     [iStage=1:NStages]
    @variable(M, e[iStage=1:NStages], Bin, base_name ="Binary Revamp")

    @variable(M, 0<= rev_vendita[iStage=1:NStages] <= max_SOH, base_name = "Vendita rev")
    @variable(M, -max_SOH<= rev_acquisto[iStage=1:NStages] <= 0, base_name = "Acquisto rev")

    #VARIABLES FOR DISCRETIZATION of Stored Energy

    @variable(M, x[iStep=1:NSteps+1], Bin, base_name = "Binary_1")
    @variable(M, y[iStep=1:NSteps+1], Bin, base_name = "Binary_2")
    @variable(M, z[iStep=1:NSteps+1], Bin, base_name = "Binary_3")
    
    @variable(M, 0<= xy[iStep=1:NSteps+1] <= 1, base_name = "xy")
    @variable(M, 0<= xz[iStep=1:NSteps+1] <= 1, base_name = "xz")
    @variable(M, 0<= yz[iStep=1:NSteps+1] <= 1, base_name = "yz")
    @variable(M, 0 <= xyz[iStep=1:NSteps+1] <=1, base_name = "xyz")
  
    # DEFINE OBJECTIVE function - length(Battery_price) = NStages+1=21

    @objective(
      M,
      MathOptInterface.MAX_SENSE, 
      sum(Power_prices[iStep]*NHoursStep*(discharge[iStep]-charge[iStep]) for iStep=1:NSteps) -               #-sum(Battery_price[iStage]*(revamping[iStage]) for iStage=1:NStages) + 
      Battery_price_purchase[1]*(revamping[1]) -
      sum(Battery_price_purchase[iStage]*(capacity[Steps_stages[iStage]+2] + rev_acquisto[iStage]) for iStage=2:NStages) +
      sum(Battery_price_sale[iStage]*(capacity[Steps_stages[iStage]+1] - rev_vendita[iStage]) for iStage=2:NStages) +
      Battery_price_sale[NStages+1]*(capacity[end]- min_SOH) - 
      sum(fix*e[iStage] for iStage=1:NStages) + 2300    
      )
         
    # DEFINE CONSTRAINTS

    @constraint(M,energy[iStep=1:NSteps], soc[iStep] + (charge[iStep]*Eff_charge-discharge[iStep]/Eff_discharge)*NHoursStep == soc[iStep+1] )

    @constraint(M, en_bal[iStep=1:NSteps+1], min_SOC + ((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]) == soc[iStep])
  
    @constraint(M, en_square[iStep=1:NSteps+1], soc_quad[iStep] == b[1]+c[1]*x[iStep]+c[2]*y[iStep]+c[3]*xy[iStep]+c[4]*z[iStep]+c[5]*xz[iStep]+c[6]*yz[iStep]+c[7]*xyz[iStep])

    # INEQUALITIES CONSTRAINTS

    @constraint(M, xy_1[iStep=1:NSteps+1], xy[iStep] <= x[iStep])
    @constraint(M, xy_2[iStep=1:NSteps+1], xy[iStep] <= y[iStep])
    @constraint(M, xy_3[iStep=1:NSteps+1], xy[iStep] >= x[iStep]+y[iStep]-1)

    @constraint(M, xz_1[iStep=1:NSteps+1], xz[iStep] <= x[iStep])
    @constraint(M, xz_2[iStep=1:NSteps+1], xz[iStep] <= z[iStep])
    @constraint(M, xz_3[iStep=1:NSteps+1], xz[iStep] >= x[iStep]+z[iStep]-1)

    @constraint(M, zy_1[iStep=1:NSteps+1], yz[iStep] <= z[iStep])
    @constraint(M, zy_2[iStep=1:NSteps+1], yz[iStep] <= y[iStep])
    @constraint(M, zy_3[iStep=1:NSteps+1], yz[iStep] >= z[iStep]+y[iStep]-1)   
    
    @constraint(M, xyz_1[iStep=1:NSteps+1], xyz[iStep] <= z[iStep])
    @constraint(M, xyz_2[iStep=1:NSteps+1], xyz[iStep] <= y[iStep])
    @constraint(M, xyz_2[iStep=1:NSteps+1], xyz[iStep] <= x[iStep])
    @constraint(M, xyz_3[iStep=1:NSteps+1], xyz[iStep] >= z[iStep]+y[iStep]+x[iStep]-2)  

    # CONSTRAINTS ON DEGRADATION
    @constraint(M, deg_1[iStep=1:NSteps], deg[iStep] >= soc_quad[iStep]/max_SOC^2 - soc_quad[iStep+1]/max_SOC^2 + (2/max_SOC)*(soc[iStep+1]-soc[iStep]))
    @constraint(M, deg_2[iStep=1:NSteps], deg[iStep] >= soc_quad[iStep+1]/max_SOC^2 - soc_quad[iStep]/max_SOC^2 + (2/max_SOC)*(soc[iStep]-soc[iStep+1]))

    #CONSTRAINT ON REVAMPING
    @constraint(M, energy_capacity[iStage=1:NStages], capacity[Steps_stages[iStage]+2] == capacity[Steps_stages[iStage]+1]+revamping[iStage]-deg[Steps_stages[iStage]+1]*k) #
   
    @constraint(M, initial_e[iStep=1], capacity[iStep] == min_SOH)

    @constraint(M,en_cap1[iStage in 1:NStages, iStep in ((Steps_stages[iStage]+2):Steps_stages[iStage+1])], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)

    @constraint(M, stop_charge[iStage in 2:NStages, iStep in (Steps_stages[iStage]:(Steps_stages[iStage]+Steps_stop[iStage-1]))], charge[iStep] <= (1-e[iStage])*max_P)
    
    @constraint(M, stop_discharge[iStage in 2:NStages, iStep in (Steps_stages[iStage]:(Steps_stages[iStage]+Steps_stop[iStage-1]))], discharge[iStep] <= (1-e[iStage])*max_P) 

    #@constraint(M, rev_3[iStage=1:NStages], capacity[Steps_stages[iStage]+2]>= capacity[Steps_stages[iStage]+1])
    @constraint(M, rev[iStage=1:NStages], revamping[iStage] <= (max_SOH-min_SOH)*e[iStage])

    # CONSTRAINT SU VARIABILI AUSILIARIE

    @constraint(M, vendita[iStage=1], rev_vendita[iStage] == 0)
    @constraint(M, vendita_1[iStage=2:NStages], rev_vendita[iStage] >= 0)
    @constraint(M, vendita_2[iStage=2:NStages], rev_vendita[iStage] >= capacity[Steps_stages[iStage]+1]- e[iStage]*Big)

    @constraint(M, acquisto_1[iStage=2:NStages], rev_acquisto[iStage] >= -(1-e[iStage])*Big)
    @constraint(M, acquisto_2[iStage=2:NStages], rev_acquisto[iStage] >= -capacity[Steps_stages[iStage]+2] )
    @constraint(M, acquisto_3[iStage=1], rev_acquisto[iStage] == 0)

    return BuildStageProblem_3(
        M,
        soc,
        soc_quad,
        charge,
        discharge,
        deg,
        x,
        y,
        z,
        xy,
        xz,
        yz,
        xyz,
        capacity,
        revamping,
        e,
        rev_vendita,
        rev_acquisto,
      )
end

