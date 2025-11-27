# STAGE MAXIMIZATION PROBLEM FORMULATION

function BuildStageProblem(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam)       

    @unpack (MIPGap, MIPFocus, Method, Cuts, Heuristics) = SolverParameters;
  
    @unpack (NYears, NMonths, NStages, Big, NHoursStep, conv, disc) = InputParameters;     #NSteps,NHoursStage
    @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull, fix) = Battery ;         

    k = min_SOH/(2*Nfull)
    Small = 0.9025

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
    #@variable(M, u[iStep=1:NSteps+1], Bin, base_name = "Binary_4")
    
    @variable(M, 0<= w_xx[iStep=1:NSteps+1] <= 1, base_name = "xx")
    @variable(M, 0<= w_yy[iStep=1:NSteps+1] <= 1, base_name = "yy")
    @variable(M, 0<= w_zz[iStep=1:NSteps+1] <= 1, base_name = "zz")
    @variable(M, 0<= w_xy[iStep=1:NSteps+1] <= 1, base_name = "xy")
    @variable(M, 0<= w_xz[iStep=1:NSteps+1] <= 1, base_name = "xz")
    @variable(M, 0<= w_zy[iStep=1:NSteps+1] <= 1, base_name = "yz")

  #=  @variable(M, 0 <= w_uu[iStep=1:NSteps+1] <=1, base_name = "uu")
    @variable(M, 0 <= w_xu[iStep=1:NSteps+1] <=1, base_name = "xu")
    @variable(M, 0 <= w_yu[iStep=1:NSteps+1] <=1, base_name = "yu")
    @variable(M, 0 <= w_zu[iStep=1:NSteps+1] <=1, base_name = "zu")=#

  
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
   # @constraint(M, en_bal[iStep=1:NSteps+1], min_SOC + ((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]+8*u[iStep]) == soc[iStep])    
    
    @constraint(M, en_square[iStep=1:NSteps+1], soc_quad[iStep] == min_SOC^2+ 2*min_SOC*((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep])+(w_xx[iStep]+4*w_xy[iStep]+8*w_xz[iStep]+4*w_yy[iStep]+16*w_zz[iStep]+16*w_zy[iStep])*((max_SOC-min_SOC)/disc)^2)
    #@constraint(M, en_square[iStep=1:NSteps+1], soc_quad[iStep] == min_SOC^2+ 2*min_SOC*((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]+8*u[iStep])+(w_xx[iStep]+4*w_xy[iStep]+8*w_xz[iStep]+16*w_xu[iStep]+4*w_yy[iStep]+16*w_zz[iStep]+16*w_zy[iStep]+32*w_yu[iStep]+64*w_zu[iStep]+64*w_uu[iStep])*((max_SOC-min_SOC)/disc)^2)

    #@constraint(M, bin_charge[iStep=1:NSteps], charge[iStep]<= binary[iStep]*max_P)
    #@constraint(M, bin_discharge[iStep=1:NSteps], discharge[iStep] <= (1-binary[iStep])*max_P)

    # INEQUALITIES CONSTRAINTS
    @constraint(M, xx_1[iStep=1:NSteps+1], w_xx[iStep] <= x[iStep])
    @constraint(M, xx_2[iStep=1:NSteps+1], w_xx[iStep] >= 2*x[iStep]-1)

    @constraint(M, xy_1[iStep=1:NSteps+1], w_xy[iStep] <= x[iStep])
    @constraint(M, xy_2[iStep=1:NSteps+1], w_xy[iStep] <= y[iStep])
    @constraint(M, xy_3[iStep=1:NSteps+1], w_xy[iStep] >= x[iStep]+y[iStep]-1)

    @constraint(M, xz_1[iStep=1:NSteps+1], w_xz[iStep] <= x[iStep])
    @constraint(M, xz_2[iStep=1:NSteps+1], w_xz[iStep] <= z[iStep])
    @constraint(M, xz_3[iStep=1:NSteps+1], w_xz[iStep] >= x[iStep]+z[iStep]-1)

    @constraint(M, yy_1[iStep=1:NSteps+1], w_yy[iStep] <= y[iStep])
    @constraint(M, yy_2[iStep=1:NSteps+1], w_yy[iStep] >= 2*y[iStep]-1)

    @constraint(M, zz_1[iStep=1:NSteps+1], w_zz[iStep] <= z[iStep])
    @constraint(M, zz_2[iStep=1:NSteps+1], w_zz[iStep] >= 2*z[iStep]-1)

    @constraint(M, zy_1[iStep=1:NSteps+1], w_zy[iStep] <= z[iStep])
    @constraint(M, zy_2[iStep=1:NSteps+1], w_zy[iStep] <= y[iStep])
    @constraint(M, zy_3[iStep=1:NSteps+1], w_zy[iStep] >= z[iStep]+y[iStep]-1)

  #=  @constraint(M, uu_1[iStep=1:NSteps+1], w_uu[iStep] <= u[iStep])
    @constraint(M, uu_2[iStep=1:NSteps+1], w_uu[iStep] >= 2*u[iStep]-1)

    @constraint(M, xu_1[iStep=1:NSteps+1], w_xu[iStep] <= x[iStep])
    @constraint(M, xu_2[iStep=1:NSteps+1], w_xu[iStep] <= u[iStep])
    @constraint(M, xu_3[iStep=1:NSteps+1], w_xu[iStep] >= x[iStep]+u[iStep]-1)

    @constraint(M, yu_1[iStep=1:NSteps+1], w_yu[iStep] <= y[iStep])
    @constraint(M, yu_2[iStep=1:NSteps+1], w_yu[iStep] <= u[iStep])
    @constraint(M, yu_3[iStep=1:NSteps+1], w_yu[iStep] >= y[iStep]+u[iStep]-1)

    @constraint(M, zu_1[iStep=1:NSteps+1], w_zu[iStep] <= z[iStep])
    @constraint(M, zu_2[iStep=1:NSteps+1], w_zu[iStep] <= u[iStep])
    @constraint(M, zu_3[iStep=1:NSteps+1], w_zu[iStep] >= z[iStep]+u[iStep]-1) =#
    

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

    return BuildStageProblem(
        M,
        soc,
        soc_quad,
        charge,
        discharge,
        #binary,
        deg,
        x,
        y,
        z,
        #u,
        w_xx,
        w_yy,
        w_zz,
        w_xy,
        w_xz,
        w_zy,
        #w_uu,
        #w_xu,
        #w_yu,
        #w_zu,
        capacity,
        revamping,
        e,
        rev_vendita,
        rev_acquisto,
      )
end


#=@constraint(M,en_cap1[iStep=(Steps_stages[1]+2):Steps_stages[2]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)
    @constraint(M,en_cap2[iStep=(Steps_stages[2]+2):Steps_stages[3]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)
    @constraint(M,en_cap3[iStep=(Steps_stages[3]+2):Steps_stages[4]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)
    @constraint(M,en_cap4[iStep=(Steps_stages[4]+2):Steps_stages[5]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)
    @constraint(M,en_cap5[iStep=(Steps_stages[5]+2):Steps_stages[6]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)
    @constraint(M,en_cap6[iStep=(Steps_stages[6]+2):Steps_stages[7]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)
    @constraint(M,en_cap7[iStep=(Steps_stages[7]+2):Steps_stages[8]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)
    @constraint(M,en_cap8[iStep=(Steps_stages[8]+2):Steps_stages[9]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)
    @constraint(M,en_cap9[iStep=(Steps_stages[9]+2):Steps_stages[10]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)
    @constraint(M,en_cap10[iStep=(Steps_stages[10]+2):Steps_stages[11]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)
    @constraint(M,en_cap11[iStep=(Steps_stages[11]+2):Steps_stages[12]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)
    @constraint(M,en_cap12[iStep=(Steps_stages[12]+2):Steps_stages[13]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)
    @constraint(M,en_cap13[iStep=(Steps_stages[13]+2):Steps_stages[14]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)
    @constraint(M,en_cap14[iStep=(Steps_stages[14]+2):Steps_stages[15]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)
    @constraint(M,en_cap15[iStep=(Steps_stages[15]+2):Steps_stages[16]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)
    @constraint(M,en_cap16[iStep=(Steps_stages[16]+2):Steps_stages[17]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)
    @constraint(M,en_cap17[iStep=(Steps_stages[17]+2):Steps_stages[18]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)
    @constraint(M,en_cap18[iStep=(Steps_stages[18]+2):Steps_stages[19]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)
    @constraint(M,en_cap19[iStep=(Steps_stages[19]+2):Steps_stages[20]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)
    @constraint(M,en_cap20[iStep=(Steps_stages[20]+2):Steps_stages[21]], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)=#

    #=CONSTRAINTS CHARGE - downtime
    @constraint(M, stop_charge_1[iStep=Steps_stages[2]:(Steps_stages[2]+Steps_stop[1])], charge[iStep]<= (1-e[2])*max_P)
    @constraint(M, stop_charge_2[iStep=Steps_stages[3]:(Steps_stages[3]+Steps_stop[2])], charge[iStep]<= (1-e[3])*max_P)
    @constraint(M, stop_charge_3[iStep=Steps_stages[4]:(Steps_stages[4]+Steps_stop[3])], charge[iStep]<= (1-e[4])*max_P)
    @constraint(M, stop_charge_4[iStep=Steps_stages[5]:(Steps_stages[5]+Steps_stop[4])], charge[iStep]<= (1-e[5])*max_P)
    @constraint(M, stop_charge_5[iStep=Steps_stages[6]:(Steps_stages[6]+Steps_stop[5])], charge[iStep]<= (1-e[6])*max_P)
    @constraint(M, stop_charge_6[iStep=Steps_stages[7]:(Steps_stages[7]+Steps_stop[6])], charge[iStep]<= (1-e[7])*max_P)
    @constraint(M, stop_charge_7[iStep=Steps_stages[8]:(Steps_stages[8]+Steps_stop[7])], charge[iStep]<= (1-e[8])*max_P)
    @constraint(M, stop_charge_8[iStep=Steps_stages[9]:(Steps_stages[9]+Steps_stop[8])], charge[iStep]<= (1-e[9])*max_P)
    @constraint(M, stop_charge_9[iStep=Steps_stages[10]:(Steps_stages[10]+Steps_stop[9])], charge[iStep]<= (1-e[10])*max_P)
    @constraint(M, stop_charge_10[iStep=Steps_stages[11]:(Steps_stages[11]+Steps_stop[10])], charge[iStep]<= (1-e[11])*max_P)
    @constraint(M, stop_charge_11[iStep=Steps_stages[12]:(Steps_stages[12]+Steps_stop[11])], charge[iStep]<= (1-e[12])*max_P)
    @constraint(M, stop_charge_12[iStep=Steps_stages[13]:(Steps_stages[13]+Steps_stop[12])], charge[iStep]<= (1-e[13])*max_P)
    @constraint(M, stop_charge_13[iStep=Steps_stages[14]:(Steps_stages[14]+Steps_stop[13])], charge[iStep]<= (1-e[14])*max_P)
    @constraint(M, stop_charge_14[iStep=Steps_stages[15]:(Steps_stages[15]+Steps_stop[14])], charge[iStep]<= (1-e[15])*max_P)
    @constraint(M, stop_charge_15[iStep=Steps_stages[16]:(Steps_stages[16]+Steps_stop[15])], charge[iStep]<= (1-e[16])*max_P)
    @constraint(M, stop_charge_16[iStep=Steps_stages[17]:(Steps_stages[17]+Steps_stop[16])], charge[iStep]<= (1-e[17])*max_P)
    @constraint(M, stop_charge_17[iStep=Steps_stages[18]:(Steps_stages[18]+Steps_stop[17])], charge[iStep]<= (1-e[18])*max_P)
    @constraint(M, stop_charge_18[iStep=Steps_stages[19]:(Steps_stages[19]+Steps_stop[18])], charge[iStep]<= (1-e[19])*max_P)
    @constraint(M, stop_charge_19[iStep=Steps_stages[20]:(Steps_stages[20]+Steps_stop[19])], charge[iStep]<= (1-e[20])*max_P)=#
 
    #=CONSTRAINTS DISCHARGE - downtime
    @constraint(M, stop_discharge_1[iStep=Steps_stages[2]:(Steps_stages[2]+Steps_stop[1])], discharge[iStep]<= (1-e[2])*max_P)
    @constraint(M, stop_discharge_2[iStep=Steps_stages[3]:(Steps_stages[3]+Steps_stop[2])], discharge[iStep]<= (1-e[3])*max_P)
    @constraint(M, stop_discharge_3[iStep=Steps_stages[4]:(Steps_stages[4]+Steps_stop[3])], discharge[iStep]<= (1-e[4])*max_P)
    @constraint(M, stop_discharge_4[iStep=Steps_stages[5]:(Steps_stages[5]+Steps_stop[4])], discharge[iStep]<= (1-e[5])*max_P)
    @constraint(M, stop_discharge_5[iStep=Steps_stages[6]:(Steps_stages[6]+Steps_stop[5])], discharge[iStep]<= (1-e[6])*max_P)
    @constraint(M, stop_discharge_6[iStep=Steps_stages[7]:(Steps_stages[7]+Steps_stop[6])], discharge[iStep]<= (1-e[7])*max_P)
    @constraint(M, stop_discharge_7[iStep=Steps_stages[8]:(Steps_stages[8]+Steps_stop[7])], discharge[iStep]<= (1-e[8])*max_P)
    @constraint(M, stop_discharge_8[iStep=Steps_stages[9]:(Steps_stages[9]+Steps_stop[8])], discharge[iStep]<= (1-e[9])*max_P)
    @constraint(M, stop_discharge_9[iStep=Steps_stages[10]:(Steps_stages[10]+Steps_stop[9])], discharge[iStep]<= (1-e[10])*max_P)
    @constraint(M, stop_discharge_10[iStep=Steps_stages[11]:(Steps_stages[11]+Steps_stop[10])], discharge[iStep]<= (1-e[11])*max_P)
    @constraint(M, stop_discharge_11[iStep=Steps_stages[12]:(Steps_stages[12]+Steps_stop[11])], discharge[iStep]<= (1-e[12])*max_P)
    @constraint(M, stop_discharge_12[iStep=Steps_stages[13]:(Steps_stages[13]+Steps_stop[12])], discharge[iStep]<= (1-e[13])*max_P)
    @constraint(M, stop_discharge_13[iStep=Steps_stages[14]:(Steps_stages[14]+Steps_stop[13])], discharge[iStep]<= (1-e[14])*max_P)
    @constraint(M, stop_discharge_14[iStep=Steps_stages[15]:(Steps_stages[15]+Steps_stop[14])], discharge[iStep]<= (1-e[15])*max_P)
    @constraint(M, stop_discharge_15[iStep=Steps_stages[16]:(Steps_stages[16]+Steps_stop[15])], discharge[iStep]<= (1-e[16])*max_P)
    @constraint(M, stop_discharge_16[iStep=Steps_stages[17]:(Steps_stages[17]+Steps_stop[16])], discharge[iStep]<= (1-e[17])*max_P)
    @constraint(M, stop_discharge_17[iStep=Steps_stages[18]:(Steps_stages[18]+Steps_stop[17])], discharge[iStep]<= (1-e[18])*max_P)
    @constraint(M, stop_discharge_18[iStep=Steps_stages[19]:(Steps_stages[19]+Steps_stop[18])], discharge[iStep]<= (1-e[19])*max_P)
    @constraint(M, stop_discharge_19[iStep=Steps_stages[20]:(Steps_stages[20]+Steps_stop[19])], discharge[iStep]<= (1-e[20])*max_P)=#



