# STAGE MAXIMIZATION PROBLEM FORMULATION

# CASE 2

function BuildRolling(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam)       #, state_variables::states When we have 2 hydropower plants- 2 turbines

    @unpack (MIPGap, MIPFocus, Method, Cuts, Heuristics) = SolverParameters;
  
    @unpack (NYears, NMonths, NStages, Big, NHoursStep, conv, disc, Hours_rolling, Hours_saved) = InputParameters;     #NSteps,NHoursStage
    @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull, fix) = Battery ;         

    k_f = min_SOH/(2*Nfull)
    Small_f = 0.64

    M_f = Model(Gurobi.Optimizer)
    set_optimizer_attribute(M_f, "MIPGap", 0.1)

    # DEFINE VARIABLES

    @variable(M_f, min_SOC <= soc_f[iStep=1:Hours_rolling+1] <= max_SOC, base_name = "Energy_f")                # MWh   energy_Capacity NSteps
    @variable(M_f, min_SOC^2 <= soc_quad_f[iStep=1:Hours_rolling+1] <= max_SOC^2, base_name = "Square energy_f")

    @variable(M_f, min_P <= charge_f[iStep=1:Hours_rolling] <= max_P, base_name= "Charge_f")      #max_disc   0<=discharge<=1
    @variable(M_f, min_P <= discharge_f[iStep=1:Hours_rolling] <= max_P, base_name= "Discharge_f")
    
    @variable(M_f, 0 <= deg_f[iStep=1:Hours_rolling] <= Small_f, base_name = "Degradation_f")

    @variable(M_f, 0 <= capacity_f[iStep=1:Hours_rolling+1] <= max_SOH, base_name = "Energy_Capacity_f")        #energy_Capacity     [iStage=1:NStages]

    #VARIABLES FOR DISCRETIZATION of Stored Energy

    @variable(M_f, x_f[iStep=1:Hours_rolling+1], Bin, base_name = "Binary_1_f")
    @variable(M_f, y_f[iStep=1:Hours_rolling+1], Bin, base_name = "Binary_2_f")
    @variable(M_f, z_f[iStep=1:Hours_rolling+1], Bin, base_name = "Binary_3_f")
    @variable(M_f, u_f[iStep=1:Hours_rolling+1], Bin, base_name = "Binary_4_f")
    
    @variable(M_f, 0<= w_xx_f[iStep=1:Hours_rolling+1] <= 1, base_name = "xx_f")
    @variable(M_f, 0<= w_yy_f[iStep=1:Hours_rolling+1] <= 1, base_name = "yy_f")
    @variable(M_f, 0<= w_zz_f[iStep=1:Hours_rolling+1] <= 1, base_name = "zz_f")
    @variable(M_f, 0<= w_xy_f[iStep=1:Hours_rolling+1] <= 1, base_name = "xy_f")
    @variable(M_f, 0<= w_xz_f[iStep=1:Hours_rolling+1] <= 1, base_name = "xz_f")
    @variable(M_f, 0<= w_zy_f[iStep=1:Hours_rolling+1] <= 1, base_name = "yz_f")

    @variable(M_f, 0 <= w_uu_f[iStep=1:Hours_rolling+1] <=1, base_name = "uu_f")
    @variable(M_f, 0 <= w_xu_f[iStep=1:Hours_rolling+1] <=1, base_name = "xu_f")
    @variable(M_f, 0 <= w_yu_f[iStep=1:Hours_rolling+1] <=1, base_name = "yu_f")
    @variable(M_f, 0 <= w_zu_f[iStep=1:Hours_rolling+1] <=1, base_name = "zu_f")

  
    # DEFINE OBJECTIVE function 

    @objective(
      M_f,
      MathOptInterface.MAX_SENSE, 
      sum(1*NHoursStep*discharge_f[iStep]-1*NHoursStep*charge_f[iStep] for iStep=1:Hours_rolling)   #-deg_f[iStep]*k_f*100
      )
         
    # DEFINE CONSTRAINTS

    @constraint(M_f, State_initial[iStep=1], soc_f[iStep] == max_SOH)
    @constraint(M_f, energy_f[iStep=1:Hours_rolling], soc_f[iStep+1] - (charge_f[iStep]*Eff_charge-discharge_f[iStep]/Eff_discharge)*NHoursStep == soc_f[iStep] )
    @constraint(M_f, cons[iStep=1:Hours_rolling], soc_f[iStep] <= capacity_f[iStep])

    #@constraint(M, en_bal[iStep=1:NSteps+1], min_SOC + ((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]) == soc[iStep])
    @constraint(M_f, en_bal_f[iStep=1:Hours_rolling+1], min_SOC + ((max_SOC-min_SOC)/disc)*(x_f[iStep]+2*y_f[iStep]+4*z_f[iStep]+8*u_f[iStep]) == soc_f[iStep])    
    
    #@constraint(M, en_square[iStep=1:NSteps+1], soc_quad[iStep] == min_SOC^2+ 2*min_SOC*((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep])+(w_xx[iStep]+4*w_xy[iStep]+8*w_xz[iStep]+4*w_yy[iStep]+16*w_zz[iStep]+16*w_zy[iStep])*((max_SOC-min_SOC)/disc)^2)
    @constraint(M_f, en_square[iStep=1:Hours_rolling+1], soc_quad_f[iStep] == min_SOC^2+ 2*min_SOC*((max_SOC-min_SOC)/disc)*(x_f[iStep]+2*y_f[iStep]+4*z_f[iStep]+8*u_f[iStep])+(w_xx_f[iStep]+4*w_xy_f[iStep]+8*w_xz_f[iStep]+16*w_xu_f[iStep]+4*w_yy_f[iStep]+16*w_zz_f[iStep]+16*w_zy_f[iStep]+32*w_yu_f[iStep]+64*w_zu_f[iStep]+64*w_uu_f[iStep])*((max_SOC-min_SOC)/disc)^2)

    # INEQUALITIES CONSTRAINTS
    @constraint(M_f, xx_1_f[iStep=1:Hours_rolling+1], w_xx_f[iStep] <= x_f[iStep])
    @constraint(M_f, xx_2_f[iStep=1:Hours_rolling+1], w_xx_f[iStep] >= 2*x_f[iStep]-1)

    @constraint(M_f, xy_1_f[iStep=1:Hours_rolling+1], w_xy_f[iStep] <= x_f[iStep])
    @constraint(M_f, xy_2_f[iStep=1:Hours_rolling+1], w_xy_f[iStep] <= y_f[iStep])
    @constraint(M_f, xy_3_f[iStep=1:Hours_rolling+1], w_xy_f[iStep] >= x_f[iStep]+y_f[iStep]-1)

    @constraint(M_f, xz_1_f[iStep=1:Hours_rolling+1], w_xz_f[iStep] <= x_f[iStep])
    @constraint(M_f, xz_2_f[iStep=1:Hours_rolling+1], w_xz_f[iStep] <= z_f[iStep])
    @constraint(M_f, xz_3_f[iStep=1:Hours_rolling+1], w_xz_f[iStep] >= x_f[iStep]+z_f[iStep]-1)

    @constraint(M_f, yy_1_f[iStep=1:Hours_rolling+1], w_yy_f[iStep] <= y_f[iStep])
    @constraint(M_f, yy_2_f[iStep=1:Hours_rolling+1], w_yy_f[iStep] >= 2*y_f[iStep]-1)

    @constraint(M_f, zz_1_f[iStep=1:Hours_rolling+1], w_zz_f[iStep] <= z_f[iStep])
    @constraint(M_f, zz_2_f[iStep=1:Hours_rolling+1], w_zz_f[iStep] >= 2*z_f[iStep]-1)

    @constraint(M_f, zy_1_f[iStep=1:Hours_rolling+1], w_zy_f[iStep] <= z_f[iStep])
    @constraint(M_f, zy_2_f[iStep=1:Hours_rolling+1], w_zy_f[iStep] <= y_f[iStep])
    @constraint(M_f, zy_3_f[iStep=1:Hours_rolling+1], w_zy_f[iStep] >= z_f[iStep]+y_f[iStep]-1)

    @constraint(M_f, uu_1_f[iStep=1:Hours_rolling+1], w_uu_f[iStep] <= u_f[iStep])
    @constraint(M_f, uu_2_f[iStep=1:Hours_rolling+1], w_uu_f[iStep] >= 2*u_f[iStep]-1)

    @constraint(M_f, xu_1_f[iStep=1:Hours_rolling+1], w_xu_f[iStep] <= x_f[iStep])
    @constraint(M_f, xu_2_f[iStep=1:Hours_rolling+1], w_xu_f[iStep] <= u_f[iStep])
    @constraint(M_f, xu_3_f[iStep=1:Hours_rolling+1], w_xu_f[iStep] >= x_f[iStep]+u_f[iStep]-1)

    @constraint(M_f, yu_1_f[iStep=1:Hours_rolling+1], w_yu_f[iStep] <= y_f[iStep])
    @constraint(M_f, yu_2_f[iStep=1:Hours_rolling+1], w_yu_f[iStep] <= u_f[iStep])
    @constraint(M_f, yu_3_f[iStep=1:Hours_rolling+1], w_yu_f[iStep] >= y_f[iStep]+u_f[iStep]-1)

    @constraint(M_f, zu_1_f[iStep=1:Hours_rolling+1], w_zu_f[iStep] <= z_f[iStep])
    @constraint(M_f, zu_2_f[iStep=1:Hours_rolling+1], w_zu_f[iStep] <= u_f[iStep])
    @constraint(M_f, zu_3_f[iStep=1:Hours_rolling+1], w_zu_f[iStep] >= z_f[iStep]+u_f[iStep]-1)
    

    # CONSTRAINTS ON DEGRADATION
    @constraint(M_f, deg_f_1[iStep=1:Hours_rolling], deg_f[iStep] >= soc_quad_f[iStep]/max_SOC^2 - soc_quad_f[iStep+1]/max_SOC^2 + (2/max_SOC)*(soc_f[iStep+1]-soc_f[iStep]) )  # 
    @constraint(M_f, deg_f_2[iStep=1:Hours_rolling], deg_f[iStep] >= soc_quad_f[iStep+1]/max_SOC^2 - soc_quad_f[iStep]/max_SOC^2 + (2/max_SOC)*(soc_f[iStep]-soc_f[iStep+1]) )  # 

    #CONSTRAINT ON ENERGY CAPACITY
    @constraint(M_f, initial_capacity_f[iStep=1], capacity_f[iStep] == min_SOH)
    @constraint(M_f, energy_capacity_f[iStep=1:Hours_rolling], capacity_f[iStep+1] == capacity_f[iStep]-deg_f[iStep]*k_f)
   
    return Rolling_problem(
        M_f,
        charge_f,
        discharge_f,
        soc_f,
        soc_quad_f,
        deg_f,
        x_f,
        y_f,
        z_f,
        u_f,
        w_xx_f,
        w_yy_f,
        w_zz_f,
        w_xy_f,
        w_xz_f,
        w_zy_f,
        w_uu_f,
        w_xu_f,
        w_yu_f,
        w_zu_f,
        capacity_f,
        State_initial,
        initial_capacity_f,
      )
end

# CASE 1


function endlife(InputParameters::InputParam, Battery::BatteryParam, ResultsOpt_3::Results_3, vec_p, Tot_steps)     # NO ROLLING PLANNING
  
  @unpack (NYears, NMonths, NStages, Big, NHoursStep, conv, disc) = InputParameters;     #NSteps,NHoursStage
  @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull, fix) = Battery ;   
  @unpack (soc,cap) = ResultsOpt_3;      

  k_e = min_SOH/(2*Nfull)
  Small_e = 0.64

  M_e = Model(Gurobi.Optimizer)
  set_optimizer_attribute(M_e, "MIPGap", 0.01)

  @variable(M_e, min_SOC <= soc_e[iStep=1:Tot_steps+1] <= max_SOC, base_name = "Energy_e")                # MWh   energy_Capacity NSteps
  @variable(M_e, min_SOC^2 <= soc_quad_e[iStep=1:Tot_steps+1] <= max_SOC^2, base_name = "Square energy_e")

  @variable(M_e, min_P <= charge_e[iStep=1:Tot_steps] <= max_P, base_name= "Charge_e")      #max_disc   0<=discharge<=1
  @variable(M_e, min_P <= discharge_e[iStep=1:Tot_steps] <= max_P, base_name= "Discharge_e")
  
  @variable(M_e, 0 <= deg_e[iStep=1:Tot_steps] <= Small_e, base_name = "Degradation_e")

  @variable(M_e, min_SOH <= capacity_e[iStep=1:Tot_steps+1] <= max_SOH, base_name = "Energy_Capacity_e")        #energy_Capacity     [iStage=1:NStages]

  #VARIABLES FOR DISCRETIZATION of Stored Energy

  @variable(M_e, x_e[iStep=1:Tot_steps+1], Bin, base_name = "Binary_1_end")
  @variable(M_e, y_e[iStep=1:Tot_steps+1], Bin, base_name = "Binary_2_end")
  @variable(M_e, z_e[iStep=1:Tot_steps+1], Bin, base_name = "Binary_3_end")
  @variable(M_e, u_e[iStep=1:Tot_steps+1], Bin, base_name = "Binary_4_end")
  
  @variable(M_e, 0<= w_xx_e[iStep=1:Tot_steps+1] <= 1, base_name = "xx_end")
  @variable(M_e, 0<= w_yy_e[iStep=1:Tot_steps+1] <= 1, base_name = "yy_end")
  @variable(M_e, 0<= w_zz_e[iStep=1:Tot_steps+1] <= 1, base_name = "zz_end")
  @variable(M_e, 0<= w_xy_e[iStep=1:Tot_steps+1] <= 1, base_name = "xy_end")
  @variable(M_e, 0<= w_xz_e[iStep=1:Tot_steps+1] <= 1, base_name = "xz_end")
  @variable(M_e, 0<= w_zy_e[iStep=1:Tot_steps+1] <= 1, base_name = "yz_end")

  @variable(M_e, 0 <= w_uu_e[iStep=1:Tot_steps+1] <=1, base_name = "uu_end")
  @variable(M_e, 0 <= w_xu_e[iStep=1:Tot_steps+1] <=1, base_name = "xu_end")
  @variable(M_e, 0 <= w_yu_e[iStep=1:Tot_steps+1] <=1, base_name = "yu_end")
  @variable(M_e, 0 <= w_zu_e[iStep=1:Tot_steps+1] <=1, base_name = "zu_end")

  @objective(
    M_e,
    MathOptInterface.MAX_SENSE, 
    sum(vec_p[iStep]*NHoursStep*discharge_e[iStep]-vec_p[iStep]*NHoursStep*charge_e[iStep] for iStep=1:Tot_steps)   #-deg_f[iStep]*k_f*100
    )

  @constraint(M_e, State_initial_e[iStep=1], soc_e[iStep] == soc_e[end])
  @constraint(M_e, energy_e[iStep=1:Tot_steps], soc_e[iStep+1] - (charge_e[iStep]*Eff_charge-discharge_e[iStep]/Eff_discharge)*NHoursStep == soc_e[iStep] )
  @constraint(M_e, cons_e[iStep=1:Tot_steps], soc_e[iStep] <= capacity_e[iStep])

  #@constraint(M, en_bal[iStep=1:NSteps+1], min_SOC + ((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]) == soc[iStep])
  @constraint(M_e, en_bal_e[iStep=1:Tot_steps+1], min_SOC + ((max_SOC-min_SOC)/disc)*(x_e[iStep]+2*y_e[iStep]+4*z_e[iStep]+8*u_e[iStep]) == soc_e[iStep])    
  
  #@constraint(M, en_square[iStep=1:NSteps+1], soc_quad[iStep] == min_SOC^2+ 2*min_SOC*((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep])+(w_xx[iStep]+4*w_xy[iStep]+8*w_xz[iStep]+4*w_yy[iStep]+16*w_zz[iStep]+16*w_zy[iStep])*((max_SOC-min_SOC)/disc)^2)
  @constraint(M_e, en_square_e[iStep=1:Tot_steps+1], soc_quad_e[iStep] == min_SOC^2+ 2*min_SOC*((max_SOC-min_SOC)/disc)*(x_e[iStep]+2*y_e[iStep]+4*z_e[iStep]+8*u_e[iStep])+(w_xx_e[iStep]+4*w_xy_e[iStep]+8*w_xz_e[iStep]+16*w_xu_e[iStep]+4*w_yy_e[iStep]+16*w_zz_e[iStep]+16*w_zy_e[iStep]+32*w_yu_e[iStep]+64*w_zu_e[iStep]+64*w_uu_e[iStep])*((max_SOC-min_SOC)/disc)^2)

  # INEQUALITIES CONSTRAINTS
  @constraint(M_e, xx_1_e[iStep=1:Tot_steps+1], w_xx_e[iStep] <= x_e[iStep])
  @constraint(M_e, xx_2_e[iStep=1:Tot_steps+1], w_xx_e[iStep] >= 2*x_e[iStep]-1)

  @constraint(M_e, xy_1_e[iStep=1:Tot_steps+1], w_xy_e[iStep] <= x_e[iStep])
  @constraint(M_e, xy_2_e[iStep=1:Tot_steps+1], w_xy_e[iStep] <= y_e[iStep])
  @constraint(M_e, xy_3_e[iStep=1:Tot_steps+1], w_xy_e[iStep] >= x_e[iStep]+y_e[iStep]-1)

  @constraint(M_e, xz_1_e[iStep=1:Tot_steps+1], w_xz_e[iStep] <= x_e[iStep])
  @constraint(M_e, xz_2_e[iStep=1:Tot_steps+1], w_xz_e[iStep] <= z_e[iStep])
  @constraint(M_e, xz_3_e[iStep=1:Tot_steps+1], w_xz_e[iStep] >= x_e[iStep]+z_e[iStep]-1)

  @constraint(M_e, yy_1_e[iStep=1:Tot_steps+1], w_yy_e[iStep] <= y_e[iStep])
  @constraint(M_e, yy_2_e[iStep=1:Tot_steps+1], w_yy_e[iStep] >= 2*y_e[iStep]-1)

  @constraint(M_e, zz_1_e[iStep=1:Tot_steps+1], w_zz_e[iStep] <= z_e[iStep])
  @constraint(M_e, zz_2_e[iStep=1:Tot_steps+1], w_zz_e[iStep] >= 2*z_e[iStep]-1)

  @constraint(M_e, zy_1_e[iStep=1:Tot_steps+1], w_zy_e[iStep] <= z_e[iStep])
  @constraint(M_e, zy_2_e[iStep=1:Tot_steps+1], w_zy_e[iStep] <= y_e[iStep])
  @constraint(M_e, zy_3_e[iStep=1:Tot_steps+1], w_zy_e[iStep] >= z_e[iStep]+y_e[iStep]-1)

  @constraint(M_e, uu_1_e[iStep=1:Tot_steps+1], w_uu_e[iStep] <= u_e[iStep])
  @constraint(M_e, uu_2_e[iStep=1:Tot_steps+1], w_uu_e[iStep] >= 2*u_e[iStep]-1)

  @constraint(M_e, xu_1_e[iStep=1:Tot_steps+1], w_xu_e[iStep] <= x_e[iStep])
  @constraint(M_e, xu_2_e[iStep=1:Tot_steps+1], w_xu_e[iStep] <= u_e[iStep])
  @constraint(M_e, xu_3_e[iStep=1:Tot_steps+1], w_xu_e[iStep] >= x_e[iStep]+u_e[iStep]-1)

  @constraint(M_e, yu_1_e[iStep=1:Tot_steps+1], w_yu_e[iStep] <= y_e[iStep])
  @constraint(M_e, yu_2_e[iStep=1:Tot_steps+1], w_yu_e[iStep] <= u_e[iStep])
  @constraint(M_e, yu_3_e[iStep=1:Tot_steps+1], w_yu_e[iStep] >= y_e[iStep]+u_e[iStep]-1)

  @constraint(M_e, zu_1_e[iStep=1:Tot_steps+1], w_zu_e[iStep] <= z_e[iStep])
  @constraint(M_e, zu_2_e[iStep=1:Tot_steps+1], w_zu_e[iStep] <= u_e[iStep])
  @constraint(M_e, zu_3_e[iStep=1:Tot_steps+1], w_zu_e[iStep] >= z_e[iStep]+u_e[iStep]-1)
  

  # CONSTRAINTS ON DEGRADATION
  @constraint(M_e, deg_e_1[iStep=1:Tot_steps], deg_e[iStep] >= soc_quad_e[iStep]/max_SOC^2 - soc_quad_e[iStep+1]/max_SOC^2 + (2/max_SOC)*(soc_e[iStep+1]-soc_e[iStep]) )  # 
  @constraint(M_e, deg_e_2[iStep=1:Tot_steps], deg_e[iStep] >= soc_quad_e[iStep+1]/max_SOC^2 - soc_quad_e[iStep]/max_SOC^2 + (2/max_SOC)*(soc_e[iStep]-soc_e[iStep+1]) )  # 

  #CONSTRAINT ON ENERGY CAPACITY
  @constraint(M_e, initial_capacity_e[iStep=1], capacity_e[iStep] == cap[end])
  @constraint(M_e, energy_capacity_e[iStep=1:Tot_steps], capacity_e[iStep+1] == capacity_e[iStep]-deg_e[iStep]*k_e)


return Case1(
        M_e,
        charge_e,
        discharge_e,
        soc_e,
        soc_quad_e,
        deg_e,
        x_e,
        y_e,
        z_e,
        u_e,
        w_xx_e,
        w_yy_e,
        w_zz_e,
        w_xy_e,
        w_xz_e,
        w_zy_e,
        w_uu_e,
        w_xu_e,
        w_yu_e,
        w_zu_e,
        capacity_e,
        State_initial_e,
        initial_capacity_e,
)

end

# CASE 3


function New_deg_formulation(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam)       

  @unpack (MIPGap, MIPFocus, Method, Cuts, Heuristics) = SolverParameters;

  @unpack (NYears, NMonths, NStages, Big, NHoursStep, conv, disc, Hours_rolling, Hours_saved) = InputParameters;     #NSteps,NHoursStage
  @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull, fix) = Battery ;         

  Small_new = 0.64

  M_new = Model(Gurobi.Optimizer)
  set_optimizer_attribute(M_new, "MIPGap", 0.1)

  # DEFINE VARIABLES

  @variable(M_new, min_SOC <= soc_new[iStep=1:Hours_rolling+1] <= max_SOC, base_name = "Energy_new")                # MWh   energy_Capacity NSteps
  @variable(M_new, min_SOC^2 <= soc_quad_new[iStep=1:Hours_rolling+1] <= max_SOC^2, base_name = "Square energy_new")

  @variable(M_new, min_P <= charge_new[iStep=1:Hours_rolling] <= max_P, base_name= "Charge_new")      #max_disc   0<=discharge<=1
  @variable(M_new, min_P <= discharge_new[iStep=1:Hours_rolling] <= max_P, base_name= "Discharge_new")
  
  @variable(M_new, 0 <= deg_new[iStep=1:Hours_rolling] <= Small_new, base_name = "Degradation_new")

  @variable(M_new, 0 <= capacity_new[iStep=1:Hours_rolling+1] <= max_SOH, base_name = "Energy_Capacity_new")        #energy_Capacity     [iStage=1:NStages]

  #VARIABLES FOR DISCRETIZATION of Stored Energy

  @variable(M_new, x_new[iStep=1:Hours_rolling+1], Bin, base_name = "Binary_1_new")
  @variable(M_new, y_new[iStep=1:Hours_rolling+1], Bin, base_name = "Binary_2_new")
  @variable(M_new, z_new[iStep=1:Hours_rolling+1], Bin, base_name = "Binary_3_new")
  @variable(M_new, u_new[iStep=1:Hours_rolling+1], Bin, base_name = "Binary_4_new")
  
  @variable(M_new, 0<= w_xx_new[iStep=1:Hours_rolling+1] <= 1, base_name = "xx_new")
  @variable(M_new, 0<= w_yy_new[iStep=1:Hours_rolling+1] <= 1, base_name = "yy_new")
  @variable(M_new, 0<= w_zz_new[iStep=1:Hours_rolling+1] <= 1, base_name = "zz_new")
  @variable(M_new, 0<= w_xy_new[iStep=1:Hours_rolling+1] <= 1, base_name = "xy_new")
  @variable(M_new, 0<= w_xz_new[iStep=1:Hours_rolling+1] <= 1, base_name = "xz_new")
  @variable(M_new, 0<= w_zy_new[iStep=1:Hours_rolling+1] <= 1, base_name = "yz_new")

  @variable(M_new, 0 <= w_uu_new[iStep=1:Hours_rolling+1] <=1, base_name = "uu_new")
  @variable(M_new, 0 <= w_xu_new[iStep=1:Hours_rolling+1] <=1, base_name = "xu_new")
  @variable(M_new, 0 <= w_yu_new[iStep=1:Hours_rolling+1] <=1, base_name = "yu_new")
  @variable(M_new, 0 <= w_zu_new[iStep=1:Hours_rolling+1] <=1, base_name = "zu_new")


  # DEFINE OBJECTIVE function 

  @objective(
    M_f,
    MathOptInterface.MAX_SENSE, 
    sum(1*NHoursStep*discharge_new[iStep]-1*NHoursStep*charge_new[iStep] for iStep=1:Hours_rolling)   #-deg_f[iStep]*k_f*100
    )
       
  # DEFINE CONSTRAINTS

  @constraint(M_new, State_initial_new[iStep=1], soc_new[iStep] == max_SOH)
  @constraint(M_new, energy_new[iStep=1:Hours_rolling], soc_new[iStep+1] - (charge_new[iStep]*Eff_charge-discharge_new[iStep]/Eff_discharge)*NHoursStep == soc_new[iStep] )
  @constraint(M_new, cons_new[iStep=1:Hours_rolling], soc_new[iStep] <= capacity_new[iStep])

  #@constraint(M, en_bal[iStep=1:NSteps+1], min_SOC + ((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]) == soc[iStep])
  @constraint(M_new, en_bal_new[iStep=1:Hours_rolling+1], min_SOC + ((max_SOC-min_SOC)/disc)*(x_new[iStep]+2*y_new[iStep]+4*z_new[iStep]+8*u_new[iStep]) == soc_new[iStep])    
  
  #@constraint(M, en_square[iStep=1:NSteps+1], soc_quad[iStep] == min_SOC^2+ 2*min_SOC*((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep])+(w_xx[iStep]+4*w_xy[iStep]+8*w_xz[iStep]+4*w_yy[iStep]+16*w_zz[iStep]+16*w_zy[iStep])*((max_SOC-min_SOC)/disc)^2)
  @constraint(M_new, en_square_new[iStep=1:Hours_rolling+1], soc_quad_new[iStep] == min_SOC^2+ 2*min_SOC*((max_SOC-min_SOC)/disc)*(x_new[iStep]+2*y_new[iStep]+4*z_new[iStep]+8*u_new[iStep])+(w_xx_new[iStep]+4*w_xy_new[iStep]+8*w_xz_new[iStep]+16*w_xu_new[iStep]+4*w_y_new[iStep]+16*w_zz_new[iStep]+16*w_zy_new[iStep]+32*w_yu_new[iStep]+64*w_zu_new[iStep]+64*w_uu_new[iStep])*((max_SOC-min_SOC)/disc)^2)

  # INEQUALITIES CONSTRAINTS
  @constraint(M_new, xx_1_new[iStep=1:Hours_rolling+1], w_xx_new[iStep] <= x_new[iStep])
  @constraint(M_new, xx_2_new[iStep=1:Hours_rolling+1], w_xx_new[iStep] >= 2*x_new[iStep]-1)

  @constraint(M_new, xy_1_new[iStep=1:Hours_rolling+1], w_xy_new[iStep] <= x_new[iStep])
  @constraint(M_new, xy_2_new[iStep=1:Hours_rolling+1], w_xy_new[iStep] <= y_new[iStep])
  @constraint(M_new, xy_3_new[iStep=1:Hours_rolling+1], w_xy_new[iStep] >= x_new[iStep]+y_new[iStep]-1)

  @constraint(M_new, xz_1_new[iStep=1:Hours_rolling+1], w_xz_new[iStep] <= x_new[iStep])
  @constraint(M_new, xz_2_new[iStep=1:Hours_rolling+1], w_xz_new[iStep] <= z_new[iStep])
  @constraint(M_new, xz_3_new[iStep=1:Hours_rolling+1], w_xz_new[iStep] >= x_new[iStep]+z_new[iStep]-1)

  @constraint(M_new, yy_1_new[iStep=1:Hours_rolling+1], w_yy_new[iStep] <= y_new[iStep])
  @constraint(M_new, yy_2_new[iStep=1:Hours_rolling+1], w_yy_new[iStep] >= 2*y_new[iStep]-1)

  @constraint(M_new, zz_1_new[iStep=1:Hours_rolling+1], w_zz_new[iStep] <= z_new[iStep])
  @constraint(M_new, zz_2_new[iStep=1:Hours_rolling+1], w_zz_new[iStep] >= 2*z_new[iStep]-1)

  @constraint(M_new, zy_1_new[iStep=1:Hours_rolling+1], w_zy_new[iStep] <= z_new[iStep])
  @constraint(M_new, zy_2_new[iStep=1:Hours_rolling+1], w_zy_new[iStep] <= y_new[iStep])
  @constraint(M_new, zy_3_new[iStep=1:Hours_rolling+1], w_zy_new[iStep] >= z_new[iStep]+y_new[iStep]-1)

  @constraint(M_new, uu_1_new[iStep=1:Hours_rolling+1], w_uu_new[iStep] <= u_new[iStep])
  @constraint(M_new, uu_2_new[iStep=1:Hours_rolling+1], w_uu_new[iStep] >= 2*u_new[iStep]-1)

  @constraint(M_new, xu_1_new[iStep=1:Hours_rolling+1], w_xu_new[iStep] <= x_new[iStep])
  @constraint(M_new, xu_2_new[iStep=1:Hours_rolling+1], w_xu_new[iStep] <= u_new[iStep])
  @constraint(M_new, xu_3_new[iStep=1:Hours_rolling+1], w_xu_new[iStep] >= x_new[iStep]+u_new[iStep]-1)

  @constraint(M_new, yu_1_new[iStep=1:Hours_rolling+1], w_yu_new[iStep] <= y_new[iStep])
  @constraint(M_new, yu_2_new[iStep=1:Hours_rolling+1], w_yu_new[iStep] <= u_new[iStep])
  @constraint(M_new, yu_3_new[iStep=1:Hours_rolling+1], w_yu_new[iStep] >= y_new[iStep]+u_new[iStep]-1)

  @constraint(M_new, zu_1_new[iStep=1:Hours_rolling+1], w_zu_new[iStep] <= z_new[iStep])
  @constraint(M_new, zu_2_new[iStep=1:Hours_rolling+1], w_zu_new[iStep] <= u_new[iStep])
  @constraint(M_new, zu_3_new[iStep=1:Hours_rolling+1], w_zu_new[iStep] >= z_new[iStep]+u_new[iStep]-1)
  

  # CONSTRAINTS ON DEGRADATION
  @constraint(M_new, deg_new_1[iStep=1:Hours_rolling], deg_new[iStep] >= (min_SOH)/(2*Nfull*max_SOC^2)*soc_quad_new[iStep] - (min_SOH)/(2*Nfull*max_SOC^2)*soc_quad_new[iStep+1] + (min_SOH/(max_SOC*Nfull))*soc_new[iStep+1]-(min_SOH/(max_SOC*Nfull))*soc_new[iStep] )  # 
  @constraint(M_new, deg_new_2[iStep=1:Hours_rolling], deg_new[iStep] >= (min_SOH)/(2*Nfull*max_SOC^2)*soc_quad_new[iStep+1] - (min_SOH)/(2*Nfull*max_SOC^2)*soc_quad_new[iStep] + (min_SOH/(max_SOC*Nfull))*soc_new[iStep]-(min_SOH/(max_SOC*Nfull))*soc_new[iStep+1] )  # 

  #CONSTRAINT ON ENERGY CAPACITY
  @constraint(M_new, initial_capacity_new[iStep=1], capacity_new[iStep] == min_SOH)
  @constraint(M_new, energy_capacity_new[iStep=1:Hours_rolling], capacity_new[iStep+1] == capacity_new[iStep]-deg_new[iStep])
 
  return Case3(
      M_new,
      charge_new,
      discharge_new,
      soc_new,
      soc_quad_new,
      deg_new,
      x_new,
      y_new,
      z_new,
      u_new,
      w_xx_new,
      w_yy_new,
      w_zz_new,
      w_xy_new,
      w_xz_new,
      w_zy_new,
      w_uu_new,
      w_xu_new,
      w_yu_new,
      w_zu_new,
      capacity_new,
      State_initial_new,
      initial_capacity_new,
    )
end
