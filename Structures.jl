# STRUCTURES USED IN THE PROBLEM

# Input data
#-----------------------------------------------

# Input parameters 
@with_kw struct InputParam{F<:Float64,I<:Int}
    NYears::F                                     # Number of years
    NMonths::I
    NStages::I                                    # Number of stages of N months in the problem FORMULATION-- calcolato come NYears/NMonths*12
    NHoursStep::F                                 # Number of hours in each time step 
    Hours_rolling::I
    Hours_saved::I
    Big::F                                        # A big number
    conv::F                                       # A small number for degradation convergence
    bin::I                                       # Discretization points
end

# Battery's characteristics
@with_kw struct BatteryParam{F<:Float64,I<:Int}
    min_SOC::F                                     # Battery's minimum energy storage capacity
    max_SOC::F                                     # Batter's maximum energy storage capacity
    Eff_charge::F                                  # Battery's efficiency for charging operation
    Eff_discharge::F                               # Battery's efficiency for discharging operation
    min_P::F
    max_P::F
    max_SOH::F                                     # Maximum SOH that can be achieved because of volume issues
    min_SOH::F                                     # Minimum SOH to be respected by contract
    Nfull::I                                       # Maximum number of full cycles for DoD=100%     
    fix::F                                         # fixed costs for battery replacement                                   
    cost::F
end
  
# solver parameters
@with_kw struct SolverParam{F<:Float64,I<:Int}
    MIPGap::F 
    MIPFocus::I
    Method::F
    Cuts::F
    Heuristics::F
end
  
# Indirizzi cartelle
@with_kw struct caseData{S<:String}
    DataPath::S
    InputPath::S
    ResultPath::S
    CaseName::S
end

# runMode Parameters
@with_kw mutable struct runModeParam{B<:Bool}

    # Solver settings
    solveMIP::B     #If using SOS2

    batterySystemFromFile::B 

    #runMode self defined reading of input 
    setInputParameters::B             #from .in file
 
    excel_savings::B 
    endlife::B
    endlife_rolling::B
    rolling_newDeg::B

end

# Optimization problem
struct BuildStageProblem_3
    M::Any
    soc::Any
    soc_quad::Any
    charge::Any 
    discharge::Any
    deg::Any
    x::Any
    y::Any
    z::Any
    xy::Any
    xz::Any
    zy::Any
    xyz::Any
    capacity::Any
    revamping::Any
    e::Any
    rev_vendita::Any
    rev_acquisto::Any
end

struct Results_3
    objective::Any
    #revenues_per_stage::Any
    gain_stage::Any
    cost_rev::Any
    deg_stage::Any
    soc::Any
    charge::Any
    discharge::Any
    #bin::Any
    deg::Any
    soc_quad::Any
    x::Any
    y::Any
    z::Any
    xy::Any
    xz::Any
    yz::Any
    xyz::Any
    rev::Any
    cap::Any
    e::Any
    rev_vendita::Any
    rev_acquisto::Any
end

# Optimization problem
struct BuildStageProblem_4
    M::Any
    soc::Any
    soc_quad::Any
    charge::Any 
    discharge::Any
    deg::Any
    x::Any
    y::Any
    z::Any
    u::Any
    xy::Any
    xz::Any
    xu::Any
    zy::Any
    yu::Any
    zu::Any
    xyzu::Any
    capacity::Any
    revamping::Any
    e::Any
    rev_vendita::Any
    rev_acquisto::Any
end

struct Results_4
    objective::Any
    #revenues_per_stage::Any
    gain_stage::Any
    cost_rev::Any
    deg_stage::Any
    soc::Any
    charge::Any
    discharge::Any
    #bin::Any
    deg::Any
    soc_quad::Any
    x::Any
    y::Any
    z::Any
    u::Any
    xy::Any
    xz::Any
    xu::Any
    yz::Any
    yu::Any
    zu::Any
    xyzu::Any
    rev::Any
    cap::Any
    e::Any
    rev_vendita::Any
    rev_acquisto::Any
end

# RESULTS ROLLING PLANNING (CASO 2)
struct ResultsEndLife       #Raccolta dati e aggiornamento parametri input
    M_end::Any
    Ch_end::Any
    Dis_end::Any
    soc_end::Any
    socq_end::Any
    deg_end::Any
    cap_end::Any
    x_end::Any
    y_end::Any
    z_end::Any
    u_end::Any
    w_xx_end::Any
    w_yy_end::Any
    w_zz_end::Any
    w_uu_end::Any
    w_xy_end::Any
    w_xz_end::Any
    w_zy_end::Any
    w_xu_end::Any
    w_yu_end::Any
    w_zu_end::Any
    new_vec::Any
end

struct Rolling_problem
    M_f::Any
    charge_f::Any
    discharge_f::Any
    soc_f::Any
    soc_quad_f::Any
    deg_f::Any
    x_f::Any
    y_f::Any
    z_f::Any
    u_f::Any
    w_xx_f::Any
    w_yy_f::Any
    w_zz_f::Any
    w_xy_f::Any
    w_xz_f::Any
    w_zy_f::Any
    w_uu_f::Any
    w_xu_f::Any
    w_yu_f::Any
    w_zu_f::Any
    capacity_f::Any
    State_initial::Any
    initial_capacity_f::Any
end

# RESULTS CASE 1 - NO ROLLING PLANNING

struct Case1
    M_e::Any
    charge_e::Any
    discharge_e::Any
    soc_e::Any
    soc_quad_e::Any
    deg_e::Any
    x_e::Any
    y_e::Any
    z_e::Any
    u_e::Any
    w_xx_e::Any
    w_yy_e::Any
    w_zz_e::Any
    w_xy_e::Any
    w_xz_e::Any
    w_zy_e::Any
    w_uu_e::Any
    w_xu_e::Any
    w_yu_e::Any
    w_zu_e::Any
    capacity_e::Any
    State_initial_e::Any
    initial_capacity_e::Any
end

struct ResultsCase1
    soc_e1::Any
    charge_e1::Any
    discharge_e1::Any
    deg_e1::Any
    soc_quad_e1::Any
    x_e1::Any
    y_e1::Any
    z_e1::Any
    u_e1::Any
    w_xx_e1::Any
    w_yy_e1::Any
    w_zz_e1::Any
    w_uu_e1::Any
    w_xy_e1::Any
    w_xz_e1::Any
    w_zy_e1::Any
    w_xu_e1::Any
    w_yu_e1::Any
    w_zu_e1::Any
    cap_e1::Any
end

# CASO 3

struct Case3
    M_new::Any
    charge_new::Any
    discharge_new::Any
    soc_new::Any
    soc_quad_new::Any
    deg_new::Any
    x_new::Any
    y_new::Any
    z_new::Any
    u_new::Any
    w_xx_new::Any
    w_yy_new::Any
    w_zz_new::Any
    w_xy_new::Any
    w_xz_new::Any
    w_zy_new::Any
    w_uu_new::Any
    w_xu_new::Any
    w_yu_new::Any
    w_zu_new::Any
    capacity_new::Any
    State_initial_new::Any
    initial_capacity_new::Any
end

struct ResultsNewDeg
    problem_new_deg::Any
    charge_new_deg::Any
    discharge_new_deg::Any
    state_new_deg::Any
    stateQuad_new_deg::Any
    degradation_new_deg::Any
    capacity_new_deg::Any
    x_new_deg::Any
    y_new_deg::Any
    z_new_deg::Any
    u_new_deg::Any
    w_xx_new_deg::Any
    w_yy_new_deg::Any
    w_zz_new_deg::Any
    w_uu_new_deg::Any
    w_xy_new_deg::Any
    w_xz_new_deg::Any
    w_zy_new_deg::Any
    w_xu_new_deg::Any
    w_yu_new_deg::Any
    w_zu_new_deg::Any
    new_new_deg::Any
end