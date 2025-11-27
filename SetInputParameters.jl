# General read and formatting functions
#--------------------------------------

function read_from_file_to_dict!(f, vars::Dict{Symbol,Bool})
  for line in eachline(f)
    if !isempty(strip(line)) && !startswith(strip(line), "#")
      var, val = strip.(split(line, "="))
      try
        vars[Symbol(strip(var))] = parse(Bool, val)
      catch
      end

    end
  end
  return vars
end

function read_from_file_to_dict!(f, vars::Dict{Symbol,Float64})
  for line in eachline(f)
    if !isempty(strip(line)) && !startswith(strip(line), "#")
      var, val = strip.(split(line, "="))
      try
        vars[Symbol(strip(var))] = parse(Float64, val)
      catch
      end

    end
  end
  return vars
end

function read_from_file_to_dict!(f, vars::Dict{Symbol,String})
  for line in eachline(f)
    if !isempty(strip(line)) && !startswith(strip(line), "#")
      var, val = strip.(split(line, "="))
      try
        vars[Symbol(strip(var))] = val
      catch
      end
    end
  end
  return vars
end

function read_from_file_to_dict!(f, vars::Dict{Symbol,Any})
  for line in eachline(f)
    if !isempty(strip(line)) && !startswith(strip(line), "#")
      var, val = strip.(split(line, "="))
      try
        vars[Symbol(strip(var))] = val
      catch
      end
    end
  end
  return vars
end

function read_from_file_to_dict!(f, vars::Dict{Symbol,Number})
  # Use Dict{Symbol, Number} to allow float and int values
  # NB: all vaues read as float from file initially
  for line in eachline(f)
    if !isempty(strip(line)) && !startswith(strip(line), "#")
      var, val = strip.(split(line, "="))
      try
        vars[Symbol(strip(var))] = parse(Float64, val)
      catch
      end
    end
  end
  return vars
end

function read_type_to_dict(file, type::DataType)
  f = open(file)
  paramDict = Dict{Symbol,type}()
  paramDict = read_from_file_to_dict!(f, paramDict)

  return paramDict
end

function set_integers!(paramDict, integers)
  for (key, value) in paramDict
    if key in integers
      paramDict[key] = parse(Int64,value)
    end
  end
  return paramDict
end

function set_floats!(paramDict, float)
  for (key, value) in paramDict
    if key in float
      paramDict[key] = parse(Float64,value)
    end
  end
  return paramDict
end

function set_vector!(paramDict, vector)
  value = paramDict[vector]
  items = split(value," ")
  NSeg = zeros(Int64,2)
  NSeg[1] = parse(Int64,items[1])
  NSeg[2] = parse(Int64,items[2])
  return NSeg
end

function read_csv(file, path)
  try
    println("Reading from: ", file)
    data = Matrix{Float64}(CSV.read(joinpath(path, file), DataFrame, header = false)) # €/MWh
    return data
  catch e
    println("Can't read ", file)
    throw(error())
  end
end

# Functions to read and prep from predefined files
#-----------------------------------------

# parameters
function read_parameters_from_config_file(file = "configParameters.in")

  #paramDict = read_type_to_dict(file,Number)
  paramDict = read_type_to_dict(file, Any)
  
  integers =[:NMonths :disc :NSteps :Hours_rolling :Hours_saved]  #NYears
  paramDict = set_integers!(paramDict, integers)

  floats =[:NYears :Big :NHoursStep :conv]
  paramDict = set_floats!(paramDict,floats)


  # scrivere codice "if" fino a quando il resto tra NYears*12/Nmonths non sia nullo - chiedere di cambiare dati
  paramDict[:NStages] = Int(paramDict[:NYears]*12/paramDict[:NMonths])
  #paramDict[:NSteps] = Int(paramDict[:NYears]*8760/paramDict[:NHoursStep])    #8760
  #paramDict[:NHoursStage] = Int(paramDict[:NSteps]/paramDict[:NStages])

  inputData = InputParam(;paramDict...)

  println("Parameters to be used:", paramDict)

  return inputData
end

function read_solverParameters_from_file(file = "solverParameters.in")
  try
    println("Reading solver parameters from config file.")
    #paramDict = read_type_to_dict(file, Number)
    paramDict = read_type_to_dict(file, Any)

    integers = [:MIPFocus]
    #integers = [:CPX_PARAM_SCRIND :CPX_PARAM_PREIND :CPX_PARAM_TILIM :CPX_PARAM_THREADS]
    paramDict = set_integers!(paramDict, integers)

    floats = [:MIPGap :Method :Cuts :Heuristics]
    paramDict = set_floats!(paramDict,floats)

    inputData = SolverParam(; paramDict...)

    return inputData
  catch e
    println("Can't set solver parameters.")
    throw(error())
  end
end


function set_parameters(runMode::runModeParam, case::caseData)
  try
    if runMode.setInputParameters
      println("Reading parameters from configuration file.")
      InputParameters = read_parameters_from_config_file("configParameters.in")
    else
      println("Reading inputparameters from JLD file.")
      InputFromFile = read_input_from_jld(case.InputPath, case.InputCase, "InputParam.jld")
      InputParameters = InputFromFile["InputParameters"]
    end

    return InputParameters
  catch e
    println("Can't set input parameters.")
    throw(error())
  end
end


function set_solverParameters()
  solverParameters = read_solverParameters_from_file()
  return solverParameters
end

# MODE SETTING
function read_runMode_file(file = "runMode.in")

  runModeDict = read_type_to_dict(file, Bool)
  runMode = runModeParam(; runModeDict...)

  println(runMode)

  return runMode
end


# INDIRIZZI CARTELLE - CASE SETTINGS
function set_runCase(file = "runCase.in")

  runCaseDict = read_type_to_dict(file, String)
  case = caseData(; runCaseDict...)

  println(caseData)

  return case
end

# BATTERY PARAMETERS SETTING
function set_battery_system(runMode::runModeParam, case::caseData)
  try
    if runMode.batterySystemFromFile
      println("Reading batterys characteristics from input files")
      Battery = read_Battery_from_file("BatteryCharacteristics.in")
    else
      println("Reading battery from JLD file.")
      InputFromFile = read_input_from_jld(case.InputPath, case.InputCase, "BatteryCharacteristics.jld")
      Battery = InputFromFile["BatteryCharacteristics"]
    end
    return Battery
  catch e
    println("Can't set battery parameters.")
    throw(error())
  end
end

function read_Battery_from_file(file = "BatteryCharacteristics.in")

  paramDict = read_type_to_dict(file, Any)
  println("Parameters to be used:", paramDict)
  integers =[:Nfull]
  paramDict = set_integers!(paramDict,integers)

  floats =[:min_SOC :max_SOC :Eff_charge :Eff_discharge :min_P :max_P :max_SOH :min_SOH :fix :cost]
  paramDict = set_floats!(paramDict,floats)

  Battery = BatteryParam(;paramDict...)

  return Battery
end

function set_price(Battery_price_purchase,cost)
  Battery_price_sale = []
  n=length(Battery_price_purchase)
  for i=1:n
    a=Battery_price_purchase[i]*cost   #0.95
    push!(Battery_price_sale,a)
  end
  return Battery_price_sale
end


function read_input_from_jld(InputPath, InputCase, filname)
  inputFile = InputCase * filname
  try
    InputFromFile = load(joinpath(InputPath, inputFile))
    return InputFromFile
  catch e
    println(string("Could not read input: ", inputFile))
  end
end

# set case run name for saving results

function set_run_name(case, ResultPath, NSteps)

  n=NSteps
  param = string("NSteps_",n)

  Finalpath = joinpath(ResultPath, splitdir(case.DataPath)[end])
  Finalpath = joinpath(Finalpath, param)

  runName = string(date, case.CaseName)

  println(string("Run name: ", runName))
  println(string("Result directory: ", joinpath(Finalpath, runName)))

  return joinpath(Finalpath, runName)
end


function calculate_coefficients(min_SOC,max_SOC,bin)
  disc = 2^bin-1
  a= zeros(disc+1);
  b= zeros(disc+1);

  step = (max_SOC-min_SOC)/disc
  a[1]= min_SOC
  b[1] = min_SOC^2
  c=0

  for i =1:disc
    a[i+1]= a[i]+step
    b[i+1] = a[i+1]^2
  end

  if bin == 3 
    include("ProblemForm.jl")
    include("solveOptimizationAlgorithm.jl")
    c = zeros(7);
    c[1]= b[2]-b[1]
    c[2]= b[3]-b[1]
    c[3]= b[4]-b[3]-b[2]+b[1]
    c[4]= b[5]-b[1]
    c[5]= b[6]-b[2]-b[5]+b[1]
    c[6]= b[7]-b[3]-b[5]+b[1]
    c[7]= 0 #b[8]-b[4]-b[6]-b[7]+b[2]+b[3]+b[5]-b[1]
  elseif bin ==4
    include("ProblemForm_4_bin.jl")
    include("solveOptimizationAlgorithm_4_bin.jl")
    c = zeros(11)
    c[1]= b[2]-b[1]
    c[2]= b[3]-b[1]
    c[3]= b[5]-b[1]
    c[4]= b[9]-b[1]
    c[5]= b[4]-b[2]-b[3]+b[1]
    c[6]= b[6]-b[2]-b[5]+b[1]
    c[7]= b[10]-b[9]-b[2]+b[1]
    c[8]= b[7]-c[3]-c[2]-b[1]
    c[9]= b[11]-b[1]-c[4]-c[2]
    c[10]= b[13]-b[1]-c[4]-c[3]
    c[11]= 0
  else
    println("Change number of binary variables")
    throw(error())
  end

  return (a,b,c,disc)
end




