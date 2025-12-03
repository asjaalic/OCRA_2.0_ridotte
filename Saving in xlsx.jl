# EXCEL SAVINGS
#using DataFrames
#using XLSX

function data_saving_3(InputParameters::InputParam,ResultsOpt_3::Results_3)

    @unpack (NYears, NMonths, NStages, Big, NHoursStep, disc) = InputParameters;       #NSteps,NHoursStage
    
   #@unpack (charge,discharge, soc, revenues_per_stage, x, y, z, w_xx, w_yy, w_zz, w_xy, w_xz, w_zy) rev_vendita, rev_acquisto = ResultsOpt;  
   @unpack (charge,discharge,rev,cap,e, soc,soc_quad,deg,deg_stage,gain_stage, cost_rev,rev_vendita, rev_acquisto,x, y, z, xy, xz, yz, xyz) = ResultsOpt_3;
   @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull,fix,cost ) = Battery ; 

    hour=string(now())
    a=replace(hour,':'=> '-')

    nameF= "Test 3 variabili 02.12 con TuneParameters"
    nameFile="Summary " 

    folder = "$nameF"
    mkdir(folder)
    cd(folder)
    main=pwd()

    general = DataFrame()
    #general_1 = DataFrame()
    battery_costs= DataFrame()

    capacity=zeros(NStages)
    capacity_finale = zeros(NStages)

    for iStage=1:NStages
        capacity[iStage]=cap[(Steps_stages[iStage]+2)]
    end

    for iStage=2:NStages
        capacity_finale[iStage-1] = cap[Steps_stages[iStage]+1]
    end
        capacity_finale[end] = cap[end]
    
    general[!,"Stage"] = 1:1:NStages
    general[!,"Initial Capacity MWh"] = capacity[:]
    general[!,"Final Capacity MWh"] = capacity_finale[:]
    general[!,"Revamping MWh"] = rev[:]
    general[!, "Binary revamp"] = e[:]
    general[!,"Degradation"] = deg_stage[:]
    #general[!,"Net_Revenues"] = revenues_per_stage[:]
    general[!,"Gain charge/discharge"] = gain_stage[:]
    general[!,"Cost revamping"] = cost_rev[:]
    general[!, "Aux acquisto"] = rev_acquisto[:]
    general[!, "Aux vendita"] = rev_vendita[:]

    battery_costs[!,"Purchase costs €/MWh"] = Battery_price_purchase[1:NStages+1]    #Battery_price_purchase
    battery_costs[!,"Sale costs €/MWh"] = Battery_price_sale[1:NStages+1] 

    XLSX.writetable("$nameFile.xlsx", overwrite=true,                                       #$nameFile
    results_stages = (collect(DataFrames.eachcol(general)),DataFrames.names(general)),
    #results_1 = (collect(DataFrames.eachcol(general_1)),DataFrames.names(general_1)),
    costs = (collect(DataFrames.eachcol(battery_costs)),DataFrames.names(battery_costs)),
    )

    for iStage=1:NStages
        steps = DataFrame()

        steps[!,"Step"] = (Steps_stages[iStage]+1):(Steps_stages[iStage+1])
        steps[!, "Energy_prices €/MWh"] = Power_prices[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Energy capacity MWh"] = cap[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "SOC MWh"] = soc[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Charge MW"] = charge[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Discharge MW"] = discharge[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "SOC_quad MWh"] = soc_quad[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Deg MWh"] = deg[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "X"] = x[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Y"] = y[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Z"] = z[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "XY"] = xy[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "XZ"] = xz[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "YZ"] = yz[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "XYZ"] = xyz[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]


        XLSX.writetable("$iStage stage $a.xlsx", overwrite=true,                                       #$nameFile
        results_steps = (collect(DataFrames.eachcol(steps)),DataFrames.names(steps)),
        )

    end

    cd(main)             # ritorno nella cartella di salvataggio dati


    return println("Saved data in xlsx")
end


function data_saving_4(InputParameters::InputParam,ResultsOpt_4::Results_4)

    @unpack (NYears, NMonths, NStages, Big, NHoursStep, bin) = InputParameters;       #NSteps,NHoursStage
    
   #@unpack (charge,discharge, soc, revenues_per_stage, x, y, z, w_xx, w_yy, w_zz, w_xy, w_xz, w_zy) rev_vendita, rev_acquisto = ResultsOpt;  
   @unpack (charge,discharge,rev,cap,e, soc,soc_quad,deg,deg_stage,gain_stage, cost_rev,rev_vendita, rev_acquisto,x, y, z,u, xy, xz, xu, yz, yu, zu, xyzu) = ResultsOpt_4;
   @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull,fix,cost ) = Battery ; 

    hour=string(now())
    a=replace(hour,':'=> '-')

    nameF= "Test 4 bin"
    nameFile="Summary " 

    folder = "$nameF"
    mkdir(folder)
    cd(folder)
    main=pwd()

    general = DataFrame()
    #general_1 = DataFrame()
    battery_costs= DataFrame()

    capacity=zeros(NStages)
    capacity_finale = zeros(NStages)

    for iStage=1:NStages
        capacity[iStage]=cap[(Steps_stages[iStage]+2)]
    end

    for iStage=2:NStages
        capacity_finale[iStage-1] = cap[Steps_stages[iStage]+1]
    end
        capacity_finale[end] = cap[end]
    
    general[!,"Stage"] = 1:1:NStages
    general[!,"Initial Capacity MWh"] = capacity[:]
    general[!,"Final Capacity MWh"] = capacity_finale[:]
    general[!,"Revamping MWh"] = rev[:]
    general[!, "Binary revamp"] = e[:]
    general[!,"Degradation"] = deg_stage[:]
    #general[!,"Net_Revenues"] = revenues_per_stage[:]
    general[!,"Gain charge/discharge"] = gain_stage[:]
    general[!,"Cost revamping"] = cost_rev[:]
    general[!, "Aux acquisto"] = rev_acquisto[:]
    general[!, "Aux vendita"] = rev_vendita[:]

    battery_costs[!,"Purchase costs €/MWh"] = Battery_price_purchase[1:NStages+1]    #Battery_price_purchase
    battery_costs[!,"Sale costs €/MWh"] = Battery_price_sale[1:NStages+1] 

    XLSX.writetable("$nameFile.xlsx", overwrite=true,                                       #$nameFile
    results_stages = (collect(DataFrames.eachcol(general)),DataFrames.names(general)),
    #results_1 = (collect(DataFrames.eachcol(general_1)),DataFrames.names(general_1)),
    costs = (collect(DataFrames.eachcol(battery_costs)),DataFrames.names(battery_costs)),
    )

    for iStage=1:NStages
        steps = DataFrame()

        steps[!,"Step"] = (Steps_stages[iStage]+1):(Steps_stages[iStage+1])
        steps[!, "Energy_prices €/MWh"] = Power_prices[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Energy capacity MWh"] = cap[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "SOC MWh"] = soc[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Charge MW"] = charge[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Discharge MW"] = discharge[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "SOC_quad MWh"] = soc_quad[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Deg MWh"] = deg[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "X"] = x[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Y"] = y[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Z"] = z[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "U"] = z[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "XY"] = xy[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "XZ"] = xz[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "XU"] = xu[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "YZ"] = yz[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "YU"] = yu[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "ZU"] = zu[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "XYZU"] = xyzu[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]

        XLSX.writetable("$iStage stage $a.xlsx", overwrite=true,                                       #$nameFile
        results_steps = (collect(DataFrames.eachcol(steps)),DataFrames.names(steps)),
        )

    end

    cd(main)             # ritorno nella cartella di salvataggio dati


    return println("Saved data in xlsx")
end


function data_saving_rolling(InputParameters::InputParam, ResultsEnd::ResultsEndLife)
    
    @unpack (NYears, NMonths, NStages, Big, NHoursStep, disc, Hours_rolling,Hours_saved) = InputParameters;       #NSteps,NHoursStage
    @unpack (Ch_end, Dis_end, soc_end,socq_end, deg_end, cap_end, new_vec) = ResultsEnd;
  
    general_r = DataFrame()
    general_r[!,"Energy Capacity MWh"] = cap_end[:]
    general_r[!,"Stored Energy MWh"] = soc_end[:]
    general_r[!,"Stored Energy Quad MWh"] = socq_end[:]

    general_s = DataFrame()
    #general_s[!,"Price €/MWh"] = new_vec[:]
    general_s[!,"Charge MW"] = Ch_end[:]
    general_s[!,"Discharge MW"] = Dis_end[:]
    general_s[!,"Degradation MWh"] = deg_end[:]

    general_p = DataFrame()
    general_p[!," Price €/MWh"] = new_vec[:]

    XLSX.writetable("Results Rolling con costo .xlsx", overwrite=true,                                       #$nameFile
        results_1 = (collect(DataFrames.eachcol(general_r)),DataFrames.names(general_r)),
        results_2 = (collect(DataFrames.eachcol(general_s)),DataFrames.names(general_s)),
        results_3 = (collect(DataFrames.eachcol(general_p)),DataFrames.names(general_p)),
        )

        return println("Saved data in xlsx")

end

function data_endlife(InputParameters::InputParam, ResEnd::ResultsCase1)

    @unpack (NYears, NMonths, NStages, Big, NHoursStep, disc, Hours_rolling,Hours_saved) = InputParameters;       #NSteps,NHoursStage
    @unpack (soc_e1, charge_e1, discharge_e1, deg_e1, cap_e1) = ResEnd;

    general_a = DataFrame()
    general_a[!,"Energy Capacity MWh"] = cap_e1[:]
    general_a[!,"Stored Energy MWh"] = soc_e1[:]

    general_b = DataFrame()
    #general_s[!,"Price €/MWh"] = new_vec[:]
    general_b[!,"Charge MW"] = charge_e1[:]
    general_b[!,"Discharge MW"] = discharge_e1[:]
    general_b[!,"Degradation MWh"] = deg_e1[:]

    general_c = DataFrame()
    general_c[!," Price €/MWh"] = new_vec[:]

    XLSX.writetable("Results Case 1 .xlsx", overwrite=true,                                       #$nameFile
        results_1 = (collect(DataFrames.eachcol(general_a)),DataFrames.names(general_a)),
        results_2 = (collect(DataFrames.eachcol(general_b)),DataFrames.names(general_b)),
        results_3 = (collect(DataFrames.eachcol(general_c)),DataFrames.names(general_c)),
        )

        return println("Saved data in xlsx Case 1")



end



function saving_new_deg(InputParameters::InputParam, NewDeg::ResultsNewDeg)
    
    @unpack (NYears, NMonths, NStages, Big, NHoursStep, disc, Hours_rolling,Hours_saved) = InputParameters;       #NSteps,NHoursStage
    @unpack (charge_new_deg, discharge_new_deg, state_new_deg, degradation_new_deg, capacity_new_deg, prices) = NewDeg;
  
    gen_r = DataFrame()
    gen_r[!,"Energy Capacity MWh"] = capacity_new_deg[:]
    gen_r[!,"Stored Energy MWh"] = state_new_deg[:]
    
    gen_s = DataFrame()
    #general_s[!,"Price €/MWh"] = new_vec[:]
    gen_s[!,"Charge MW"] = charge_new_deg[:]
    gen_s[!,"Discharge MW"] = discharge_new_deg[:]
    gen_s[!,"Degradation MWh"] = degradation_new_deg[:]

    gen_p = DataFrame()
    gen_p[!," Price €/MWh"] = prices[:]

    XLSX.writetable("Results Rolling con costo .xlsx", overwrite=true,                                       #$nameFile
        res_1 = (collect(DataFrames.eachcol(gen_r)),DataFrames.names(gen_r)),
        res_2 = (collect(DataFrames.eachcol(gen_s)),DataFrames.names(gen_s)),
        res_3 = (collect(DataFrames.eachcol(gen_p)),DataFrames.names(gen_p)),
        )

        return println("Saved data in xlsx Case 3")

end







