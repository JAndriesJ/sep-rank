module batch
srcDir  = dirname(@__FILE__)*"\\";

include(srcDir*"Model\\R_sep_Model.jl")
include(srcDir*"Model\\C_sep_Model.jl")
include(srcDir*"Model\\Utils_Model.jl")
include(srcDir*"sep_Compute.jl")
using .R_sep_Model
using .C_sep_Model
using .Utils_Model
using .sep_Compute


using JuMP
using CSV
using DataFrames


export batch_model
default_sdir = dirname(dirname(@__FILE__))*"\\assets\\bounds\\"
""" Produces models for a set of constraints."""
function batch_model(t,ρdic,df;sdir = default_sdir,con_list= ["S1G","S2G","S3G"]) # ,
    for ex in keys(ρdic)
        ρ = ρdic[ex]
        d = filter(:Example => ==(ex),df).Size[1]
        for con ∈  con_list, isR ∈ ['R','C'] # , 
            (!filter(:Example => ==(ex),df).isReal[1] && (isR == 'R')) ? continue : 0
            println("Currently $isR-modeling example $ex with constraints $con")
            sep_mod = (isR == 'R') ? R_sep_Model.Modelξₜˢᵉᵖ(ρ,d,t;con_list=con) : C_sep_Model.Modelξₜˢᵉᵖ(ρ,d,t;con_list=con,noBlock = false)

            Utils_Model.export_model(sep_mod, sdir*"$(ex)_$(con)_$(isR)_$t.dat-s")
        end
    end
end

""" Reads all .dat-s files and computes bounds storing result summary in a .csv file"""
function batch_Computeξₜˢᵉᵖ(boundsDir = default_sdir)
    datsFiles = [file for file in readdir(boundsDir,join = true) if contains(file,".dat-s")]
    file_loc = boundsDir*"Summary.csv"
    touch(file_loc)
    open(file_loc,"a") do io
        write(io, "|Example|Model|Constraint|Primal|Dual|obj_val|\n")
    end
    for file in datsFiles
        name = split(basename(file),".")[1]
        ex   = split(name,"_")[1]
        con  = split(name,"_")[2]
        num  = split(name,"_")[3]
        # try
        mod  = Utils_Model.read_model(file)
        omod = sep_Compute.Computeξₜˢᵉᵖ(mod)
        # catch
        #     open(file_loc,"a") do io; write(io, "|$ex|$con|err.|err.|err.| \n"); end
        #     continue
        # end
        pstat   = JuMP.primal_status(omod)
        dstat   = JuMP.dual_status(omod)
        ov_temp = round(JuMP.objective_value(omod),digits=3)
        ov = (string(pstat) == string(dstat) == "FEASIBLE_POINT") ?
             ov_temp :
             (contains(string(pstat)*string(dstat),"CER") ? "*" : NaN )

        open(file_loc,"a") do io
            write(io, "|$ex|$num|$con|$pstat|$dstat|$ov|\n")
        end
    end
end

"""Postprocesses the output of batch_Computeξₜˢᵉᵖ into a nices tabular structure"""
function unstack_constraints(df,bDir = default_sdir)
    temp_df  = CSV.read("C:\\Users\\andries\\all-my-codes\\sep-rank\\assets\\bounds\\"*"Summary.csv", DataFrame)
    temp2_df = select(innerjoin(df, temp_df, on = :Example),[:Example,:isReal,:isSeparable,:Model,:Size,:Bi_rank,:Constraint,:obj_val])
    temp4_df = unstack(temp2_df, :Model, :obj_val)

    temp4R_df = select(temp4_df,[:Example,:isReal,:isSeparable,:Size,:Bi_rank,:Constraint,:R])
    temp4Rus_df = unstack(temp4R_df, :Constraint, :R)
    rename!(temp4Rus_df,Dict(:S1G=>:RS1G, :S2G => :RS2G, :S3G => :RS3G));

    temp4C_df = select(temp4_df,[:Example,:isReal,:isSeparable,:Size,:Bi_rank,:Constraint,:C])
    temp4Cus_df = unstack(temp4C_df, :Constraint, :C)
    rename!(temp4Cus_df,Dict(:S1G=>:CS1G, :S2G => :CS2G, :S3G => :CS3G));

    temp5_df = innerjoin(temp4Cus_df,select(temp4Rus_df,[:Example ,:RS1G, :RS2G, :RS3G]), on = :Example)

    CSV.write(bDir*"SummaryUnstacked.csv", temp5_df, delim="&")
    final_df = DataFrames.select(temp5_df,[:Example,:Size,:Bi_rank,:CS1G, :CS2G, :CS3G, :RS1G, :RS2G, :RS3G, :isSeparable])
    rename!(final_df,Dict(:Example=>:ρ, :Size => :d₁d₂));
    CSV.write(bDir*"t=__.csv", final_df, delim="&")
    return final_df
end

end
