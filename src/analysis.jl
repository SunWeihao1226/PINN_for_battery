using Plots

include("ode_spm.jl")
@time include("pinn_spm.jl")

# Setting plot parameters
plt_size = (800,600)
x_ticks = ts
legends = [:bottomright, :topright]

# calculating errors w.r.t. dt
err_c_sp = []
err_c_sn = []
for i in (1:iter)
    append!(err_c_sp, abs.(pred_c_sp_dt[i] - sol_c_sp_dt[i]))
    append!(err_c_sn, abs.(pred_c_sn_dt[i] - sol_c_sn_dt[i]))
end 


# Function to create the plot
function find_max(data)
    max = 0
    for i in (1:iter)
        if (maximum(data[i]) >= max)
            max = maximum(data[i])
        end
    end
    return max
end

function find_min(data)
    min = minimum(data[1])
    for i in (1:iter)
        if (minimum(data[i]) <= min)
            min = minimum(data[i])
        end
    end
    return min
end

function drawPlt(data, title)
    upper = find_max(data)
    lower = find_min(data)
    temp_plt = plot(size=plt_size, title = title, xlabel = "dt", xticks=x_ticks, ylabel = "Concentration")
    ani =  @animate for i in (1:iter)
        lab = string("dr= ", string((i-1)/(iter-1)))
        plot(x_axis, data[i], size=plt_size, title = title, xlabel = "dt", ylabel = "Concentration", ylims=(lower,upper), xticks=x_ticks, legends=false)
    end
    return ani
end


# Create plots for PINN predictions
# Posivive and Negative Concentrations vs dt 

plt_c_sp = drawPlt(pred_c_sp_dt, "Positive Concentration by PINN vs. dt")
plt_c_sn = drawPlt(pred_c_sn_dt, "Negative Concentration by PINN vs. dt")
gif(plt_c_sp, "../plots/c_sp_predict.gif", fps=5)
gif(plt_c_sn, "../plots/c_sn_predict.gif", fps=5)

# Create Plots for MethidOfLine predictions

plt_c_sp_mol = drawPlt(sol_c_sp_dt, "Positive Concentration by MoL vs. dt")
plt_c_sn_mol = drawPlt(sol_c_sn_dt, "Negative Concentration by MoL vs. dt")
gif(plt_c_sp_mol, "../plots/c_sp_sol.gif", fps=5)
gif(plt_c_sn_mol, "../plots/c_sn_sol.gif", fps=5)

# Create error heat plots

err_csp_plt = plot(ts, rs, err_c_sp, size=(800,800), linetype=:contourf, title = "c_sp Error", xlabel = "dt", ylabel = "dr", xticks=x_ticks, yticks = x_ticks)
err_csn_plt = plot(ts, rs, err_c_sn, size=(800,800), linetype=:contourf, title = "c_sn Error", xlabel = "dt", ylabel = "dr", xticks=x_ticks, yticks = x_ticks)
savefig(err_csp_plt,"../plots/c_sp_error.png")
savefig(err_csn_plt,"../plots/c_sn_error.png")


# Save plots
# savefig(plt_c_sp,"../plots/c_sp_predict.")
# savefig(plt_c_sn,"../plots/c_sn_predict.png")
# savefig(plt_c_sp_mol,"../plots/c_sp_sol.png")
# savefig(plt_c_sn_mol,"../plots/c_sn_sol.png")
# savefig(err_csp_plt,"../plots/c_sp_error.png")
# savefig(err_csn_plt,"../plots/c_sn_error.png")