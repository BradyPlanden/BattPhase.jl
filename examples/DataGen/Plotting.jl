using CSV, DataFrames, Plots
function ImportData(PredLoc, GraphLoc)
    df_pred = CSV.read(PredLoc, DataFrame)
    df_graph = CSV.read(GraphLoc, DataFrame)
    return df_pred, df_graph
end


df_pred, df_graph = ImportData("92183031418602_DataOut_Pred_010","92183031418602_DataOut_Graph_010")
NN = MM = 43
νg = Array{Float64}(undef,NN,MM,size(df_graph,2))
νp = copy(νg)
for i ∈ 1:size(df_graph,2)
    νg[:,:,i] = reshape(df_graph[:,i],(MM,NN))
    νp[:,:,i] = reshape(df_pred[:,i],(MM,NN))
end

νt = copy(νg)
νi = copy(νp)

#1D to 2D
for k ∈ 1:size(νg,3)
    for j ∈ 1:size(νt,2)
        if isodd(j) == true
            νt[:,j,k] .= reverse(νt[:,j,k])
            νi[:,j,k] .= reverse(νi[:,j,k])
        else
            νt[:,j,k] .= νt[:,j,k]
            νi[:,j,k] .= νi[:,j,k]
        end
    end
end
plotly()
l = @layout [a b c{0.33w}; d e f{0.33w}]
p1 = plot(heatmap(νi[:,:,1],clims=(0, 1), color = :davos), colorbar=false)
p2 = plot(heatmap(νi[:,:,10],clims=(0, 1), color = :davos), colorbar=false)
p3 = plot(heatmap(νi[:,:,20],clims=(0, 1), color = :davos))

p4 = plot(heatmap(νt[:,:,1],clims=(0, 1), color = :davos), xlabel = "Timestep 25")
p5 = plot(heatmap(νt[:,:,10],clims=(0, 1), color = :davos), xlabel = "Timestep 35" )
p6 = plot(heatmap(νt[:,:,20],clims=(0, 1), color = :davos), xlabel = "Timestep 45")
plot(p1,p2,p3,p4,p5,p6, color = :davos, layout=l, size=(1280,720))

