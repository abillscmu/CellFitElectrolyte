using LinearAlgebra,Statistics, CSV, DataFrames, MATLABPlots, StatsBase

function movingaverage(vec,M)
    N = length(vec)
    newvec = zeros(N-M)
    side = Int(M/2)
    for m in 1:N-M
    newvec[m] = mean(vec[m:m+2*side])
    end
    return newvec
end

function movingstdev(vec,M)
    N = length(vec)
    newvec = zeros(N-M)
    side = Int(M/2)
    for m in 1:N-M
    newvec[m] = std(vec[m:m+2*side])
    end
    return newvec
end

function movingmax(vec,M)
    N = length(vec)
    newvec = zeros(N-M)
    x = similar(newvec)
    side = Int(M/2)
    for m in 1:N-M
    newvec[m] = maximum(vec[m:m+2*side])
    x[m] = mean(collect(m:m+2*side))
    end
    return newvec,x
end

function movingmin(vec,M)
    N = length(vec)
    newvec = zeros(N-M)
    side = Int(M/2)
    for m in 1:N-M
    newvec[m] = minimum(vec[m:m+2*side])
    end
    return newvec
end

df = CSV.read("results/VAH02_1204.csv",DataFrame)
cell = 2
filter!(row->row.cell==cell,df)
sort!(df,"cycle")

bottom = 2
top = 300

s = "εₛ⁻"

color1 = [230,159,0]./255
color2 = [86,180,233]./255
color3 = [0, 158, 115]./255
color7 = [240, 228, 66]./255
color5 = [0, 114, 178]./255
color6 = [213, 94, 0]./255
color4 = [204, 121, 167]./255
color8 = [0, 0, 0]./255

mycolors = [color1,color2,color3,color4,color5,color6,color7,color8]
i=5
averaging=20

param = "Root Mean Squared Error[V]"
VAHNUM = lpad(cell,2,"0")
thistitle = "VAH$VAHNUM"

figure(1)
clf()

subplot(2,1,1)
scatter(collect(df.cycle[2:end]),collect(df[2:end,s]))
hold_on()
title(thistitle)
R_0_movingaverage = movingaverage(df[2:end,s],averaging)
R_0_95CI = 2 .*movingstdev(df[2:end,s],averaging)
R_0_MAX,x = movingmax(df[2:end,s],averaging)
R_0_MIN = movingmin(df[2:end,s],averaging)
l = length(R_0_movingaverage)
lower = averaging/2
#x = collect(range(start=lower,length=l,step=1))
#plot(x,R_0_movingaverage,options=Dict("LineWidth"=>3,"Color"=>mycolors[i]))
ylabel("RMSE [V]")
xlabel("Cycle Number")
box_on()
grid_on()
setgca(Dict("FontSize"=>16,"FontName"=>"Helvetica"))
#plot(x,R_0_movingaverage.-R_0_95CI,options=Dict("LineWidth"=>1,"Color"=>mycolors[2]))
#plot(x,R_0_movingaverage.+R_0_95CI,options=Dict("LineWidth"=>1,"Color"=>mycolors[2]))
#stairs(x,R_0_MAX,options=Dict("LineWidth"=>1,"Color"=>mycolors[3]))
#stairs(x,R_0_MIN,options=Dict("LineWidth"=>1,"Color"=>mycolors[3]))

ac = autocor(df[!,s])
subplot(2,1,2)
scatter(collect(range(0,stop=length(ac)-1)),ac,125 .*ones(length(ac)),"filled")
hold_on()
plot(collect(range(0,stop=length(ac)+1)), (2.807/sqrt(length(df[!,s]))).*ones(length(ac)+2),"k",options=Dict("LineWidth"=>3))
plot(collect(range(0,stop=length(ac)+1)), -(2.807/sqrt(length(df[!,s]))).*ones(length(ac)+2),"k",options=Dict("LineWidth"=>3))
setgca(Dict("FontSize"=>16,"FontName"=>"Helvetica"))
ylabel("Autocorrelation")
xlabel("Lag")
grid_on()
box_on()

#saveas("$(thistitle)_autocor.eps",formattype="epsc")





