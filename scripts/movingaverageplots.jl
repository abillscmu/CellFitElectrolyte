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

df = CSV.read("data/outputs/outputs0728_interpolative.csv",DataFrame)
sort!(df,"cycle")

bottom = 2
top = 300

s = "θₛ⁺"

color1 = [230,159,0]./255
color2 = [86,180,233]./255
color3 = [0, 158, 115]./255
color7 = [240, 228, 66]./255
color5 = [0, 114, 178]./255
color6 = [213, 94, 0]./255
color4 = [204, 121, 167]./255
color8 = [0, 0, 0]./255

mycolors = [color1,color2,color3,color4,color5,color6,color7,color8]
i=1
averaging=20

figure(1)
clf()
subplot(2,1,1)
scatter(collect(df.cycle[2:end]),collect(df[2:end,s]),"\'filled\'",options=Dict("MarkerFaceAlpha"=>0.4,"MarkerFaceColor"=>mycolors[i],"MarkerEdgeColor"=>mycolors[i],"MarkerEdgeAlpha"=>0,"HandleVisibility"=>"off"))
hold_on()
R_0_movingaverage = movingaverage(df[2:end,s],averaging)
R_0_95CI = 2 .*movingstdev(df[2:end,s],averaging)
R_0_MAX,x = movingmax(df[2:end,s],averaging)
R_0_MIN = movingmin(df[2:end,s],averaging)
l = length(R_0_movingaverage)
lower = averaging/2
#x = collect(range(start=lower,length=l,step=1))
plot(x,R_0_movingaverage,options=Dict("LineWidth"=>3,"Color"=>mycolors[i]))
#plot(x,R_0_movingaverage.-R_0_95CI,options=Dict("LineWidth"=>1,"Color"=>mycolors[2]))
#plot(x,R_0_movingaverage.+R_0_95CI,options=Dict("LineWidth"=>1,"Color"=>mycolors[2]))
stairs(x,R_0_MAX,options=Dict("LineWidth"=>1,"Color"=>mycolors[3]))
#stairs(x,R_0_MIN,options=Dict("LineWidth"=>1,"Color"=>mycolors[3]))
grid_on()
setgca(Dict("FontSize"=>24,"FontName"=>"Helvetica"))
subplot(2,1,2)
histogram(df[2:end,s],options=Dict("NumBins"=>50,"FaceColor"=>mycolors[i],"EdgeColor"=>mycolors[i],"FaceAlpha"=>0.5,"EdgeAlpha"=>0.5))
setgca(Dict("FontSize"=>24,"FontName"=>"Helvetica"))


