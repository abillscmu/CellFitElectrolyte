


VAH = "VAH02"

arr = readdir("data/cycle_individual_data/")

thisvec = zeros(625)


for a in arr
    if occursin(VAH,a)
        vec = split(a,['.','_'])
        local cycle = parse(Int,vec[2])
        ndf = CSV.read("data/cycle_individual_data/$(a)",DataFrame)
        filter!(row->row.Ns>=4,ndf)
        try
            thisvec[cycle+1] = ndf.QDischargemAh[end]
        catch
            @warn "problem on cycle $(cycle)"
        end
    end
end


figure(1)
clf()
plot(thisvec,"o",options=Dict("color"=>mycolors[1],"linewidth"=>3))
grid_on()
xlabel("Cycle")
ylabel("Voltage at End of Rest")
setgca(Dict("FontName"=>"Helvetica","FontSize"=>24))
