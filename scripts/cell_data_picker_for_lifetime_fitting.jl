using DataFrames, CSV, PythonPlot, Turing, JLD2, ProgressMeter

files = readdir("results/0302/")
x = map( s -> s[1:5], files) 
cells = unique(x)
files_per_cell = 10

files_to_use = Dict()


for cell in cells
    these_cells = filter(s-> s[1:5] == cell, files)
    these_cycs = map(s -> parse(Int, split(s, '_')[2]), these_cells)
    maximum_cell = maximum(these_cycs)
    minimum_cell = minimum(these_cycs)
    if minimum_cell == 0
        if 2 in these_cycs
            minimum_cell = 2
        else
            error("shit")
        end
    end
    arr =Int.(round.(collect(range(minimum_cell, stop=maximum_cell, length=files_per_cell))))
    successarr = Array{Bool,1}(undef, files_per_cell)
    for (i,a) in enumerate(arr)
        df = CSV.read("data/cycle_individual_data/$(cell)_$a.csv", DataFrame)
        ttest = ((df.times[end] - df.times[1]) < 20000.0)
        ctest = try
            d = load("results/0302/$(cell)_$(attempt)_HMC.jld2")
            chain = d["chain"]
            ctest = length(chain[:ω].data[:,1]) == 1000
            ctest
        catch
            false
        end
        ttest = ctest & ttest
        successarr[i] = false
        if ((a in these_cycs) & ttest)
            continue
            successarr[i] = true
        else
            evencount=1
            oddcount=1
            totalcount=0
            notdone = true
            while notdone
                if iseven(totalcount)
                    attempt = a + evencount
                    if attempt in these_cycs
                        ttest=try
                            df = CSV.read("data/cycle_individual_data/$(cell)_$attempt.csv", DataFrame)
                            ttest = ((df.times[end] - df.times[1]) < 20000.0)
                            ttest
                        catch
                            false
                        end
                        if ttest
                            ctest = try
                                d = load("results/0302/$(cell)_$(attempt)_HMC.jld2")
                                chain = d["chain"]
                                ctest = length(chain[:ω].data[:,1]) == 1000
                                ctest
                            catch
                                false
                            end
                            if ctest
                                notdone=false
                                successarr[i] = true
                                arr[i] = attempt
                            end
                        end
                    end
                    evencount+=1
                else
                    attempt = a - oddcount
                    if attempt in these_cycs
                        ttest=try
                            df = CSV.read("data/cycle_individual_data/$(cell)_$attempt.csv", DataFrame)
                            ttest = ((df.times[end] - df.times[1]) < 20000.0)
                            ttest
                        catch
                            false
                        end
                        if ttest
                            ctest = try
                                d = load("results/0302/$(cell)_$(attempt)_HMC.jld2")
                                chain = d["chain"]
                                ctest = length(chain[:ω].data[:,1]) == 1000
                                ctest
                            catch
                                false
                            end
                            if ctest
                                notdone=false
                                successarr[i] = true
                                arr[i] = attempt
                            end
                        end
                    end
                    oddcount+=1
                end
                totalcount+=1
                if totalcount > 50
                    successarr[i] = false
                    @warn "no suitable match found for $cell $a"
                    notdone = false
                end
            end
        end
    end
    for (i,a) in enumerate(arr)
        ## SPECIAL CASES 
        if (cell == "VAH05") & (a == 1547)
            arr[i] = 1546
            a = 1546
        end
        if (cell == "VAH07") & (a == 275)
            arr[i] = 273
            a = 273
        end
        if (cell == "VAH13") & (a == 810)
            arr[i] = 809
            a = 809
        end
        if (cell == "VAH16") & (a == 494)
            arr[i] = 485
            a = 485
        end
        if (cell == "VAH25") & (a == 432)
            arr[i] = 430
            a = 430
        end
        if (cell == "VAH25") & (a == 493)
            arr[i] = 492
            a = 492
        end
        if (cell == "VAH26") & (a == 1161)
            arr[i] = 1158
            a = 1158
        end
        if !(successarr[i])
            continue
        end

        df = CSV.read("data/cycle_individual_data/$(cell)_$a.csv", DataFrame)
        fig,ax=subplots()
        ax.plot(df.times,df.EcellV)
        fig.savefig("cycleplots/$(cell)_$a.png")
        @info "$cell $a $(df.times[end] - df.times[1])"
        

    end
    arr = unique(arr)
    if length(arr) != 10
        error("shid")
    end
    

    



    files_to_use[cell] = arr
end

params = (:ω, :εₑ⁺, :εₑ⁻, :frac_sol_am_neg, :frac_sol_am_pos)
@showprogress for cell in keys(files_to_use)
    if cell in ["VAH01", "VAH02", "VAH05", "VAH06", "VAH07", "VAH09", "VAH10", "VAH11", "VAH12"]
        continue
    end
    for p in params
        fig, ax = subplots()
        for c in files_to_use[cell]
            d = load("results/0302/$(cell)_$(c)_HMC.jld2")
            chain = d["chain"]
            datas = chain[p].data[:,1]
            ax.scatter(c*ones(1000), datas)
        end
        fig.savefig("chain_plots/$(cell)_$(p).pdf", bbox_inches = "tight")
    end
end



@save "files_to_use.jld2" files_to_use

