using DataFrames, CSV, PythonPlot

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

@save "files_to_use.jld2" files_to_use

