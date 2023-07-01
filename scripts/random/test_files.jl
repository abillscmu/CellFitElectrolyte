FOLDERNAME = "results/0302/"


files = readdir(FOLDERNAME)


good_files = []
bad_files = []
bad_results = []


@showprogress for file in files
    cell_to_load = "$(FOLDERNAME)$(file)"
    CELL = split(file, "_")[1]
    VAH = file[1:end-9]
    chain = try
        d = load(cell_to_load)
        chain = d["chain"]
        if chain[:lp][end] > 0
            push!(good_files, VAH)
        else
            push!(bad_results, VAH)
        end
    catch
        if isfile(cell_to_load)
            push!(bad_files, VAH)
        end
        continue
    end
end

io = open("results/good.txt", "w")
for f in good_files
    println(io, f)
end

close(io)
