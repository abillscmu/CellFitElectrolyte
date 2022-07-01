#=
annealOptions.jl
- Author: Leif Fredericks
- Date: 2019-03-13
- Part of the simulated annealing adapted from the matlab code written by Joachim Vandekerckhove
https://www.mathworks.com/matlabcentral/fileexchange/10548-general-simulated-annealing-algorithm
=#
# These are the default values in simAnneal.jl, matching the defaults of the source material
#
#     InitTemp ::Float64    = 1.
#     MaxConsRej ::Int64    = 1000
#     MaxSuccess::Int64     = 20
#     MaxTries::Int64       = 300
#     StopTemp::Float64     = 1e-8
#     StopVal::Float64      = -Inf
#     Verbosity::Int64      = 1
# end
struct annealOptions
    InitTemp ::Float64  # Initial annealing temperature
    MaxConsRej ::Int64  # Maximum number of rejected new solutions without an acceptance
    MaxSuccess::Int64   # Maximum number of acceptances at a given temperature: length of Markov chain
    MaxTries::Int64     # Maximum total new solutions that can be tried, acceptance or rejection, at a given temperature
    StopTemp::Float64   # Temperature that indicates the end of simAnneal when reached
    StopVal::Float64    # Value of loss() objective cost function that is sufficient to stop fitting
    Verbosity::Int64    # Reporting to command line: (0) Suppress output, (1) Final report,
                        #                            (2) Report at each new temeperature or when new best is found
end
####################################################################################################
# Function call with minimum, a loss function and parent array
####################################################################################################
function anneal(loss::Function, parent,λ::Float64)
    # Initialize Upper Bound to Infinity
    ub = Array{Float64, 1}(undef, length(parent))
            fill!(ub,Inf)
    # Initialize Lower Bound to negative Infintity
    lb = Array{Float64, 1}(undef, length(parent))
            fill!(lb,-Inf)
    # Search function for new set of parameters
    function newsol(input, ub, lb)
        # (x.+(randperm(length(x))==length(x))*randn()/100)
        # ie the one from the original code, pick a random parameter and add random normal /100 to it
        Lp      =   length(input)                   # Length of parameter vector
        output  =   Array{Float64, 1}(undef, Lp)    # Initialize output vector
        idx     =   rand(1:Lp)                      # Pick a random parameter index
        # Fill each index of output array
        for i in 1:Lp
            # Change the value of the previously picked random index
            if i == idx
                # Suggest output value by adding random normal scaled by 1/100
                proposedNew = input[i]+(randn().*λ)
                # Check to keep proposal within bounds
                if proposedNew < lb[i]
                    # Set output to lower bound if proposed is too low
                    output[i] = lb[i]
                elseif proposedNew > ub[i]
                    # Set output to upper bound if proposed is too high
                    output[i] = ub[i]
                else
                    # Set output to proposal if within the bounds
                    output[i] = proposedNew
                end
            else
                # Set output to input for all indices besides the randomly selected one
                output[i] = input[i]
            end
        end

        # Calculate loss of new parameters
        newenergy   =   loss(output)

        return output, newenergy
    end
    # Reduce the annealing temperature
    function cool(Tin::Float64)
        # Set new temperature to 80% old temperature
        return 0.8*Tin
    end
    # Set type to annealing process from source material
    type = "Kirkpatrick"
    # Set options to those used in source material
    options::annealOptions = annealOptions( 1., 1000, 2000, 3000, 1e-8, -Inf, 2)
    # Run Main Function with all inputs now set

    return anneal(loss, parent, ub, lb, newsol, cool, type, options)
end


function anneal_step(loss::Function, parent, λ; options::annealOptions = annealOptions( 1., 1000, 2000, 3000, 1e-8, -Inf, 1))
    # println(options)
    # Initialize Upper Bound to Infinity
    ub = Array{Float64, 1}(undef, length(parent))
            fill!(ub,Inf)
    # Initialize Lower Bound to negative Infintity
    lb = Array{Float64, 1}(undef, length(parent))
            fill!(lb,-Inf)
    # Search function for new set of parameters
    function newsol(input, ub, lb)
        # (x.+(randperm(length(x))==length(x))*randn()/100)
        # ie the one from the original code, pick a random parameter and add random normal /100 to it
        Lp      =   length(input)                   # Length of parameter vector
        output  =   Array{Float64, 1}(undef, Lp)    # Initialize output vector
        idx     =   rand(1:Lp)                      # Pick a random parameter index
        # Fill each index of output array
        for i in 1:Lp
            # Change the value of the previously picked random index
            if i == idx
                # Suggest output value by adding random normal scaled by 1/100
                proposedNew = input[i]+(randn().*λ[i])
                # Check to keep proposal within bounds
                if proposedNew < lb[i]
                    # Set output to lower bound if proposed is too low
                    output[i] = lb[i]
                elseif proposedNew > ub[i]
                    # Set output to upper bound if proposed is too high
                    output[i] = ub[i]
                else
                    # Set output to proposal if within the bounds
                    output[i] = proposedNew
                end
            else
                # Set output to input for all indices besides the randomly selected one
                output[i] = input[i]
            end
        end

        # Calculate loss of new parameters
        newenergy   =   loss(output)

        return output, newenergy
    end
    # Reduce the annealing temperature
    function cool(Tin::Float64)
        # Set new temperature to 80% old temperature
        return 0.8*Tin
    end
    # Set type to annealing process from source material
    type = "Kirkpatrick"
    # Set options to those used in source material
    # options::annealOptions = annealOptions( 1., 1000, 2000, 3000, 1e-8, -Inf, 1)
    # Run Main Function with all inputs now set

    return anneal(loss, parent, ub, lb, newsol, cool, type, options)
end
####################################################################################################
# Function call with upper and lower bounds and logarithmic 
####################################################################################################
function anneal_log(loss::Function, parent, ub, lb; λ = 1/100., options::annealOptions = annealOptions( 1., 1000, 2000, 3000, 1e-8, -Inf, 1))
    # Search function for new set of parameters
    function newsol(input, ub, lb)
        # (x.+(randperm(length(x))==length(x))*randn()/100)
        # ie the one from the original code, pick a random parameter and add random normal /100 to it
        Lp      =   length(input)                   # Length of parameter vector
        output  =   Array{Float64, 1}(undef, Lp)    # Initialize output vector
        idx     =   rand(1:Lp)                      # Pick a random parameter index
        # Fill each index of output array
        for i in 1:Lp
            # Change the value of the previously picked random index
            if i == idx
                if ub[i] > 0.0 && lb[i] > 0.0 && input[i] > 0.0
                    # Suggest output value by adding random normal scaled 
                    logRange    =   log10(ub[i] / lb[i]) 
                    logAddend   =   logRange * randn() * λ
                    exponential =   log10(input[i]) + logAddend
                    proposedNew =   10^(exponential)

                elseif ub[i] < 0.0 && lb[i] < 0.0 && input[i] < 0.0
                    logRange    =   log10(ub[i] / lb[i]) 
                    logAddend   =   logRange * randn() * λ
                    exponential =   log10(-input[i]) + logAddend
                    proposedNew =   -10^(exponential)

                else # If the range spans positive and negative (or includes zero), logarithmic doesn't make sense
                    proposedNew = input[i]+(ub[i].-lb[i]).*randn() * λ
                end


                #println(i)

                #println(proposedNew)
                # Check to keep proposal within bounds
                if proposedNew < lb[i]
                    # Set output to lower bound if proposed is too low
                    output[i] = lb[i]

                elseif proposedNew > ub[i]
                    # Set output to upper bound if proposed is too high
                    output[i] = ub[i]
                else
                    # Set output to proposal if within the bounds
                    output[i] = proposedNew
                end
            else
                # Set output to input for all indices besides the randomly selected one
                output[i] = input[i]
            end
        end

        # Calculate loss of new parameters
        newenergy   =   loss(output)

        return output, newenergy
    end
    # Reduce the annealing temperature
    function cool(Tin::Float64)
        # Set new temperature to 80% old temperature
        return 0.8*Tin
    end
    # Set type to annealing process from source material
    type = "Kirkpatrick"
    # Set options to those used in source material
    # options::annealOptions = annealOptions( 1., 1000, 20, 300, 1e-8, -Inf, 1)
    # Run Main Function with all inputs now set
    return anneal(loss, parent, ub, lb, newsol, cool, type, options)
end
####################################################################################################
# Function call above additionally with upper and lower bounds
####################################################################################################
function anneal(loss::Function, parent, ub, lb; λ::Float64=1/100, options::annealOptions = annealOptions( 1., 1000, 2000, 3000, 1e-8, -Inf, 2))
    # Search function for new set of parameters
    function newsol(input, ub, lb)
        # (x.+(randperm(length(x))==length(x))*randn()/100)
        # ie the one from the original code, pick a random parameter and add random normal /100 to it
        Lp      =   length(input)                   # Length of parameter vector
        output  =   Array{Float64, 1}(undef, Lp)    # Initialize output vector
        idx     =   rand(1:Lp)                      # Pick a random parameter index
        # Fill each index of output array
        for i in 1:Lp
            # Change the value of the previously picked random index
            if i == idx
                # Suggest output value by adding random normal scaled by 1/100
                proposedNew = input[i]+(ub[i].-lb[i]).*randn()*λ
                #println(i)

                #println(proposedNew)
                # Check to keep proposal within bounds
                if proposedNew < lb[i]
                    # Set output to lower bound if proposed is too low
                    output[i] = lb[i]

                elseif proposedNew > ub[i]
                    # Set output to upper bound if proposed is too high
                    output[i] = ub[i]
                else
                    # Set output to proposal if within the bounds
                    output[i] = proposedNew
                end
            else
                # Set output to input for all indices besides the randomly selected one
                output[i] = input[i]
            end
        end

        # Calculate loss of new parameters
        newenergy   =   loss(output)

        return output, newenergy
    end
    # Reduce the annealing temperature
    function cool(Tin::Float64)
        # Set new temperature to 80% old temperature
        return 0.8*Tin
    end
    # Set type to annealing process from source material
    type = "Kirkpatrick"
    # Set options to those used in source material
    # options::annealOptions = annealOptions( 1., 1000, 20, 300, 1e-8, -Inf, 1)
    # Run Main Function with all inputs now set
    return anneal(loss, parent, ub, lb, newsol, cool, type, options)
end
####################################################################################################
# Function call above additionally with a custom new solution function
####################################################################################################
function anneal(loss::Function, parent, ub, lb, newsol::Function)
    # Reduce the annealing temperature
    function cool(Tin::Float64)
        # Set new temperature to 80% old temperature
        return 0.8*Tin
    end
    # Set type to annealing process from source material
    type = "Kirkpatrick"
    # Set options to those used in source material
    options::annealOptions = annealOptions( 1., 1000, 20, 300, 1e-8, -Inf, 2)
    # Run Main Function with all inputs now set
    return anneal(loss, parent, ub, lb, newsol, cool, type, options)
end
####################################################################################################
# Function call above additionally with a custom cooling
####################################################################################################
function anneal(loss::Function, parent, ub, lb, newsol::Function,
    cool::Function)
    # Set type to annealing process from source material
    type = "Kirkpatrick"
    # Set options to those used in source material
    options::annealOptions = annealOptions( 1., 1000, 20, 300, 1e-8, -Inf, 2)
    # Run Main Function with all inputs now set
    return anneal(loss, parent, ub, lb, newsol, cool, type, options)
end
####################################################################################################
# Function call above additionally with a specified variation
####################################################################################################
function anneal(loss::Function, parent, ub, lb, newsol::Function,
    cool::Function, type::String)
    # Set options to those used in source material
    options::annealOptions = annealOptions( 1., 1000, 20, 300, 1e-8, -Inf, 2)
    # Run Main Function with all inputs now set
    return anneal(loss, parent, ub, lb, newsol, cool, type, options)
end
####################################################################################################
# Main Function: call above additionally with custom options
####################################################################################################
function anneal(loss::Function, parent, ub, lb, newsol::Function,
    cool::Function, type::String, options)

    ####################################################################################################
    # Description from source material
    ####################################################################################################
    # ANNEAL  Minimizes a function with the method of simulated annealing
    # (Kirkpatrick et al., 1983)
    #
    #  ANNEAL takes three input parameters, in this order:
    #
    #  LOSS is a function handle [anonymous function or inline()] with a loss()
    #  function(), which may be of any type(), and needn't be continuous. It does
    #  however, need to return a single value.
    #
    #  PARENT is a vector with initial guess parameters. You must input an
    #  initial guess.
    #
    #  OPTIONS is a structure with settings for the simulated annealing. If no
    #  OPTIONS structure is provided, ANNEAL uses a default structure. OPTIONS
    #  can contain any or all of the following fields (missing fields are
    #  filled with default values):
    #
    #       Verbosity: Controls output to the screen.
    #                  0 suppresses all output
    #                  1 gives final report only [default]
    #                  2 gives temperature changes and final report
    #       Generator: Generates a new solution from an old one.
    #                  Any function handle that takes a solution as input and
    #                  gives a valid solution (i.e. some point in the solution
    #                  space) as output.
    #                  The default function generates a row vector which()
    #                  slightly differs from the input vector in one element:
    #                  @(x) (x+(randperm(length(x))==length(x))*randn()/100)
    #                  Other examples of possible solution generators:
    #                  @(x) (rand(3,1)): Picks a random point in the unit cube
    #                  @(x) (ceil([9 5].*rand(2,1))): Picks a point in a 9-by-5
    #                                                 discrete grid()
    #                  Note that if you use the default generator, ANNEAL only
    #                  works on row vectors. For loss functions that operate on
    #                  column vectors, use this generator instead of the
    #                  default:
    #                  @(x) (x[:]"+(randperm(length(x))==length(x))*randn()/100)"
    #        InitTemp: The initial temperature, can be any positive number.
    #                  Default is 1.
    #        StopTemp: Temperature at which to stop(), can be any positive number
    #                  smaller than InitTemp.
    #                  Default is 1e-8.
    #         StopVal: Value at which to stop immediately, can be any output of
    #                  LOSS that is sufficiently low for you.
    #                  Default is() -Inf().
    #       CoolSched: Generates a new temperature from the previous one.
    #                  Any function handle that takes a scalar as input and
    #                  returns a smaller but positive scalar as output.
    #                  Default is() @(T) (.8*T)
    #      MaxConsRej: Maximum number of consecutive rejections, can be any()
    #                  positive number.
    #                  Default is 1000.
    #        MaxTries: Maximum number of tries within one temperature, can be
    #                  any positive number.
    #                  Default is 300.
    #      MaxSuccess: Maximum number of successes within one temperature, can
    #                  be any positive number.
    #                  Default is 20.
    #
    #
    #  Usage:
    #     [MINIMUM,FVAL] = ANNEAL[LOSS,NEWSOL,[OPTIONS]]
    #          MINIMUM is the solution which generated the smallest encountered
    #          value when input into LOSS.
    #          FVAL is the value of the LOSS function evaluated at MINIMUM.
    #     OPTIONS = ANNEAL[]
    #          OPTIONS is the default options structure.
    #
    #
    #  Example:
    #     The so-called "six-hump camelback" function has several local minima
    #     in the range -3<=x<=3 and -2<=y<=2. It has two global minima, namely
    #     f[-0.0898,0.7126] = f[0.0898,-0.7126] = -1.0316. We can define and
    #     minimise it as follows:
    #          camel = (x,y) ->(4-2.1*x.^2+x.^4/3).*x.^2+x.*y+4*(y.^2-1).*y.^2
    #          loss = (p) ->camel(p[1],p[2])
    #          [x f] = ANNEAL[loss(),[0 0]]
    #     We get output:
    #               Initial temperature:     	1
    #               Final temperature:       	3.21388e-007
    #               Consecutive rejections:  	1027
    #               Number of function calls:	6220
    #               Total final loss():        	-1.03163
    #               x =
    #                  -0.0899    0.7127
    #               f =
    #                  -1.0316
    #     Which reasonably approximates the analytical global minimum (note
    #     that due to randomness, your results will likely not be exactly the
    #     same).
    #  Reference:
    #    Kirkpatrick, S., Gelatt, C.D., & Vecchi, M.P. (1983). Optimization by
    #    Simulated Annealing. _Science, 220_, 671-680.
    #   joachim.vandekerckhove@psy.kuleuven.be
    #   $Revision: v5 $  $Date: 2006/04/26 12:54:04 $ Explanation from Joachim
    ####################################################################################################
    # End description from source material
    ####################################################################################################

    ####################################################################################################
    # Set up annealing conditions and behavior
    ####################################################################################################

    # Settings from options
    Tinit                   = options.InitTemp  # Initial annealing temperature
    minT                    = options.StopTemp  # Temperature that indicates the end of simAnneal when reached
    minF                    = options.StopVal   # Value of loss() objective cost function that is sufficient to stop fitting
    max_consec_rejections   = options.MaxConsRej# Maximum number of rejected new solutions without an acceptance
    max_try                 = options.MaxTries  # Maximum total new solutions that can be tried, acceptance or rejection, at a given temperature
    max_success             = options.MaxSuccess# Maximum number of acceptances at a given temperature: length of Markiv chain
    report                  = options.Verbosity # Reporting to command line: (0) Suppress output, (1) Final report,
                                                #                            (2) Report at each new temeperature or when new best is found

    # Initialize annealing properties
    itry ::Int64    = 0             # Set number of attempted solutions at initial temperature to zero
    success ::Int64 = 0             # Set number of accepted solutions at initial temperature to zero
    finished ::Bool = false         # Initialize algorithm as not finished
    consec ::Int64  = 0             # Set consecutive failures to zero
    T               = Tinit         # Assign initial temperature
    initenergy      = loss(parent)  # Calculate loss energy associated with input parent
    oldenergy       = initenergy    # Initialize the previous energy for comparison to the initial
    total::Int64    = 0             # Initialize total attempted solutions for entire annealing algorithm

    bestEnergy      = initenergy        # Initialize lowest energy found over course of annealing
    bestParams      = deepcopy(parent)  # Initialize the best parameter set corresponding to the lowest energy found

    println(bestParams)
    println(initenergy)
    
    # report = 2
    # If requested, print initial temperature, energy, and parameter array
    if report==2
         println("T = ", T, ", loss = ", oldenergy)
         println(parent)
         flush(stdout)
    end

    ####################################################################################################
    # Run annealing algorithm until an end condition is reached
    ####################################################################################################

    # Loop until end condition indicated
    while !finished

        # itry counts the number of searches at EACH TEMP, and gets reset with each cooling
        itry        +=  1

        # Parent is the inital guess, current is the current parameter set, parent gets updated later
        current     =   parent


        ####################################################################################################
        # Check to see if the algorithm is finished with the current temperature
        ####################################################################################################

        # You have either hit a predefined attempt limit, or max successes WITHIN ONE TEMPERATURE
        if itry >= max_try || success >= max_success

            # Only way to end the loop: cooled all the way or too many rejections (new changes not accepted)
            if T < minT || consec >= max_consec_rejections

                # Indicated algorithm completion
                finished    =   true

                # Keeps track of all attempts (multiple within each temperature)
                total       =   total + itry
                break

            # Drop temperature and try again if not finished with algorithm
            else
                # Decrease T according to cooling schedule function
                T = cool(T)

                # Report new temperature, energy, and parameter array if requested
                if report==2
                    println("T = ", T, ", loss = ", oldenergy)
                    println(current)
                    flush(stdout)
                end

                # Update counters
                total           = total + itry      # Keeps track of all attempts (multiple within each temperature)
                itry    = 1                 # Reset number of tries cause we at a new temp
                success = 1                 # Reset number of successive successes (lol) cause new temp
            end
        end

        ####################################################################################################
        # New attempt: generate new solution (parameter array) and process according to its energy
        ####################################################################################################

        # Perturb the system based on the generator function: newsol (if newsol doesn't create new array,
        # parent will already be = to newparam and rejecting the change doesn't work
        newparam, newenergy = newsol(current, ub, lb)

        # Determine cost of those new parameters in the defined cost function: loss
        acceptanceEnergy = oldenergy-newenergy

        # If we hit some predefined option for "good enough", call it and return
        if newenergy < minF
            parent = newparam       # Set current parameter array to new values
            oldenergy = newenergy   # Update current (old) energy to new value
            break
        end

        # If we get a lower cost than the previous parameters, accept new parameters
        if oldenergy-newenergy > 1e-6

            parent      = newparam  # Add the new parameters to the chain and set them as parent (will now begin outer while loop again)
            oldenergy   = newenergy # Associated energy with the new parent
            success     += 1        # Length of success chain increased
            consec  = 0     # Reset consecutive failures to 0 cuase this is a success


            # Keep track of global best energy found
            # Copy over best so far here, if we just assign a pointer it may get lost, so you need to deepcopy
            if oldenergy < bestEnergy
                bestEnergy  =   oldenergy           # Update best overall energy
                bestParams  =   deepcopy(parent)    # Update best overall parameter array

                # Print new best parameters and associated loss if requested
                if report == 2
                    println("vvv here is best at ", bestEnergy, " vvv")
                    println(bestParams)
                    flush(stdout)
                    flush(stdout)
                end
            end


        # If energy is higher, use acceptance criteria to decide whether or not to accept new value
        # CONSECUTIVE REJECTIONS IS NOT RESET: IT IS ONLY RESET WHEN A LOWER ENERGY IS FOUND
        else

            # According to type of simulated anneal variation, and current temperature, potentially accept
            if acceptance(acceptanceEnergy, T, type)
        #
                parent      =   newparam    # Update current parameter array
                oldenergy   =   newenergy   # Update current (old) energy
                success     =   success+1   # Add one to success count

            # If new, higher energy is rejected: don't update parameter array and energy
            else
                consec += 1                 # Add one to consecutive rejection array
            end

        end
    end

    ####################################################################################################
    # Report and return the best energy and corresponding parameter array
    ####################################################################################################
    minimum = bestParams
    fval = bestEnergy

    # Print comment about any difference between the final and best parameter arrays
    if report > 0 && oldenergy-bestEnergy > 1e-6 # in other words, bestEnergy < oldenergy
        # minimum = bestParams
        # fval = bestEnergy
        println("You have a different minimum than where the search ended, hmmmmmm?")
        println("bestEnergy = ", bestEnergy, " at ",  bestParams)
        println("end energy = ", oldenergy, " at ",  parent)
        flush(stdout)
    # If there is no difference, sreturn the final array and energy
    # Changing to return the best, not sure why originally did otherwise
    # else
        
        # minimum = parent
        # fval = oldenergy
    end

    # If requested, report the end conditions of the annealing algorithm
    if report > 0
        println("Initial temperature: ", Tinit)
        println("Final temperature: ", T)
        println("Consecutive rejections: ", consec)
        println("Number of function calls: ", total)
        println("Total final loss(): ", fval)
        println("End parameter: ", minimum)
        flush(stdout)
    end

    # Return best energy and corresponding parameter array
    return [minimum,fval]

end

####################################################################################################
# Function that accepts new higher energy with probability weighted by energy gain and temperature
####################################################################################################
function acceptance(acceptanceEnergy::Float64, T::Float64,
                    type::String)

    # Stand-in for Boltzmann constant
    k::Float64 = 1.

    # Allow for type variations
    if type == "Kirkpatrick"
        # Probability of acceptance based on energy difference and annealing temperature
        if rand() < exp( (acceptanceEnergy)/(k*T) )
            return true
        else
            return false
        end
    else
        error("no valid type defined")
    end
end
