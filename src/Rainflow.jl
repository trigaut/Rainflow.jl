module Rainflow
import Base.show

export sort_peaks, find_boundary_vals, count_cycles

immutable Cycle  # This is the information stored for each cycle found
    count::Float64
    range::Float64
    mean::Float64  # value
    v_s::Float64   # value start
    t_s::Float64   # time start
    v_e::Float64   # value end
    t_e::Float64   # time end
end

type Cycles_bounds #
    min_mean::Float64
    max_mean::Float64
    max_range::Float64
end

""" This function sort out points where the slope is changing sign."""
function sort_peaks(signal::AbstractArray{Float64,1}, dt=[0.:length(signal)-1.])
    slope = diff(signal)
    extremas = slope[1:end-1].*slope[2:end].<-1e-6
    # Determines if the point is local extremum
    is_extremum = [true, extremas..., true]
    return signal[is_extremum] , dt[is_extremum]
end


function cycle(count::Float64, v_s::Float64, t_s::Float64, v_e::Float64, t_e::Float64)
    Cycle(count, abs(v_s-v_e), (v_s+v_e)/2, v_s, t_s, v_e, t_e)
end

""" Count the cycles from the data. """
function count_cycles_from_ext(ext_in::Array{Float64,1},t::Array{Float64,1})
    ext = copy(ext_in) # Makes a copy because there is going to be sorted in the vectors
    time = copy(t)
    i = 1
    j = 2
    cycles = Cycle[]
    @inbounds begin
    while length(ext)>(i+1)
        Y = abs(ext[i+1]-ext[i])
        X = abs(ext[j+1]-ext[j])
        if X>=Y
            if i == 1 # This case counts a half and cycle deletes the poit that is counted
                push!(cycles,cycle(0.5 ,ext[i], time[i], ext[i+1],time[i+1]))
                shift!(ext) # Removes the first entrance in ext and time
                shift!(time)
            else # This case counts one cycle and deletes the poit that is counted
                push!(cycles,cycle(1. ,ext[i], time[i], ext[i+1],time[i+1]))
                splice!(ext,i+1)  # Removes the i and i+1 entrance in ext and time
                splice!(ext,i)
                splice!(time,i+1)
                splice!(time,i)
            end
            i = 1
            j = 2
        else
            i += 1
            j += 1
        end
    end
    for i=1:length(ext)-1 # This counts the rest of the points that have not been counted as a half cycle
        push!(cycles,cycle(0.5 ,ext[i], time[i], ext[i+1],time[i+1]))
    end
    end
    return cycles
end

function count_cycles(signal::Array{Float64,1}, dt = 1:1.:length(signal))
    ext, t = sort_peaks(signal, dt)
    count_cycles_from_ext(ext, t)
end

""" Find the minimum and maximum mean value and maximum range from a vector of cycles. """
function find_boundary_vals(cycles::Array{Cycle,1})
    bounds = Cycles_bounds(Inf, -Inf, -Inf)
    for cycle in cycles
        cycle.mean > bounds.max_mean && bounds.max_mean = cycle.mean
        cycle.mean < bounds.min_mean && bounds.min_mean = cycle.mean
        cycle.range > bounds.max_range && bounds.max_range = cycle.range
    end
    bounds
end


end 
