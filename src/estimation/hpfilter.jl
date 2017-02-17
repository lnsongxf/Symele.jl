"""
    hpfilter(y,phi::Number)
HP filter the data in y with a phi as given. 
Returns two arrays: trend,cycle
"""
function hpfilter(y,phi::Number)
    T = size(y,1)
    HP = zeros(Float64,T,T)
    if !isa(y,Array)
        # Probably a data array. Need to convert because A_ldiv_B doesn't exist for DataArrays
        y = convert(Array{eltype(y)},y)
    end
    # End points need different treatment
    HP[1,1] = 1.0+phi
    HP[1,2] = -2.0*phi
    HP[1,3] = phi
    HP[2,1] = -2.0*phi
    HP[2,2] = 1.0+5.0*phi
    HP[2,3] = -4.0*phi
    HP[2,4] = phi
    HP[T-1,T-3] = phi
    HP[T-1,T-2] = -4.0*phi
    HP[T-1,T-1] = 1.0+5.0*phi
    HP[T-1,T] = -2.0*phi
    HP[T,T-2] = phi
    HP[T,T-1] = -2.0*phi
    HP[T,T] = 1.0+phi

    for i in 3:T-2
        HP[i,i-2] = phi
        HP[i,i-1] = -4.0*phi
        HP[i,i] = 1.0+6.0*phi
        HP[i,i+1] = -4.0*phi
        HP[i,i+2] = phi
    end

    trend = HP\y
    cycle = y-trend
    return trend,cycle
end
