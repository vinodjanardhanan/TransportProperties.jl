module TransportProperties

using LinearAlgebra, StaticArrays
using LightXML, Printf

using RxnHelperUtils
using IdealGas

include("Constants.jl")

export viscosity, D_km!, thermal_coductivity, D_ij, D_ii!, create_transport_data, transport_properties

#defne the collision integral Omega11
omega11 = @SMatrix[
    4.008e+000 -1.0380e+00  5.9659e+00 -2.9977e+00  4.9812e-01;
    3.1300e+00 -4.6244e-01  2.2622e+00 -8.4707e-01  1.0961e-01;
    2.6490e+00 -3.6366e-01  1.4543e+00 -4.7309e-01  5.2330e-02;
    2.3140e+00 -2.6379e-01  1.0326e+00 -2.9782e-01  2.7717e-02;
    2.0660e+00 -1.9750e-01  8.0301e-01 -2.1612e-01  1.7991e-02;
    1.8770e+00 -2.5097e-01  9.7452e-01 -4.5890e-01  9.1504e-02;
    1.7290e+00 -1.1931e-01  5.5652e-01 -1.4021e-01  1.0704e-02;
    1.6122e+00 -8.8902e-02  4.6422e-01 -1.0763e-01  6.7026e-03;
    1.5170e+00 -6.6921e-02  3.9820e-01 -8.7569e-02  4.7731e-03;
    1.4400e+00 -5.6468e-02  3.5399e-01 -7.6371e-02  4.0889e-03;
    1.3204e+00 -4.0079e-02  2.8619e-01 -5.9986e-02  3.2544e-03;
    1.2340e+00 -3.4407e-02  2.4399e-01 -5.2429e-02  3.4858e-03;
    1.1680e+00 -1.4465e-02  1.7825e-01 -2.5006e-02 -7.2435e-04;
    1.1166e+00 -1.3643e-02  1.5968e-01 -2.5965e-02  6.4166e-04;
    1.0750e+00 -6.3727e-03  1.2820e-01 -1.4746e-02 -9.7281e-04;
    1.0006e+00 -9.1095e-03  1.0503e-01 -1.6931e-02  6.9808e-04;
    9.5000e-01 -4.6060e-03  7.7117e-02 -1.0881e-02  3.7158e-04;
    9.1310e-01 -1.9581e-03  5.6961e-02 -5.2681e-03 -3.0468e-04;
    8.8450e-01 -2.9286e-03  5.0573e-02 -7.3139e-03  4.8588e-04;
    8.4280e-01 -1.4158e-03  3.3049e-02 -2.7314e-03 -1.2165e-04;
    8.1300e-01 -1.3882e-03  2.4887e-02 -2.0834e-03 -2.3635e-05;
    7.8980e-01 -1.9075e-05  1.7423e-02 -4.1225e-04 -2.2042e-04;
    7.7110e-01  2.0888e-05  1.3897e-02 -2.4013e-04 -1.7556e-04;
    7.5550e-01 -4.5415e-05  1.1903e-02 -6.0156e-04 -3.7943e-05;
    7.4220e-01  5.4903e-04  8.6391e-03  2.7416e-04 -1.5004e-04;
    7.2022e-01 -2.7033e-04  7.7439e-03 -6.1756e-04  4.3613e-05;
    7.0250e-01  8.3428e-04  3.8508e-03  8.2584e-04 -2.0902e-04;
    6.8776e-01  3.9894e-05  3.9815e-03  2.0106e-04 -8.6618e-05;
    6.7510e-01  3.8786e-05  3.3483e-03  5.6978e-05 -3.8795e-05;
    6.6400e-01  3.2925e-04  2.3460e-03  3.2880e-04 -9.2696e-05;
    6.4140e-01 -1.6252e-04  2.1257e-03 -4.1127e-05 -1.6239e-05;
    6.2350e-01  3.8876e-05  1.4220e-03 -2.7050e-05 -1.7042e-06;
    6.0882e-01 -1.2964e-04  1.5607e-03 -4.0253e-04  8.8076e-05;
    5.9640e-01  2.6913e-04  8.8982e-05  6.1738e-04 -1.4247e-04;
    5.7630e-01 -1.8613e-04  9.0887e-04 -2.7653e-04  7.1825e-05;
    5.4150e-01 -5.7840e-05  4.9071e-04 -2.3284e-04  5.6760e-05;
    5.1800e-01  1.5100e-04  7.2383e-04 -5.6918e-04  1.2030e-04 
]

# Collision integral Omega22
omega22  = @SMatrix [
    4.1000e+00 -4.9400e-01  5.1705e+00 -2.4579e+00  3.8700e-01;
    3.2630e+00 -4.8580e-01  2.4444e+00 -9.0334e-01  1.1276e-01;
    2.8400e+00 -4.7159e-01  1.5619e+00 -4.5765e-01  4.0936e-02;
    2.5310e+00 -3.6736e-01  1.0905e+00 -2.6097e-01  1.3879e-02;
    2.2840e+00 -2.8365e-01  8.4186e-01 -1.8024e-01  5.8124e-03;
    2.0840e+00 -2.1794e-01  6.9371e-01 -1.4468e-01  4.3636e-03;
    1.9220e+00 -2.0893e-01  6.9031e-01 -1.8193e-01  1.4363e-02;
    1.7902e+00 -1.2932e-01  5.1435e-01 -1.0742e-01  3.8897e-03;
    1.6820e+00 -1.0196e-01  4.5847e-01 -9.7322e-02  4.1072e-03;
    1.5930e+00 -8.5797e-02  4.2360e-01 -9.6560e-02  5.8448e-03;
    1.4550e+00 -4.8607e-02  3.2642e-01 -6.4005e-02  1.9493e-03;
    1.3550e+00 -3.2972e-02  2.7273e-01 -5.1800e-02  1.5361e-03;
    1.2800e+00 -2.4970e-02  2.3581e-01 -4.5617e-02  1.8995e-03;
    1.2220e+00 -1.3736e-02  1.9366e-01 -3.2066e-02  3.8064e-04;
    1.1760e+00 -1.2000e-02  1.6956e-01 -2.7242e-02  3.7995e-04;
    1.0933e+00 -6.3451e-03  1.2403e-01 -1.7897e-02  1.2523e-04;
    1.0390e+00 -2.1846e-03  9.0757e-02 -1.0363e-02 -2.7932e-04;
    9.9960e-01 -1.8453e-03  7.2675e-02 -8.3026e-03  1.1244e-05;
    9.6990e-01 -2.0468e-03  5.9254e-02 -5.6424e-03 -2.8151e-04;
    9.2680e-01 -1.5662e-03  4.1725e-02 -3.9558e-03 -4.8331e-05;
    8.9620e-01 -5.9689e-04  2.9341e-02 -1.8027e-03 -1.5170e-04;
    8.7270e-01  1.2046e-04  2.1443e-02 -6.0037e-04 -2.1344e-04;
    8.5380e-01  2.3847e-04  1.6539e-02 -7.9724e-05 -2.1286e-04;
    8.3790e-01 -2.8887e-04  1.4856e-02 -8.8946e-04 -3.7426e-06;
    8.2430e-01  5.6032e-04  1.0655e-02  3.2629e-04 -1.8930e-04;
    8.0180e-01  5.2094e-04  7.2131e-03  6.8005e-04 -2.1524e-04;
    7.8360e-01 -1.0792e-04  6.8273e-03 -4.4591e-04  5.5406e-05;
    7.6830e-01  5.7900e-04  3.4869e-03  1.1370e-03 -2.8172e-04;
    7.5520e-01  1.0191e-04  3.7993e-03  8.7340e-05 -2.1863e-05;
    7.4360e-01  1.1106e-04  3.1512e-03  1.9598e-04 -6.2322e-05;
    7.1982e-01  5.9318e-05  2.1699e-03  9.4738e-05 -2.8537e-05;
    7.0100e-01  1.2696e-04  1.3406e-03  2.2449e-04 -5.1240e-05;
    6.8545e-01  5.2502e-04 -3.2893e-05  9.4361e-04 -2.1258e-04;
    6.7230e-01  3.0299e-05  1.4509e-03 -5.6073e-04  1.4627e-04;
    6.5100e-01 -3.5478e-05  6.2081e-04  4.7129e-05 -2.3793e-05;
    6.1400e-01 -5.1551e-05  1.6994e-03 -1.2009e-03  2.2986e-04;
    5.8870e-01  4.7810e-04  3.2386e-03 -2.7109e-03  5.3094e-04
]

# Reduced temperature tstart

tstar = @SVector[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
                2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0,
                25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0 
]

#=
Structure for defining the transport properties of pure species as read from transport.dat file
    While storing all parameters are stored in SI units
=#
struct TransportData 
    name::String            # Name of the species
    geom::Int64             # monoatomic, linear, non-linear
    ϵ::Float64              # LJ potential well depth (J)
    σ::Float64              # LJ collision diameter (m)
    μ::Float64              # dipole moment in Debye (m^3/2 J^1/2)
    α::Float64              # polarizability in cubic angstroms (m^3)
    Zrot::Float64           # rotational relaxation collision number
end


"""
viscosity(sp_trd::Array{TransportData}, T::Float64, molwt::Array{Float64})
Function for the calculation of pure species viscosity in Kg/m-s
#   Usage
    viscosity(sp_trd, T, molwt)
-   sp_trd : Array of species transport data
-   T : Temperature
-   molwt : molecular weights array
"""
function viscosity(sp_trd::Array{TransportData}, T::Float64, molwt::Array{Float64})
    #calculate the reduced temperature
    eps_reciprocal = [1/i.ϵ for i in sp_trd]    
    red_Temp = kB*T*eps_reciprocal
    
    #reduced dipole moment        
    δ = [td.μ^2/(2*td.ϵ*td.σ^3) for td in sp_trd]
    Ω22 = get_omega(red_Temp,δ,key="Omega22")
    
    #Calculate pure species viscosity
    vector_term = [sqrt(molwt[i])/(sp_trd[i].σ^2*Ω22[i]) for i in eachindex(sp_trd)]
    return (5.0*(sqrt(π*kB*T/Na))/(16.0*π))*vector_term    
end

"""
viscosity(sp_trd::Array{TransportData}, T::Float64, molwt::Array{Float64}, mole_fracs::Array{Float64})
Calculation of mixture viscosity in Kg/m-s
#   Usage
    viscosity(sp_trd, T, molwt, mole_fracs)
-   sp_trd: Species transport data
-   T: Temperature
-   molwt : molecular weight vector 
-   mole_fracs : Mole fractions vector 
"""
function viscosity(sp_trd::Array{TransportData}, T::Float64, molwt::Array{Float64}, mole_fracs::Array{Float64})
    μ = viscosity(sp_trd, T, molwt)
    n = length(mole_fracs)
    Φ_kj = zeros(n)
    X = zeros(n)
    for k in 1:n        
        for j in 1:n
            ml_wt_r = molwt[k]/molwt[j]
            visc_r = μ[k]/μ[j]
            # Φ_kj[j] = (1/sqrt(8))*(1+(molwt[k]/molwt[j]))^(-0.5) * (1+ sqrt(μ[k]/μ[j]) * (molwt[j]/molwt[k])^(0.25) )^2
            Φ_kj[j] = (1.0/sqrt(8.0))*(1+ml_wt_r)^(-0.5)*(1+ visc_r^(0.5) *(1.0/ml_wt_r)^(0.25))^2             
        end
        denom = sum(mole_fracs .* Φ_kj)
        X[k] = μ[k]*mole_fracs[k]/denom
        
    end
    return sum(X)
end

"""
Function for calculating diffusion coefficient of a species in a mixture
# Usage 
    D_km!(Dkm, sp_trd, T, p, molwt, molefracs)
-   Dkm : Array for storing the mixture diffusion coefficients
-   sp_trd : Array of the type TransportData
-   T : temperature in K 
-   p : Pressure in Pa 
-   molwt : species molecular weights 
-   molefracs : species mole fractions 
"""
function D_km!(Dkm::Array{Float64}, sp_trd::Array{TransportData}, T::Float64, p::Float64, molwt::Array{Float64}, molefracs::Array{Float64})
    massfracs = similar(molefracs)
    molefrac_to_massfrac!(massfracs,molefracs,molwt)
    #Get the binary diffusion coefficients
    bdc = D_ij(sp_trd,T,p,molwt)    
    
    for j in 1:length(molefracs)
        denom = 0.0                
        for k in 1:length(molefracs)
            if j!=k 
                denom += molefracs[k]/bdc[j,k]
            end
        end        
        Dkm[j] = (1-massfracs[j])/denom        
    end    
end


"""
A function to calculate the mixture diffusion coefficients given the binary diffusion coefficients and the mole fractions
# Usage:  
D_km!(Dkm, D_ij, molefracs, molwt )
-   Dkm : Array to store the mixture diffusion coefficients (size N)
-   D_ij : N X N Matrix of binary diffusion coefficients
-   molefracs : Mole fractions 
-   molwt : species molecular weights 
"""
function D_km!(Dkm::Array{Float64}, D_ij, molefracs, molwt)
    massfracs = similar(molefracs)
    molefrac_to_massfrac!(massfracs,molefracs,molwt)
    for j in 1:length(molefracs)
        denom = 0.0                
        for k in 1:length(molefracs)
            if j!=k 
                denom += molefracs[k]/D_ij[j,k]
            end
        end        
        Dkm[j] = (1-massfracs[j])/denom        
    end    
    return Dkm
end

"""
Calculates the thermal conductivity of the mixture 
# Usage 
thermal_coductivity(sp_trd, T, p, thermo_obj, molefracs)
sp_trd : Species transport data
T : Temperature in K 
p : Pressure in Pa 
thermo_ob : SpeciesThermoObj 
molefracs : Species mole fractions 
"""
function thermal_coductivity(sp_trd::Array{TransportData}, T::Float64, p::Float64, thermo_obj::IdealGas.SpeciesThermoObj,molefracs::Array{Float64})
    molwt = thermo_obj.molwt    
    sdc = zeros(length(molwt))
    D_ii!(sdc, sp_trd,T,p,molwt)
    μ = viscosity(sp_trd,T,molwt)    
    λ = 0
    cp = IdealGas.cp_all(thermo_obj,T)
    cv = cp .- R
    
    cv_rot = cv_vib = 0
    #Vector of density
    ρ = (p/R/T)*molwt
    t1 = t2 = 0    
    cv_trans = 1.5R
    f_trans = f_rot = f_vib = 0
    for i in eachindex(sp_trd)
        if sp_trd[i].geom == 0             
            f_trans = 2.5             
        else
            if sp_trd[i].geom == 1            
                cv_rot = R 
                cv_vib = cv[i] - 2.5R
            else
                cv_rot = 1.5R 
                cv_vib = cv[i]-3R
            end

            val = sp_trd[i].ϵ/kB/T
            val0 = sp_trd[i].ϵ/kB/298.0

            ft(x) = 1+ 0.5*π^(1.5)*sqrt(x) + ( (π^2/4.0)+2 )*x + π^(1.5)*x^(1.5)
            Zrot = sp_trd[i].Zrot*ft(val0)/ft(val)

            f_vib = ρ[i]*sdc[i]/μ[i]
            B = Zrot + (2.0/π)*( (5cv_rot/3R) + f_vib )
            A = 2.5-f_vib
            AB = A/B            
            f_trans = 2.5*(1- (2cv_rot/π/cv_trans)*AB )
            f_rot = f_vib*(1+2AB/π)
        end        
        λ = μ[i]*(f_trans*cv_trans + f_rot*cv_rot + f_vib*cv_vib)/molwt[i]           
        t1 += molefracs[i]*λ
        t2 += molefracs[i]/λ        
    end
    return 0.5*(t1+1.0/t2)    
end


"""
Function for the calculation of binary diffusion coefficients
# Usage     
    D_ij(sp_trd, T, p, molwt)
-   sp_trd : Array of Species transport data 
-   T : Temperature K 
-   p : Pressure Pa 
-   molwt : Species molecular weight 
"""
function D_ij(sp_trd::Array{TransportData}, T::Float64, p::Float64, molwt::Array{Float64})
    n = length(molwt)

    function reduced_T_and_delta(trd_p, trd_np)
        
        ϵ_jk = σ_jk = μ_jk_sq = 0
        if trd_p.α == 0.0 && trd_np.α == 0.0 || trd_p.α != 0.0 && trd_np.α != 0.0
            ϵ_jk = sqrt(trd_p.ϵ*trd_np.ϵ)
            σ_jk= 0.5*(trd_p.σ+trd_np.σ)
            μ_jk_sq = trd_p.μ * trd_np.μ
        else
            if trd_p.α != 0
                trd_p, trd_np = trd_np, trd_p
            end
            α_n = trd_np.α/trd_np.σ^3
            μ_p = trd_p.μ/sqrt(trd_p.ϵ*trd_p.σ^3)
            ξ = 1+ 0.25 * α_n * μ_p * sqrt(trd_p.ϵ/trd_np.ϵ) 
            ϵ_jk = ξ^2 * sqrt(trd_p.ϵ * trd_np.ϵ)
            σ_jk = 0.5*(trd_np.σ + trd_p.σ)*ξ^(-1.0/6.0)
            μ_jk_sq = 0
        end
        rT = kB*T/ϵ_jk
        δ_jk = 0.5*μ_jk_sq/(ϵ_jk*σ_jk^3)
        return rT,δ_jk,σ_jk
    end

    red_Temp = zeros(n,n)
    δ = zeros(n,n)
    σ = zeros(n,n)
    for j in 1:n
        red_Temp[j,j],δ[j,j],σ[j,j] = reduced_T_and_delta(sp_trd[j],sp_trd[j])            
        for k in j+1:n              
            red_Temp[j,k],δ[j,k],σ[j,k] = reduced_T_and_delta(sp_trd[j],sp_trd[k])            
            red_Temp[k,j] = red_Temp[j,k]
            δ[k,j] = δ[j,k]           
            σ[k,j] = σ[j,k]
        end
    end

    D = zeros(n,n)
    const_term = 3.0*sqrt(2*π*Na*kB^3*T^3)/(16*p*π)
    for j in 1:n
        Ω11 = get_omega(red_Temp[j,:],δ[j,:],key="Omega11")             
        for k in j+1:n
            m_jk = molwt[j]*molwt[k]/(molwt[j]+molwt[k])
            D[j,k] = const_term * sqrt(1/m_jk)/(σ[j,k]^2*Ω11[k])
            D[k,j] = D[j,k]
        end
        # Following two lines calculate the self diffusion coefficients
        m_jj = 0.5*molwt[j]
        D[j,j] = const_term * sqrt(1/m_jj)/(σ[j,j]^2*Ω11[j])
    end
    
   return D
    
end


"""
Function to calculate self diffusion coefficients
# Usage
    D_ii!(D::Array{Float64},sp_trd::Array{TransportData}, T::Float64, p::Float64, molwt::Array{Float64})
-   D : Array to store the diffusion coefficients (size N)    
-   sp_trd : Array of species transport data 
-   T : Temperature (K) 
-   p : Pressure (Pa 
-   molwt : specoes molecular weights

"""
function D_ii!(D::Array{Float64},sp_trd::Array{TransportData}, T::Float64, p::Float64, molwt::Array{Float64})
    n = length(molwt)    
    #Calculate reduced temperature     
    eps_reciprocal = [1/i.ϵ for i in sp_trd]    
    red_Temp = kB*T*eps_reciprocal
    #Dipole moment
    δ = [td.μ^2/(2*td.ϵ*td.σ^3) for td in sp_trd]

    Ω11 = get_omega(red_Temp,δ,key="Omega11")    
    const_term = 3.0*sqrt(2*π*Na*(kB*T)^3)/(16*p*π)
    for j in 1:n
        D[j]= const_term * sqrt(1/molwt[j])/(sp_trd[j].σ^2*Ω11[j])
    end
    
end


#=
calculation of collision integral 
=#
function get_omega(red_Temp::Array{Float64}, δ::Array{Float64};kwargs...)    
    T_coordinate = Array{Int64,1}()
    for tr in red_Temp
        local T_xmin = 0
        for i in eachindex(tstar)
            if tr < tstar[i+1] && tr > tstar[i]
                T_xmin = i
                break
            end
        end
        T_xmax = min(length(tstar),T_xmin+3)
        if T_xmax == length(tstar)
            T_xmin = T_xmax-3            
        end
        push!(T_coordinate,T_xmin-1)
    end    
    
    o_xx = Matrix{Float64}(undef,length(T_coordinate),3)
    for i in eachindex(T_coordinate)        
        if δ[i] < 1e-2             
            if kwargs[Symbol("key")]=="Omega22"
                o_xx[i,:] = omega22[T_coordinate[i]:T_coordinate[i]+2,1]  
            else
                o_xx[i,:] = omega11[T_coordinate[i]:T_coordinate[i]+2,1]  
            end
        else
            if kwargs[Symbol("key")]=="Omega22"
                o_xx[i,:] = polyfit(δ[i],T_coordinate[i],key="Omega22")
            else
                o_xx[i,:] = polyfit(δ[i],T_coordinate[i],key="Omega11")
            end
        end                    
    end
    red_Temp = [ t < tstar[1] ? tstar[2] : t > 500 ? 500 : t for t in red_Temp]
    
    Ω_XX = Array{Float64,1}()
    quad_interpolation!(Ω_XX,red_Temp,T_coordinate,o_xx)
    return Ω_XX
end


function quad_interpolation!(Ω::Array{Float64}, red_Temp::Array{Float64}, tc::Array{Int64} ,o22::Matrix{Float64})     
    for i in eachindex(red_Temp)
        if red_Temp[i] < 100             
            x12 = 0.5*(tstar[tc[i]]+tstar[tc[i]+1])
            x23 = 0.5*(tstar[tc[i]+1]+tstar[tc[i]+2])
            x2x = 0.5*(tstar[tc[i]+1]+red_Temp[i])
            m12 = (o22[i,2]-o22[i,1])/(tstar[tc[i]+1]-tstar[tc[i]])
            m23 = (o22[i,3]-o22[i,2])/(tstar[tc[i]+2]-tstar[tc[i]+1])
            m2x = (m12*(x23-x2x) + m23*(x2x-x12))/(x23-x12)         
            # println(o22[i,2], "\t", m2x, "\t", red_Temp[i], "\t", tstar[tc[i]+1])   
            push!(Ω,o22[i,2] + m2x*(red_Temp[i]-tstar[tc[i]+1]))
        else
            push!(Ω,0.703 + red_Temp[i]*(-1.46e-3+ red_Temp[i]*(3.57e-6+ red_Temp[i]*(-3.43e-9))))
        end
    end
end


function polyfit(rd::Float64,tc::Int64;kwargs...)    
    oxx = Array{Float64,1}()
    n_column = size(omega22,2)
    rdvec = [ rd^i for i in range(n_column-1,step=-1,stop=0)]        
    reverse!(rdvec)
    
    if kwargs[Symbol("key")]=="Omega22"                        
        o22_times_rd = omega22[tc:tc+2,:]' .* rdvec
        oxx = [sum(o22_times_rd[:,i]) for i in 1:size(o22_times_rd,2)]
    elseif kwargs[Symbol("key")]=="Omega11"
        o11_times_rd = omega11[tc:tc+2,:]' .* rdvec
        oxx = [sum(o11_times_rd[:,i]) for i in 1:size(o11_times_rd,2)]
    end
    return oxx
end


"""
A function to read the transport.dat file and create the structure TransportData
#   Usage
    create_transport_data(gasphase, trans_file)
-   gasphase: list of gasphase species 
-   trans_file: path to the transport.dat file
"""
function create_transport_data(gasphase::Array{String} , trans_file::AbstractString)
    species_trasport_data = Array{TransportData,1}(undef,length(gasphase))
    open(trans_file) do io
        while !eof(io)
            data_string = readline(io)
            data = split(data_string)
            if in(String(strip(data[1])),gasphase)
                species_trasport_data[get_index(String(strip(data[1])),gasphase)] = get_td_struct(data_string)
            end
        end
    end
    return species_trasport_data
end

function get_td_struct(data_str::String)
    data = split(data_str)
    name = String(strip(data[1]))
    geom = parse(Int64,data[2])
    ljwd = parse(Float64,data[3])*kB   # What is given in the input file is ϵ/kB (Values are stored in J)
    ljcd = parse(Float64,data[4])*Angstroms # Values are stored in m
    dipole = parse(Float64,data[5])*Debye   # Values are stored in m^(3/2) J^(1/2)
    polarizability = parse(Float64,data[6])*Angstroms^3 # Values are in m^3    
    zrot = parse(Float64,data[7]) # Unitless
    return TransportData(name,geom,ljwd,ljcd,dipole,polarizability,zrot)
end




"""
This function is for testing the transport code
"""
function transport_properties(file_path::AbstractString, lib_dir::AbstractString)
    xmldoc = parse_file(file_path)
    xmlroot = root(xmldoc)

    #get the gasphase species
    gasphase = get_collection_from_xml(xmlroot,"gasphase")
    thermo_file = joinpath(lib_dir,"therm.dat") 
    thermo_all = IdealGas.create_thermo(gasphase,thermo_file)
    mole_fracs = get_molefraction_from_xml(xmlroot,thermo_all.molwt,gasphase)
    local T = get_value_from_xml(xmlroot,"T")
    local p = get_value_from_xml(xmlroot,"p")
    trans_file =  joinpath(lib_dir, "transport.dat") 
    sp_tr_data = create_transport_data(gasphase,trans_file )
    
    #Evaluate viscosity
    μ =  viscosity(sp_tr_data,T,thermo_all.molwt)    
    println("\nPure species viscosity:")  
    println("     Species \t viscosity(Kg/m-s)")
    println("-----------------------------------")  
    for i in eachindex(gasphase)        
        @printf("%12s \t      %.4e\n", gasphase[i], μ[i])
    end
    μ = viscosity(sp_tr_data,T,thermo_all.molwt,mole_fracs)
    @printf("Mixture viscosity: %.4e Kg/m-s\n",μ)
    
    #Evaluate diffusion coefficients
    println("\nMixture diffusion coefficients:")
    Dkm = zeros(length(gasphase))
    D_km!(Dkm,sp_tr_data,T,p,thermo_all.molwt,mole_fracs)
    println("     Species \t Diff.Coeff(m^2/s)")
    println("-----------------------------------")  
    for i in eachindex(gasphase)
        @printf("%12s \t      %.4e\n",gasphase[i], Dkm[i])
    end
    bdc =  D_ij(sp_tr_data,T,p,thermo_all.molwt)
    println("\nBinary diffusion coefficients:")
    # println("Note: self diffusion coefficients are not printed below:")
    println("-----------------------------------")  
    for j in 1:length(gasphase)
        @printf("\t%11s",gasphase[j])
    end
    println()
    for j in 1:length(gasphase)
        @printf("%12s\t",gasphase[j])
        for k in 1:length(gasphase)
            @printf("%.4e\t",bdc[j,k])
        end
        println()
    end

    #Evaluate thermal conductivity
    λ = thermal_coductivity(sp_tr_data,T,p,thermo_all,mole_fracs)
    @printf("\nThermal conductivity of mixture: %.4e (W/m-K)\n", λ)
    
    return Symbol("Success")
end

# end of module
end
