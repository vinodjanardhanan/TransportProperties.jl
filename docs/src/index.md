```@meta
CurrentModule = TransportProperties
```

# TransportProperties
TransportProperties lets you calculate the viscosity, thermal conductivity and diffusion coefficients.
The code relies on *transport.dat* for the evaluation of these properties. A sample of the *transport.dat* file is shown below.

```
AR                 0   136.500     3.330     0.000     0.000     0.000
CH4                2   141.400     3.746     0.000     2.600    13.000
CO2                1   244.000     3.763     0.000     2.650     2.100
CO                 1    98.100     3.650     0.000     1.950     1.800
H2                 1    38.000     2.920     0.000     0.790   280.000
H2O                2   572.400     2.605     1.844     0.000     4.000
O2                 1   107.400     3.458     0.000     1.600     3.800
```

The different columns present in the above table is described below.

- Column-1: Name of species
- Column-2: Geometric configuration of the species; 0- single atom, 1- linear molecule, 2-nonlinear molecule
- Column-3: Lennard-Jones potential well depth expressed as $\epsilon/k_B$ and has units K
- Column-4: Lennard-Jones collision diameter ($\sigma$) in Angstroms
- Column-5: Dipole moment ($\mu$) in Debye
- Column-6: Polarizability ($\alpha$) in cubic Angstroms
- Column-7: Rotational relaxation collision number ($Z_\mathrm{rot}$) at 298 K

Documentation for [TransportProperties](https://github.com/vinodjanardhanan/TransportProperties.jl).

## Installation
To install the package, use the following commands in the julia REPL
```julia
julia> using Pkg
julia> Pkg.add("TransportProperties")
```

## General interfaces
```@index
```

```@autodocs
Modules = [TransportProperties]
```
## Diffusion coefficient calculation
### Binary diffusion coefficients
The binary diffusion coefficient is expressed as 
```math
D_{jk} = \frac{3}{16} \frac{ \sqrt{2\pi N_Ak_B^3T^3/m_{jk}} }{p\pi \sigma_{jk}^2 \Omega^{(1,1)} }
```
Here $N_A$ is the Avogadro's number, $k_B$ is the Boltzmann constant, $p$ is the pressure, and $T$ is the temperature. The reduced molar mass of the species pair (j,k) is defined as

```math
m_{jk} = \frac{m_jm_k}{m_j+m_k}
```

The collision integral $\Omega^{(1,1)}$ is determined using the reduced temperature $T^*_{jk}$ and the reduced dipole moment $\delta^*_{jk}$

```math
T^*_{jk} = \frac{k_BT}{\epsilon_{jk}}
```

```math
\delta^*_{jk} = \frac{1}{2} \frac{\mu^2_{jk}}{\epsilon_{jk}\sigma_{jk}^3}
```

The reduced dipole moment depends on the polarizability of the interacting molecules. In the case where both molecules are 
either polar or non-polar, then it follows

```math
\epsilon_{jk} = \sqrt{\epsilon_j \epsilon_k}
```

```math
\mu_{jk}^2 = \mu_j \mu_k
```

and the reduced collision diameter $\sigma_{jk}$ is defined as

```math
\sigma_{jk} = \frac{\sigma_k+\sigma_j}{2}
```

In the case of interaction between a polar and non-polar molecule

```math
\epsilon_{jk}= \xi^2  \sqrt{\epsilon_j \epsilon_k}
```

```math
\sigma_{jk} = \frac{1}{2} (\sigma_j + \sigma_k) \xi^{-1/6}
```

```math
\mu_{jk}= 0
```

```math
\xi = 1+\frac{1}{4} \alpha_n^* \mu_p^* \sqrt{ \frac{\epsilon_p}{\epsilon_n} }
```
Note that the subscripts $p$ and $n$ represents either $j$ or $k$. If $j$ is polar species, then the subscript $p$ refers to that species

```math
\alpha_n^* = \frac{\alpha_n}{\sigma_n^3}
```

```math
\mu_p^* = \frac{\mu_p}{\sqrt{\epsilon_p\sigma_p^3}}
```
The estimation of $\Omega^{(1,1)}$ is a table look up procedure using the values of $T^*_{jk}$ and $\delta^*_{jk}$.

### Mixture diffusion coefficients
The mixture diffusion coefficients are calculated from

```math
D_{k,m} = \frac{1-Y_k}{\sum_{j\ne k}^N X_j/D_{jk}}
```
In the above equation, $Y_k$ is the mass fraction of the species $k$, and $X_j$ is the mole fraction of species $j$. $N$ is the total number of species

## Viscosity
The viscosity of the mixture is calculated using pure species viscosity. The pure species viscosity is defined as

```math
\eta_k = \frac{5}{16} \frac{\sqrt{\pi m_k k_BT/N_A}}{\pi \sigma_k^2 \Omega^{(2,2)}}
```
The Lennard-Jones collision integral $\Omega^{(2,2)}$ is estimated using a table lookup procedure and depends on the reduced temperature. 

```math
T^*_k = \frac{k_BT}{\epsilon_k}
```
and the reduced dipole moment

```math
\delta^*_k = \frac{1}{2} \frac{\mu_k^2}{\epsilon_k \sigma_k^3}
```
The viscosity of the mixture is then defined as

```math
\eta = \sum_{k=1}^{N} \frac{X_k\eta_k}{\sum{j=1}^K X_j\Phi_{kj}}
```

```math
\Phi_{kj}= \frac{1}{\sqrt{8}} \left(1+\frac{M_k}{M_j}\right)^{-1/2} \left( 1+ \left(\frac{\eta_k}{\eta_j}\right)^{1/2} \left(\frac{M_j}{M_k}\right)^{1/4}\right)^2
```

## Thermal conductivity
Similar to calculating viscosity, the thermal conductivity of a mixture is calculated from the pure species thermal conductivity. 
The pure species thermal conductivity is defined as

```math
\lambda_k = \frac{\eta_k}{M_k} \left( f_t C_{v,t} + f_r C_{v,r} + f_v C_{v,v} \right)
```

```math
f_t = \frac{5}{2}\left( 1-\frac{2}{\pi} \frac{C_{v,r}}{C_{v,t}} \frac{A}{B} \right)
```

```math
f_v = \frac{\rho D_{kk}}{\eta_k}
```


```math
f_r = f_v \left( 1+ \frac{2}{\pi} \frac{A}{B}\right)
```

```math
A = 2.5-f_v
```

```math
B = Z_{rot} + \frac{2}{\pi} \left( \frac{5}{3} \frac{C_{v,r}}{R} + f_v\right)
```
The molar heat capacity $C_v$ for rotational, vibrational or transnational mode depends on the molecule's geometry. For linear molecules

```math
\frac{C_{v,t}}{R} = \frac{3}{2}
```
```math
\frac{C_{v,r}}{R} = 1
``` 
```math
\frac{C_{v,v}}{R} = C_v  - 2.5R
```
```math
C_v = C_p -R
```

For non-linear molecules
```math
\frac{C_{v,t}}{R} = \frac{3}{2}
```
```math
\frac{C_{v,r}}{R} = \frac{3}{2}
``` 
```math
\frac{C_{v,v}}{R} = C_v  - 3R
```
```math
C_v = C_p -R
```

## Executing the code
To see all properties as a screen output
```julia
julia>using TransportProperties
julia>transport_properties("transport.xml", "lib_dir")
```

In the above call, it is assumed that the input file *transport.xml* is present in the working directory and *lib_dir* is the path to the *lib* directory relative to the current working directory. The structure of the *transport.xml* input file is shown below.
```
<?xml version="1.0" encoding="ISO-8859-1"?>
<trans>
	<gasphase>CH4 H2O H2 CO CO2 O2 N2</gasphase>
	<molefractions>CH4=0.125, H2O=0.252, CO2=0.084, N2=0.539</molefractions>
	<T>1073.15</T>
	<p>1e5</p>
</trans>
```
The meaning of the different xml elements are as follows
- <gasphase> : list of gasphase species separated by space
- <molefractions> : mole fractions of the species (instead of <molefractions>, <massfractions> may also be specified)
- <T> : temperature in K
- <p> : pressure in Pa

## Input file download
The xml input file and the *lib* directory containig other required input files may be downloaded from [here](https://github.com/vinodjanardhanan/TransportProperties.jl/tree/main/test).


## Calculation of properties
The following methods may be used to calculate the properties of pure species or mixtures 
###  Pure species voscosity
```julia
julia>using TransportProperties, IdealGas
julia>gasphase = ["CH4", "CO2", "H2O", "H2", "CO"]
julia>sp_tr_data = create_transport_data(gasphase,"transport.dat")
julia>thermo_all = IdealGas.create_thermo(gasphase,"therm.dat")
julia>mu = viscosity(sp_tr_data,T,thermo_all.molwt)
```

###  Mixture viscosity
```julia
julia>using TransportProperties, IdealGas
julia>gasphase = ["CH4", "CO2", "H2O", "H2", "CO"]
julia>sp_tr_data = create_transport_data(gasphase,"transport.dat")
julia>thermo_all = IdealGas.create_thermo(gasphase,"therm.dat")
julia>mu = viscosity(sp_tr_data,T,thermo_all.molwt,mole_fracs)    
```

###  Binary diffusion coefficients
```julia
julia>using TransportProperties, IdealGas
julia>gasphase = ["CH4", "CO2", "H2O", "H2", "CO"]
julia>sp_tr_data = create_transport_data(gasphase,"transport.dat")
julia>thermo_all = IdealGas.create_thermo(gasphase,"therm.dat")
julia>Dij = D_ij(sp_tr_data,T,p,thermo_all.molwt)
```


###  Mixture diffusion coefficients
```julia
julia>using TransportProperties, IdealGas
julia>gasphase = ["CH4", "CO2", "H2O", "H2", "CO"]
julia>sp_tr_data = create_transport_data(gasphase,"transport.dat")
julia>thermo_all = IdealGas.create_thermo(gasphase,"therm.dat")
julia>Dkm = zeros(length(gasphase))
julia>D_km!(Dkm,sp_tr_data,T,p,thermo_all.molwt,mole_fracs)        
```

or

```julia
julia>using TransportProperties, IdealGas
julia>gasphase = ["CH4", "CO2", "H2O", "H2", "CO"]
julia>sp_tr_data = create_transport_data(gasphase,"transport.dat")
julia>thermo_all = IdealGas.create_thermo(gasphase,"therm.dat")
julia>Dkm = zeros(length(gasphase))
julia>Dij = D_ij(sp_tr_data,T,p,thermo_all.molwt)
julia>D_km!(Dkm, bdc, mole_fracs, thermo_all.molwt)
```

###  Thermal conductivity
```julia
julia>using TransportProperties, IdealGas
julia>gasphase = ["CH4", "CO2", "H2O", "H2", "CO"]
julia>sp_tr_data = create_transport_data(gasphase,"transport.dat")
julia>thermo_all = IdealGas.create_thermo(gasphase,"therm.dat")
julia>tc = thermal_coductivity(sp_tr_data,T,p,thermo_all,mole_fracs)
```

In the above calls *mole_fracs* is an array of mole fractions; must be of same size as the number of 
gasphase species. *T* and *p* are respectively the temperature (K) and pressure (Pa)

## Output
The method transport_properties creates a screen output
```
Pure species viscosity:
     Species 	 viscosity(Kg/m-s)
-----------------------------------
         CH4 	      2.8992e-05
         H2O 	      3.8879e-05
          H2 	      2.0596e-05
          CO 	      4.2735e-05
         CO2 	      4.3121e-05
          O2 	      5.0182e-05
          N2 	      4.3463e-05
Mixture viscosity: 4.1079e-05 Kg/m-s

Mixture diffusion coefficients:
     Species 	 Diff.Coeff(m^2/s)
-----------------------------------
         CH4 	      2.2091e-04
         H2O 	      2.5665e-04
          H2 	      6.8730e-04
          CO 	      1.9235e-04
         CO2 	      1.4955e-04
          O2 	      1.9554e-04
          N2 	      1.8077e-04

Binary diffusion coefficients:
Note: self diffusion coefficients are not printed below:
-----------------------------------
	        CH4	        H2O	         H2	         CO	        CO2	         O2	         N2
         CH4	0.0000e+00	2.4993e-04	6.3858e-04	2.0108e-04	1.6734e-04	2.0549e-04	2.0277e-04	
         H2O	2.4993e-04	0.0000e+00	8.2047e-04	2.3675e-04	1.8888e-04	2.4327e-04	2.3919e-04	
          H2	6.3858e-04	8.2047e-04	0.0000e+00	6.6030e-04	5.8600e-04	6.9276e-04	6.6647e-04	
          CO	2.0108e-04	2.3675e-04	6.6030e-04	0.0000e+00	1.4726e-04	1.8413e-04	1.8318e-04	
         CO2	1.6734e-04	1.8888e-04	5.8600e-04	1.4726e-04	0.0000e+00	1.4794e-04	1.4849e-04	
          O2	2.0549e-04	2.4327e-04	6.9276e-04	1.8413e-04	1.4794e-04	0.0000e+00	1.8572e-04	
          N2	2.0277e-04	2.3919e-04	6.6647e-04	1.8318e-04	1.4849e-04	1.8572e-04	0.0000e+00	

Thermal conductivity of mixture: 8.1707e-02 (W/m-K)
```
