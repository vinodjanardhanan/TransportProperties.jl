using TransportProperties
using Test
using LightXML, RxnHelperUtils, IdealGas
@testset "TransportProperties.jl" begin

    if Sys.isapple() || Sys.islinux()
        lib_dir = "lib/"
    else
        lib_dir = "lib\\"
    end

    xmldoc = parse_file("transport.xml")
    xmlroot = root(xmldoc)

    #get the gasphase species
    gasphase = get_collection_from_xml(xmlroot,"gasphase")
    thermo_file = joinpath(lib_dir,"therm.dat") 
    thermo_all = IdealGas.create_thermo(gasphase,thermo_file)
    mole_fracs = get_molefraction_from_xml(xmlroot,thermo_all.molwt,gasphase)
    local T = get_value_from_xml(xmlroot,"T")
    local p = get_value_from_xml(xmlroot,"p")
    trans_file =  joinpath(lib_dir, "transport.dat") 
    sp_tr_data = create_transport_data(gasphase,trans_file)

    @testset "Testing pure species viscosity " begin
        #Evaluate viscosity
        μ =  viscosity(sp_tr_data,T,thermo_all.molwt)    
        @test 1e-6 < μ[1] < 1e-4
    end

    @testset "Testing mixture viscosity " begin                
        μ = viscosity(sp_tr_data,T,thermo_all.molwt,mole_fracs)    
        @test 1e-6 < μ < 1e-4
    end
    
    @testset "Testing mixture diffusion coefficients method-1" begin
        #Evaluate diffusion coefficients
        Dkm = zeros(length(gasphase))
        D_km!(Dkm,sp_tr_data,T,p,thermo_all.molwt,mole_fracs)        
        @test 1e-4 < Dkm[1] < 3e-4
    end

    @testset "Testing binary diffusion coefficients" begin
        bdc =  D_ij(sp_tr_data,T,p,thermo_all.molwt)
        @test 2.1000e-04 < bdc[1] < 2.2407e-04
    end

    @testset "Testing mixture diffusion coefficients method-2" begin
        bdc =  D_ij(sp_tr_data,T,p,thermo_all.molwt)
        Dkm = zeros(length(gasphase))
        D_km!(Dkm, bdc, mole_fracs, thermo_all.molwt)
        @test 1e-4 < Dkm[1] < 3e-4
    end


    @testset "Testing thermal conductivity " begin
        λ = thermal_coductivity(sp_tr_data,T,p,thermo_all,mole_fracs)
        @test ceil(λ) == 1.0
    end

    @testset "Testing all " begin
        retcode = transport_properties("transport.xml", lib_dir)
        @test Symbol("Success") == retcode
    end

end
