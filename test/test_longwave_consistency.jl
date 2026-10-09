using Test, Chion

@testset "Thin-layer Robin boundary preserves insulated equilibrium" begin
    m=BESSIModel(SnowpackGrid(1);Ntot=1,D_sh=0.,ϵ_snow=0.,ice_substrate_layers=0)
    for mass in (100.,1e-8)
        s=BESSIState(m);s.N[1]=1;s.mass[1,1]=mass;s.density[1,1]=300.
        s.temperature[1,1]=260.;s.Tsrf[1]=260.
        scratch=Chion.EnergyWorkspace(s.mass,Float64,1,1)
        e=Chion.go_energy_flux!(s,1,260.,0.,0.,0.,10800.;scratch,
            q_sw_net=0.,q_lw_down=0.,q_sh=0.,q_lh=0.)
        @test s.Tsrf[1] ≈ 260. atol=1e-10
        @test s.temperature[1,1] ≈ 260. atol=1e-10
        @test scratch.diag[1,1]==1.
        @test e.melt_energy_available==0.
    end
end

@testset "Graybody longwave equilibrium with prescribed incident radiation" begin
    gray=BESSIModel(SnowpackGrid(1);Ntot=2,ϵ_snow=.98,eps_ice=.98,ice_substrate_layers=0)
    legacy=BESSIModel(SnowpackGrid(1);seb_scheme=:bessi,Ntot=2,ϵ_snow=.98,eps_ice=.98,ice_substrate_layers=0)
    c=gray.c;dt=10800.;incident=c.σ*c.T0^4
    f=Chion.SnowpackStepForcing(c.T0,.125,0.,0.,0.,0.;
        q_sw_net=0.,has_q_sw_net=true,q_lw_down=incident,has_q_lw_down=true,
        q_sh=0.,has_q_sh=true,q_lh=0.,has_q_lh=true)
    # An isothermal blackbody radiation environment cannot melt gray ice:
    # absorbed incident LW and emitted LW must cancel, regardless of epsilon.
    bare=Chion._bare_ice_ablation_mass(c,f,dt)
    @test abs(bare.longwave_flux)<1e-10
    @test abs(bare.melt_mass)<1e-10
    old=Chion._bare_ice_ablation_mass(legacy.c,f,dt)
    @test old.longwave_flux ≈ (1-c.ϵ_snow)*incident
    @test old.melt_mass>0
    s=BESSIState(gray);s.N[1]=2;s.mass[:,1].=(100.,300.)
    s.density[:,1].=(300.,500.);s.temperature[:,1].=c.T0;s.Tsrf[1]=c.T0
    scratch=Chion.EnergyWorkspace(s.mass,Float64,2,1)
    energy=Chion.go_energy_flux!(s,1,c.T0,0.,0.,0.,dt;scratch,
        q_sw_net=0.,q_lw_down=incident,q_sh=0.,q_lh=0.)
    @test s.Tsrf[1] ≈ c.T0 atol=1e-10
    @test all(abs.(s.temperature[:,1].-c.T0).<1e-10)
    @test abs(energy.melt_energy_available)<1e-5
end
