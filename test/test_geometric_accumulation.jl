using Test, Chion

@testset "Geometric layers preserve stratification during zero-input accumulation" begin
    model = BESSIModel(SnowpackGrid(1); Ntot=8,
        near_surface_layer_max_thicknesses_m=(0.02, 0.05, 0.10, 0.30))
    state = BESSIState(model)
    state.N[1] = 5
    state.mass[1:5, 1] .= [6., 15., 30., 90., 300.]
    state.density[1:5, 1] .= 300.
    state.temperature[1:5, 1] .= [250., 255., 260., 265., 270.]
    before = copy(state.temperature)
    mass_before = copy(state.mass)
    Chion._apply_accumulation_resolved!(state.N, state.mass, state.mass_w,
        state.density, state.temperature, state.mass_base, state.smb_ice,
        state.runoff, state.Tsrf, state.albedo, 1, model.c, model.Ntot,
        model.mass_max, model.mass_split, 0., 0., 0., 86400., 250., 5.)
    Chion._remesh_near_surface_layers!(state.N, state.mass, state.mass_w,
        state.density, state.temperature,
        1, model.Ntot, model.near_surface_layer_max_thicknesses_m, model.c)
    @test state.N[1] == 5
    @test state.mass ≈ mass_before
    @test state.temperature ≈ before
end
