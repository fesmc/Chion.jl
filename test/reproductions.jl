# The named setups against pinned totals of one MAR year from the initial state (no spin-up). The
# values are those of the code that produced the paper (bitwise on the CPU); skipped where the MAR
# files are not available.
using Test
using Chion

const _MAR_GRIS = "/p/projects/ou/labs/ai/Nils/MAR_daily/mar_daily/MARv3.14.3-10km-daily-ERA5-{year}.nc"
const _MAR_AIS = "/p/projects/ou/labs/ai/Nils/MAR_daily/mar_daily_antarctica/MARv3.14.3-27.5km-daily-ERA5-{year}.nc"
_mar_file(pattern, year) = replace(pattern, "{year}" => string(year))

@testset "Reproductions" begin
    if isfile(_mar_file(_MAR_GRIS, 1980))
        sim = SimulationGreenlandMARv1(_MAR_GRIS; years = 1980:1980, spinup = NoSpinup(), verbose = false)
        run!(sim)
        totals = smbstats(sim).total_Gt
        pinned = (smb_ice = -344.867614949, melt = 549.290053937, runoff = 480.006162542,
            refreezing = 92.2942549142, sublimation = 49.161123928)
        for name in keys(pinned)
            @test totals[name] ≈ pinned[name] rtol = 1e-10
        end
    else
        @info "Skipping the Greenland reproduction: no MAR files at $(_MAR_GRIS)."
    end
    if isfile(_mar_file(_MAR_AIS, 1980))
        sim = SimulationAntarcticaMARv1(_MAR_AIS, _mar_file(_MAR_AIS, 1979); years = 1980:1980, spinup_years = 0,
            area_file = _mar_file(_MAR_AIS, 2025), verbose = false)
        run!(sim)
        totals = smbstats(sim).total_Gt
        pinned = (smb_ice = -152.00110585, melt = 267.882189556, runoff = 199.924626571,
            refreezing = 57.8244884923, sublimation = 186.438177062)
        for name in keys(pinned)
            @test totals[name] ≈ pinned[name] rtol = 1e-10
        end
        @test sim.spinup === NoSpinup()
    else
        @info "Skipping the Antarctic reproduction: no MAR files at $(_MAR_AIS)."
    end
end
