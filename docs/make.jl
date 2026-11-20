using Documenter, Chion

DocMeta.setdocmeta!(Chion, :DocTestSetup, :(using Chion); recursive=true)

makedocs(
    modules=[Chion],
    format=Documenter.HTML(edit_link=nothing),
    sitename="Chion.jl Documentation",
    checkdocs=:exports,
    warnonly=[:docs_block, :missing_docs],
    pages=[
        "Home" => "index.md",
        "Running Simulations" => "simulations.md",
        "Models And Coupling" => "models.md",
        "Model State And Step Flow" => "model_state.md",
        "Processes" => [
            "Albedo" => "processes/albedo.md",
            "Accumulation And Melt" => "processes/accumulation_ablation.md",
            "Layer Structure And Basal Transfer" => "processes/layer_structure.md",
            "Densification" => "processes/densification.md",
            "Energy Balance" => "processes/energy.md",
            "Diurnal Shortwave Cycle" => "processes/diurnal_cycle.md",
            "Percolation" => "processes/percolation.md",
            "Refreezing" => "processes/refreezing.md",
        ],
        "Parameters" => "parameters.md",
        "Verification" => [
            "Verification And Validation" => "validation.md",
            "Fresh-Snow Density Test" => "tests/fresh_snow_density.md",
            "Energy Flux Analytical Tests" => "tests/test_energy_flux_analytical.md",
        ],
    ],
)

deploydocs(
    repo = "github.com/fesmc/Chion.jl.git",
)
