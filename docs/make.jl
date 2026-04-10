import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Documenter, Chion

makedocs(
    modules=[Chion],
    format=Documenter.HTML(edit_link=nothing),
    sitename="Chion.jl Documentation",
    checkdocs=:exports,
    pages=[
        "Home" => "index.md",
        "Model State And Step Flow" => "model_state.md",
        "Processes" => [
            "Albedo" => "processes/albedo.md",
            "Accumulation And Melt" => "processes/accumulation_ablation.md",
            "Layer Structure And Basal Transfer" => "processes/layer_structure.md",
            "Densification" => "processes/densification.md",
            "Energy Balance" => "processes/energy.md",
            "Percolation" => "processes/percolation.md",
            "Refreezing" => "processes/refreezing.md",
        ],
        "Case API And Outputs" => "case_api.md",
        "Reference Utilities" => "reference.md",
        "Validation And Audit" => [
            "Formula Audit" => "validation.md",
            "Fresh-Snow Density Test" => "tests/fresh_snow_density.md",
            "Energy Flux Analytical Tests" => "tests/test_energy_flux_analytical.md",
        ],
    ],
)
