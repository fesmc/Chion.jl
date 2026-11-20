# Shared setup of the test suite. Each test file can also run on its own, from the package root:
#   julia --project=test -e 'include("test/setup.jl"); include("test/test_simulation.jl")'
using Test
using Dates
using NCDatasets
using CUDA      # activates ChionCUDAExt; the GPU tests run only on a functional device
using Chion
