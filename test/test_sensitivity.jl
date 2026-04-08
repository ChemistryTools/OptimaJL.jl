@testset "Sensitivity ∂n*/∂b and ∂n*/∂μ⁰" begin
    # Same 1-row ideal Gibbs problem as test_solver
    μ⁰ = [0.0, 1.0, 2.0]

    function G(n, p)
        μ0 = p.μ⁰
        return sum(n[i] * (μ0[i] + log(n[i])) for i in eachindex(n))
    end

    function ∇G!(grad, n, p)
        μ0 = p.μ⁰
        for i in eachindex(n)
            grad[i] = μ0[i] + log(n[i]) + 1
        end
    end

    A = ones(1, 3)
    b = [1.0]
    prob = OptimaProblem(A, b, G, ∇G!; lb = fill(1.0e-16, 3), p = (μ⁰ = μ⁰,))
    result = solve(prob, OptimaOptions(tol = 1.0e-12))
    @test result.converged

    # Build Hessian diagonal at the solution
    n = result.n
    hf = OptimaKit.gibbs_hessian_diag(n)
    # barrier at convergence: μ is very small, use eps as proxy
    μ_conv = result.iterations > 0 ? 1.0e-14 : 1.0e-14
    h = OptimaKit.hessian_diagonal(prob, n, μ_conv, hf)

    sens = sensitivity(prob, n, result.y, h, μ_conv)

    @test size(sens.∂n_∂b) == (3, 1)
    @test size(sens.∂n_∂μ0) == (3, 3)
    @test all(isfinite, sens.∂n_∂b)
    @test all(isfinite, sens.∂n_∂μ0)

    # ∂n*/∂b: increasing total moles should increase all nᵢ proportionally
    @test all(sens.∂n_∂b .> 0)
    @test sum(sens.∂n_∂b[:, 1]) ≈ 1.0 atol = 1.0e-8  # Σ ∂nᵢ/∂b = 1

    # Finite-difference check for ∂n*/∂b
    δb = 1.0e-5
    prob_pert = OptimaProblem(A, b .+ δb, G, ∇G!; lb = fill(1.0e-16, 3), p = (μ⁰ = μ⁰,))
    result_pert = solve(prob_pert, OptimaOptions(tol = 1.0e-12); u0 = result)
    ∂n_∂b_fd = (result_pert.n .- result.n) ./ δb
    @test sens.∂n_∂b[:, 1] ≈ ∂n_∂b_fd atol = 1.0e-5

    # Finite-difference check for ∂n*/∂μ⁰₁
    δμ = 1.0e-5
    μ⁰_pert = copy(μ⁰); μ⁰_pert[1] += δμ
    prob_μ = OptimaProblem(A, b, G, ∇G!; lb = fill(1.0e-16, 3), p = (μ⁰ = μ⁰_pert,))
    result_μ = solve(prob_μ, OptimaOptions(tol = 1.0e-12); u0 = result)
    ∂n_∂μ0_fd1 = (result_μ.n .- result.n) ./ δμ
    @test sens.∂n_∂μ0[:, 1] ≈ ∂n_∂μ0_fd1 atol = 1.0e-5
end
