@testset "AD compatibility (ForwardDiff)" begin
    # ── Setup ─────────────────────────────────────────────────────────────────
    # gibbs_hessian_diag is ForwardDiff-compatible
    n0 = [0.5, 0.3, 0.2]
    g_h = ForwardDiff.jacobian(n -> Optima.gibbs_hessian_diag(n), n0)
    @test all(isfinite, g_h)

    # ── kkt_residual is ForwardDiff-compatible through n ──────────────────────
    μ⁰ = [0.0, 1.0, 2.0]
    A = ones(Float64, 1, 3)
    b = ones(Float64, 1)

    function G(n, p)
        return sum(n[i] * (p.μ⁰[i] + log(n[i])) for i in eachindex(n))
    end
    function ∇G!(grad, n, p)
        for i in eachindex(n)
            grad[i] = p.μ⁰[i] + log(n[i]) + one(eltype(n))
        end
    end

    prob = OptimaProblem(A, b, G, ∇G!; lb = fill(1.0e-16, 3), p = (μ⁰ = μ⁰,))

    n_test = [0.6, 0.3, 0.1]
    y_test = [0.5]

    jac_ex = ForwardDiff.jacobian(
        n -> begin
            g = similar(n)
            prob.g!(g, n, prob.p)
            Optima.kkt_residual(prob, n, y_test, g, 1.0e-4).ex
        end,
        n_test,
    )
    @test all(isfinite, jac_ex)

    # ── hessian_diagonal is ForwardDiff-compatible ────────────────────────────
    h_jac = ForwardDiff.jacobian(
        n -> begin
            hf = Optima.gibbs_hessian_diag(n)
            Optima.hessian_diagonal(prob, n, 1.0e-4, hf)
        end,
        n_test,
    )
    @test all(isfinite, h_jac)

    # ── sensitivity matrices are finite ───────────────────────────────────────
    result = solve(prob, OptimaOptions(tol = 1.0e-12))
    @test result.converged

    n_sol = result.n
    hf_sol = Optima.gibbs_hessian_diag(n_sol)
    h_sol = Optima.hessian_diagonal(prob, n_sol, 1.0e-14, hf_sol)
    sens = sensitivity(prob, n_sol, result.y, h_sol, 1.0e-14)

    @test all(isfinite, sens.∂n_∂b)
    @test all(isfinite, sens.∂n_∂μ0)

    # ── objective function is ForwardDiff-compatible through μ⁰ ──────────────
    # Differentiate the objective value at the current point w.r.t. μ⁰
    g_obj = ForwardDiff.gradient(
        μ0 -> G(n_test, (μ⁰ = μ0,)),
        μ⁰,
    )
    @test all(isfinite, g_obj)
    @test g_obj ≈ n_test atol = 1.0e-12   # ∂G/∂μ⁰_i = n_i
end
