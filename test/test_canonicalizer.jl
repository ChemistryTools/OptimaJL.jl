@testset "Canonicalizer" begin
    # 2 elements, 4 species: H₂O, H⁺, OH⁻, H₂
    # Conservation: [H, O] — rows correspond to elements
    A = Float64[
        2  1  1  2   # H: 2 in H₂O, 1 in H⁺, 1 in OH⁻, 2 in H₂
        1  0  1  0   # O: 1 in H₂O, 0 in H⁺, 1 in OH⁻, 0 in H₂
    ]

    can = Canonicalizer(A)

    @test can.m == 2
    @test can.ns == 4
    @test can.rank_A == 2
    @test length(can.jb) == 2
    @test length(can.jn) == 2

    # Basic columns must span the row space
    @test rank(A[:, can.jb]) == 2

    # Permutation consistency
    @test sort(vcat(can.jb, can.jn)) == 1:4

    # Schur complement is symmetric positive semi-definite
    h = ones(4)
    S = OptimaKit.schur_complement(can, h)
    @test S ≈ S'
    @test all(eigvals(S) .>= -1.0e-12)

    # solve_B and solve_Bt: round-trip
    rhs = ones(2)
    x = OptimaKit.solve_B(can, rhs)
    @test A[:, can.jb] * x ≈ rhs atol = 1.0e-12

    xt = OptimaKit.solve_Bt(can, rhs)
    @test A[:, can.jb]' * xt ≈ rhs atol = 1.0e-12
end
