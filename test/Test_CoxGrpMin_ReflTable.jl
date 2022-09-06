using Test
import CoxeterGroups: QInt, create_reflection_table_coxeter, create_reflection_table_gcm, cartan_A, cartan_B, cartan_C, cartan_D, cartan_E, cartan_F4, cartan_G2

@testset "QInt operations" begin
    # Test minimal polynomials for those QInt(2m, 2) such
    # that the min poly is quadratic.
    m3 = QInt(6, 2)
    m4 = QInt(8, 2)
    m5 = QInt(10, 2)
    m6 = QInt(12, 2)

    @test QInt(0) == m3 - QInt(1)
    @test QInt(0) == m4 * m4 - QInt(2)
    @test QInt(0) == m5 * m5 - m5 - QInt(1)
    @test QInt(0) == m6 * m6 - QInt(3)
end

@testset "Minimal root table A$rank" for rank in 0:10
    gcm = cartan_A(rank)
    @test create_reflection_table_coxeter(gcm_to_coxeter_matrix(gcm)) == create_reflection_table_gcm(gcm)
end

@testset "Minimal root table B$rank" for rank in 2:10
    gcm = cartan_B(rank)
    @test create_reflection_table_coxeter(gcm_to_coxeter_matrix(gcm)) == create_reflection_table_gcm(gcm)
end

@testset "Minimal root table C$rank" for rank in 2:10
    gcm = cartan_C(rank)
    @test create_reflection_table_coxeter(gcm_to_coxeter_matrix(gcm)) == create_reflection_table_gcm(gcm)
end

@testset "Minimal root table D$rank" for rank in 2:10
    gcm = cartan_D(rank)
    @test create_reflection_table_coxeter(gcm_to_coxeter_matrix(gcm)) == create_reflection_table_gcm(gcm)
end

@testset "Minimal root table E$rank" for rank in 6:8
    gcm = cartan_E(rank)
    @test create_reflection_table_coxeter(gcm_to_coxeter_matrix(gcm)) == create_reflection_table_gcm(gcm)
end

@testset "Minimal root table F4" begin
    gcm = cartan_F4()
    @test create_reflection_table_coxeter(gcm_to_coxeter_matrix(gcm)) == create_reflection_table_gcm(gcm)
end

@testset "Minimal root table G2" begin
    gcm = cartan_G2()
    @test create_reflection_table_coxeter(gcm_to_coxeter_matrix(gcm)) == create_reflection_table_gcm(gcm)
end

# Explicit reflection tables below here are copied in from the PrintReflectionTables.m Magma script.
@testset "Minimal root table H3" begin
    coxeter_mat = [
        1 5 2
        5 1 3
        2 3 1
    ]
    refl_tab = [
        [ 0  5  3  7  2  9  4 10  6  8 13 14 11 12 15]
        [ 4  0  6  1  7  3  5  8 11 12  9 10 13 15 14]
        [ 1  6  0  8  9  2 10  4  5  7 12 11 14 13 15]
    ]
    @test create_reflection_table_coxeter(coxeter_mat) == refl_tab
end

@testset "Minimal root table H4" begin
    coxeter_mat = [
        1 5 2 2
        5 1 3 2
        2 3 1 3
        2 2 3 1
    ]
    refl_tab = [
        [ 0  6  3  4  9  2 11  8  5 13  7 16 10 18 19 12 21 14 15 23 17 26 20 28 25 22 30 24 29 27 34 32 36 31 38 33 37 35 39 40 44 42 43 41 47 46 45 50 51 48 49 53 52 55 54 57 56 58 59 60]
        [ 5  0  7  4  1  9  3 12  6 10 15  8 17 14 11 20 13 22 19 16 25 18 23 24 21 29 27 31 26 33 28 35 30 37 32 39 34 41 36 40 38 45 43 44 42 48 49 46 47 52 51 50 53 54 56 55 58 57 59 60]
        [ 1  7  0  8 10 11  2  4 13  5  6 12  9 14 17 16 15 18 21 24 19 27 28 20 25 30 22 23 32 26 31 29 35 34 33 38 40 36 42 37 45 39 46 47 41 43 44 48 49 50 51 54 55 52 53 56 57 59 58 60]
        [ 1  2  8  0  5  6 12  3  9 14 16  7 18 10 20 11 22 13 23 15 26 17 19 27 29 21 24 30 25 28 33 32 31 36 35 34 39 38 37 43 41 46 40 44 48 42 50 45 52 47 53 49 51 54 55 56 57 58 60 59]
    ]
    @test create_reflection_table_coxeter(coxeter_mat) == refl_tab
end

@testset "Minimal root table 353" begin
    coxeter_mat = [
        [1 3 2 2]
        [3 1 5 2]
        [2 5 1 3]
        [2 2 3 1]
    ]
    refl_tab = [
        [ 0  5  3  4  2  9 12  8  6 16 15  7 18 20 11 10 23 13 25 14 27 28 17 29 19  0 21 22 24  0 31  0]
        [ 5  0  7  4  1 10  3 13 14  6 17 12  8  9 22 20 11 18 26 16 21 15 28 30  0 19 31 23  0 24 27 32]
        [ 1  6  0  8  9  2 10  4  5  7 11 16 19 21 15 12 24 25 13 27 14  0 29 17 18 26 20  0 23 32 31 30]
        [ 1  2  8  0  5 11 13  3 15 17  6 18  7 22  9 23 10 12 24 28  0 14 16 19 29 30  0 20 25 26  0 32]
    ]
    @test create_reflection_table_coxeter(coxeter_mat) == refl_tab
end
