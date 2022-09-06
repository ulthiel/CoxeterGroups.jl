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

@testset "Minimal root table rand" begin
    coxeter_mat = [
        [1 3 4 2 2 2 2]
        [3 1 3 3 2 2 2]
        [4 3 1 9 2 2 2]
        [2 3 9 1 8 3 2]
        [2 2 2 8 1 3 2]
        [2 2 2 3 3 1 5]
        [2 2 2 2 2 5 1]
    ]
    refl_tab = [
        [0 8 12 4 5 6 7 2 9 0 21 3 24 0 15 16 17 18 19 20 11 38 39 13 0 41 0 0 0 30 31 48 33 34 35 36 37 22 23 55 26 0 0 0 0 46 47 32 61 50 51 52 53 65 40 66 0 0 59 60 49 62 63 73 54 56 74 75 0 78 71 72 64 67 68 84 83 70 87 80 81 89 77 76 92 91 79 95 82 97 86 85 100 99 88 104 90 105 94 93 109 110 108 96 98 115 116 103 101 102 121 122 119 114 106 107 127 128 113 120 111 112 130 131 0 126 117 118 129 123 124 137 138 0 135 136 132 133 0 142 141 140 0 145 144 0 148 147 0 0 0 0 0 0 0 0]
        [8 0 10 11 5 6 7 1 0 3 4 0 0 0 22 23 32 18 19 20 21 15 16 0 0 0 0 0 0 0 40 17 49 34 35 54 37 38 39 31 0 0 0 0 0 0 56 48 33 0 51 64 53 36 55 47 0 0 70 68 61 62 63 52 65 66 76 60 0 59 79 72 73 82 75 67 85 78 71 88 81 74 90 89 77 93 87 80 84 83 98 97 86 102 95 96 92 91 107 105 101 94 112 114 100 106 99 118 120 116 111 103 125 104 126 110 117 108 0 109 129 128 123 134 113 115 135 122 121 136 0 139 133 124 127 130 0 141 132 143 138 0 140 146 0 144 150 0 149 147 151 152 153 154 156 155]
        [9 10 0 14 5 6 7 0 1 2 0 12 25 4 28 29 0 18 19 20 0 0 0 0 13 43 44 15 16 0 45 0 0 34 35 0 37 0 0 0 0 57 26 27 31 0 58 0 0 0 51 0 53 0 0 0 42 47 0 69 0 62 63 0 0 0 0 0 60 0 0 72 0 0 0 0 0 0 0 0 81 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [1 11 13 0 17 16 7 21 0 0 2 24 3 27 30 6 5 0 31 36 8 0 23 12 42 26 14 0 0 15 19 32 50 0 52 20 0 0 39 40 41 25 0 57 0 46 59 48 0 33 0 35 0 54 55 67 44 0 47 71 0 0 0 64 65 74 56 77 0 76 60 0 73 66 83 70 68 84 85 80 0 82 75 78 79 86 92 96 89 90 91 87 101 94 104 88 97 106 99 109 93 111 113 95 115 98 117 119 100 121 102 124 103 114 105 127 107 0 108 126 110 131 132 112 134 120 116 0 135 137 122 123 140 125 129 0 130 142 139 133 0 138 143 144 145 149 147 148 146 151 150 152 153 155 154 156]
        [1 2 3 15 0 18 7 8 9 10 22 12 0 28 4 0 33 6 34 37 38 11 0 0 0 0 0 14 0 46 0 49 17 19 53 0 20 21 0 0 0 0 0 0 0 30 0 61 32 50 62 0 35 0 0 0 0 0 0 0 48 51 72 0 0 0 0 0 0 0 0 63 0 0 0 0 0 0 0 0 81 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [1 2 3 16 18 0 20 8 9 10 23 12 26 29 0 4 0 5 35 7 39 0 11 41 43 13 0 0 14 0 47 0 0 51 19 36 37 0 21 56 24 0 25 0 58 0 31 0 0 0 34 59 62 54 66 40 0 45 52 60 0 53 63 70 65 55 67 68 69 64 80 81 78 74 75 76 86 73 88 71 72 82 91 84 93 77 95 79 89 98 83 100 85 103 87 101 105 90 108 92 96 112 94 109 97 106 118 99 104 122 123 102 113 120 115 128 0 107 119 114 130 110 111 132 125 126 0 116 136 121 137 124 133 139 0 129 131 138 134 144 141 145 146 140 142 143 147 148 149 150 152 151 154 153 155 156]
        [1 2 3 4 5 19 0 8 9 10 11 12 13 14 15 31 17 34 6 35 21 22 40 24 25 0 27 28 45 30 16 32 33 18 20 52 53 38 55 23 0 42 0 44 29 46 60 48 49 50 63 36 37 64 39 68 57 69 71 47 61 72 51 54 73 75 77 56 58 79 59 62 65 83 66 85 67 87 70 80 81 90 74 92 76 94 78 88 97 82 99 84 102 86 95 96 89 107 91 110 111 93 103 104 116 117 98 108 121 100 101 112 113 114 127 105 106 118 119 129 109 122 133 124 125 135 115 128 120 138 131 140 123 134 126 141 142 130 143 132 136 137 139 147 148 150 144 145 151 146 149 153 152 154 155 156]
    ]
    @test create_reflection_table_coxeter(coxeter_mat) == refl_tab
end
