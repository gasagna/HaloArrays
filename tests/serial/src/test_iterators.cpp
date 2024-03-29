#include "DArrays.hpp"
#include "catch.hpp"
#include <iostream>

TEST_CASE("Testing linear range", "[LinRange]") {
    int i = 0;
    
    SECTION("case 1") {
        auto r = DArrays::Iterators::LinRange(5);
        std::array<int, 5> exact = {0, 1, 2, 3, 4};
        for (auto val : r ) {
            REQUIRE( val == exact[i++] );
        }
        REQUIRE(i == 5);
    }

    SECTION("case 2") {
        auto r = DArrays::Iterators::LinRange(0, 5);
        std::array<int, 5> exact = {0, 1, 2, 3, 4};
        for (auto val : r) {
            REQUIRE(val == exact[i++]);
        }
        REQUIRE(i == 5);
    }

    SECTION("case 3") {
        auto r = DArrays::Iterators::LinRange(0, 5, 1);
        std::array<int, 5> exact = {0, 1, 2, 3, 4};
        for (auto val : r) {
            REQUIRE(val == exact[i++]);
        }
        REQUIRE(i == 5);
    }

    SECTION("case 4") {
        auto r = DArrays::Iterators::LinRange(0, 6, 2);
        std::array<int, 3> exact = {0, 2, 4};
        for (auto val : r) {
            REQUIRE(val == exact[i++]);
        }
        REQUIRE(i == 3);
    }

    SECTION("case 5") {
        auto r = DArrays::Iterators::LinRange(0, 7, 2);
        std::array<int, 4> exact = {0, 2, 4, 6};
        for (auto val : r) {
            REQUIRE(val == exact[i++]);
        }
        REQUIRE(i == 4);
    }

    SECTION("case 6") {
        auto r = DArrays::Iterators::LinRange(-1, 2, 1);
        std::array<int, 3> exact = {-1, 0, 1};
        for (auto val : r) {
            REQUIRE(val == exact[i++]);
        }
        REQUIRE(i == 3);
    }

    SECTION("case 7") {
        auto r = DArrays::Iterators::LinRange(1, -2, -1);
        std::array<int, 3> exact = {1, 0, -1};
        for (auto val : r) {
            REQUIRE(val == exact[i++]);
        }
        REQUIRE(i == 3);
    }
    
    SECTION("case 8") {
        auto r = DArrays::Iterators::LinRange(-2, 0, 1);
        std::array<int, 2> exact = {-2, -1};
        for (auto val : r) {
            REQUIRE(val == exact[i++]);
        }
        REQUIRE(i == 2);
    }
}

TEST_CASE("Testing index range", "[IndexRange]") {
    int i = 0;
    
    SECTION("values") {
        SECTION("case 1d - array") {
            std::array<std::array<int, 1>, 5> exact = {{{0}, {1}, {2}, {3}, {4}}};
            std::array<int, 1> size = {5};
            for (auto val : DArrays::Iterators::IndexRange<1>(size)) {
                REQUIRE(val == exact[i++]);
            }
            REQUIRE(i == 5);
        }

        SECTION("case 1d - integers") {
            std::array<std::array<int, 1>, 5> exact = {{{0}, {1}, {2}, {3}, {4}}};
            for (auto val : DArrays::Iterators::IndexRange<1>(5)) {
                REQUIRE(val == exact[i++]);
            }
            REQUIRE(i == 5);
        }

        SECTION("case 2d - array") {
            SECTION("test 1") {           
                std::array<std::array<int, 2>, 6> exact = {{{0, 0}, {1, 0}, {0, 1}, {1, 1}, {0, 2}, {1, 2}}};
                std::array<int, 2> size = {2, 3};
                for (auto val : DArrays::Iterators::IndexRange<2>(size)) {
                    REQUIRE(val == exact[i++]);
                }
                REQUIRE(i == 6);
            }

            SECTION("test 2") {           
                auto b = DArrays::Iterators::IndexRange<2>(3, 4).begin();
                for (int j = 0; j < 4; j++) {
                    for (int i = 0; i < 3; i++) {
                        REQUIRE( (*b)[0] == i );
                        REQUIRE( (*b)[1] == j );
                        b += 1;
                    }
                }
            }
        }

        SECTION("case 2d - integers") {
            std::array<std::array<int, 2>, 6> exact = {{{0, 0}, {1, 0}, {0, 1}, {1, 1}, {0, 2}, {1, 2}}};
            for (auto val : DArrays::Iterators::IndexRange<2>(2, 3)) {
                REQUIRE(val == exact[i++]);
            }
            REQUIRE(i == 6);
        }

        SECTION("case 3d - array") {
            SECTION("test 1") {               
                std::array<std::array<int, 3>, 18> exact = {{{0, 0, 0}, {1, 0, 0}, {2, 0, 0},
                                                             {0, 1, 0}, {1, 1, 0}, {2, 1, 0}, 
                                                             {0, 2, 0}, {1, 2, 0}, {2, 2, 0}, 
                                                             {0, 0, 1}, {1, 0, 1}, {2, 0, 1},
                                                             {0, 1, 1}, {1, 1, 1}, {2, 1, 1}, 
                                                             {0, 2, 1}, {1, 2, 1}, {2, 2, 1}}};
                for (auto val : DArrays::Iterators::IndexRange<3>(3, 3, 2)) {
                    REQUIRE(val == exact[i++]);
                }
                REQUIRE(i == 18);
            }

            SECTION("test 2") {           
                auto b = DArrays::Iterators::IndexRange<3>(3, 4, 5).begin();
                for (int k = 0; k < 5; k++) {
                    for (int j = 0; j < 4; j++) {
                        for (int i = 0; i < 3; i++) {
                            REQUIRE( (*b)[0] == i );
                            REQUIRE( (*b)[1] == j );
                            REQUIRE( (*b)[2] == k );
                            b += 1;
                        }
                    }
                }
            }
        }
    }

    SECTION("random access iterator interface") {
        DArrays::Iterators::IndexRange<3> rng(3, 4, 5);
        auto b = rng.begin();

        SECTION("in place addition/subtraction") {
            SECTION("case - 1") {
                b += 1;
                std::array<int, 3> expected = {1, 0, 0};
                REQUIRE( *b == expected );
            }

            SECTION("case - 2") {
                b += 4;
                std::array<int, 3> expected = {1, 1, 0};
                REQUIRE( *b == expected );
            }

            SECTION("case - 3") {
                b += 12;
                std::array<int, 3> expected = {0, 0, 1};
                REQUIRE( *b == expected );
            }
            
            SECTION("case - 4") {
                b += 15;
                std::array<int, 3> expected = {0, 1, 1};
                REQUIRE( *b == expected );
            }

            SECTION("case - 5") {
                b += 3;
                b -= 3;
                std::array<int, 3> expected = {0, 0, 0};
                REQUIRE( *b == expected );
            }

            SECTION("case - 6") {
                b += 1;
                b += 1;
                b += 1;
                b += 1;
                std::array<int, 3> expected = {1, 1, 0};
                REQUIRE( *b == expected );
            }
            
            SECTION("case - 7") {
                b += 2;
                b += 2;
                b += 2;
                b += 2;
                b += 2;
                std::array<int, 3> expected = {1, 3, 0};
                REQUIRE( *b == expected );
            }
            
            SECTION("case - 8") {
                b += 12;
                std::array<int, 3> expected = {0, 0, 1};
                REQUIRE( *b == expected );

                b += 12;
                expected = {0, 0, 2};
                REQUIRE( *b == expected );
            }
        }

        SECTION("comparison") {
            auto c = rng.begin(); c+= 1;
            REQUIRE( (c == b) == false );
            REQUIRE( (c >  b) == true );
            REQUIRE( (c <  b) == false );
            REQUIRE( (c != b) == true );
        }
    }
}

