#include "DArrays.hpp"
#include <catch.hpp>
#include <iostream>
#include <array>
#include <map>

// import all
using namespace DArrays;

TEST_CASE("haloregionspec", "test_1") {

    HaloRegionSpec<1> s1(Boundary::RIGHT);
    REQUIRE( s1.hash(HaloIntent::SEND) ==  -4  );
    REQUIRE( s1.hash(HaloIntent::RECV) ==  +4  );
    REQUIRE( s1[0] == Boundary::RIGHT );
    REQUIRE_THROWS( s1[-1] );
    REQUIRE_THROWS( s1[ 1] );

    HaloRegionSpec<1> os1 = opposite(s1);
    REQUIRE( os1[0] == Boundary::LEFT );
    REQUIRE( os1.hash(HaloIntent::SEND) ==  -1 );

    HaloRegionSpec<2> s2(Boundary::LEFT, Boundary::WILDCARD);
    REQUIRE( s2.hash(HaloIntent::SEND) ==  -81  );
    REQUIRE( s2.hash(HaloIntent::RECV) ==  +81  );
    REQUIRE( s2[0] == Boundary::LEFT );
    REQUIRE( s2[1] == Boundary::WILDCARD );
    REQUIRE_THROWS( s2[-1] );
    REQUIRE_THROWS( s2[ 2] );    

    HaloRegionSpec<2> os2 = opposite(s2);
    REQUIRE( os2[0] == Boundary::RIGHT );
    REQUIRE( os2[1] == Boundary::WILDCARD );
    REQUIRE( os2.hash(HaloIntent::SEND) ==  -84 );

    HaloRegionSpec<3> s3(Boundary::LEFT, Boundary::CENTER, Boundary::WILDCARD);
    REQUIRE( s3.hash(HaloIntent::SEND) ==  -821  );
    REQUIRE( s3.hash(HaloIntent::RECV) ==  +821  );
    REQUIRE_THROWS( s3[-1] );
    REQUIRE_THROWS( s3[ 3] );    

    HaloRegionSpec<3> os3 = opposite(s3);
    REQUIRE( os3[0] == Boundary::RIGHT );
    REQUIRE( os3[1] == Boundary::CENTER );
    REQUIRE( os3[2] == Boundary::WILDCARD );
    REQUIRE( os3.hash(HaloIntent::SEND) ==  -824 );

}