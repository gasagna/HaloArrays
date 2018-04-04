#pragma once
#include <array>
#include <tuple>

////////////////////////////////////////////////////////
//                     HALO REGION                    //
////////////////////////////////////////////////////////

namespace DArrays {

// tags for the the left, center and right boundaries
enum class HaloRegionTag: int {LEFT = 1, CENTER = 2, RIGHT = 3};

// tags to specificy whether we are pointing at the region
// that is supposed to be sent or the region that is supposed
// to be received
enum class HaloIntent : int {SEND = 0, RECV = 1};

// specification of the boundary region is an NDIMS array of HaloTags
template<size_t NDIMS>
using HaloRegion = std::array<HaloRegionTag, NDIMS>;

// THESE VECTORS SHOULD BE CONSTRUCTED PROGRAMMATICALLY WITH A CARTESIAN PRODUCT ITERATOR!

// 1D arrangement
// L X R                                                                         // message_tag
static const std::array<HaloRegion<1>, 2> _bnds_1d = {{{HaloRegionTag::LEFT},    // 10
                                                       {HaloRegionTag::RIGHT}}}; // 30

// 2D arrangement
// LL LC LR
// CL XX CR
// RL RC RR
static const std::array<HaloRegion<2>, 8> _bnds_2d = {{{HaloRegionTag::LEFT,   HaloRegionTag::LEFT},    // 10 + 100
                                                       {HaloRegionTag::LEFT,   HaloRegionTag::CENTER},  // 10 + 200
                                                       {HaloRegionTag::LEFT,   HaloRegionTag::RIGHT},   // 10 + 300
                                                       {HaloRegionTag::CENTER, HaloRegionTag::LEFT},    // 20 + 100
                                                       {HaloRegionTag::CENTER, HaloRegionTag::RIGHT},   // 20 + 300
                                                       {HaloRegionTag::RIGHT,  HaloRegionTag::LEFT},    // 30 + 100
                                                       {HaloRegionTag::RIGHT,  HaloRegionTag::CENTER},  // 30 + 200
                                                       {HaloRegionTag::RIGHT,  HaloRegionTag::RIGHT}}}; // 30 + 300 

// 3D arrangement
// LLL LCL LRL  // LLC LCC LRC  // LLR LCR LRR
// CLL CCL CRL  // CLC XXX CRC  // CLR CCR CRR
// RLL RCL RRL  // RLC RCC RRC  // RLR RCR RRR
static const std::array<HaloRegion<3>, 26> _bnds_3d = {{{HaloRegionTag::LEFT,   HaloRegionTag::LEFT,   HaloRegionTag::LEFT},    //  10 + 100 + 1000
                                                        {HaloRegionTag::CENTER, HaloRegionTag::LEFT,   HaloRegionTag::LEFT},    //  10 + 200 + 1000                            
                                                        {HaloRegionTag::RIGHT,  HaloRegionTag::LEFT,   HaloRegionTag::LEFT},    //  10 + 300 + 1000                            
                                                        {HaloRegionTag::LEFT,   HaloRegionTag::CENTER, HaloRegionTag::LEFT},    //  20 + 100 + 1000                           
                                                        {HaloRegionTag::CENTER, HaloRegionTag::CENTER, HaloRegionTag::LEFT},    //  20 + 200 + 1000                           
                                                        {HaloRegionTag::RIGHT,  HaloRegionTag::CENTER, HaloRegionTag::LEFT},    //  20 + 300 + 1000                           
                                                        {HaloRegionTag::LEFT,   HaloRegionTag::RIGHT,  HaloRegionTag::LEFT},    //  30 + 100 + 1000                            
                                                        {HaloRegionTag::CENTER, HaloRegionTag::RIGHT,  HaloRegionTag::LEFT},    //  30 + 200 + 1000                            
                                                        {HaloRegionTag::RIGHT,  HaloRegionTag::RIGHT,  HaloRegionTag::LEFT},    //  30 + 300 + 1000                            
                                                        {HaloRegionTag::LEFT,   HaloRegionTag::LEFT,   HaloRegionTag::CENTER},  //  10 + 100 + 2000                          
                                                        {HaloRegionTag::CENTER, HaloRegionTag::LEFT,   HaloRegionTag::CENTER},  //  10 + 200 + 2000                          
                                                        {HaloRegionTag::RIGHT,  HaloRegionTag::LEFT,   HaloRegionTag::CENTER},  //  10 + 300 + 2000                          
                                                        {HaloRegionTag::LEFT,   HaloRegionTag::CENTER, HaloRegionTag::CENTER},  //  20 + 100 + 2000                          
                                                        {HaloRegionTag::RIGHT,  HaloRegionTag::CENTER, HaloRegionTag::CENTER},  //  20 + 200 + 2000                          
                                                        {HaloRegionTag::LEFT,   HaloRegionTag::RIGHT,  HaloRegionTag::CENTER},  //  30 + 100 + 2000                          
                                                        {HaloRegionTag::CENTER, HaloRegionTag::RIGHT,  HaloRegionTag::CENTER},  //  30 + 200 + 2000                          
                                                        {HaloRegionTag::RIGHT,  HaloRegionTag::RIGHT,  HaloRegionTag::CENTER},  //  30 + 300 + 2000                          
                                                        {HaloRegionTag::LEFT,   HaloRegionTag::LEFT,   HaloRegionTag::RIGHT},   //  10 + 100 + 3000                           
                                                        {HaloRegionTag::CENTER, HaloRegionTag::LEFT,   HaloRegionTag::RIGHT},   //  10 + 200 + 3000                           
                                                        {HaloRegionTag::RIGHT,  HaloRegionTag::LEFT,   HaloRegionTag::RIGHT},   //  10 + 300 + 3000                           
                                                        {HaloRegionTag::LEFT,   HaloRegionTag::CENTER, HaloRegionTag::RIGHT},   //  20 + 100 + 3000                           
                                                        {HaloRegionTag::CENTER, HaloRegionTag::CENTER, HaloRegionTag::RIGHT},   //  20 + 200 + 3000                           
                                                        {HaloRegionTag::RIGHT,  HaloRegionTag::CENTER, HaloRegionTag::RIGHT},   //  20 + 300 + 3000                           
                                                        {HaloRegionTag::LEFT,   HaloRegionTag::RIGHT,  HaloRegionTag::RIGHT},   //  30 + 100 + 3000                           
                                                        {HaloRegionTag::CENTER, HaloRegionTag::RIGHT,  HaloRegionTag::RIGHT},   //  30 + 200 + 3000                           
                                                        {HaloRegionTag::RIGHT,  HaloRegionTag::RIGHT,  HaloRegionTag::RIGHT}}}; //  30 + 300 + 3000                           


// tuple to collect all bits together
static const auto _halo_regions = std::make_tuple(0, _bnds_1d, _bnds_2d, _bnds_3d);

// ===================================================================== //
// OBTAIN AN ARRAY OF HALO REGIONS TO ITERATE OVER
template <size_t NDIMS>
auto HaloRegions() -> decltype(std::get<NDIMS>(_halo_regions)) {
    return std::get<NDIMS>(_halo_regions);
}

}