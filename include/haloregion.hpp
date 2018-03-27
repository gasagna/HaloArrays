#pragma once
#include <array>
#include <numeric>
#include <initializer_list>

////////////////////////////////////////////////////////
//                     HALO REGION                    //
////////////////////////////////////////////////////////

namespace DArrays {

// tags for the the left, center and right boundaries
enum class HaloTag: int {LEFT = 1, CENTER = 2, RIGHT = 3};

// tags to specificy whether we are pointing at the region
// that is supposed to be sent or the region that is supposed
// to be received
enum class HaloIntent : int {SEND = 0, RECV = 1};

// specification of the boundary region is an NDIMS array of HaloTags
template<size_t NDIMS>
using HaloRegion = std::array<HaloTag, NDIMS>;

// THESE VECTORS SHOULD BE CONSTRUCTED PROGRAMMATICALLY WITH A CARTESIAN PRODUCT ITERATOR!

// 1D arrangement
// L X R                                                                 // message_tag
static const std::array<HaloRegion<1>, 2> _bnds_1d = {{HaloTag::LEFT},   // 10
                                                      {HaloTag::RIGHT}}; // 30

// 2D arrangement
// LL LC LR
// CL XX CR
// RL RC RR
static const std::array<HaloRegion<2>, 8> _bnds_2d = {{HaloTag::LEFT,   HaloTag::LEFT},   // 10 + 100
                                                      {HaloTag::LEFT,   HaloTag::CENTER}, // 10 + 200
                                                      {HaloTag::LEFT,   HaloTag::RIGHT},  // 10 + 300
                                                      {HaloTag::CENTER, HaloTag::LEFT},   // 20 + 100
                                                      {HaloTag::CENTER, HaloTag::RIGHT},  // 20 + 300
                                                      {HaloTag::RIGHT,  HaloTag::LEFT},   // 30 + 100
                                                      {HaloTag::RIGHT,  HaloTag::CENTER}, // 30 + 200
                                                      {HaloTag::RIGHT,  HaloTag::RIGHT}}; // 30 + 300 

// 3D arrangement
// LLL LCL LRL  // LLC LCC LRC  // LLR LCR LRR
// CLL CCL CRL  // CLC XXX CRC  // CLR CCR CRR
// RLL RCL RRL  // RLC RCC RRC  // RLR RCR RRR
static const std::array<HaloRegion<3>, 26> _bnds_3d = {{HaloTag::LEFT,   HaloTag::LEFT,   HaloTag::LEFT},   //  10 + 100 + 1000
                                                       {HaloTag::LEFT,   HaloTag::CENTER, HaloTag::LEFT},   //  10 + 200 + 1000                            
                                                       {HaloTag::LEFT,   HaloTag::RIGHT,  HaloTag::LEFT},   //  10 + 300 + 1000                            
                                                       {HaloTag::CENTER, HaloTag::LEFT,   HaloTag::LEFT},   //  20 + 100 + 1000                           
                                                       {HaloTag::CENTER, HaloTag::CENTER, HaloTag::LEFT},   //  20 + 200 + 1000                           
                                                       {HaloTag::CENTER, HaloTag::RIGHT,  HaloTag::LEFT},   //  20 + 300 + 1000                           
                                                       {HaloTag::RIGHT,  HaloTag::LEFT,   HaloTag::LEFT},   //  30 + 100 + 1000                            
                                                       {HaloTag::RIGHT,  HaloTag::CENTER, HaloTag::LEFT},   //  30 + 200 + 1000                            
                                                       {HaloTag::RIGHT,  HaloTag::RIGHT,  HaloTag::LEFT},   //  30 + 300 + 1000                            
                                                       {HaloTag::LEFT,   HaloTag::LEFT,   HaloTag::CENTER}, //  10 + 100 + 2000                          
                                                       {HaloTag::LEFT,   HaloTag::CENTER, HaloTag::CENTER}, //  10 + 200 + 2000                          
                                                       {HaloTag::LEFT,   HaloTag::RIGHT,  HaloTag::CENTER}, //  10 + 300 + 2000                          
                                                       {HaloTag::CENTER, HaloTag::LEFT,   HaloTag::CENTER}, //  20 + 100 + 2000                          
                                                       {HaloTag::CENTER, HaloTag::RIGHT,  HaloTag::CENTER}, //  20 + 200 + 2000                          
                                                       {HaloTag::RIGHT,  HaloTag::LEFT,   HaloTag::CENTER}, //  30 + 100 + 2000                          
                                                       {HaloTag::RIGHT,  HaloTag::CENTER, HaloTag::CENTER}, //  30 + 200 + 2000                          
                                                       {HaloTag::RIGHT,  HaloTag::RIGHT,  HaloTag::CENTER}, //  30 + 300 + 2000                          
                                                       {HaloTag::LEFT,   HaloTag::LEFT,   HaloTag::RIGHT},  //  10 + 100 + 3000                           
                                                       {HaloTag::LEFT,   HaloTag::CENTER, HaloTag::RIGHT},  //  10 + 200 + 3000                           
                                                       {HaloTag::LEFT,   HaloTag::RIGHT,  HaloTag::RIGHT},  //  10 + 300 + 3000                           
                                                       {HaloTag::CENTER, HaloTag::LEFT,   HaloTag::RIGHT},  //  20 + 100 + 3000                           
                                                       {HaloTag::CENTER, HaloTag::CENTER, HaloTag::RIGHT},  //  20 + 200 + 3000                           
                                                       {HaloTag::CENTER, HaloTag::RIGHT,  HaloTag::RIGHT},  //  20 + 300 + 3000                           
                                                       {HaloTag::RIGHT,  HaloTag::LEFT,   HaloTag::RIGHT},  //  30 + 100 + 3000                           
                                                       {HaloTag::RIGHT,  HaloTag::CENTER, HaloTag::RIGHT},  //  30 + 200 + 3000                           
                                                       {HaloTag::RIGHT,  HaloTag::RIGHT,  HaloTag::RIGHT}}; //  30 + 300 + 3000                           


// tuple to collect all bits together
static const auto _halo_regions = std::make_tuple(0, _bnds_1d, _bnds_2d, _bnds3);

// ===================================================================== //
// OBTAIN AN ARRAY OF HALO REGIONS TO ITERATE OVER
template <size_t NDIMS>
auto HaloRegions() -> decltype(std::get<NDIMS>(_halo_regions)) {
    return std::get<NDIMS>(_halo_regions);
}