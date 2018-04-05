#pragma once
#include <array>

////////////////////////////////////////////////////////
//                 BOUNDARY REGION                    //
////////////////////////////////////////////////////////

namespace DArrays {

// tags for the the left, center and right boundaries
enum class BoundaryTag: int {LEFT = 1, CENTER = 2, RIGHT = 3};

// tags to specificy whether we are pointing at the region
// that is supposed to be sent or the region that is supposed
// to be received
enum class BoundaryIntent : int {SEND = 0, RECV = 1};

// specification of the boundary region is an NDIMS array of BoundaryTags
template<size_t NDIMS>
using Boundary = std::array<BoundaryTag, NDIMS>;

// THESE VECTORS SHOULD BE CONSTRUCTED PROGRAMMATICALLY WITH A CARTESIAN PRODUCT ITERATOR!

// 1D arrangement
// L X R                                                                     // message_tag
static const std::array<Boundary<1>, 2> _bnds_1d = {{{BoundaryTag::LEFT},    // 10
                                                     {BoundaryTag::RIGHT}}}; // 30

// 2D arrangement
// LL LC LR
// CL XX CR
// RL RC RR
static const std::array<Boundary<2>, 8> _bnds_2d = {{{BoundaryTag::LEFT,   BoundaryTag::LEFT},    // 10 + 100
                                                     {BoundaryTag::LEFT,   BoundaryTag::CENTER},  // 10 + 200
                                                     {BoundaryTag::LEFT,   BoundaryTag::RIGHT},   // 10 + 300
                                                     {BoundaryTag::CENTER, BoundaryTag::LEFT},    // 20 + 100
                                                     {BoundaryTag::CENTER, BoundaryTag::RIGHT},   // 20 + 300
                                                     {BoundaryTag::RIGHT,  BoundaryTag::LEFT},    // 30 + 100
                                                     {BoundaryTag::RIGHT,  BoundaryTag::CENTER},  // 30 + 200
                                                     {BoundaryTag::RIGHT,  BoundaryTag::RIGHT}}}; // 30 + 300 

// 3D arrangement
// LLL LCL LRL  // LLC LCC LRC  // LLR LCR LRR
// CLL CCL CRL  // CLC XXX CRC  // CLR CCR CRR
// RLL RCL RRL  // RLC RCC RRC  // RLR RCR RRR
static const std::array<Boundary<3>, 26> _bnds_3d = {{{BoundaryTag::LEFT,   BoundaryTag::LEFT,   BoundaryTag::LEFT},    //  10 + 100 + 1000
                                                      {BoundaryTag::CENTER, BoundaryTag::LEFT,   BoundaryTag::LEFT},    //  10 + 200 + 1000                            
                                                      {BoundaryTag::RIGHT,  BoundaryTag::LEFT,   BoundaryTag::LEFT},    //  10 + 300 + 1000                            
                                                      {BoundaryTag::LEFT,   BoundaryTag::CENTER, BoundaryTag::LEFT},    //  20 + 100 + 1000                           
                                                      {BoundaryTag::CENTER, BoundaryTag::CENTER, BoundaryTag::LEFT},    //  20 + 200 + 1000                           
                                                      {BoundaryTag::RIGHT,  BoundaryTag::CENTER, BoundaryTag::LEFT},    //  20 + 300 + 1000                           
                                                      {BoundaryTag::LEFT,   BoundaryTag::RIGHT,  BoundaryTag::LEFT},    //  30 + 100 + 1000                            
                                                      {BoundaryTag::CENTER, BoundaryTag::RIGHT,  BoundaryTag::LEFT},    //  30 + 200 + 1000                            
                                                      {BoundaryTag::RIGHT,  BoundaryTag::RIGHT,  BoundaryTag::LEFT},    //  30 + 300 + 1000                            
                                                      {BoundaryTag::LEFT,   BoundaryTag::LEFT,   BoundaryTag::CENTER},  //  10 + 100 + 2000                          
                                                      {BoundaryTag::CENTER, BoundaryTag::LEFT,   BoundaryTag::CENTER},  //  10 + 200 + 2000                          
                                                      {BoundaryTag::RIGHT,  BoundaryTag::LEFT,   BoundaryTag::CENTER},  //  10 + 300 + 2000                          
                                                      {BoundaryTag::LEFT,   BoundaryTag::CENTER, BoundaryTag::CENTER},  //  20 + 100 + 2000                          
                                                      {BoundaryTag::RIGHT,  BoundaryTag::CENTER, BoundaryTag::CENTER},  //  20 + 200 + 2000                          
                                                      {BoundaryTag::LEFT,   BoundaryTag::RIGHT,  BoundaryTag::CENTER},  //  30 + 100 + 2000                          
                                                      {BoundaryTag::CENTER, BoundaryTag::RIGHT,  BoundaryTag::CENTER},  //  30 + 200 + 2000                          
                                                      {BoundaryTag::RIGHT,  BoundaryTag::RIGHT,  BoundaryTag::CENTER},  //  30 + 300 + 2000                          
                                                      {BoundaryTag::LEFT,   BoundaryTag::LEFT,   BoundaryTag::RIGHT},   //  10 + 100 + 3000                           
                                                      {BoundaryTag::CENTER, BoundaryTag::LEFT,   BoundaryTag::RIGHT},   //  10 + 200 + 3000                           
                                                      {BoundaryTag::RIGHT,  BoundaryTag::LEFT,   BoundaryTag::RIGHT},   //  10 + 300 + 3000                           
                                                      {BoundaryTag::LEFT,   BoundaryTag::CENTER, BoundaryTag::RIGHT},   //  20 + 100 + 3000                           
                                                      {BoundaryTag::CENTER, BoundaryTag::CENTER, BoundaryTag::RIGHT},   //  20 + 200 + 3000                           
                                                      {BoundaryTag::RIGHT,  BoundaryTag::CENTER, BoundaryTag::RIGHT},   //  20 + 300 + 3000                           
                                                      {BoundaryTag::LEFT,   BoundaryTag::RIGHT,  BoundaryTag::RIGHT},   //  30 + 100 + 3000                           
                                                      {BoundaryTag::CENTER, BoundaryTag::RIGHT,  BoundaryTag::RIGHT},   //  30 + 200 + 3000                           
                                                      {BoundaryTag::RIGHT,  BoundaryTag::RIGHT,  BoundaryTag::RIGHT}}}; //  30 + 300 + 3000                           


// tuple to collect all bits together
static const auto _boundaries = std::make_tuple(0, _bnds_1d, _bnds_2d, _bnds_3d);

// ===================================================================== //
// OBTAIN AN ARRAY OF BOUNDARY REGIONS TO ITERATE OVER
template <size_t NDIMS>
auto AllBoundaries() -> decltype(std::get<NDIMS>(_boundaries)) {
    return std::get<NDIMS>(_boundaries);
}

}