/*!
\file field_ops.hh
\brief Operations on Fields types
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_FIELD_OPS_HH
#define SPECTRUM_FIELD_OPS_HH

#include "../compiletime_lists.hh"
#include <utility>


/*
 * These operations are intended to work on the Fields type
 * with the Field type as member types in the typelist.
 * The implementation is based on guidelines in:
 * Vandevoorde, Josuttis, Gregor. C++ Templates: The Complete Guide. 2nd edition. Addison-Wesley, 2018.
 *
 */


namespace Spectrum {

template <typename FConfig, typename... Ts>
class Fields;

}


namespace Spectrum::FieldOps {



namespace Impl {

// This establishes the depth-1 signature of the impl type
template<typename Fields>
struct Front_impl;

// This provides the template system (TS) with a rule for syntax pattern-matching (depth 2)
template<typename FConfig, typename T, typename... Ts>
struct Front_impl<Fields<FConfig, T, Ts...>> {
   using Front = T;
};

}

// Now the TS can parse the impl type as an expression, and we can extract the 'front' of the type
template<typename Fields>
using Front = Impl::Front_impl<Fields>::Front;


/*
 * We can repeat this {depth-1, depth-2, definition} pattern as many times as we like.
 *
 */


namespace Impl {

template<typename Fields>
class Pop_impl;

template<typename FConfig, typename T, typename... Ts>
struct Pop_impl<Fields<FConfig, T, Ts...>> {
   using Pop = Fields<FConfig, Ts...>;
};

}

template<typename Fields>
using Pop = typename Impl::Pop_impl<Fields>::Pop;






namespace Impl {

template<typename Fields, typename U>
class Push_impl;

template<typename U, typename FConfig, typename T, typename... Ts>
struct Push_impl<Fields<FConfig, T, Ts...>, U> {
   using Push = Fields<FConfig, U, T, Ts...>;
};

}

template <typename Fields, typename U>
using Push = typename Impl::Push_impl<Fields, U>::Push;









namespace Impl {

template<typename Fields, int N>
struct Nth_impl : public Nth_impl<Pop<Fields>, N - 1> {
};

template<typename Fields>
struct Nth_impl<Fields, 0> : public Front_impl<Fields> {
};

}

template <typename Fields, int N>
using Nth = typename Impl::Nth_impl<Fields, N>::Front;







namespace Impl {

template<typename Fields>
struct Size_impl;

template<typename FConfig, typename... Ts>
struct Size_impl<Fields<FConfig, Ts...>> {
   static constexpr int Size = sizeof...(Ts);
};

}

template <typename Fields>
constexpr int Size = Impl::Size_impl<Fields>::Size;









namespace Impl {

template<typename Fields, typename T>
struct Found_impl {
private:
   static constexpr bool FrontFound = std::same_as<Front<Fields>, T>;
   static constexpr bool RemainderFound = Found_impl<Pop<Fields>, T>::Found;
public:
   static constexpr bool Found = FrontFound || RemainderFound;
};

template<typename FConfig, typename T>
struct Found_impl<Fields<FConfig>, T> {
   static constexpr bool Found = false;
};

}

template <typename Fields, typename T>
constexpr bool Found = Impl::Found_impl<Fields, T>::Found;









namespace Impl {

template<typename Fields1, typename Fields2>
struct Concat_impl;

template<typename FConfig1, typename FConfig2, typename... Ts, typename... Us>
struct Concat_impl<Fields<FConfig1, Ts...>, Fields<FConfig2, Us...>> {
   using Concat = Fields<FConfig1, Ts..., Us...>;
};

}

/*
 * Caveat: Information will be lost, namely the information in the second FConfig.
 *
 */
template <typename Fields1, typename Fields2>
using Concat = typename Impl::Concat_impl<Fields1, Fields2>::Concat;














namespace Impl {

template<typename Fields, typename U>
class Push_Back_impl;

template<typename FConfig, typename... Ts, typename U>
struct Push_Back_impl<Fields<FConfig, Ts...>, U> {
   using Push_Back = Concat<Fields<FConfig, Ts...>, Fields<FConfig, U>>;
};

}

template <typename Fields, typename U>
using Push_Back = typename Impl::Push_Back_impl<Fields, U>::Push_Back;










namespace Impl {

template <typename Fields>
struct Reverse_impl {
private:
   using Front = Front<Fields>;
   using Remainder = Pop<Fields>;
   using RemainderReverse = typename Reverse_impl<Remainder>::Reverse;
public:
   using Reverse = Push_Back<RemainderReverse, Front>;
};

template <typename FConfig>
struct Reverse_impl<Fields<FConfig>> {
   using Reverse = Fields<FConfig>;
};

}

template <typename Fields>
using Reverse = typename Impl::Reverse_impl<Fields>::Reverse;










//--------------------------------------------------------------------
// Set-yielding operations


/*
 * All operations after this are guaranteed to produce 'Sets'
 * as defined by the following operation.
 * Caveat: this introduces inefficiency if operations are deeply nested.
 *
 */


namespace Impl {

template<typename Fields, int N>
struct Set_impl {
private:
   using Front = Front<Fields>;
   using Remainder = Pop<Fields>;
   using RemainderSet = typename Set_impl<Remainder, N - 1>::Set;
public:
   using Set = Cond<Found<Remainder, Front>, RemainderSet, Push_Back<RemainderSet, Front>>;
};

template<typename Fields>
struct Set_impl<Fields, 0> {
   using Set = Fields;
};

template<typename Fields>
using SetReverse = typename Set_impl<Fields, Size<Fields>>::Set;

}

template <typename Fields>
using Set = Reverse<Impl::SetReverse<Fields>>;








template <typename Fields1, typename Fields2>
using Union = Set<Concat<Fields1, Fields2>>;














namespace Impl {

template<typename Fields1, typename Fields2, int N>
struct Difference_impl {
private:
   using Front = Front<Fields1>;
   using Remainder = Pop<Fields1>;
   using RemainderDifference = typename Difference_impl<Remainder, Fields2, N-1>::Difference;
public:
   using Difference = Cond<Found<Remainder, Front> || Found<Fields2, Front>, RemainderDifference, Push_Back<RemainderDifference, Front>>;
};

template<typename Fields1, typename Fields2>
struct Difference_impl<Fields1, Fields2, 0> {
   using Difference = Fields1;
};


}

template <typename Fields1, typename Fields2>
using Difference = Reverse<typename Impl::Difference_impl<Fields1, Fields2, Size<Fields1>>::Difference>;












namespace Impl {

template<typename Fields1, typename Fields2, int N>
struct Intersection_impl {
private:
   using Front = Front<Fields1>;
   using Remainder = Pop<Fields1>;
   using RemainderIntersection = typename Intersection_impl<Remainder, Fields2, N-1>::Intersection;
public:
   using Intersection = Cond<Found<Remainder, Front> || !Found<Fields2, Front>, RemainderIntersection, Push_Back<RemainderIntersection, Front>>;
};

template<typename Fields1, typename Fields2>
struct Intersection_impl<Fields1, Fields2, 0> {
   using Intersection = Fields1;
};


}

template <typename Fields1, typename Fields2>
using Intersection = Reverse<typename Impl::Intersection_impl<Fields1, Fields2, Size<Fields1>>::Intersection>;




}


#endif
