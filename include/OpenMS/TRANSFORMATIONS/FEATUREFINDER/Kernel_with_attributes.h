#ifndef CGAL_KERNEL_WITH_ATTRIBUTES_H
#define CGAL_KERNEL_WITH_ATTRIBUTES_H

#include <CGAL/basic.h>
#include <CGAL/Kernel/Type_equality_wrapper.h>

#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_void.hpp>

CGAL_BEGIN_NAMESPACE

template< class K, class K_base, class Attributes >
struct Kernel_with_attributes_base
  : public K_base::template Base< K >::Type
{
    typedef typename K_base::template Base< K >::Type Kernel_base;
    template < typename K2 >
    struct Base {
        typedef Kernel_with_attributes_base<K2, K_base, Attributes> Type;
    };

#define CGAL_KERNEL_WITH_ATTRIBUTES_CONSTRUCTOR(Z,N,TEXT)            \
        template< BOOST_PP_ENUM_PARAMS(N, typename A) >              \
        TEXT ## _base( BOOST_PP_ENUM_BINARY_PARAMS(N, const A, &a) ) \
          : Kernel_base::TEXT( BOOST_PP_ENUM_PARAMS(N, a) )          \
        {}

#define CGAL_Kernel_obj(X)                                                           \
    template< bool DummyParam >                                                      \
    struct X ## _base : public Kernel_base::X {                                      \
        typename Attributes::X ## _attribute attribute;                              \
        BOOST_PP_REPEAT_FROM_TO( 1, 15, CGAL_KERNEL_WITH_ATTRIBUTES_CONSTRUCTOR, X)  \
        X ## _base() {}                                                              \
    };                                                                               \
                                                                                     \
    typedef                                                                          \
        typename ::boost::mpl::if_c<                                                 \
            ::boost::is_void< typename Attributes::X ## _attribute >::value,         \
            typename Kernel_base::X,                                                 \
            X ## _base<true>                                                         \
        >::type X;

#include "my_kernel_objects.h"

#undef CGAL_KERNEL_WITH_ATTRIBUTES_CONSTRUCTOR
};


template< class K, class Attributes >
struct Kernel_with_attributes
  : public Type_equality_wrapper<
                Kernel_with_attributes_base< Kernel_with_attributes<K,
                                                                    Attributes>,
                                            K,
                                            Attributes >,
                Kernel_with_attributes<K,
                                        Attributes>
            >
{};

struct Kernel_with_attributes_default_attributes {
#define CGAL_Kernel_obj(X) typedef void X ## _attribute;
#include "my_kernel_objects.h"
};

template< typename Attribute >
struct Kernel_with_attributes_uniform_attributes {
#define CGAL_Kernel_obj(X) typedef Attribute X ## _attribute;
#include "my_kernel_objects.h"  
};

CGAL_END_NAMESPACE

#endif // CGAL_KERNEL_WITH_ATTRIBUTES_H

