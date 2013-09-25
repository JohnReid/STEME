/** Copyright John Reid 2012
 *
 * \file
 * \brief Default traits for indexes.
 */


#ifndef STEME_JR_13JUL2011_SPEC_DEFAULT_H_
#define STEME_JR_13JUL2011_SPEC_DEFAULT_H_

#include <steme/defs.h>
#include <steme/seqan_types.h>
#include <steme/background_model_on_the_fly.h>

#include <seqan/basic.h>

namespace steme {




/// STEME specification for genome sized indices.
struct default_spec { };







} // namespace steme




//
// Specialise various types to reduce memory requirements for index. If you don't use steme::default_spec
// you will need to do repeat this for your spec (see genome_spec for an example).
//
namespace seqan
{

    /**
     * Specialise the SA fibre.
     */
    template<>
    struct Fibre< steme::index_meta< steme::default_spec >::type, FibreSA >
    {
        typedef String<
                        Pair<
                             steme::index_sizes_meta< steme::default_spec >::seq_index_t,
                             steme::index_sizes_meta< steme::default_spec >::seq_offset_t,
                             SEQAN_PACK
                            >,
                        DefaultIndexStringSpec< steme::index_meta< steme::default_spec >::type >::Type
                       > Type;

        // Use a mmapped string
        // typedef String< Pair<unsigned char, unsigned int, SEQAN_PACK>, MMap<> >   Type;
    };


    /**
     * Specialise the LCP fibre.
     */
    template <>
    struct Fibre< steme::index_meta< steme::default_spec >::type, FibreLcp >
    {
        typedef String<
                        steme::index_sizes_meta< steme::default_spec >::text_offset_t,
                        DefaultIndexStringSpec< steme::index_meta< steme::default_spec >::type>::Type
                       > Type;
    };


    /**
     * Specialise the Childtab fibre.
     */
    template <>
    struct Fibre< steme::index_meta< steme::default_spec >::type, FibreChildtab >
    {
        typedef String<
                        steme::index_sizes_meta< steme::default_spec >::text_offset_t,
                        DefaultIndexStringSpec< steme::index_meta< steme::default_spec >::type>::Type
                       > Type;
    };
}



#endif //STEME_JR_13JUL2011_SPEC_DEFAULT_H_
