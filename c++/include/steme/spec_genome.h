/** Copyright John Reid 2012
 *
 * \file
 * \brief Specialises some traits for genome sized indexes.
 */


#ifndef STEME_JR_11JUL2011_GENOME_H_
#define STEME_JR_11JUL2011_GENOME_H_

#include <steme/defs.h>
#include <steme/seqan_types.h>
#include <steme/background_model_on_the_fly.h>

namespace steme {




/// STEME specification for genome sized indices.
struct genome_spec { };





/// Specialisation for genome indexes.
template<>
struct index_sizes_meta< genome_spec > {
    typedef unsigned char seq_index_t;   ///< The type used to index sequences. Our index will work on up to numeric_limits< seq_index_t >::max() sequences.
    typedef unsigned int seq_offset_t;  ///< The type used to index positions in sequences. Our index will work on sequences up to numeric_limits< seq_offset_t >::max() base pairs long.
    typedef unsigned int text_offset_t; ///< The type used to index positions in the whole text. Our index will work on sequences of total length up numeric_limits< text_offset_t >::max() base pairs long.
};



/// Specialisation for genome indexes.
template<>
struct background_model_meta< genome_spec > {
    typedef on_the_fly_bg_model< genome_spec > type; ///< Background model type.
};

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
    struct Fibre< steme::index_meta< steme::genome_spec >::type, FibreSA >
    {
        typedef String<
                        Pair<
                             steme::index_sizes_meta< steme::genome_spec >::seq_index_t,
                             steme::index_sizes_meta< steme::genome_spec >::seq_offset_t,
                             SEQAN_PACK
                            >,
                        DefaultIndexStringSpec< steme::index_meta< steme::genome_spec >::type >::Type
                       > Type;

        // Use a mmapped string
        // typedef String< Pair<unsigned char, unsigned int, SEQAN_PACK>, MMap<> >   Type;
    };


    /**
     * Specialise the LCP fibre.
     */
    template <>
    struct Fibre< steme::index_meta< steme::genome_spec >::type, FibreLcp >
    {
        typedef String<
                        steme::index_sizes_meta< steme::genome_spec >::text_offset_t,
                        DefaultIndexStringSpec< steme::index_meta< steme::genome_spec >::type>::Type
                       > Type;
    };


    /**
     * Specialise the Childtab fibre.
     */
    template <>
    struct Fibre< steme::index_meta< steme::genome_spec >::type, FibreChildtab >
    {
        typedef String<
                        steme::index_sizes_meta< steme::genome_spec >::text_offset_t,
                        DefaultIndexStringSpec< steme::index_meta< steme::genome_spec >::type>::Type
                       > Type;
    };
}



#endif //STEME_JR_11JUL2011_GENOME_H_
