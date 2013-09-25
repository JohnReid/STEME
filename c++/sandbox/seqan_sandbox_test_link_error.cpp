/*
Getting link error:

bin/gcc-4.6/release/sandbox/seqan_sandbox_test_link_error.o: In function `seqan::Pipe<seqan::Pipe<seqan::Pool<seqan::Pair<unsigned long, unsigned long, seqan::Tag<seqan::Compressed_> >, seqan::Mappe
    rSpec<seqan::MapperConfigSize<seqan::_skew7NMapSliced<seqan::Pair<unsigned long, unsigned long, seqan::Tag<seqan::Compressed_> >, unsigned long>, unsigned long, seqan::File<seqan::Async<void> >
    > > >, seqan::Filter<seqan::filterI2<seqan::Pair<unsigned long, unsigned long, seqan::Tag<seqan::Compressed_> >, unsigned long> > >, seqan::LarssonSadakane>::Pipe(seqan::Pipe<seqan::Pool<seqan::
    Pair<unsigned long, unsigned long, seqan::Tag<seqan::Compressed_> >, seqan::MapperSpec<seqan::MapperConfigSize<seqan::_skew7NMapSliced<seqan::Pair<unsigned long, unsigned long, seqan::Tag<seqan:
    :Compressed_> >, unsigned long>, unsigned long, seqan::File<seqan::Async<void> > > > >, seqan::Filter<seqan::filterI2<seqan::Pair<unsigned long, unsigned long, seqan::Tag<seqan::Compressed_> >,
    unsigned long> > >&)':
seqan_sandbox_test_link_error.cpp:(
        .text._ZN5seqan4PipeINS0_INS_4PoolINS_4PairImmNS_3TagINS_11Compressed_EEEEENS_10MapperSpecINS_16MapperConfigSizeINS_16_skew7NMapSlicedIS6_mEEmNS_4FileINS_5AsyncIvEEEEEEEEEENS_6FilterINS_8filt
        erI2IS6_mEEEEEENS_15LarssonSadakaneEEC2ERSM_[
        _ZN5seqan4PipeINS0_INS_4PoolINS_4PairImmNS_3TagINS_11Compressed_EEEEENS_10MapperSpecINS_16MapperConfigSizeINS_16_skew7NMapSlicedIS6_mEEmNS_4FileINS_5AsyncIvEEEEEEEEEENS_6FilterINS_8filterI2IS
        6_mEEEEEENS_15LarssonSadakaneEEC5ERSM_]+0x20a): undefined reference to
   `void seqan::createSuffixArray<seqan::String<unsigned long, seqan::Alloc<void> >, seqan::String<unsigned long, seqan::Alloc<void> > >(
        seqan::String<unsigned long, seqan::Alloc<void> > &, seqan::String<unsigned long, seqan::Alloc<void> > &, seqan::LarssonSadakane const &, unsigned int
    )'
collect2: ld returned 1 exit status
 */

#include <seqan/index.h>
#include <seqan/sequence.h>

using namespace seqan;

typedef String< Dna5 >                                      string_t;     /**< A string of Dna5. */
typedef seqan::StringSet< string_t >                        string_set_t; /**< StringSet type. */
typedef Index< string_set_t >                               index_t;      /**< The index over our strings. */
typedef seqan::Iterator< index_t, seqan::TopDown<> >::Type  top_down_it;  /**< A iterator over the index type. */

int
main( int argc, char * argv[] ) {
    string_set_t sequences;
    index_t index( sequences );
    top_down_it it( index );
    return 0;
}
