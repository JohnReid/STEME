#include <stdio.h>


/**
 * Macro from http://dbp-consulting.com/tutorials/SuppressingGCCWarnings.html to
 * control GCC warnings
 */
#if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 402
# define GCC_DIAG_STR(s) #s
# define GCC_DIAG_JOINSTR(x,y) GCC_DIAG_STR(x ## y)
# define GCC_DIAG_DO_PRAGMA(x) _Pragma (#x)
# define GCC_DIAG_PRAGMA(x) GCC_DIAG_DO_PRAGMA(GCC diagnostic x)
# if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406
#  define GCC_DIAG_OFF(x) GCC_DIAG_PRAGMA(push) \
    GCC_DIAG_PRAGMA(ignored GCC_DIAG_JOINSTR(-W,x))
#  define GCC_DIAG_ON(x) GCC_DIAG_PRAGMA(pop)
# else
#  define GCC_DIAG_OFF(x) GCC_DIAG_PRAGMA(ignored GCC_DIAG_JOINSTR(-W,x))
#  define GCC_DIAG_ON(x)  GCC_DIAG_PRAGMA(warning GCC_DIAG_JOINSTR(-W,x))
# endif
#else
# define GCC_DIAG_OFF(x)
# define GCC_DIAG_ON(x)
#endif

/**
 * It seems that gcc 4.6 and 4.7 name this warning differently
*/
#if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 407
# define GCC_UNINIT_WARNING uninitialized
#else
# define GCC_UNINIT_WARNING uninitialized
#endif


GCC_DIAG_OFF(GCC_UNINIT_WARNING)
int main()
{
    int foo;
    printf("I am a number: %d \n", foo);
    return 0;
}
GCC_DIAG_ON(GCC_UNINIT_WARNING)
