// from
// http://askubuntu.com/questions/36211/undefined-reference-error-dl-stack-flags-with-gcc-and-pthreads

#include <pthread.h>

static void *
thread_start(void *arg)
{
    return 0;
}

int
main(int argc, char **argv)
{
    pthread_t thread_id = 0;
    void *result = NULL;

    pthread_create(&thread_id, NULL, &thread_start, NULL);
    pthread_join(thread_id, &result);

    return 0;
}
