#ifndef CUBISM_H
#define CUBISM_H

#ifdef _USE_NUMA_
#warning  _USE_NUMA_ is deprecated, use CUBISM_USE_NUMA instead.
#define CUBISM_USE_NUMA
#endif

#ifdef _ON_FERMI_
#warning  _ON_FERMI_ is deprecated, use CUBISM_ON_FERMI instead.
#define CUBISM_ON_FERMI
#endif

#ifdef _ALIGNBYTES_
#warning _ALIGNBYTES_ is deprecated, use CUBISM_ALIGNMENT instead.
#define CUBISM_ALIGNMENT _ALIGNBYTES_
#elif !defined(CUBISM_ALIGNMENT)
#define CUBISM_ALIGNMENT  // If you get duplicate definition, put all Cubism
                          // includes after the main header include.
#endif

#ifndef CUBISM_NAMESPACE_BEGIN
#define CUBISM_NAMESPACE_BEGIN namespace cubism {
#endif

#ifndef CUBISM_NAMESPACE_END
#define CUBISM_NAMESPACE_END   }  // namespace cubism
#endif

#endif
