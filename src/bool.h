#ifndef BOOL_INCLUDED
#define BOOL_INCLUDED

#if __STDC_VERSION__ < 202311L
/* typedef enum{false,true} bool; */
typedef unsigned char bool;
#endif

#define false 0
#define true 1



#endif
