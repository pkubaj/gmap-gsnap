/* $Id: f78e41fe0b4fb22c5e3c132969387f7eb0722c5e $ */
#ifndef PATH_LEARN_INCLUDED
#define PATH_LEARN_INCLUDED


#include "path.h"
#include "pathpair.h"
#include "univcoord.h"
#include "uintlist.h"


extern void
Path_learn_defect_rate (Path_T this, unsigned long long *total_mismatches,
			unsigned long long *total_querylength);

extern void
Path_learn_introns (Path_T this, Univcoordlist_T *donor_startpoints, Univcoordlist_T *donor_partners,
		    Univcoordlist_T *acceptor_startpoints, Univcoordlist_T *acceptor_partners,
		    Univcoordlist_T *antidonor_startpoints, Univcoordlist_T *antidonor_partners,
		    Univcoordlist_T *antiacceptor_startpoints, Univcoordlist_T *antiacceptor_partners);

extern void
Path_learn_indels (Path_T this, Univcoordtable_T indel_table);

extern void
Pathpair_learn_insertlengths (Pathpair_T pathpair, Uintlist_T *insertlengths);

extern void
Pathpair_analyze_insertlengths (int *expected_pairlength, int *pairlength_deviation,
				Uintlist_T insertlengths);

#endif


