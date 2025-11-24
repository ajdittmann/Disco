#ifndef DISCO_PROFILER_H
#define DISCO_PROFILER_H

#include "disco.h"

void start_clock( struct domain * theDomain );
int count_cells( struct domain * theDomain );
void generate_log( struct domain * theDomain );

void prof_init(struct profiler *prof);
void prof_tick(struct profiler *prof, int label);
void prof_tock(struct profiler *prof, int label);

#endif
