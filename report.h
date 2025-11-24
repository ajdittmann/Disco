#ifndef DISCO_REPORT_H
#define DISCO_REPORT_H

#include "disco.h"

//These are defined in each report setup file
void setReportParams(struct domain *theDomain);

int num_shared_reports();
int num_distributed_aux_reports();
int num_distributed_integral_reports();

void get_shared_reports(double *Q, struct domain *theDomain);
void get_distributed_aux_reports(double *Q, struct domain *theDomain);
void get_distributed_integral_reports(const double *x, const double *prim,
                                      double *Q, struct domain *theDomain);

//These are in report.c
void initializeReport(struct domain *theDomain);
void report(struct domain *theDomain);

#endif
