void
setgval (SNP ** xsnps, int nrows, Indiv ** indivmarkers, int numindivs,
         int *xindex, int *xtypes, int ncols);
void
unsetgval ();
int
getgval (int row, int col, double *val);
int
getggval (int indindx, int col, double *val);

void
set_ind_mask ();

size_t
get_nrows ();
size_t
get_ncols ();

void
kjg_geno_get_normalized_row (const size_t snp_index, double* y);
size_t
kjg_geno_get_normalized_rows (const size_t i, const size_t r, double* Y);
