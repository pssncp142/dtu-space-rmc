double vec_norm(double*);

int vec_unit(double*);
int vec_unit_wc(double*,double*);

double vec_dotp(double*,double*);

int vec_cross(double*,double*);
int vec_cross_wc(double*,double*,double*);

int vec_scap(double,double*);
int vec_scap_wc(double*,double,double*);

int vec_proj(double*,double*);
int vec_proj_wc(double*,double*,double*);

int vec_ort(double*,double*);
int vec_ort_wc(double*,double*,double*);

int vec_csum(double*,double*,double,double);
int vec_csum_wc(double*,double*,double*,double,double);

int vec_add(double*,double*);
int vec_add_wc(double*,double*,double*);

int vec_subt(double*,double*);
int vec_subt_wc(double*,double*,double*);

double vec_angle(double*,double*,double*);

int vec_rotate(double*,double,double*);
int vec_rotate_wc(double*,double*,double,double*);
