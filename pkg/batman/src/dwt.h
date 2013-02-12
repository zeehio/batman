// written by Dr William Astle
#include <vector>

using namespace std;

typedef vector<double>::iterator db_vit;

void conv_dec(const db_vit &v1start_vit, const db_vit &v1end_vit, int sign, const db_vit &v2start_vit, const db_vit &v2end_vit, int dec, db_vit &out_vit);

void dwt(vector<double> &x, unsigned int lev, vector<double> &h, vector<double> &wx);
void idwt(vector<double>& wx, vector<int>& levsize, vector<double>& h, vector<double>& x);

void my_for_fath(vector<double>& fath_vec, vector<double>& h,  vector<double> &fath_dn_vec);
void my_for_moth(vector<double>& fath_vec, vector<double>& h,  db_vit &moth_dn_vit);
size_t dwt_size(size_t xsize, unsigned int lev, size_t hsize);
size_t dwt_nlev(size_t xsize, size_t hsize);
void dydadup(vector<double>& vin, vector<double>& vout);
void my_rev_fath(vector<double>& fath_vec, vector<double>& moth_vec, vector<double>& h_vec, const db_vit &fath_up_vit);
void dwt_levvec(int n, int lev, int hsize, vector<int>& levsize);

#define SYM6 0.015404109327027373,0.0034907120842174702,-0.11799011114819057,-0.048311742585632998,0.49105594192674662,0.787641141030194,0.3379294217276218,-0.072637522786462516,-0.021060292512300564,0.044724901770665779,0.0017677118642428036,-0.007800708325034148
