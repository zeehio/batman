// written by Dr. Jie Hao, Dr William Astle
#include "dwt.h"
#include "myheader.h"

// convolve two vectors and decimate 
// could make more abstraction by templating element type
void conv_dec(const db_vit &v1start_vit, const db_vit &v1end_vit, int sign, const db_vit &v2start_vit, const db_vit &v2end_vit, int dec, db_vit &out_vit)
{  
	
    signed int it1, it2;
    signed int size1=v1end_vit-v1start_vit;
    signed int size2=v2end_vit-v2start_vit;
    signed int lb, ub;
    double buf;
    
    for(it1=0;it1<size1+size2-1;it1++)
    {
        if((it1%2)||(dec==0))
        {
            lb= (0<it1-size2+1) ? it1-size2+1 : 0; 
            ub= (size1-1<it1) ? size1-1 : it1;
            *out_vit=0; 
            for(it2=0;it2<ub-lb+1;it2++)
            {
                  //which direction
                if(sign%2==1)
                    buf=*(v1end_vit-1-(it2+lb))*(*(v2start_vit+it1-it2-lb));
                else
                    buf=*(v1start_vit+it2+lb)*(*(v2start_vit+it1-it2-lb));
                  
                if((sign<2)&&((it2+lb)%2==0))
                    *out_vit-=buf; 
                else
                    *out_vit+=buf; 
            }
            ++out_vit;
        }
    }
}


void my_for_fath(vector<double>& fath_vec, vector<double>& h, vector<double> &fath_dn_vec)
{
  db_vit fath_dn_vit=fath_dn_vec.begin();
  conv_dec(h.begin(), h.end(), 2, fath_vec.begin(), fath_vec.end(), 1,  fath_dn_vit);
}


void my_for_moth(vector<double>& fath_vec, vector<double>& h_vec, db_vit &moth_dn_vit)
{
  conv_dec(h_vec.begin(), h_vec.end(), 1, fath_vec.begin(), fath_vec.end(),1, moth_dn_vit);
}

void my_rev_fath(vector<double>& fath_vec, vector<double>& moth_vec, vector<double>& h_vec,  vector<double>& fath_up_vec)
{
  unsigned int len=fath_vec.size();  
  unsigned int lenup=fath_up_vec.size();  
  vector<double> dyfath_vec(2*len-1,0.0);
  vector<double> dymoth_vec(2*len-1,0.0);
  
  dydadup(fath_vec,dyfath_vec);
  dydadup(moth_vec,dymoth_vec);
  
  vector<double> fconv(h_vec.size()+2*len-2,0.0);
  vector<double> mconv(h_vec.size()+2*len-2,0.0);  
 
  db_vit fconv_vit=fconv.begin(); 
  db_vit mconv_vit=mconv.begin();
  
  vector<double> sconv(h_vec.size()+2*len-1,0.0);
  conv_dec(h_vec.begin(), h_vec.end(), 3, dyfath_vec.begin(), dyfath_vec.end(), 0, fconv_vit);
  conv_dec(h_vec.begin(), h_vec.end(), 0, dymoth_vec.begin(), dymoth_vec.end(), 0, mconv_vit);
  VVdif(&fconv,&mconv, &sconv);
  
  for(unsigned int it=0;it<lenup;it++)
    fath_up_vec[it]=*(sconv.begin()+h_vec.size()-2+it);  
}

void dwt(vector<double>& x, unsigned int lev, vector<double>& h_vec, vector<double>& wx)
{
  vector<double>* pfath, pnewfath;
  pfath= new vector<double>(x);
  
  size_t lfath=pfath->size();
  unsigned int offset=0;
  db_vit wx_it;
  unsigned it;
  //  *(levsize.end()-1)=x.size();
  for(it=0;it<lev;it++)
  {  
      vector<double>* pnewfath=new vector<double>((h_vec.size()+lfath-2+lfath%2)/2);
      lfath=pnewfath->size();
      offset+=lfath;
      //      *(levsize.end()-it-2)=(int)lfath;
      wx_it=wx.end()-offset;
      my_for_moth(*pfath, h_vec, wx_it);
      my_for_fath(*pfath, h_vec, *pnewfath);
      delete pfath; 
      pfath=new vector<double> (*pnewfath);
      delete pnewfath;
  }
  for(it=0;it<lfath;it++)
    *(wx.begin()+it)=(*pfath)[it]; 
  //  *(levsize.end()-lev-2)=(int)lfath;
  delete pfath;
}

void idwt(vector<double>& wx, vector<int>& levsize, vector<double>& h_vec, vector<double>& x)
{
 
  vector<double>* pfath_vec, *pnewfath_vec;
  pfath_vec=new vector<double>(wx.begin(),wx.begin()+levsize[0]);
  unsigned int it;
  unsigned int offset=levsize[0];
  unsigned int nlevs=levsize.size()-2;
 
  for(it=0;it<nlevs;it++)
  {   
      vector<double> moth_vec(wx.begin()+offset, wx.begin()+offset+levsize[1+it]);
      // for debug
      offset+=levsize[1+it];
      pnewfath_vec=new vector<double>(levsize[2+it]);//2*pfath_vec->size()-h_vec.size()+2);
      my_rev_fath(*pfath_vec, moth_vec, h_vec,  *pnewfath_vec);
      delete pfath_vec; 
      pfath_vec=new vector<double> (*pnewfath_vec);
      delete pnewfath_vec; 
  } 

  x=*pfath_vec;
  delete pfath_vec;
}

void dwt_levvec(int n, int lev, int hsize, vector<int>& levsize)
{
  *(levsize.end()-1)=n;
  int it;
  for(it=0;it<lev;it++)
    *(levsize.end()-it-2)=(int)(hsize+*(levsize.end()-it-1)-2+*(levsize.end()-it-1)%2)/2;
  
  *(levsize.end()-lev-2)=(int)*(levsize.end()-it-1);
}


size_t dwt_nlev(size_t xsize, size_t hsize)
{
  // only works for certain classes of wavelet
  double val=log(((double)xsize)/((double)(hsize-1)))/log(2.0);
  return val> 0.0 ? (size_t) floor(val) : 0 ;
}


size_t dwt_size(size_t xsize, unsigned int lev, size_t hsize)
{
  size_t lfath=xsize;
  size_t size=0;
  unsigned int buf;
  for(unsigned int it=0;it<lev;it++)
  {  
      buf=(hsize+lfath-2+lfath%2)/2;
      size+=buf;
      lfath=buf;
  }
  size+=lfath;
  return size;
}


void dydadup(vector<double>& vin, vector<double>& vout)
{
   for(size_t it=0;it<vout.size();it++)
   { 
   if(it%2)
    vout[it]=0.0;
   else	
	vout[it]=vin[it/2];
   }
}

