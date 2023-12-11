//
// Created by lidon on 2023/12/11.
//
#include "decimate.h"
////输入的是列向量
void decimate(const MatrixXf &x,int r,MatrixXf &y){
    int nd=x.rows();
    int nout=int(ceil(float(nd)/float(r)));
    MatrixXf b(31,1);
    b<<-0.000647960381752297,	-0.00144397905293930,	-0.00270225796524958,	-0.00444074474768837,	-0.00619148453717438,	-0.00695915901448257,
            -0.00537067608563326,	2.39068676949218e-18,	0.0102068181832813,	0.0255224469502992,	0.0451695881610100,	0.0672889149763633,	0.0891803923257266,
            0.107780638659668,	0.120271315654830,	0.124672293747482,	0.120271315654830,	0.107780638659668,	0.0891803923257266,	0.0672889149763633,
            0.0451695881610100,	0.0255224469502992,	0.0102068181832813,	2.39068676949218e-18,	-0.00537067608563326,	-0.00695915901448257,	-0.00619148453717438,
            -0.00444074474768837,	-0.00270225796524958,	-0.00144397905293930,	-0.000647960381752297;
    int nfilt=31;
    MatrixXf q=2*x(0,0)-x(seq(nfilt,1,fix<-1>),0).array();
    MatrixXf zi=MatrixXf::Zero(30,1);
    MatrixXf a(1,1);
    a<<1;
    MatrixXf z(0,1);
    filter(b,a,q,z,zi);
    MatrixXf zf=MatrixXf::Zero(30,1);
    MatrixXf c=filter(b,a,x,zi,zf);
    MatrixXf itemp1=2*x(nd-1,0)-x(seq(nd-2,nd-2*nfilt-1,fix<-1>),0).array();
    MatrixXf itemp=filter1(b,a,itemp1,zf);
    cout<<itemp<<endl;
    int len=int(ceil(float(nd-16)/float(r)));
    MatrixXi list(len,1);
    for(int i=0;i<len;i++){
        list(i,0)=16+i*r;
    }
    MatrixXf c1=c(list.col(0).array()-1,0);
    int nlen=nout-len;
    int nbeg=r-(nd-list(len-1,0));
    y=MatrixXf ::Zero(nout,1);
    y.topRows(len)=c1;
    y.bottomRows(nlen)=itemp(seq(nbeg-1,nbeg+nlen*r-2,r),0);
}



