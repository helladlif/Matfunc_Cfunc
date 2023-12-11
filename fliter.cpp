#include "fliter.h"


#define FILTER_LEN 10
///b,a,x,z均为列向量


MatrixXf filter(const MatrixXf &b,const MatrixXf &a,const MatrixXf &x,const MatrixXf &z ,MatrixXf &zf){
    int lb=b.rows();
    int la=a.rows();
    int len=int(max(b.rows(),a.rows()))-1;
    MatrixXf y=MatrixXf ::Zero(x.rows()+len,1);
    MatrixXf a1=a;
    MatrixXf b1=b;
    if(a(0,0)!=1){
        a1=a/a(0,0);
        b1=b/a(0,0);
    }
    for(int i=0;i<y.rows();i++){
        for(int k=0;k<lb;k++){
            if(i-k>=0&&i-k<x.rows()){
                y(i,0)+=b1(k,0)*x(i-k,0);}
        }
        for(int k=1;k<la;k++){
            if(i-k>=0&&i-k<x.rows())
                y(i,0)-=a1(k,0)*y(i-k,0);
        }
        if(z.size()!=0){
            if(i<z.rows())
                y(i,0)+=z(i,0);
        }
    }

    zf=y.bottomRows(len);
    y.conservativeResize(x.rows(),1);
    return y;
}



MatrixXf filter1(const MatrixXf &b,const MatrixXf &a,const MatrixXf &x,const MatrixXf &z){
    int lb=b.rows();
    int la=a.rows();
    MatrixXf y=x;


    y.setZero();
    MatrixXf a1=a;
    MatrixXf b1=b;
    if(a(0,0)!=1){
        a1=a/a(0,0);
        b1=b/a(0,0);
    }
    for(int i=0;i<x.rows();i++){
        for(int k=0;k<lb;k++){
            if(i-k>=0){
                y(i,0)+=b1(k,0)*x(i-k,0);}
        }
        for(int k=1;k<la;k++){
            if(i-k>=0)
                y(i,0)-=a1(k,0)*y(i-k,0);
        }
        if(z.size()!=0){
            if(i<z.rows())
                y(i,0)+=z(i,0);
        }

    }
    return y;
}

///b,a,x均为列向量
MatrixXf filtfilt(const MatrixXf &b,const MatrixXf &a,const MatrixXf &x){
    int npts=x.rows();
//    cout<<npts<<endl;
    int l=b.cols();
//    cout<<"l===="<<l<<endl;
    MatrixXf b1=b/a(0,0);
    MatrixXf a1=a/a(0,0);
    int nb=b.rows();
    int na=a.rows();
    int nfilt=max(nb,na);
    int nfact = max(1,3*(nfilt-1)); ///length of edge transients
    MatrixXf yout=x;
    if(nb<nfilt){
        b1.conservativeResize(nfilt,1);
        b1.bottomRows(nfilt-nb).setZero();
    } else if (na<nfilt){
        a1.conservativeResize(nfilt,1);
        a1.bottomRows(nfilt-na).setZero();
    }
    MatrixXf zi;
    if(nfilt>1){
//        MatrixXf z_t1=MatrixXf ::Identity(nfilt-1,nfilt-1);
///zi = ( eye(nfilt-1) - [-a1(2:nfilt,1), [eye(nfilt-2);zeros(1,nfilt-2)]] ) \
//     z_t3                           z_t2               z_t1
/// ( b1(2:nfilt,1) - b1(1,1)*a1(2:nfilt,1) );
///                   z_t4
        MatrixXf z_t1=MatrixXf::Zero(nfilt-1,nfilt-2);
        z_t1.topRows(nfilt-2)=MatrixXf ::Identity(nfilt-2,nfilt-2);
        MatrixXf z_t2(nfilt-1,nfilt-1);
        z_t2.leftCols(1)=a1(seq(1,nfilt-1,fix<1>),0)*(-1.0);
        z_t2.rightCols(nfilt-2)=z_t1;
        MatrixXf z_t3=MatrixXf ::Identity(nfilt-1,nfilt-1).array()-z_t2.array();

        MatrixXf z_t4(nfilt-1,1);
        z_t4=b1(seq(1,nfilt-1,fix<1>),0)-b1(0,0)*a1(seq(1,nfilt-1,fix<1>),0);
        zi=z_t3.inverse()*z_t4;
    }else{
        zi(0,1);////////?
    }

    for(int i=0;i<l;i++){
        MatrixXf ytemp(nfact*2+npts,1);
        ytemp.topRows(nfact)=2*yout(0,0)-yout(seq(nfact,1,fix<-1>),0).array();
        ytemp.bottomRows(nfact)=2*yout(npts-1,0)-yout(seq(npts-2,npts-nfact-1,fix<-1>),0).array();
        ytemp.middleRows(nfact,npts)=yout;
//        cout<<"ans1::"<<endl;
//        cout<<ytemp<<endl;

//        cout<<"ans2::"<<endl;
//        cout<<b1<<endl;
//        cout<<"ans3::"<<endl;
//        cout<<a1<<endl;
//        cout<<"ans4::"<<endl;
//        cout<<zi<<endl;

        MatrixXf ytemp1=filter1(b1.col(i),a1.col(i),ytemp.col(0),zi.col(i)*ytemp(0,0));
//        cout<<"ans5::"<<endl;
//        cout<<ytemp1<<endl;

        MatrixXf ytemp2=ytemp1(seq(nfact*2+npts-1,0,fix<-1>),0);
//        cout<<"ans6::"<<endl;
//        cout<<ytemp2<<endl;
        MatrixXf ytemp3=filter1(b1.col(i),a1.col(i),ytemp2.col(0),zi.col(i)*ytemp2(0,0));
//        cout<<"ans7::"<<endl;
//        cout<<ytemp3<<endl;

        yout=ytemp3(seq(nfact+npts-1,nfact,fix<-1>),0);
    }
    return yout;
}


//else{
//            for(int i=0;i<l;i++){
//                MatrixXf xt=-1*yout(seq(nfact,1,fix<-1>),0).array()+2*yout(0,0);
//                MatrixXf zo=filter(b1.col(i),a1.col(i),xt.col(0),zi.col(i)*xt(0,0));/////
//                MatrixXf yc2=filter(b1.col(i),a1.col(i),xt.col(0),zo);
//                xt=-1*yout(seq(last-1,last-nfact,fix<-1>),0).array()+2*yout(last,0);
//                MatrixXf yc3=filter(b1.col(i),a1.col(i),xt.col(0),zo);
//                zo=filter(b1.col(i),a1.col(i),yc3(seq(last,0,fix<-1>),0),zi.col(i)*yc3(last,0));
//                MatrixXf yc5=filter(b1.col(i),a1.col(i),yc2(seq(last,0,fix<-1>),0),zo);
//                yout=yc5(seq(last,0,fix<-1>),0);
//            }
//    }

