
#include "pwelch.h"

fftw_complex *matrix2fftw(MatrixXf a){
    fftw_complex *temp = (fftw_complex *)malloc(a.size() * sizeof(fftw_complex));
    for(int i=0;i<a.rows();i++){
        temp[i][0]=a(i,0);
        temp[i][1]=0.0;
    }
    return temp;
}

void pwelch(MatrixXf x,int fs ,Ref<MatrixXf> Px,Ref<MatrixXf> Fv){
    int Nx=x.rows();
    int a=floor(Nx/4.5);
    int b=a/2;
    int d=ceil(log2(a));
    int e=pow(2,d);
    int c=std::max(256,e);
    int n=floor(c/2)+1;
    for(int i=0;i<n;i++){
        Fv(i,0)=float(i)*float(fs)/float(c);
    }
    int k=floor(float(Nx-b)/(a-b));
    int Lmin=a-b;
//    cout<<a<<" "<<b<<" "<<c<<" "<<n<<" "<<k<<" "<<Lmin<<" "<<endl;
    MatrixXi xstart(k,1);
    for(int i=0;i<k;i++){
        xstart(i,0)=1+i*Lmin;
    }
    MatrixXi xend=xstart.array()+a-1;
//    cout<<"start::"<<endl;
//    cout<<xstart<<endl;
//    cout<<"xend::"<<endl;
//    cout<<xend<<endl;
    MatrixXf in(a,1);
    for(int i=0;i<a;i++){
        in(i,0)=2*Pi*i/float(a-1);
    }
    MatrixXf win=0.54-0.46*(in.array().cos());
//    cout<<"win::"<<endl;
//    cout<<win<<endl;
    float u=win.array().square().mean()*a;
    MatrixXf X=MatrixXf ::Zero(c,1);
    fftw_complex* xin = fftw_alloc_complex(c);
    fftw_complex* yout = fftw_alloc_complex(c);
    for(int i=1;i<=k;i++){

        MatrixXf xt=MatrixXf ::Zero(c,1);
        xt.topRows(a)=x.middleRows(xstart(i-1,0)-1,xend(i-1,0)-xstart(i-1,0)+1).array()*win.array();
//        cout<<"xin::"<<endl;
//        cout<<xt.rows()<<xt.cols()<<endl;
        xin=matrix2fftw(xt);
//        for(int i=0;i<a;i++){
//            cout<<"xin::"<<endl;
//            cout<<xin[i][0]<<" "<<xin[i][1]<<endl;
//        }
        fftw_plan plan = fftw_plan_dft_1d(c, xin, yout, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(plan);
//        cout<<"i::"<<i<<endl;
//        for(int i=c-10;i<c;i++){
//            cout<<"yout::"<<endl;
//            cout<<yout[i][0]<<" "<<yout[i][1]<<endl;
//        }


//        MatrixXcf yo=fftw2matrix(yout,c);
        MatrixXf Pxx(c,1);
        for(int i=0;i<c;i++){
            Pxx(i,0)=(yout[i][0]*yout[i][0]+yout[i][1]*yout[i][1])/(float(k)*float(u));
        }

//        Pxx=Pxx.array()/k;
        X=X.array()+Pxx.array();
//        cout<<"X::"<<endl;
//        cout<<X<<endl;
        fftw_destroy_plan(plan);


    }
    fftw_free(yout);
//    fftw_free(xin);
    MatrixXi sel(c/2+1,1);
    for(int i=0;i<sel.size();i++){
        sel(i,0)=i;
    }
    MatrixXf Px_unscaled=X(sel.col(0),0);
    Px=2*Px_unscaled;
    Px.row(0)=Px_unscaled.row(0);
    Px.row(Px.rows()-1)=Px_unscaled.row(Px.rows()-1);
    Px=Px/fs;
}