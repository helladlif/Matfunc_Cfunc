
#include "interpft.h"
//引入第三方库fftw

//MatrixXf fftw2matrix1(fftw_complex *data,int len){
//    MatrixXf a(len,2);
//    for(int i=0;i<len;i++){
//        a(i,0)=data[i][0];
//        a(i,1)=data[i][1];
//    }
//    return a;
//}

void fft_shift(fftw_complex *data, int size,int nyquist,int m,fftw_complex *temp) {


    for(int i=0;i<size;i++){
        temp[i][0]=0.0;
        temp[i][1]=0.0;
    }
    for (int i = 0; i < nyquist; ++i) {
        temp[i][0] = data[i][0];
        temp[i][1] = data[i][1];
    }
    for (int i = size-m+nyquist; i < size; ++i) {
        temp[i][0] = data[nyquist+i-size+m-nyquist][0];
        temp[i][1] = data[nyquist+i-size+m-nyquist][1];
    }
}

MatrixXf interpft(fftw_complex *x, int n,int N) {
    // Compute the FFT of the original signal
    fftw_complex* y = fftw_alloc_complex(n);
    fftw_plan plan = fftw_plan_dft_1d(n, x, y, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    int nyqst=ceil((n+1)/2.0);

    // Insert zeros to increase the length of the FFT result
    fftw_complex* y_interp = fftw_alloc_complex(N);

    fftw_complex *y_ex = fftw_alloc_complex(N);
    fft_shift(y,N,nyqst,n,y_ex);
    if(n-int(floor(float (n)/float(2)))==0){
        y_ex[nyqst-1][0]=y_ex[nyqst-1][0]/2.0;
        y_ex[nyqst-1][1]=y_ex[nyqst-1][1]/2.0;
        y_ex[nyqst+N-n-1][0]=y_ex[nyqst-1][0];
        y_ex[nyqst+N-n-1][1]=y_ex[nyqst-1][1];
    }

    // Compute the inverse FFT to get the interpolated signal
    fftw_plan inv_plan = fftw_plan_dft_1d(N, y_ex, y_interp, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(inv_plan);
    MatrixXf y_out=MatrixXf ::Zero(N,1);
    for (int i = 0; i < N; ++i) {
        y_out(i,0) = y_interp[i][0]/float(N);
    }
    fftw_free(y_interp);
    y_out=y_out*float(N)/float(n);
    fftw_destroy_plan(plan);
    fftw_destroy_plan(inv_plan);
    fftw_free(y);
    fftw_free(y_ex);
    return y_out;
}

///输入为列向量
void interp(const MatrixXf &x,int r,MatrixXf &y){
    int len=x.rows();
    int rl=r*len;
    int n=4;
    int rn=r*n;
    MatrixXf b(65,1);
    b<<4.33317300404410e-12,
            -0.00243121300213272,-0.00455931939103976,-0.00608673948086279,-0.00677751383073431,-0.00648585232767677,
            -0.00517775986060116,-0.00294360862771454,-1.87526641740302e-11,0.0135716402609144,0.0257844009638513,
            0.0349066333082616,0.0394577742311000,0.0383799884227652,0.0311866090741486,0.0180758645307547,4.38222685193681e-11,
            -0.0449928600232999,-0.0877010984414075,-0.122184827314076,-0.142658093428576,-0.143967522815231,-0.122046528396825,
            -0.0743154081440588,-6.94981651553863e-11,0.136968069748474,0.291005785229400,0.452218214986458,0.609836360661724,
            0.752957912301643,0.871305419836420,0.955956266335061,1.00000000008053,0.955956266335061,0.871305419836420,
            0.752957912301643,0.609836360661724,0.452218214986458,0.291005785229400,0.136968069748474,-6.94981651553863e-11,
            -0.0743154081440588,-0.122046528396825,-0.143967522815231,-0.142658093428576,-0.122184827314076,-0.0877010984414075,
            -0.0449928600232999,4.38222685193681e-11,0.0180758645307547,0.0311866090741486,0.0383799884227652,0.0394577742311000,
            0.0349066333082616,0.0257844009638513,0.0135716402609144,-1.87526641740302e-11,-0.00294360862771454,-0.00517775986060116,
            -0.00648585232767677,-0.00677751383073431,-0.00608673948086279,-0.00455931939103976,-0.00243121300213272,4.33317300404410e-12;
    MatrixXf y_tmp=MatrixXf ::Zero(rl,1);
    y_tmp(seq(0,rl-1,r),0)=x;
    MatrixXf od=MatrixXf ::Zero(2*rn,1);
    od(seq(0,2*rn-1,r),0)=2*x(0,0)-x(seq(2*n,1,fix<-1>),0).array();
    MatrixXf a(1,1);
    a<<1;
    MatrixXf z(0,1);
    MatrixXf zi=MatrixXf ::Zero(b.rows()-1,1);
    filter(b,a,od,z,zi);
    MatrixXf zf=MatrixXf ::Zero(b.rows()-1,1);
    y_tmp=filter(b,a,y_tmp,zi,zf);
    y=y_tmp;
    y(seq(0,(len-n)*r-1,1),0)=y_tmp(seq(rn,rl-1,1),0);
    od=MatrixXf ::Zero(2*rn,1);
    od(seq(0,2*rn-1,r),0)=2*x(len-1,0)-x(seq(len-2,len-2*n-1,fix<-1>),0).array();
    od=filter1(b,a,od,zf);
    y(seq(rl-rn,rl-1,1),0)=od.topRows(rn);
}


