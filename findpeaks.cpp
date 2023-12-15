//
// Created by lidon on 2023/12/12.
//

#include "findpeaks.h"

/**
 *
 * @param y  输入为列向量
 * @param th  设置峰值阈值
 * @param peakInterval  设置峰峰间隔
 * @param Loc 返回值 返回找到的峰值位置
 * @return null
 * @author lidongfang
 * @update 2023/12/12
 */

void getAllPeaks(const MatrixXf &y,MatrixXi &iPk){
    int len=int (y.rows())+2;
    MatrixXf yTemp=MatrixXf ::Zero(len,1);
    typedef std::numeric_limits<float> Info;
    yTemp(0,0)=Info::quiet_NaN();
    yTemp(len-1,0)=Info::quiet_NaN();
    yTemp.middleRows(1,y.rows())=y;
    MatrixXi iTemp=MatrixXi ::Zero(len,1);
    for(int i=0;i<len;i++){
        iTemp(i,0)=i;
    }
    MatrixXi yFinite=(!yTemp.array().isNaN()).cast<int>();
//    cout<<"yFinite::"<<endl;
//    cout<<yFinite<<endl;
    MatrixXi uneq=(yTemp.topRows(len-1).array()!=yTemp.bottomRows(len-1).array()).cast<int>();
//    cout<<"uneq::"<<endl;
//    cout<<uneq<<endl;
    MatrixXi yF_tmp=yFinite.topRows(len-1)+yFinite.bottomRows(len-1);
    MatrixXi iNeq_tmp=find_loc((uneq.array()*yF_tmp.array()).cast<float>());//////////
    MatrixXi iNeq=MatrixXi ::Zero(iNeq_tmp.rows()+1,1);
    iNeq.bottomRows(iNeq_tmp.rows())=iNeq_tmp.array()+1;
//    cout<<"iNeq::"<<endl;
//    cout<<iNeq<<endl;
    iTemp=iTemp(iNeq.col(0),0);
//    cout<<"iTemp::"<<endl;
//    cout<<iTemp;
    MatrixXf ytt=yTemp(iTemp.col(0),0);
    MatrixXf yt=diff(ytt);
    MatrixXf s=yt.cwiseSign().cast<float>();
    MatrixXi iMax=1+find_loc(-1*diff(s)).array();
    cout<<"iMax  "<<iMax<<endl;
//    MatrixXi iAny=1+find_loc((s.topRows(s.rows()-1).array()!=s.bottomRows(s.rows()-1).array()).cast<int>()).array();
    iPk=iTemp(iMax.col(0),0).array()-1;
    cout<<"iPk::"<<endl;
    cout<<iPk<<endl;
}


void removePeaksBelowMinPeakHeight(const MatrixXf &y,MatrixXi &iPk,float th){
    MatrixXi tmp=iPk;
    if(tmp.size()!=0){
        MatrixXf yt=y(tmp.col(0),0);
        MatrixXf ytt=(yt.array()>th).cast<float>();
        iPk=tmp(find_loc(ytt).col(0),0);
    }
}

void removePeaksBelowThreshold(const MatrixXf &y,MatrixXi &iPk){
    cout<<"new iPk"<<endl;
    cout<<iPk<<endl;
    MatrixXf base=MatrixXf ::Zero(iPk.rows(),1);
    for(int i=0;i<iPk.rows();i++){
        if(y(iPk(i,0)-1,0)>y(iPk(i,0)+1,0))
            base(i,0)=y(iPk(i,0)-1,0);
        else
            base(i,0)=y(iPk(i,0)+1,0);
    }
    iPk=iPk(find_loc((y(iPk.col(0),0).array()>=base.array()).cast<float>()).col(0),0);
}

void orderPeaks(const MatrixXf &y,const MatrixXi &iPk,const MatrixXi &b,MatrixXi &idx){
    if(b.size()==0){
        return;
    }
    MatrixXi s=b;
    MatrixXf x=y(iPk.col(0),0);
    quicksort(x,s,0,int(iPk.rows())-1,1);
    idx=b(s.col(0),0);
}

void findPeaksSeparatedByMoreThanMinPeakDistance(const MatrixXf &y,const MatrixXi &x,const MatrixXi &iPk,int pd,MatrixXi &idx){
    int len=iPk.rows();
    MatrixXi b(len,1);
    for(int i=0;i<len;i++){
        b(i,0)=i;
    }
    if(iPk.size()==0||pd==0){
        orderPeaks(y,iPk,b,idx);
    }
    MatrixXf pks=y(iPk.col(0),0);
    MatrixXi locs_tmp=x(iPk.col(0),0);
    MatrixXi sortIdx=b;
    quicksort(pks,sortIdx,0,len-1,0);
    MatrixXi locs=locs_tmp(sortIdx.col(0),0);
    MatrixXi idelete=MatrixXi ::Zero(len,1);
    for(int i=0;i<len;i++){
        if(idelete(i,0)==0){
            idelete=(idelete.array()+((locs.array()>=locs(i,0)-pd)&&(locs.array()<=locs(i,0)+pd)).array().cast<int>()).cast<bool>().cast<int>();
            idelete(i,0)=0;
        }
    }
    idelete=(!idelete.array().cast<bool>()).cast<int>();
    MatrixXf idx_tmp=sortIdx(find_loc(idelete.cast<float>()).col(0),0).cast<float>();
    quicksort(idx_tmp,b.topRows(idx_tmp.rows()),0,int(idx_tmp.rows())-1,1);
    idx=idx_tmp.cast<int>();
}

void findPeaks(const MatrixXf &y,float th,int peakInterval,MatrixXi &Loc){
    int len=y.rows();
    int maxN=len;
    MatrixXi x=MatrixXi ::Zero(len,1);
    for(int i=0;i<len;i++){
        x(i,0)=i;
    }
    MatrixXi iPk(1,1);
    iPk<<1;
    getAllPeaks(y,iPk);
    cout<<"1"<<endl;
    removePeaksBelowMinPeakHeight(y,iPk,th);
    cout<<"removePeaksBelowMinPeakHeight:::"<<endl;
    cout<<iPk<<endl;
    removePeaksBelowThreshold(y,iPk);
    cout<<"3"<<endl;
    if(iPk.size()==0){
        Loc.resize(0,0);
    }else {
        MatrixXi idx(1, 1);
        idx << 1;
        findPeaksSeparatedByMoreThanMinPeakDistance(y, x, iPk, peakInterval, idx);
        cout << 4 << endl;
        MatrixXi ipk1 = iPk(idx.col(0), 0);
        Loc = x(ipk1.col(0), 0);
    }
}