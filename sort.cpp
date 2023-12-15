#include "sort.h"
//int getpivot_num(int *arr1, int left, int right) {   //每个循环的条件都是left<right
//    double pivot1 = arr1[left];
//    //把最左侧的值赋给支点
//    while (left < right) {
//        while (arr1[right] >= pivot1 && left < right)   right--;
//        arr1[left] = arr1[right];
//        while (arr1[left] < pivot1&&left < right)   left++;
//        arr1[right] = arr1[left];
//    }
//    arr1[right] = pivot1;// 这里left==right 所以都可以 下同
//    return left;
//}
//
//void quicksort_num(int *arr1,int left, int right)
//{
//    if (left < right) {
//        int pivot = getpivot_num(arr1, left, right);//找中间位置左右分割
//        quicksort_num(arr1, left, pivot - 1);//处理左边，这里是一个递归的过程
//        quicksort_num(arr1,  pivot + 1, right);//处理右边 ，这里是一个递归的过程
//    }
//}
////arr1为列向量，left = 0 ，right = length-1
//////仅对列向量 type=1 升序，type=0 降序
int getpivot(Ref<MatrixXf> arr1,Ref<MatrixXi> arr2,int left, int right,int type) {//每个循环的条件都是left<right
    float pivot1 = arr1(left,0);
    int pivot2 = arr2(left,0);
    //把最左侧的值赋给支点
    if(type ==1) {
        while (left < right) {
            while (arr1(right, 0) >= pivot1 && left < right) right--;
            arr1(left, 0) = arr1(right, 0);
            arr2(left, 0) = arr2(right, 0);
            while (arr1(left, 0) < pivot1 && left < right) left++;
            arr1(right, 0) = arr1(left, 0);
            arr2(right, 0) = arr2(left, 0);
        }
        arr1(right, 0) = pivot1;// 这里left==right 所以都可以 下同
        arr2(right, 0) = pivot2;
        return left;
    }else if(type ==0){
        while (left < right) {
            while (arr1(right, 0) <= pivot1 && left < right) right--;
            arr1(left, 0) = arr1(right, 0);
            arr2(left, 0) = arr2(right, 0);
            while (arr1(left, 0) > pivot1 && left < right) left++;
            arr1(right, 0) = arr1(left, 0);
            arr2(right, 0) = arr2(left, 0);
        }
        arr1(right, 0) = pivot1;// 这里left==right 所以都可以 下同
        arr2(right, 0) = pivot2;
        return left;
    }else{
        exit(-1);
        printf("type must be 1 or 0 !!!\n");
    }

}

//arr1是重排后的数组，arr2是记录的位置
//int getpivot(double *arr1,int *arr2,int left, int right) {   //每个循环的条件都是left<right
//    double pivot1 = arr1[left];
//    int pivot2 = arr2[left];
//    //把最左侧的值赋给支点
//    while (left < right) {
//        while (arr1[right] >= pivot1 && left < right)   right--;
//        arr1[left] = arr1[right];
//        arr2[left] = arr2[right];
//        while (arr1[left] < pivot1&&left < right)   left++;
//        arr1[right] = arr1[left];
//        arr2[right] = arr2[left];
//    }
//    arr1[right] = pivot1;// 这里left==right 所以都可以 下同
//    arr2[right] = pivot2;
//    return left;
//}

void quicksort(Ref<MatrixXf> arr1,Ref<MatrixXi> arr2,int left, int right,int type)
{
//    double start = clock();
    if (left < right) {
        int pivot = getpivot(arr1,arr2,left, right,type);//找中间位置左右分割
        quicksort(arr1,arr2,left, pivot - 1,type);//处理左边，这里是一个递归的过程
        quicksort(arr1,arr2,pivot + 1, right,type);//处理右边 ，这里是一个递归的过程
    }
//    double finish = clock();
//    double duration = (double)(finish - start) / CLOCKS_PER_SEC;
//    printf("本次quick排序用了%f seconds\n", duration);
}

//void sort1(MatrixXf &arr,MatrixXi &index) {
//    int left = 0, right = arr.size()-1;
//    //double pivot = rms_all[left];
//    quicksort(arr,index,left, right);
//
//}//返回包含排列后得数据和数据对应位置的链表