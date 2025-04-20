---
title: Gaussian Elimination(A=LU)
date: 2025-03-29 14:30
categories: Numerical Linear Algebra
mathjex: true
---

# LU分解
``` python
def gaussian_elimination(A):
    if len(A) != len(A[0]):
        raise ValueError(f"A is not a square matrix!")
    for i in range(len(A)):
        for j in range(i+1, len(A)):
            A[j][i] /= A[i][i]
            for k in range(i+1, len(A)):
                A[j][k] -= A[i][k] * A[j][i]
    return A
```
运行这段代码会完成对**A=LU**的分解，下三角阵L（除主对角线外）的元素储存在A的对应位置，上三角阵U同理。

# 列主元Gauss
这是一段改进的Gauss分解，通过选取列主元来保证L中元素绝对值小于等于1：
``` python
def col_major_gaussian(A):
    index = [i for i in range(len(A))]
    if len(A) != len(A[0]):
        raise ValueError(f"A is not a square matrix!")
    for col_num in range(len(A)):
        major = abs(A[col_num][col_num])
        index[col_num] = col_num
        for row_num in range(col_num+1, len(A)):
            if abs(A[row_num][col_num]) > major:
                major = abs(A[row_num][col_num])
                index[col_num] = row_num
        if index[col_num] != col_num:
            B = [0 for _ in range(col_num, len(A))]
            for k in range(col_num, len(A)):
                B[k-col_num] = A[col_num][k]
                A[col_num][k] = A[index[col_num]][k]
                A[index[col_num]][k] = B[k-col_num]
        for row_num in range(col_num+1, len(A)):
            A[row_num][col_num] /= A[col_num][col_num]
            for k in range(col_num+1, len(A)):
                A[row_num][k] -= A[col_num][k] * A[row_num][col_num]
    for col_num in range(1, len(A)-1):
        if index[col_num] != col_num:
                    B = [0 for _ in range(col_num)]
                    for k in range(col_num):
                        B[k] = A[col_num][k]
                        A[col_num][k] = A[index[col_num]][k]
                        A[index[col_num]][k] = B[k]    
    return A, index       
```

