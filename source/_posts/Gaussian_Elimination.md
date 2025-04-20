---
title: Cholesky Factorization
date: 2025-03-29 14:50:02
categories: Numerical Linear Algebra
mathjex: true
---

# Cholesky分解
对于有着特殊性质的矩阵，在这里我们考虑**对称正定阵**，我们可以设计一些奇妙的算法来对他进行分解
```python
def cholesky_factorization(A):              # Cholesky分解，使用直接法减少计算量，得到的下三角矩阵L进行了本地化存储
    if len(A) != len(A[0]):
        raise ValueError(f"A is not a square matrix!")
    for col_num in range(len(A[0])):
        A[col_num][col_num] = math.sqrt(A[col_num][col_num])    # 对角线元素
        for row_num in range(col_num+1, len(A)):    
            A[row_num][col_num] /= A[col_num][col_num]          # 该列的下三角部分，计算后直接原地储存
        for col_num_2 in range(col_num+1, len(A[0])):           # 修改右下角剩余矩阵元素的值，提前减去求和中的一项，减少重复计算
            for row_num_2 in range(col_num, len(A)):
                A[row_num_2][col_num_2] -= A[row_num_2][col_num]*A[col_num_2][col_num] 
    return A
```

# LDLt分解
注意到上面的算法在每次迭代时需对主元进行开方操作，而在矩阵的条件数较大且矩阵本身规模较大时，可能会出现
由于计算过程中的舍入误差，导致中间出现了负数列主元，从而导致算法失败。

为此，我们可以对Cholesky分解进行改进：
```python
def LDLt(A):                # 改进后的Cholesky分解，L中除主对角线外元素存储在初始矩阵的对应位置，对角线位置用于存放D中的元素
    if len(A) != len(A[0]):
        raise ValueError(f"A is not a square matrix!")
    for col_num in range(len(A[0])):
        for col_num_2 in range(col_num):
            s = A[col_num_2][col_num_2]*A[col_num][col_num_2]
            A[col_num][col_num] -= A[col_num][col_num_2]*s
            for row_num in range(col_num+1, len(A)):
                A[row_num][col_num] -= A[row_num][col_num_2]*s
        for row_num in range(col_num+1, len(A)):
            A[row_num][col_num] /= A[col_num][col_num]
    return A
```

test post update