# NumericalCalculation — 数值计算方法 C++ 实现

## 目录

- [1. 项目概述](#1-项目概述)
- [2. 架构设计](#2-架构设计)
- [3. 模块详解](#3-模块详解)
  - [3.1 公共工具模块 (common_fd)](#31-公共工具模块-common_fd)
  - [3.2 文件日志模块 (filelog)](#32-文件日志模块-filelog)
  - [3.3 表达式解析与科学计算器 (formula)](#33-表达式解析与科学计算器-formula)
  - [3.4 线性方程组基类 (linear_equations)](#34-线性方程组基类-linear_equations)
  - [3.5 直接法求解线性方程组 (directmethod)](#35-直接法求解线性方程组-directmethod)
  - [3.6 迭代法求解线性方程组 (iterationmethod)](#36-迭代法求解线性方程组-iterationmethod)
  - [3.7 插值法 (interpolation)](#37-插值法-interpolation)
  - [3.8 函数最优逼近 (optimal_approx)](#38-函数最优逼近-optimal_approx)
  - [3.9 数值积分与数值微分 (int_diff)](#39-数值积分与数值微分-int_diff)
  - [3.10 非线性方程/方程组 (nonlinear_equations)](#310-非线性方程方程组-nonlinear_equations)
  - [3.11 矩阵特征值与特征向量 (eigen_val_vec)](#311-矩阵特征值与特征向量-eigen_val_vec)
- [4. 编译与运行](#4-编译与运行)
- [5. 使用说明](#5-使用说明)
- [6. 类继承关系图](#6-类继承关系图)

---

## 1. 项目概述

本项目是一套**数值计算方法**的 C++ 控制台程序，涵盖了大学《计算方法》课程的核心内容。程序提供交互式菜单，支持中英文双语界面，允许用户输入数据并选择不同的数值算法进行求解，最终将结果输出到屏幕和控制台，并可输出 MATLAB `.m` 脚本或 Plotly.js `.html` 网页用于可视化。

### 功能总览

| 编号 | 功能 | 包含算法数 |
|------|------|-----------|
| 0 | 科学计算器 | — |
| 1 | 直接法求解线性方程组 | 8 |
| 2 | 迭代法求解线性方程组 | 6 |
| 3 | 插值法 | 3 |
| 4 | 函数最优逼近 | 2 |
| 5 | 数值积分与数值微分 | 2 |
| 6 | 非线性方程（组） | 9 |
| 7 | 矩阵特征值与特征向量 | 2 |

---

## 2. 架构设计

### 文件组织

```
NumericalCalculation/
├── main.cpp                  # 主入口，菜单调度
├── common_fd.h / .cpp        # 公共工具函数（数学、字符串、矩阵/向量范数）
├── filelog.h / .cpp          # 文件日志输出（支持 .txt / .m / .html 格式）
├── formula.h / .cpp          # 表达式解析器 & 科学计算器
├── linear_equations.h / .cpp # 线性方程组基类 Ax_b
├── directmethod.h / .cpp     # 直接法（8种）
├── iterationmethod.h / .cpp  # 迭代法（6种）
├── interpolation.h / .cpp    # 插值法（3种）
├── optimal_approx.h / .cpp   # 最优逼近（2种）
├── int_diff.h / .cpp         # 数值积分与微分（2种）
├── nonlinear_equations.h/.cpp# 非线性方程（5种单方程 + 4种方程组）
├── eigen_val_vec.h / .cpp    # 特征值与特征向量（2种）
├── testdata.txt              # 测试数据样例
├── Makefile.win              # Windows Makefile (Dev-C++)
└── README.md                 # 简要说明
```

### 设计模式

- **策略模式**：每个算法模块使用基类定义 `calc()` / `out_result()` 虚函数接口，具体算法子类实现不同策略
- **工厂模式**：`main.cpp` 中通过 `switch` 语句根据用户选择实例化对应的算法子类
- **运算符重载**：`filelog` 类重载 `<<` 实现链式输出；`formula` 类重载 `*` 实现函数乘法

### 语言切换

通过 `common_fd.h` 中的 `#define CHINESE_VERSION` 宏控制：
- **定义时**：界面显示中文
- **注释掉时**：界面显示英文

---

## 3. 模块详解

### 3.1 公共工具模块 (common_fd)

**文件**: `common_fd.h`, `common_fd.cpp`

提供项目各处使用的公共函数和模板：

| 函数/模板 | 说明 |
|-----------|------|
| `sign(double x)` | 符号函数（返回 1.0 / 0.0 / -1.0） |
| `in_int()` | 从控制台安全读取整数，输入非法时抛出异常 |
| `calc_fac(k)` | 阶乘函数（int / double 重载） |
| `vecnorm1 / vecnorm2 / vecnorminf` | 向量的 1-范数、2-范数、无穷范数 |
| `matnorm1 / matnorminf` | 矩阵的 1-范数、无穷范数 |
| `mat_inverse(A)` | 矩阵求逆（Gauss-Jordan 消元） |
| `mat_multi(A, B)` | 矩阵乘法 |
| `mat_multi_vec(A, a)` | 矩阵-向量乘法 |
| `int_pow<T>(a, x)` | 模板：整数次幂 |
| `StringToNum<T>(str)` | 模板：字符串转数值 |
| `sqr(x) / cube(x)` | 宏：平方 / 立方 |

---

### 3.2 文件日志模块 (filelog)

**文件**: `filelog.h`, `filelog.cpp`

封装文件输出操作，支持 `.txt`（普通文本）、`.m`（MATLAB 脚本）和 `.html`（Plotly.js 可视化网页）三种格式。

**关键特性**：
- 自动检测文件是否存在，询问覆盖或追加
- 通过 `set_precision(eps)` 根据误差限自动设置输出精度
- 通过 `set_mat_size()` / `set_array_len()` 配合重载的 `<<` 运算符输出矩阵和数组
- 15 位有效数字输出精度

**使用示例**：
```cpp
filelog fl;
fl.init("result.txt");         // 打开文件
fl.set_mat_size(3, 3) << A;   // 输出 3x3 矩阵
fl.set_array_len(3) << b;     // 输出长度为 3 的向量
```

---

### 3.3 表达式解析与科学计算器 (formula)

**文件**: `formula.h`, `formula.cpp`

这是项目中最复杂的模块，实现了一个完整的数学表达式解析器和求值引擎。

**架构**：
```
输入字符串 → rf_str() 预处理 → trans_rpn() 逆波兰转换 → f(x) 求值
```

**支持的运算**：
- 基本运算：`+` `-` `*` `/` `^`
- 函数：`sqrt` `exp` `ln` `log` `sin` `cos` `tan` `arcsin` `arccos` `arctan`
- 常量：`e` `pi`
- 绝对值：`|x|`
- 科学计数法：`9.425E+3`
- 隐式乘法：`2x` 自动识别为 `2*x`

**关键函数**：

| 函数 | 说明 |
|------|------|
| `rf_str(str)` | 预处理：清除空格、补全括号、处理 `log` → `ln/log(a)`、处理科学计数法、插入隐式乘号、分词 |
| `trans_rpn(str)` | 将中缀表达式转为逆波兰表示（调度场算法） |
| `f(x)` | 给定 x 值，基于逆波兰序列求值 |
| `f_xnum(x)` | 多变量版本，用于方程组 |
| `matlab_format()` | 转换为 MATLAB 兼容格式 |

**`formulae` 类**：管理多变量函数组，用于非线性方程组的输入和求值。通过 `recog_xnum()` 自动识别 `x1, x2, ...` 等变量。

---

### 3.4 线性方程组基类 (linear_equations)

**文件**: `linear_equations.h`, `linear_equations.cpp`

抽象基类 `Ax_b`，封装了线性方程组 `Ax = b` 的通用数据和方法：

| 成员 | 说明 |
|------|------|
| `m, n` | 矩阵行数、列数 |
| `A, x, b` | 系数矩阵、解向量、右端向量 |
| `A_init()` | 分配内存 |
| `enter_Ab()` | 从控制台输入 A 和 b |
| `exchange_row()` | 交换两行（选主元时使用） |
| `check_symmetry()` | 检查矩阵是否对称 |
| `check_tridiagonal()` | 检查矩阵是否三对角 |
| `construct_sym_tri()` | 构造对称三对角矩阵（用于快速测试） |
| `construct_tri()` | 用三条对角线构造三对角矩阵 |

---

### 3.5 直接法求解线性方程组 (directmethod)

**文件**: `directmethod.h`, `directmethod.cpp`

继承自 `Ax_b`，所有直接法求解器的抽象基类为 `Direct_method`。

| 类名 | 算法 | 适用条件 |
|------|------|---------|
| `Gauss` | 高斯消去法 | 通用（要求主元非零） |
| `CP_Gauss` | 列主元高斯消去法 | 通用（更稳定） |
| `Doolittle` | Doolittle 分解（LU分解） | 方阵，顺序主子式非零 |
| `Cholesky` | Cholesky 分解（平方根法） | 对称正定矩阵 |
| `Improved_sqrt` | 改进平方根法（LDL^T） | 对称矩阵 |
| `Chasing` | 追赶法 | 三对角矩阵 |
| `Givens` | Givens 变换（QR分解） | 通用（含超定方程组最小二乘） |
| `Householder` | Householder 变换（QR分解） | 通用（含超定方程组最小二乘） |

**算法流程**（以 Cholesky 为例）：
1. 检查对称性 → 2. 执行 Cholesky 分解 A = G G^T → 3. 解 G y = b → 4. 解 G^T x = y

**输出**：解向量 x，中间矩阵（L/U/G 等），结果保存至 `Direct_method.txt`。

---

### 3.6 迭代法求解线性方程组 (iterationmethod)

**文件**: `iterationmethod.h`, `iterationmethod.cpp`

继承自 `Ax_b`，抽象基类为 `Iteration_method`。

| 类名 | 算法 | 特点 |
|------|------|------|
| `Jacobi` | Jacobi 迭代 | 简单，收敛慢 |
| `Gauss_Seidel` | Gauss-Seidel 迭代 | 比 Jacobi 快 |
| `SOR` | 逐次超松弛迭代 | 需指定松弛因子 ω |
| `steepest_descent` | 最速下降法 | 要求对称正定 |
| `conjugate_gradient` | 共轭梯度法 | 要求对称正定，收敛快 |
| `GMRES` | 广义极小残量法 | 适用一般矩阵，需指定子空间维度 m |

**公共参数**：
- `eps`：迭代误差限
- `maxcounter`：最大迭代次数
- `itrerr[]`：记录每步迭代误差

**输出**：解向量 x、迭代次数、误差表，自动生成 `Iteration_method.m` 绘制误差收敛曲线。

---

### 3.7 插值法 (interpolation)

**文件**: `interpolation.h`, `interpolation.cpp`

抽象基类为 `Interpolation`，支持从连续函数采样或离散列表输入数据。

| 类名 | 算法 | 说明 |
|------|------|------|
| `Newton_Ip` | Newton 插值 | 支持动态增加插值节点 |
| `Hermite_Ip` | Hermite 插值 | 支持指定各阶导数值（Newton 型） |
| `cube_spline` | 三次样条插值 | 支持三种边界条件 |

**Hermite 插值说明**：
输入时需指定每个点的导数阶数。若有重复节点，程序按导数阶数升序排列。差商表在重节点处使用 Taylor 展开替代。

**三次样条边界条件**：
1. 给定两端二阶导数
2. 给定两端一阶导数
3. 周期边界条件

**输出**：插值多项式（分段）、自动生成 `Interpolation.m` 绘制插值曲线。

---

### 3.8 函数最优逼近 (optimal_approx)

**文件**: `optimal_approx.h`, `optimal_approx.cpp`

抽象基类为 `optimal_approx`。

| 类名 | 方法 | 说明 |
|------|------|------|
| `sqr_approx` | 最佳平方逼近 | 支持连续函数和列表函数；基函数可选递推正交多项式或自定义 |
| `uni_approx` | 近似最佳一致逼近 | 两种子方法：Chebyshev 插值多项式 / 截断 Chebyshev 级数 |

**最佳平方逼近流程**：
1. 选择基函数（递推生成正交多项式或自定义）
2. 计算 Gram 矩阵元素（连续函数用 Romberg 积分，列表函数用求和）
3. 解法方程组得到系数
4. 计算均方误差

**截断 Chebyshev 级数法**：
1. 变量代换 x → cos θ，将 f(x) 变为 f(cos θ)
2. 计算 Chebyshev 系数（Romberg 积分）
3. 截断后反代换回 x 的多项式

**输出**：逼近多项式、均方误差，自动生成 `Optimal_approx.m` 绘图。

---

### 3.9 数值积分与数值微分 (int_diff)

**文件**: `int_diff.h`, `int_diff.cpp`

抽象基类为 `Int_Diff`。

| 类名 | 算法 | 说明 |
|------|------|------|
| `Romberg` | Romberg 积分 | 梯形公式 + Richardson 外推 |
| `DerivExtra` | 外推求导法 | 中心差分 + Richardson 外推，支持任意阶导数 |

**Romberg 积分流程**：
```
T_1 → T_2 → S_1 → T_4 → S_2 → C_1 → T_8 → S_4 → C_2 → R_1 → ...
```
每步判断 |R_k - R_{k-1}| < eps。

**DerivExtra**：
- 一阶导数：中心差分 + 外推
- 高阶导数：递归调用低阶导数 + 中心差分 + 外推
- 支持多变量偏导数计算（`deriv(xx, xnum, j)` 对第 j 个变量求偏导）

**输出**：积分/导数值、误差表，自动生成 `Integration.m` 绘制误差曲线。

---

### 3.10 非线性方程/方程组 (nonlinear_equations)

**文件**: `nonlinear_equations.h`, `nonlinear_equations.cpp`

分为两大类：单变量方程 `nl_eq` 和多变量方程组 `nl_eqs`。

**单变量方程方法**：

| 类名 | 算法 | 说明 |
|------|------|------|
| `eq_Simple` | 简单迭代法 | x_{k+1} = g(x_k) |
| `eq_Newton` | Newton 法 | 使用 DerivExtra 自动求导 |
| `eq_Secant` | 弦割法 | 两点迭代，避免求导 |
| `eq_Relaxation` | 松弛加速法 | 在简单迭代基础上引入松弛因子 |
| `eq_Aitken` | Aitken 加速法 (Steffensen) | 加速简单迭代的收敛 |

**多变量方程组方法**：

| 类名 | 算法 | 说明 |
|------|------|------|
| `eqs_Simple` | 简单迭代法 | x_{k+1} = G(x_k) |
| `eqs_Newton` | Newton 法 | 使用 DerivExtra 计算 Jacobi 矩阵，Householder 法解线性方程组 |
| `eqs_Secant` | 离散 Newton 法 | 差商代替 Jacobi 矩阵 |
| `eqs_Broyden` | Broyden 秩1方法 | 拟 Newton 法，避免重复计算 Jacobi 逆 |

**寻根策略**：
- 输入区间 `[a, b]` 时，自动以 0.01*(b-a) 步长扫描，用二分法隔离所有根
- 输入初值时，从初值开始迭代

---

### 3.11 矩阵特征值与特征向量 (eigen_val_vec)

**文件**: `eigen_val_vec.h`, `eigen_val_vec.cpp`

抽象基类为 `Eigen_val_vec`。

| 类名 | 算法 | 说明 |
|------|------|------|
| `Power_method` | 幂法 | 求按模最大的特征值，支持原点平移加速 |
| `Inverse_power` | 反幂法 | 给定近似特征值，求更精确的特征值和对应特征向量 |

**幂法流程**：
1. 取初始向量 z^(0)
2. y^(k) = A z^(k-1)
3. m_k = max(y^(k))
4. z^(k) = y^(k) / m_k
5. 判断 |m_k - m_{k-1}| < eps 或 ||z^(k) - z^(k-1)|| < eps

**反幂法**：利用 LU 分解求解 (A - pI) y = z，加速收敛到最接近 p 的特征值。

**输出**：特征值 λ、特征向量、迭代次数；结果保存至 `Eigen_val_vec.txt`。

---

## 4. 编译与运行

### 环境要求

- **编译器**：支持 C++11 的 g++（开发环境为 TDM-GCC 4.9.2）
- **编码**：所有源文件为 UTF-8
- **操作系统**：Windows

### 编译方式

若 g++ 已在 PATH 中：

```bash
g++ -std=c++11 -o NumericalCalculation src/*.cpp && rm -f src/*.o
```

若不在 PATH 中，设置 `MINGW` 变量指向 MinGW 安装目录后编译：

```bash
MINGW="/path/to/MinGW64"
$MINGW/bin/g++.exe -std=c++11 -o NumericalCalculation src/*.cpp \
    -I$MINGW/include \
    -I$MINGW/x86_64-w64-mingw32/include \
    -I$MINGW/lib/gcc/x86_64-w64-mingw32/4.9.2/include \
    -I$MINGW/lib/gcc/x86_64-w64-mingw32/4.9.2/include/c++ \
    && rm -f src/*.o
```

> **注意**：部分旧版 g++（如 TDM-GCC 4.9.2）不会自动搜索 include 目录，必须通过 `-I` 显式指定。
```

**方式三：使用 Makefile**

```bash
mingw32-make -f Makefile.win
```

### 切换语言

编辑 `common_fd.h`：
- 保留 `#define CHINESE_VERSION` → 中文界面
- 注释掉 `//#define CHINESE_VERSION` → 英文界面

### 运行

```bash
./NumericalCalculation
```

按照菜单提示选择功能编号，输入对应数据即可。

---

## 5. 使用说明

### 5.0 科学计算器（选项 0）

输入表达式后直接求值。支持 x 变量时可反复代入不同 x 值计算。输入 `000` 退出。

```
示例：
f(x) = x^2 + sin(x)
  x  = 3.14
f(3.14) = 9.8596...
```

### 5.1 直接法 & 5.2 迭代法（选项 1-2）

输入矩阵维数 → 选择是否构造对称三对角矩阵 → 输入矩阵元素 → 选择算法 → 得到解。

### 5.3 插值法（选项 3）

可选择连续函数或离散列表：
- 连续函数：输入表达式，自动采样
- 离散列表：直接输入 (x, y) 点对

### 5.4 最优逼近（选项 4）

- 最佳平方逼近：支持自定义基函数或自动生成正交多项式
- 一致逼近：Chebyshev 插值或截断级数

### 5.5 数值积分/微分（选项 5）

- 积分：输入被积函数和区间，自动 Romberg 积分
- 微分：输入阶数和求导点，外推法求导

### 5.6 非线性方程（选项 6）

- 单方程：输入 f(x) 和初值/区间
- 方程组：输入方程个数和各方程表达式，以及初值向量

### 5.7 特征值（选项 7）

输入矩阵 → 选择幂法或反幂法 → 得到特征值和特征向量。

---

## 6. 类继承关系图

```
Ax_b (linear_equations)
├── Direct_method (directmethod)
│   ├── Gauss
│   ├── CP_Gauss
│   ├── Doolittle
│   ├── Cholesky
│   ├── Improved_sqrt
│   ├── Chasing
│   ├── Givens
│   └── Householder
│
└── Iteration_method (iterationmethod)
    ├── Jacobi
    ├── Gauss_Seidel
    ├── SOR
    ├── steepest_descent
    ├── conjugate_gradient
    └── GMRES

Interpolation (interpolation)
├── Newton_Ip
│   └── Hermite_Ip
└── cube_spline

optimal_approx (optimal_approx)
├── sqr_approx
└── uni_approx

Int_Diff (int_diff)
├── Romberg
└── DerivExtra

nl_eq (nonlinear_equations)
├── eq_Simple
│   ├── eq_Relaxation
│   └── eq_Aitken
├── eq_Newton
└── eq_Secant

nl_eqs (nonlinear_equations)
├── eqs_Simple
├── eqs_Newton
├── eqs_Secant
└── eqs_Broyden

Eigen_val_vec (eigen_val_vec)
├── Power_method
└── Inverse_power

formula       — 表达式解析器
formulae      — 多变量函数组
filelog       — 文件日志输出
```

---

## 附录：数据文件说明

### testdata.txt

存放测试用矩阵数据，是一组 4x4 线性方程组 Ax=b 的 A 和 b：
```
A =                      b =
2.52   0.95   1.25  -0.85    1.38
0.39   1.69  -0.45   0.49   -0.34
0.55  -1.25   1.96  -0.98    0.67
0.23  -1.15  -0.45   2.31    1.52
```

### 自动生成的文件

| 文件 | 生成模块 |
|------|---------|
| `Direct_method.txt` | 直接法求解结果 |
| `Iteration_method.txt` | 迭代法求解结果 |
| `Iteration_method.m/.html` | 迭代误差收敛曲线 |
| `Interpolation.txt` | 插值结果 |
| `Interpolation.m/.html` | 插值函数图像 |
| `Optimal_approx.txt` | 逼近结果 |
| `Optimal_approx.m/.html` | 逼近函数图像 |
| `Int_Diff.txt` | 积分/微分结果 |
| `Integration.m/.html` | 误差收敛曲线 |
| `Nonlinear_equations.txt` | 非线性方程求解结果 |
| `Eigen_val_vec.txt` | 特征值/特征向量结果 |
