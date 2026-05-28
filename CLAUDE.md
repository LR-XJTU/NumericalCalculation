# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 项目概述

数值计算方法的 C++ 控制台程序，涵盖计算方法课程核心内容：线性方程组（直接法/迭代法）、插值、逼近、数值积分微分、非线性方程（组）、特征值。

## 构建

使用 TDM-GCC 4.9.2（随 Dev-C++ 安装在 `D:/Dev-Cpp/MinGW64/`），不在 PATH 中，且该版本不会自动搜索 include 目录。

编译时使用以下命令（`$MINGW` 指向 MinGW64 安装根目录）：

```bash
MINGW="D:/Dev-Cpp/MinGW64"
$MINGW/bin/g++.exe -std=c++11 -o NumericalCalculation src/*.cpp \
    -I$MINGW/include \
    -I$MINGW/x86_64-w64-mingw32/include \
    -I$MINGW/lib/gcc/x86_64-w64-mingw32/4.9.2/include \
    -I$MINGW/lib/gcc/x86_64-w64-mingw32/4.9.2/include/c++ \
    && rm -f src/*.o
```

## 编码

**所有 `.cpp` 和 `.h` 文件是 UTF-8 编码**。Dev-C++ 5.11 编辑器不支持 UTF-8，中文在 Dev-C++ 中显示为乱码属正常现象——**编辑请用 Cursor 或 VS Code，Dev-C++ 仅用于编译**。

编译时 `-finput-charset=UTF-8` 可加可不加（不用 `-fexec-charset` 时 g++ 直接透传字节，不影响结果）。`main.cpp` 启动时通过 `system("chcp 65001 > nul")` 将控制台切换为 UTF-8。

## 双语机制

`src/common_fd.h` 中的 `#define CHINESE_VERSION` 控制全局语言切换：

- **定义** → 中文界面
- **注释** → 英文界面

所有 12 个 `.cpp` 文件的用户交互文本均通过 `#ifdef CHINESE_VERSION` / `#else` / `#endif` 包裹。修改时注意——只有用户交互文本（`cout << "..."` 和 `fl << "..."` 中的描述文字）需要双语包裹，文件名、数学符号、MATLAB 生成代码、纯变量数据输出不需要。

## Web 可视化

`src/common_fd.h` 中的 `#define WEB_VISUALIZATION` 控制可视化输出格式：

- **定义** → 生成 `.html` 文件（Plotly.js，浏览器打开）
- **注释** → 生成 `.m` 文件（MATLAB 脚本）

与 `CHINESE_VERSION` 同理，切换后必须 `make clean && make`。覆盖模块：迭代法、插值（牛顿/样条）、数值积分、最佳逼近（平方/一致）。

## 架构

### 类继承体系

所有模块采用**策略模式**：基类定义 `calc()` / `out_result()` 虚函数接口，子类实现具体算法。`main.cpp` 通过 `switch` 创建子类，通过基类指针调用并 `delete`。

```
Ax_b (linear_equations)
├── Direct_method (directmethod) — 8种直接法
│   ├── Gauss, CP_Gauss, Doolittle, Cholesky
│   ├── Improved_sqrt, Chasing, Givens, Householder
└── Iteration_method (iterationmethod) — 6种迭代法
    ├── Jacobi, Gauss_Seidel, SOR, steepest_descent
    ├── conjugate_gradient, GMRES

Interpolation — Newton_Ip → Hermite_Ip, cube_spline

optimal_approx — sqr_approx, uni_approx

Int_Diff — Romberg, DerivExtra

nl_eq  — eq_Simple → eq_Relaxation, eq_Aitken, eq_Newton, eq_Secant
nl_eqs — eqs_Simple, eqs_Newton, eqs_Secant, eqs_Broyden

Eigen_val_vec — Power_method, Inverse_power
```

**重要**：所有基类已有 `virtual ~XXX() = default;`，子类析构函数能正确调用。新增基类时务必加虚析构。

### 核心工具类

- **`formula`** — 数学表达式解析器（逆波兰求值）。流程：`rf_str()` 预处理 → `trans_rpn()` 调度场算法 → `f(x)` 求值。支持隐式乘法、科学计数法（仅大写 `E`，小写 `e` 是自然对数底数）、`log` 底数转换。`formulae` 是多变量版本，变量名为 `x1, x2, ...`。
- **`filelog`** — 文件日志，重载 `<<` 输出矩阵/向量（接受 `const std::vector<double>&` 和 `const std::vector<std::vector<double>>&`），支持 `.txt`、`.m`（MATLAB 脚本）、`.html`（Web 可视化）。配合 `set_mat_size()` / `set_array_len()` 使用。
- **`Ax_b`** — 线性方程组基类，封装矩阵/向量（`std::vector<std::vector<double>> A`、`std::vector<double> x`、`std::vector<double> b`）、输入/输出、对称性/三对角检测、行交换（O(1) `std::swap`）等。

### 内存管理

全部使用 `std::vector` / `std::vector<std::vector<double>>`，**无 `new[]`/`delete[]`**。RAII 自动管理，无需手动释放。所有仅做 `delete[]` 的析构函数已删除，编译器自动生成。

### 函数约定

- 矩阵参数用 `const std::vector<std::vector<double>>&`
- 向量参数用 `const std::vector<double>&` 或 `std::vector<double>&`（输出参数）
- 范数函数 `vecnorm*`、`matnorm*` 只接受 vector 版本（无原始指针重载）
- `Ax_b::get_x(std::vector<double>&)` 直接赋值（`xx = x`），无需传大小

## 已知问题

- **`bool no_init` 构造函数**（`directmethod.cpp` 8 个类 + `interpolation.cpp` Newton_Ip）：`else ClassName();` 创建临时对象而非委托构造，`no_init=false` 路径不工作。当前所有调用都传 `true`，暂时潜伏。不要新增 `false` 调用。
- **`cube_spline::calc()`** `bdcd_flag == 2` 时 `M_solve[-1]` 越界，需对照算法修正。
- **`exchange_diag_no_0()`** 负数索引检查在访问之后，极罕见触发。

## 注意事项

- 改了 `src/common_fd.h` 中的宏定义后必须重新编译所有文件
- 自动生成的输出文件（`.txt`、`.m`、`.html`）在项目根目录，已在 `.gitignore` 中
- `testdata.txt` 是 4×4 线性方程组的测试数据（A 矩阵 + b 向量），用于快速验证线性方程组求解
- 构建 `resultstr`（用于 `formula` 解析的表达式字符串）时，`stringstream` 需加 `uppercase`，确保科学计数法输出大写 `E`
- 所有 `.cpp` 和 `.h` 使用 `std::vector`，**不要引入 `new[]`/`delete[]`**。需要动态数组时用 `vector::resize()`
