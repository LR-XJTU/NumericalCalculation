# Bug 修复记录

> 记录使用Claude code协助编程期间发现并修复的所有问题。标注了"原始代码"还是"重构引入"。

## 一、原始代码的 Bug（已修复）

### B1. `optimal_approx::inner_product()` — `delete [] fy, gy;` 仅释放了一个数组

- **文件**: `src/optimal_approx.cpp`，`inner_product(formula, formula, double*, int, double*)`
- **严重度**: 高
- **原因**: `delete [] fy, gy;` 中的逗号是**逗号运算符**，解析为 `(delete[] fy), gy;`，`gy` 从未被 `delete[]`
- **修复**: 拆成两行独立的 `delete[] fy; delete[] gy;`

### B2. 9 个基类缺少虚析构函数 — 派生类析构函数永不执行

- **文件**: `src/linear_equations.h`（`Ax_b`）、`src/directmethod.h`（`Direct_method`）、`src/iterationmethod.h`（`Iteration_method`）、`src/interpolation.h`（`Interpolation`）、`src/optimal_approx.h`（`optimal_approx`）、`src/int_diff.h`（`Int_Diff`）、`src/nonlinear_equations.h`（`nl_eq`、`nl_eqs`）、`src/eigen_val_vec.h`（`Eigen_val_vec`）
- **严重度**: 高
- **原因**: `main.cpp` 中全部 8 处 `delete` 都是通过基类指针进行。C++ 标准规定基类无虚析构时行为未定义；实际编译器跳过派生类析构，导致 `Givens::Q`、`Householder::Q/d/alpha`、所有迭代法子类的 `itrerr[]`、插值类的全部数据数组等资源永久泄漏
- **修复**: 每个基类添加 `virtual ~XXX() = default;`

### B3. `cube_spline::~cube_spline()` — 析构函数中误调 `init()`

- **文件**: `src/interpolation.cpp`
- **严重度**: 高（因 B2 虚析构缺失，析构函数本身不会被调用，故未实际触发）
- **原因**: 析构函数第一行 `init()` 会创建 `filelog` 对象（打开 `Interpolation.txt`）并调用 `in_int()` 等待用户输入。显然是复制构造函数时保留了第一行
- **修复**: 删除析构函数中的 `init()` 调用

### B4. `formula::get_dfy()` — 返回了错误的成员变量

- **文件**: `src/formula.cpp`
- **严重度**: 中
- **原因**: 函数体 `return dfx;`，应为 `return dfy;`。`operator=` 中唯一调用处 `dfy = fx.get_dfy();` 实际把对方的 `dfx`（x 坐标值）赋给了 `dfy`（y 值）
- **修复**: `return dfx;` → `return dfy;`

### B5. 科学计数法输出小写 `e` 被 `formula` 误解析为欧拉数

- **文件**: `src/optimal_approx.cpp`、`src/interpolation.cpp`、`src/filelog.cpp`、`src/main.cpp`
- **严重度**: 中
- **原因**: `polytostr()` 用 `stringstream` 输出浮点数时默认产生小写 `e`（如 `7.0689e-007`），这些表达式字符串传给 `formula` 解析时，`rf_str()` 将小写 `e` 识别为自然对数底数 2.71828，求值完全错误。修复分两步：① `polytostr()` 的 `stringstream` 加 `uppercase`，确保生成的公式字符串用大写 `E`（optimal_approx.cpp 5 处、interpolation.cpp 2 处）；② `filelog::init()` 和 `main()` 加 `uppercase`，确保文件输出（.txt）和控制台输出的科学计数法也统一用大写 `E`

### B6. `cube_spline::calc()` — `M_solve[-1]` 越界访问

- **文件**: `src/interpolation.cpp`，case 2 分支
- **严重度**: 高
- **原因**: 三种边界条件中，case 1 和 case 3 的 `M_solve` 存储了 M[1..]（偏移 1），case 2 存储了全部 M[0..np1-1]（偏移 0）。但还原代码从 case 1 复制过来，`M_solve[i-1]` 没改。
- **修复**: `M_solve[i - 1]` → `M_solve[i]`

### B7. `exchange_diag_no_0()` — 负数索引越界 + 外层死循环

- **文件**: `src/iterationmethod.cpp`
- **严重度**: 高
- **原因**: ①内层 while 先访问 `A[i+temp][i]` 再检查 `i+temp < 0`，temp 为负时越界；②外层 while(1) 找到候选行但该行不可用时没有 `temp--`，卡死在同一条
- **修复**: ① `while (i+temp>=0 && A[i+temp][i]==0.0)` 先检查再访问；② 拆开 `||` 为两个独立 `if`，中间加 `temp--`

### B8. `conjugate_gradient::calc()` — 忽略用户设置的 `maxcounter`

- **文件**: `src/iterationmethod.cpp`
- **严重度**: 中
- **原因**: 循环条件硬编码 `while (err > eps && itcounter < 5 * n)`，用户通过 `set_max()` 设置的值被无视
- **修复**: `5 * n` → `maxcounter`

### B9. 10 个类的 `bool` 构造函数 — `else ClassName()` 创建临时对象

- **文件**: `src/directmethod.cpp`（Gauss, CP_Gauss, Doolittle, Cholesky, Improved_sqrt, Chasing, Givens, Householder）、`src/interpolation.cpp`（Newton_Ip）、`src/int_diff.cpp`（DerivExtra）
- **严重度**: 中（潜伏，所有调用方都传 `true`）
- **原因**: `else ClassName();` 不是委托构造，是在栈上创建临时对象然后销毁。`false` 路径完全不工作
- **修复**: 删除死代码分支。构造体简化为 `{}`（Doolittle 保留 `print_flag = false`）

### B10. `Householder(bool)` — `print_flag` 未初始化

- **文件**: `src/directmethod.cpp`
- **严重度**: 中（潜伏，该构造函数无人调用）
- **原因**: 构造函数直接返回，`print_flag` 是垃圾值。原代码 `if(no_init) return;` 同样没赋值
- **修复**: 加 `print_flag = false;`

### B11. `int_diff.cpp` — 两处 `double* e` 内存泄漏

- **文件**: `src/int_diff.cpp`，`Romberg::out_result()` 和 `DerivExtra::out_result()`
- **严重度**: 中
- **原因**: `double* e = new double[err.size()]` 从未被 `delete[]`
- **修复**: 删除临时数组，直接使用 `err`（已是 `vector<double>`）

### B12. 缺少 `#include <algorithm>`

- **文件**: `src/directmethod.cpp`、`src/linear_equations.cpp`
- **严重度**: 低（当前编译器恰好通过传递包含编译）
- **原因**: 使用了 `std::swap` 和 `std::min` 但没有包含对应头文件
- **修复**: 添加 `#include <algorithm>`

### B13. 全项目 `-Wsign-compare` 警告（约 25 处）

- **文件**: `formula.cpp`、`nonlinear_equations.cpp`、`int_diff.cpp` 等
- **严重度**: 低
- **原因**: 循环变量 `int i` 与 `.size()`（返回 `size_t`）比较
- **修复**: 循环变量改为 `size_t`，不能用 `size_t` 的地方加显式转换

### B14. `-Wunused-parameter` 警告

- **文件**: `src/int_diff.cpp`，`Int_Diff(bool)` 构造函数
- **严重度**: 低
- **原因**: 参数 `no_init` 声明但未使用
- **修复**: 省略参数名 → `Int_Diff(bool)`

---

## 二、重构本身引入的问题（已修复）

**无。** 所有 B1-B9 都是原始代码中的问题，重构没有引入新的 bug。

重构中有一处被我过度修复后撤回：

- `conjugate_gradient` 参数构造器的 `fl.init()` — 加上了又撤回。该构造器是内部求解器模式，不需要文件输出，原设计正确。

---

## 三、重构主要变更（非 Bug 修复）


| 变更                 | 详情                                                                    |
| ------------------ | --------------------------------------------------------------------- |
| 动态数组 → vector      | `double*`/`double`** 成员全部改为 `vector<double>`/`vector<vector<double>>` |
| 删除手动析构             | ~11 个仅做 `delete[]` 的析构函数移除，RAII 自动管理                                  |
| `exchange_row` 优化  | 逐元素循环 → `std::swap`（O(n) → O(1)）                                      |
| `add_point` 简化     | 34 行手动重分配 → 6 行 `push_back`/`resize`                                  |
| `resize_itrerr` 简化 | 7 行手动复制 → 1 行 `itrerr.resize(itcounter)`                              |
| 函数签名简化             | `(double* arr, int size)` → `(const vector<double>& arr)`             |
| 原始指针重载清理           | `vecnorm*`、`matnorm*`、`filelog operator<<` 的指针版本删除                    |
| CLAUDE.md 更新       | 内存管理、函数约定、已知问题、编译规范                                                   |


