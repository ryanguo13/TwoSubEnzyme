# ========================= convenience.jl =========================
# 通用 convenience kinetics（Liebermeister & Klipp 2006）
# 适用：任意计量学的可逆酶促反应
# ==================================================================

using ModelingToolkit, Symbolics, Latexify, Plots
import Symbolics: value


export convenience_rate, h_act, h_inh

"""
    convenience_rate(species_left, species_right, stoich_left, stoich_right,
                     Km_left, Km_right, kcat_plus, kcat_minus, E_tot)

返回 convenience kinetics 的符号速率表达式。

参数
------
- species_left / right : 元组/向量，物种符号
- stoich_left / right  : 对应计量系数（整数 ≥ 0）
- Km_left / right      : 对应半饱和常数
- kcat_plus / minus    : 宏观正/逆催化常数
- E_tot                : 酶总浓度

示例
------
v = convenience_rate((S,T), (P,Q), (1,1), (1,1),
                     (Ks,Kt), (Kp,Kq), k⁺, k⁻, Etot)
"""
function convenience_rate(species_left, species_right,
                          stoich_left, stoich_right,
                          Km_left, Km_right,
                          kcat_plus, kcat_minus, E_tot)

    # 归一化浓度
    x_left  = [s / k for (s, k) in zip(species_left,  Km_left)]
    x_right = [s / k for (s, k) in zip(species_right, Km_right)]

    # 分子
    numerator = kcat_plus * prod(x_left  .^ stoich_left) -
                kcat_minus * prod(x_right .^ stoich_right)

    # 分母：多项式 (1 + x + ... + x^α)
    poly_left  = [sum([xi^k for k in 0:α]) for (xi, α) in zip(x_left,  stoich_left)]
    poly_right = [sum([xi^l for l in 0:β]) for (xi, β) in zip(x_right, stoich_right)]

    denominator = prod(poly_left) * prod(poly_right) - 1

    return E_tot * numerator / denominator
end

# ---------- 可选：激活 / 抑制因子 ----------
h_act(d, Ka) = 1 + d / Ka
h_inh(d, Ki) = Ki / (Ki + d)

# ==================================================================
# 示例 1：单底物 ↔ 单产物  A ⇌ B
# ==================================================================
@variables t A(t) B(t) Ka Kb k⁺ k⁻ Etot
v1 = convenience_rate((A,), (B,), (1,), (1,), (Ka,), (Kb,), k⁺, k⁻, Etot)

# ==================================================================
# 示例 2：双底物 ↔ 双产物  S + T ⇌ P + Q
# ==================================================================
@variables S(t) T(t) P(t) Q(t) Ks Kt Kp Kq
v2 = convenience_rate((S,T), (P,Q), (1,1), (1,1),
                      (Ks,Kt), (Kp,Kq), k⁺, k⁻, Etot)

# ==================================================================
# 示例 3：2A + B ⇌ 3C
# ==================================================================
@variables A(t) B(t) C(t) Ka Kb Kc
v3 = convenience_rate((A,B), (C,), (2,1), (3,),
                      (Ka,Kb), (Kc,), k⁺, k⁻, Etot)

# ==================================================================
# 一键 LaTeX 输出
# ==================================================================
println("示例 2 速率表达式:")
println(latexify(v2))

# ==================================================================
# 示例：带激活/抑制的双底物反应
# ==================================================================
@variables Act(t) Inh(t) Ka_act Ki_inh
v2_mod = v2 * h_act(Act, Ka_act) * h_inh(Inh, Ki_inh)
println("\n带激活/抑制:")
println(latexify(v2_mod))




# 设定变量的数值范围
S_vals = range(0, 10, length=100)
T_vals = range(0, 10, length=100)

# 固定参数
P_val = 0.0
Q_val = 0.0
Ks_val = 1.0
Kt_val = 1.0
Kp_val = 1.0
Kq_val = 1.0
k⁺_val = 1.0
k⁻_val = 1.0
Etot_val = 1.0

# 构建网格
S_grid = repeat(S_vals', length(T_vals), 1)
T_grid = repeat(T_vals, 1, length(S_vals))

# 计算 v2 在网格上的数值
v2_vals = [Symbolics.value(substitute(v2, Dict(
    S=>Sg, T=>Tg, P=>P_val, Q=>Q_val,
    Ks=>Ks_val, Kt=>Kt_val, Kp=>Kp_val, Kq=>Kq_val,
    k⁺=>k⁺_val, k⁻=>k⁻_val, Etot=>Etot_val
))) for (Sg, Tg) in zip(S_grid[:], T_grid[:])]

v2_mat = reshape(v2_vals, length(T_vals), length(S_vals))

# 绘制 3D 曲面图
surface(S_vals, T_vals, v2_mat,
    xlabel="S", ylabel="T", zlabel="v₂",
    title="Convenience rate v₂ as a function of S and T",
    color=:viridis,
    size=(800, 600)
)
savefig("convenience_rate_3d.svg")