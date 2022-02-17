<center>
吴子博
</center>
<center>
2020年2月15日
</center>
***
#### Zhu-Nakamura公式与约束分子动力学方法笔记
## 1.非绝热动力学

1. **Ehrenfest 动力学**(Mean-Field方法)：电子在单一有效势能面运动，布局$|a_i| ^2$决定其权重大小
2. **Surface hopping**(系间窜越) 原子核在单一势能面上演化，电子可以发生在不同势能面上的跃迁，从而使得动力学过程在不同的绝热态之间转化。

## 2.轨迹系间窜越(TSH)

轨迹系间窜越是将量子力学与非绝热动力学模拟相结合的半经典方法。模拟过程中，每一时刻轨迹都在单一势能面上演化。但在特殊点，例如势能面的**圆锥相交点(conical intersection, CI)**、避免相交点，轨迹可以从一个势能面跃迁到另一个势能面。跃迁几率取决与电子态之间的耦合。跃迁一般发生在电子态能差较小的地方。

### 2.1 Fewest Switching方法（考虑电子跃迁的分子动力学）

轨线系间窜越方法FS允许跃迁在任何时间发生。体系在电子态之间的随机转化取决于电子态的几率布局随时间的演化，在给定的时刻，任意一条模拟轨迹都是在**单一势能面上运行**，而不在平均有效的势能面上运行。

### 2.2 Landau-Zener方法及其改进的Zhu-Nakamura(ZN)方法

LZ方法主要用于描述非绝热动力学过程在两个透热态之间的跃迁几率
$$
 P_{LZ} =\exp \left[ -\frac{\pi}{4ab}  \right]  \tag{1} 
$$
其中$a^2,b^2$分别为有效耦合参数和有效碰撞能参数 
$$
a^2 = \frac{\hbar ^2}{2\mu}\frac{\sqrt{|F_2 F_1|}|F_2-F_1|}{(2V_{12})^3}
$$
$$
b^2 = (E_t- E_x)\frac{|F_2-F_1|}{\sqrt{F_2F_1}{(2V_{12})}}
$$
其中$F_1$和$F_2$分别为两个透热面上的力，$V_{12}$是非绝热耦合项，$\mu$是双原子的约化质量，$E_x$是交叉点的能量，$E_t$为非绝热跃迁耦合量方向上的动能与势能总和。
朱超原和Nakamura改进了LZ计算跃迁几率：
$$
P_{ZN} =\exp \left[ -\frac{\pi}{4\sqrt{a^2}} \sqrt{\frac{2}{b^2 + \sqrt{|b^4 \pm 1|}}} \right]  \tag {2}
$$
ZN方法与其他TSH方法的主要区别在于如何求非绝热解跃迁几率，在最初的ZN方法中，对于非绝热跃迁几率的求解涉及到了两个绝热态的能量，梯度和跃迁方向。跃迁方向可以近似通过计算**两个态Hessian矩阵差**或**非绝热耦合因子**确定。在ZN方法中，跃迁方向可以直接通过**数值差分**方法得到。
在新的ZN方法中，$a^2,b^2$被重新定义:
$$
a^2 = \frac{\hbar}{2m}\frac{(F_2-F_1)^2}{8V_X^3} \\
 b^2 = (E-E_x)\frac{1}{2V_X}
$$
$F=\sqrt{|F_1F_2|},F_1F_2 \not ={0}$,其中$F_j(j=1,2),E_x和V_x$分别是透热面势能斜率,势能面交叉点的能量和透热耦合量 。$E_x和V_x$可以重新定义成：
$$
E_x =\frac{U_+(q^2)+U_{-}(q^2)}{2} \\
V_x = \frac{U_+(q^2)-U_{-}(q^2)}{2}
$$
其中$U_+(q^2)和U_-(q^2)$分别为较高态的势能和较低态的势能在避免相交点处的值。
对于交点处有N个原子的体系，
$$
\frac{(F_2-F_1)^2}{m} = \sum_{i=1}^{3N} \frac{(F_i^2-F_i^1)^2}{m_i}  \\
 F_i^1(\mathbf{q}) = \frac{1}{q_i^3 -q_i^1}\left[   \frac{\partial U_-}{\partial q_i^3}(q_i -q_i^1)  -\frac{\partial U_+}{\partial q_i^1
 }(q_i -q_i^3)\right] \\
 F_i^2(\mathbf{q}) = \frac{1}{q_i^3 -q_i^1}\left[   \frac{\partial U_+}{\partial q_i^3}(q_i -q_i^1)  -\frac{\partial U_-}{\partial q_i^1
 }(q_i -q_i^3)\right]
$$
其中$\frac{\partial U+}{\partial q_i^1}和 \frac{\partial U_+}{\partial q_i^3}$ 为原子1和3 在构型$q^1和q^3$处较高的电子态密度，上述公式可以进一步写成：
$$
\frac{1}{m}(F_2-F_1)^2 = \sum_{i=1}^{3N}\frac{1}{m_i}[F_i^2(\mathbf{q}^2) -F_i^1(\mathbf{q}^2)][F_i^2(\mathbf{q}^2) -F_i^1(\mathbf{q}^2)]
$$
定义一个新向量 ***S*** 与$\frac{1}{\sqrt{m}}(F_2-F_1)$有相同方向的单位向量：
$$
\mathbf{s_i} = \frac{[F_i^2(\mathbf{q}^2) -F_i^1(\mathbf{q}^2)]\frac{1}{\sqrt{m}}}{\sqrt{(F_2-F_1)^2 \frac{1}{m}}}
$$
因此 ZN方法中每一个原子的跃迁方向可以定义为$\mathbf{n_i}=\frac{\mathbf{s_i}}{|\mathbf{s_i|}}$。 在轨迹模拟中，整个过程应该保持能量守恒即：
$$
U_+(\mathbf{q}^2)+\sum_{i=1}^{N}\frac{\mathbf{P}_i^2(+)}{2m_i} = U_-(\mathbf{q}^2)+\sum_{i=1}^{N}\frac{\mathbf{P}_i^2(-)}{2m_i}
$$
在电子跃迁过程中，沿跃迁方向的分量$\mathbf{P}_{i\parallel}$重新分配，而垂直于跃迁方向的分量$\mathbf{P}_{i\perp}$保持不变，其中满足$\mathbf{P}_{i\parallel} =(\mathbf{P}_i \cdot \mathbf{n}_i)\mathbf{n}_i ,\mathbf{P}_i = \mathbf{P}_{i\parallel}+\mathbf{P}_{i\perp}$和$\mathbf{P}_i = \mathbf{P}_{i\parallel}^2+\mathbf{P}_{i\perp}^2$。
对于从较高电子态到较低电子态的跃迁有:
$$
U_+(\mathbf{q}^2)-U_-(\mathbf{q}^2)+\sum_{i=1}^{N}\frac{\mathbf{P}_{i\parallel}^2(+)}{2m_i} = \sum_{i=1}^{N}\frac{\mathbf{P}_{i\parallel}^2(-)}{2m_i}
$$
其中动力分量$\mathbf{P}_{i\parallel},\mathbf{P}_{i\perp}$满足:
$$
\mathbf{P}_{i\parallel}(-) = k\mathbf{P}_{i\parallel}(+)=k(\mathbf{P}_i(+)\cdot\mathbf{n}_i)\mathbf{n}_i\\ 
\mathbf{P}_{i\perp}(-) =\mathbf{P}_{i\perp}(+)=\mathbf{P}_i(+)-(\mathbf{P}_i(+)\cdot\mathbf{n}_i)\mathbf{n}_i
$$
其中增量参数 *k*满足：
$$
k=\sqrt{1 + \frac{U_+(\mathbf{q})^2- U_-(\mathbf{q})^2}{\sum_{i=1}^N\frac{\mathbf{P}_{i\parallel}^2(+)}{2m_i}}} > 1
$$
因此:
$$
\mathbf{P}_i(-) = \mathbf{P}_{i\parallel}(-)+\mathbf{P}_{i\parallel}(+)=\mathbf{P}_i(+)+(k-1)(\mathbf{P}_i(+)\cdot\mathbf{n}_i)\mathbf{n}_i \tag{3}
$$
同理，从较低电子态跃迁到较高电子态后有:
$$
k=\sqrt{1-\frac{U_+(\mathbf{q})^2- U_-(\mathbf{q})^2}{\sum_{i=1}^N\frac{\mathbf{P}_{i\parallel}^2(-)}{2m_i}}} <1
$$
因此动量$\mathbf{P}_i(+)$可以表示为:
$$
\mathbf{P}_i(+) = \mathbf{P}_{i\parallel}(-)+\mathbf{P}_{i\parallel}(+)=\mathbf{P}_i(-)+(k-1)(\mathbf{P}_i(-)\cdot\mathbf{n}_i)\mathbf{n}_i \tag{4}
$$
采用新的ZN方法模拟动力学过程中，电子跃迁可以只发生在所涉及到的两个电子态的绝热能量间隔的局域最小点。可以根据 **公式(2)** 计算非绝热跃迁几率，当此几率在[0,1]之间时跃迁发生。其动量$\mathbf{P}_i$可以根据 **公式(3)** 和 **公式(4)** 重新分配。其后轨迹停留在跃迁后的绝热势能面直到另一个跃迁发生或动力学模拟结束.

***