# 二维 Elfouhaily 海面局部卷浪

本模块在二维 Elfouhaily 随机海面上自动选择具备破碎几何特征的波峰，并对局部结构网格施加有限支撑的参数化卷曲。它不使用一维高斯波横向复制，也不把卷曲结果重新压缩成单值高度场。

## 卷浪位置判定

对轻度平滑后的海面，沿主传播方向计算一阶和二阶导数。候选点必须同时满足：

```text
z_s >= q_0.88(z_s)
abs(z_u) <= 0.32 std(z_u)
z_uu <= -q_0.70(max(-z_uu,0))
```

评分函数联合考虑高程、负曲率和前向肩部坡度：

```text
Q = 0.50 Q_height + 0.35 Q_curvature + 0.15 Q_forward_slope.
```

程序在联合候选中选择最高评分事件，然后逐个波峰向截面跟踪真实脊线 `u_c(v)`。这套条件用于放置 breaking-eligible event，是静态几何代理，不是普适的流体破碎起始判据。

## 局部卷曲

卷曲前先通过宽过渡窗轻微改变周围非卷曲表面，使后肩、波峰和前坡连续演化。卷曲核心不再只修改高度，而是以原网格传播向材料坐标 `q=u-u_c(v)` 构造多值参数剖面。它包含：

1. 后肩到上唇的平滑前移，且唇顶略低于原始最高点；
2. 半椭圆参数鼻端，使表面先向前、再向下并向后折返；
3. 低位下支继续向前延伸，随后平滑回接原 Elfouhaily 前坡。

鼻端水平半径、前伸量和下落量均由检测到的局部波高 `Hlocal` 缩放。它保留 `Ideal_Curl_Wave_Echo.m` 的“材料点越过波峰后在传播向-垂向平面形成回转支路”思想，但删除固定大支点深度和单角旋转。沿波峰方向使用紧支撑余弦窗，因此卷曲不会横向复制成柱状面。卷曲结果保存为 `(X,Y,Z)` 参数曲面，三角面连接沿用原结构网格。

几何翻卷由传播方向映射雅可比判断：

```text
J_u = d(u_final)/du < 0.
```

## 默认运行结果

固定随机种子 `20260825` 的结果为：

- 选中波峰：`(9.100, 26.800, 0.230) m`；
- 估计局部波高：`0.1225 m`；
- 最大前向位移：`0.1890 m`，即 `1.543 Hlocal`；
- 最大向下位移：`0.0741 m`，即 `0.605 Hlocal`；
- 最大向上位移：`0.0047 m`，即 `0.038 Hlocal`；
- 最小传播雅可比：`-0.7447`；
- `J_u<0` 的翻卷网格点：`93` 个，且全部位于波峰前侧；
- 紧支撑区域外位移：数值零。

## 文件

- `default_elfouhaily_ideal_curl_config.m`：海面、判定和卷曲参数；
- `generate_elfouhaily_ideal_curl_surface.m`：海面生成、候选判定、脊线跟踪和参数曲面卷曲；
- `run_elfouhaily_ideal_curl_demo.m`：运行断言并分别绘制五张图。

## 输出图

- `01_breaking_location_conditions.png`：高度、近零坡度、负曲率以及联合候选与最终选址；
- `02_elfouhaily_selected_crest.png`：原始海面和选中波峰；
- `03_elfouhaily_local_curled_surface.png`：完整卷曲海面；
- `04_local_curl_center_section.png`：卷曲前后中心参数剖面；
- `05_local_curl_closeup.png`：按材料坐标着色的局部前卷三维侧视图；
- `elfouhaily_local_curl.mat`：可选完整结构体数据，设置 `cfg.output.saveMatFile=true` 时生成。

## 运行

```matlab
run_elfouhaily_ideal_curl_demo
```
