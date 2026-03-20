这是一个Python脚本，用于自动化处理qPCR数据与结果可视化。

本脚本数据来源于软件`Bio-Rad CFX Maestro`，导出格式为csv。
导出步骤为：Export -> Custom Export。
设置如下图所示：

<div style="text-align: center;">
<img width="464" height="509" alt="73fc05193c176d5aea19b7894d0119a0" src="https://github.com/user-attachments/assets/a8f6a6e8-2ae8-4563-8642-ee906569745f" />
</div>

在本文件中，Target为目的片段，Sample为样本名。由于本数据来源为虫体实验，因此仅计算 $2^{-\Delta Cq}$ 即可。
本文件数据为一块96孔板的数据，除去内参基因片段，还可扩增2个目的基因片段。每个基因有3个处理，每个处理有3个技术重复。

本数据处理的基本思路为：
1. 计算 $\Delta Cq$
2. 计算 $2^{-\Delta Cq}$
3. 计算对照组的 $2^{-\Delta Cq}$ 的均值，并将对照组、处理组1、处理组2归一化。

本脚本的基本流程为：
1. 数据预处理：导入数据，筛除后续分析所需的数据（`Target`，`Sample`和`Cq`）
2. 数据处理：如上所示
3. 结果导出：导出为2个csv文件，可使用其他软件（如`GraphPad Prism`）进行可视化处理
4. 可视化：进行统计学检验并进行可视化

- Ver260320：该版本为原始版本，能够基本实现功能，但注释尚不完善，且部分代码需要手动修改，且仅进行一个Target的统计学检验和可视化。
