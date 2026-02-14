# 🚀 AFDM Communication System Simulation

![MATLAB](https://img.shields.io/badge/MATLAB-R2023b%2B-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg) 本项目是一个基于 MATLAB 的 **仿射频分复用 (AFDM, Affine Frequency Division Multiplexing)** 通信系统链路级仿真平台。AFDM 是一种针对高移动性、双色散（时变多径）信道设计的新型波形，能够通过调整 Chirp 参数获得全分集增益。

本工程采用了高度封装的面向对象设计 (`classdef`)，集成了从信道编码、调制、信道估计到高级均衡算法的完整物理层收发链路，特别侧重于在深衰落信道下的性能评估。

## ✨ 核心特性 (Key Features)

* **面向对象架构**: 核心逻辑封装于 `AfdmSystem` 类中，参数配置集中、易于扩展和维护。
* **双色散信道模型 (LTV Channel)**: 支持自定义多径时延和多普勒频移的线性时变瑞利衰落信道。
* **先进信道估计**: 包含基于导频的信道估计，支持分数阶多普勒精搜索 (Fractional Doppler Search) 与动态加窗去噪。
* **高性能均衡算法对比**:
  * **MMSE (最小均方误差)**: 稳健的线性均衡器。
  * **Weighted MRC-DFE (加权最大比合并判决反馈均衡器)**: 针对双色散信道优化的非线性均衡器，有效抑制符号间干扰。

## 🛠️ 环境依赖 (Prerequisites)

* **MATLAB R2023b 或更高版本**
* **Communications Toolbox**

## 📁 项目结构 (Project Structure)

```text
├── AfdmSystem.m          % AFDM 核心收发机类 (包含编码、调制、估计、均衡)
├── LtvChannel.m          % 物理层线性时变 (LTV) 信道生成函数
├── main.m                % 仿真主脚本 (用于配置参数并绘制 BER vs SNR 曲线)
└── README.md             % 项目说明文档
