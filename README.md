# ALARM: Average Likelihood for Attack-Resilient Multi-Object Filtering


This repository contains an implementation of the **Average Likelihood for Attack-Resilient Multi-Object (ALARM) filtering** method, as proposed in:

> V. Ghorbani, A. K. Gostar, Z. Tari, N. Sohrabi, A. Ghorbani & R. Hoseinnezhad,  
> “Robust Filtering for Multi-Object Tracking Against Stealthy Measurement-Oriented Adversarial Attacks,” *Authorea Preprints*, 2025.  


The ALARM filter is built on top of the RFS Tracking Toolbox provided by Ba-Tuong Vo:  
https://ba-tuong.vo-au.com/codes.html

---

## Features

- Implements the ALARM filtering algorithm for robust multi-object tracking under adversarial measurement attacks  
- Built upon the well-known RFS Tracking Toolbox 
- MATLAB-based scripts with reproducible results  

## Usage

1. **Open MATLAB** and change directory to the project root.  
2. **Ensure required directories are on the MATLAB path**:
   ```matlab
   addpath('_common');
   addpath('_network');
3. **Run the main.m script**:
