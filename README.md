# Low-rank components and Graph Estimation (LGE)
This code implements the low-rank components and graph estimation (LGE) algorithm described in the paper: 

Liu, R., Nejati, H., & Cheung, N. M. (2018). Joint Estimation of Low-Rank Components and Connectivity Graph in High-Dimensional Graph Signals: Application to Brain Imaging. 

[arXiv preprint: https://arxiv.org/abs/1801.02303]

## Acknowledgement
We reuse some code for “Learning Laplacian Matrix in Smooth Graph Signal Representation “, originally written by Xiaowei Dong, for a function in our implementation. 
We reuse some function in the GSP toolbox for initial graph generation. 
We also use CVX toolbox for optimisation.

## Usage
You need to have GSP toolbox (https://epfl-lts2.github.io/gspbox-html/) and CVX toolbox (http://cvxr.com/cvx/download/) installed. Other functions are already included. Once you have everything in place, you should be able to run “Demo.m” and get results of LGE running on synthetic data.

## Contact
If you have any questions, please contact: rui_liu@mymail.sutd.edu.sg.
Thank you for your interests in our work.
